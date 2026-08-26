#!/usr/bin/env python3
"""
Aegis — Tumor pipeline (CNV-refined epithelial malignant calling)
For a cancer type: select + download tumor scRNA samples, classify each cell
into a compartment (immune / stromal / epithelial), run infercnvpy CNV using
immune cells as the diploid reference, and label as MALIGNANT only the
epithelial cells that sit in aneuploid (high-CNV) clusters. TME is retained but
flagged, so downstream aggregation can use malignant-only signal.

Usage:
  python tumor_pipeline.py --cancer crc --n 4 --download --dry-run
  python tumor_pipeline.py --cancer crc --n 4 --download
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, logging
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp
import infercnvpy as cnv

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("aegis.tumor")

AEGIS = Path(__file__).resolve().parent.parent
PROJECT = AEGIS.parent
META = PROJECT / "data" / "scbasecount" / "sample_metadata.parquet"
POS = AEGIS / "refs" / "gene_positions.parquet"
RAW = AEGIS / "raw_tumor"
ANNOT = AEGIS / "annot_tumor"

CANCERS = {
    "breast":   r"breast (cancer|carcinoma|adenocarc)|\bBRCA\b|mammary carc|triple.?negative|\bTNBC\b",
    "breast_her2": r"(breast|mammary|brca).{0,45}her2\s*\+|her2\s*\+.{0,45}(breast|mammary|cancer)|erbb2[- ]?amplif",
    "luad":     r"lung adenocarc|\bLUAD\b|non-small cell|\bNSCLC\b|lung carcinoma",
    "crc":      r"colorectal|colon (adeno)?carc|rectal (adeno)?carc|\bCRC\b|colon cancer|colonic adeno",
    "hcc":      r"hepatocellular|\bHCC\b|\bLIHC\b",
    "pdac":     r"pancreatic (ductal )?adeno|\bPDAC\b|pancreatic cancer|pancreatic carcinoma",
    "ovarian":  r"ovarian (cancer|carcinoma)|\bHGSOC\b|high.?grade serous|ovarian adeno",
    "gastric":  r"gastric (cancer|carcinoma|adeno)|stomach (cancer|adeno)",
}
CANCER_LABEL = {"breast":"breast tumor","breast_her2":"breast HER2+ tumor","luad":"lung adeno tumor",
    "crc":"colorectal tumor","hcc":"HCC tumor","pdac":"pancreatic tumor","ovarian":"ovarian tumor",
    "gastric":"gastric tumor"}
# genes with strong single-nucleus dropout (membrane/cytoplasmic) — prefer single-cell tumors
BAD = (r"organoid|cell line|cultured|ipsc|-derived|\bpdx\b|xenograft|spheroid|sorted|"
       r"epcam\+|cd45|flow.?sort|facs|nucle(i|us)")

import re
# EPITHELIAL checked FIRST so tumor epithelial subtypes win; word boundaries on
# short immune tokens (t cell / b cell) so they don't match goblet/tuft "…t cell".
EPI_RX = re.compile(
    r"epitheli|carcinoma|tumou?r|malignan|hepatocyte|cholangiocyte|enterocyte|goblet|"
    r"colonocyte|ductal|acinar|alveolar|\bat1\b|\bat2\b|club cell|basal cell|secretory|"
    r"luminal|keratinocyte|squamous|urothelial|foveolar|paneth|tuft|enteroendocrine|"
    r"ciliated|ionocyte|serous cell|mucous|chief cell|parietal cell|hillock|"
    r"pneumocyte|glandular|adenocarc", re.I)
IMMUNE_RX = re.compile(
    r"\bt cell|\bb cell|\bnk\b|\bnkt\b|natural killer|myeloid|macrophage|monocyte|dendritic|"
    r"\bdc\b|mast cell|plasma cell|neutrophil|granulocyte|lymphocyte|microglia|kupffer|"
    r"langerhans|\bcd4|\bcd8|treg|regulatory t|\bilc\b|megakaryocyte|erythro|hematopoietic|"
    r"leukocyte|thymocyte", re.I)
STROMAL_RX = re.compile(
    r"fibroblast|endothelial|pericyte|smooth muscle|mesenchym|\bstroma|myofibro|stellate|"
    r"mesothelial|adipocyte|\bvsmc\b|schwann", re.I)


def compartment(ct):
    c = str(ct).lower()
    if EPI_RX.search(c): return "epithelial"
    if IMMUNE_RX.search(c): return "immune"
    if STROMAL_RX.search(c): return "stromal"
    return "other"


def select(cancer, n, lo=3000, hi=120000):
    df = pd.read_parquet(META)
    df = df[df["organism"].fillna("").str.contains("Homo sapiens|human", case=False, na=False)]
    df["disease"] = df["disease"].fillna(""); df["tissue"] = df["tissue"].fillna("")
    m = df[df["disease"].str.contains(CANCERS[cancer], case=False, regex=True, na=False)].copy()
    m = m[~m["tissue"].str.contains(BAD, case=False, regex=True, na=False)]
    m = m[~m["disease"].str.contains(BAD, case=False, regex=True, na=False)]
    if "cell_line" in m.columns:
        cl = m["cell_line"].fillna("").astype(str).str.strip().str.lower()
        m = m[cl.isin({"","nan","none","na","unsure","unknown","not applicable","other"})]
    # prefer single-cell over single-nucleus (snRNA drops membrane antigens e.g. CEACAM5)
    if "cell_prep" in m.columns:
        cp = m["cell_prep"].fillna("").astype(str).str.lower()
        sc_only = m[cp.str.contains("cell") & ~cp.str.contains("nucl")]
        if len(sc_only) >= 3:
            m = sc_only
    m = m[(m["obs_count"] >= lo) & (m["obs_count"] <= hi)]
    m = m.sort_values("obs_count", ascending=False)
    key = "entrez_id" if "entrez_id" in m.columns else "srx_accession"
    return m.drop_duplicates(subset=[key]).head(n)


def load_positions():
    p = pd.read_parquet(POS)
    by_ens = p.dropna(subset=["ensembl_id"]).drop_duplicates("ensembl_id").set_index("ensembl_id")
    by_sym = p.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")
    return by_ens, by_sym


def annotate_sample(path, by_ens, by_sym):
    a = sc.read_h5ad(path)
    ens = a.var["ensembl_id"].astype(str).values if "ensembl_id" in a.var.columns else a.var_names.astype(str).values
    if "gene_symbols" in a.var.columns and str(a.var_names[0]).upper().startswith("ENSG"):
        a.var_names = a.var["gene_symbols"].astype(str).values
    a.var_names_make_unique()
    # QC
    a.var["mt"] = a.var_names.str.upper().str.startswith(("MT-","MT."))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(a, min_genes=200)
    a = a[a.obs["n_genes_by_counts"] < 9000, :].copy()
    a = a[a.obs["pct_counts_mt"] < 30, :].copy()
    sc.pp.filter_genes(a, min_cells=10)
    if a.n_obs < 300:
        return None, "too few cells"
    # cell types
    ctcol = "cell_type_scbasecount" if "cell_type_scbasecount" in a.obs.columns else "cell_type"
    if ctcol not in a.obs.columns:
        return None, "no cell_type"
    cts = a.obs[ctcol].astype(str).str.strip()
    if (cts.isin(["","nan","None"]).mean()) > 0.5:
        return None, "unlabeled (>50%)"
    a.obs["cell_type"] = cts.values
    a.obs["compartment"] = [compartment(x) for x in cts.values]
    comp = a.obs["compartment"].value_counts().to_dict()
    n_imm = comp.get("immune", 0); n_epi = comp.get("epithelial", 0)
    if n_epi < 50:
        return None, f"few epithelial ({n_epi})"
    # store raw + lognorm
    a.layers["raw_counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    a.layers["log_normalized"] = a.X.copy()

    # genomic positions on var (map by ensembl id first, else symbol)
    var = a.var
    ens_ser = pd.Series(ens[:var.shape[0]], index=var.index) if len(ens) >= var.shape[0] else None
    chrom = pd.Series(index=var.index, dtype=object); st = pd.Series(index=var.index, dtype=float); en = pd.Series(index=var.index, dtype=float)
    if ens_ser is not None:
        hit = ens_ser.isin(by_ens.index)
        chrom.loc[hit] = by_ens.loc[ens_ser[hit].values, "chromosome"].values
        st.loc[hit] = by_ens.loc[ens_ser[hit].values, "start"].values
        en.loc[hit] = by_ens.loc[ens_ser[hit].values, "end"].values
    miss = chrom.isna()
    sym = pd.Series(var.index, index=var.index)
    hit2 = miss & sym.isin(by_sym.index)
    chrom.loc[hit2] = by_sym.loc[sym[hit2].values, "chromosome"].values
    st.loc[hit2] = by_sym.loc[sym[hit2].values, "start"].values
    en.loc[hit2] = by_sym.loc[sym[hit2].values, "end"].values
    a.var["chromosome"] = chrom.values; a.var["start"] = st.values; a.var["end"] = en.values
    n_pos = int(chrom.notna().sum())

    epi_mask = a.obs["compartment"].values == "epithelial"
    fallback = n_imm < 50 or n_pos < 3000
    if not fallback:
        try:
            cnv.tl.infercnv(a, reference_key="compartment", reference_cat=["immune"],
                            window_size=100, step=10)
            cnv.tl.pca(a, n_comps=20)
            cnv.pp.neighbors(a)
            cnv.tl.leiden(a, resolution=1.0)
            cnv.tl.cnv_score(a)  # obs['cnv_score'] = per-cnv_leiden mean CNV burden
            cl = a.obs["cnv_leiden"].astype(str).values
            comp = a.obs["compartment"].values
            cdf = pd.DataFrame({"cl": cl, "comp": comp, "score": a.obs["cnv_score"].values})
            g = cdf.groupby("cl")
            frac_epi = g["comp"].apply(lambda s: (s == "epithelial").mean())
            frac_imm = g["comp"].apply(lambda s: (s == "immune").mean())
            cscore = g["score"].mean()
            # baseline CNV burden = highest among immune-dominated (normal) clusters
            base = float(cscore[frac_imm > 0.5].max()) if (frac_imm > 0.5).any() else float(cscore.median())
            malig_cl = set(cscore.index[(frac_epi > 0.5) & (cscore > base)])
            cand = epi_mask & pd.Series(cl).isin(malig_cl).values
            if cand.sum() >= 0.15 * max(epi_mask.sum(), 1):
                a.obs["malignant"] = cand; method = "cnv-refined"
            else:  # CNV didn't cleanly separate a tumor subclone → keep epithelial compartment
                a.obs["malignant"] = epi_mask; method = "epithelial(cnv-inconclusive)"
        except Exception as e:
            log.warning(f"    CNV failed ({type(e).__name__}: {e}); epithelial fallback")
            a.obs["malignant"] = epi_mask; method = f"epithelial(cnv-error:{type(e).__name__})"
    else:
        a.obs["malignant"] = epi_mask
        method = f"epithelial(imm={n_imm},pos={n_pos})"

    a.obs["cell_type_final"] = np.where(a.obs["malignant"].values, "malignant",
                                        a.obs["cell_type"].values)
    return a, method


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cancer", required=True, choices=list(CANCERS))
    ap.add_argument("--n", type=int, default=5)
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    sel = select(args.cancer, args.n)
    print(f"\n=== {args.cancer} — {len(sel)} samples ===")
    for _, r in sel.iterrows():
        print(f"  {r['srx_accession']:>14} {int(r['obs_count']):>7} cells  {str(r['disease'])[:44]}")
    if args.dry_run:
        print("[dry run]"); return

    if not POS.exists():
        raise SystemExit("Run build_gene_positions.py first.")
    by_ens, by_sym = load_positions()
    rawd = RAW / args.cancer; rawd.mkdir(parents=True, exist_ok=True)
    annd = ANNOT / args.cancer; annd.mkdir(parents=True, exist_ok=True)

    if args.download:
        import gcsfs; fs = gcsfs.GCSFileSystem(token="anon")
        for _, r in sel.iterrows():
            out = rawd / f"{r['srx_accession']}.h5ad"
            if out.exists(): continue
            try:
                fs.get(r["file_path"], str(out)); log.info(f"  ↓ {out.name} ({out.stat().st_size/1e6:.0f} MB)")
            except Exception as e:
                log.warning(f"  FAILED {r['srx_accession']}: {e}")

    for f in sorted(rawd.glob("*.h5ad")):
        out = annd / f.name
        if out.exists():
            log.info(f"  {f.stem}: already done"); continue
        try:
            a, method = annotate_sample(f, by_ens, by_sym)
        except Exception as e:
            log.warning(f"  {f.stem}: FAILED ({type(e).__name__}: {e})"); continue
        if a is None:
            log.warning(f"  {f.stem}: skipped ({method})"); continue
        nmal = int(a.obs["malignant"].sum()); nepi = int((a.obs["compartment"]=="epithelial").sum())
        log.info(f"  {f.stem}: {a.n_obs} cells | malignant {nmal}/{nepi} epi | {method}")
        a.write_h5ad(out)
    log.info(f"Done {args.cancer}. Annotated → {annd}")


if __name__ == "__main__":
    main()
