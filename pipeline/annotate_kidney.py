#!/usr/bin/env python3
"""
Aegis — Kidney annotator
scBaseCount kidney samples ship UNLABELED and there is no CellTypist kidney
model, so we annotate by canonical marker signatures at the cluster level:
Leiden clusters -> assign each cluster the kidney cell type whose marker
signature scores highest. Writes labeled h5ads (with raw_counts +
log_normalized layers + obs.cell_type) that build_summary can consume directly.

Usage:
  python annotate_kidney.py --download --max 14        # fetch more + annotate all
  python annotate_kidney.py                            # annotate whatever is in raw/kidney
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, logging
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("aegis.kidney")

AEGIS = Path(__file__).resolve().parent.parent
PROJECT = AEGIS.parent
META = PROJECT / "data" / "scbasecount" / "sample_metadata.parquet"
RAW = AEGIS / "raw" / "kidney"
ANNOT = AEGIS / "annot" / "kidney"

# Canonical human kidney marker signatures (cluster-level assignment by argmax)
SIGS = {
    "proximal tubule": ["LRP2","CUBN","SLC34A1","SLC5A12","SLC13A3","GATM","MIOX",
                        "SLC22A6","SLC22A8","ANPEP","PCK1","SLC5A2","SLC3A1"],
    "thick ascending limb": ["UMOD","SLC12A1","CLDN16","CASR","KCNJ1","ENOX1"],
    "distal convoluted tubule": ["SLC12A3","TRPM6","CALB1","KLHL3","PVALB","TRPV5"],
    "principal cell (CD)": ["AQP2","AQP3","SCNN1G","SCNN1B","GATA3","FXYD4","HSD11B2"],
    "intercalated cell": ["SLC4A1","ATP6V0D2","ATP6V1B1","SLC26A4","FOXI1","KIT","DMRT2"],
    "podocyte": ["NPHS1","NPHS2","PODXL","PTPRO","WT1","PLA2R1","SYNPO"],
    "parietal epithelial cell": ["CLDN1","PAX8","CFH","VCAM1","ALDH1A2"],
    "endothelial cell": ["PECAM1","FLT1","EMCN","KDR","CD34","EHD3","PLVAP","CLDN5"],
    "mesangial / vSMC": ["PDGFRB","ACTA2","MYH11","REN","ITGA8","PIEZO2","NOTCH3"],
    "fibroblast": ["PDGFRA","COL1A1","COL1A2","DCN","LUM","C7","MEG3"],
    "T cell": ["CD3D","CD3E","CD3G","IL7R","CD8A","CD2","CCL5"],
    "NK cell": ["NKG7","GNLY","KLRD1","NCAM1","KLRF1"],
    "myeloid": ["LYZ","CD68","C1QA","C1QB","CD14","FCGR3A","ITGAM","AIF1"],
    "B / plasma cell": ["MS4A1","CD79A","CD79B","IGHG1","MZB1","JCHAIN"],
    "mast cell": ["TPSAB1","CPA3","MS4A2","KIT"],
}

KIDNEY_RX = r"\bkidney\b|\brenal\b|nephron|glomerul|renal cortex|renal medulla"
BAD = (r"adrenal|organoid|cell line|cultured|ipsc|-derived|sorted|enrich|proximal tubule cell|"
       r"tumou?r|carcinoma|cancer|transplant|nephrectomy|rejection|allograft|"
       r"lupus|nephropathy|nephritis|glomerulonephritis|fetal|embryo")


def fetch(maxn):
    import gcsfs
    df = pd.read_parquet(META)
    df = df[df["organism"].fillna("").str.contains("Homo sapiens|human", case=False, na=False)]
    df["tissue"] = df["tissue"].fillna(""); df["disease"] = df["disease"].fillna("")
    m = df[df["tissue"].str.contains(KIDNEY_RX, case=False, regex=True, na=False)].copy()
    m = m[~m["tissue"].str.contains(BAD, case=False, regex=True, na=False)]
    m = m[~m["disease"].str.contains(BAD, case=False, regex=True, na=False)]
    if "cell_line" in m.columns:
        cl = m["cell_line"].fillna("").astype(str).str.strip().str.lower()
        m = m[cl.isin({"","nan","none","na","unsure","unknown","not applicable","other"})]
    m = m[(m["obs_count"] >= 3000) & (m["obs_count"] <= 90000)]
    m = m.sort_values("obs_count", ascending=False)
    key = "entrez_id" if "entrez_id" in m.columns else "srx_accession"
    m = m.drop_duplicates(subset=[key]).head(maxn)
    fs = gcsfs.GCSFileSystem(token="anon"); RAW.mkdir(parents=True, exist_ok=True)
    for _, r in m.iterrows():
        out = RAW / f"{r['srx_accession']}.h5ad"
        if out.exists(): continue
        try:
            fs.get(r["file_path"], str(out))
            log.info(f"  ↓ {out.name} ({out.stat().st_size/1e6:.0f} MB, {int(r['obs_count'])} cells)")
        except Exception as e:
            log.warning(f"  FAILED {r['srx_accession']}: {e}")


def annotate_one(path):
    a = sc.read_h5ad(path)
    if "gene_symbols" in a.var.columns and str(a.var_names[0]).upper().startswith("ENSG"):
        a.var_names = a.var["gene_symbols"].astype(str).values
    a.var_names_make_unique()
    a.var["mt"] = a.var_names.str.upper().str.startswith(("MT-","MT."))
    sc.pp.calculate_qc_metrics(a, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)
    sc.pp.filter_cells(a, min_genes=200)
    a = a[a.obs["n_genes_by_counts"] < 8000, :].copy()
    a = a[a.obs["pct_counts_mt"] < 50, :].copy()   # kidney PT is mito-rich → lenient
    sc.pp.filter_genes(a, min_cells=10)
    if a.n_obs < 200:
        return None
    a.layers["raw_counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4); sc.pp.log1p(a)
    a.layers["log_normalized"] = a.X.copy()
    sc.pp.highly_variable_genes(a, n_top_genes=2000, flavor="seurat_v3", layer="raw_counts")
    ah = a[:, a.var["highly_variable"]].copy(); sc.pp.scale(ah, max_value=10)
    sc.tl.pca(ah, n_comps=40); a.obsm["X_pca"] = ah.obsm["X_pca"]
    sc.pp.neighbors(a, n_pcs=30); sc.tl.leiden(a, resolution=1.0)
    # score each signature (genes present)
    scores = {}
    for name, genes in SIGS.items():
        g = [x for x in genes if x in a.var_names]
        if len(g) < 3:
            scores[name] = np.full(a.n_obs, -np.inf); continue
        sc.tl.score_genes(a, g, score_name="_s")
        scores[name] = a.obs["_s"].values.copy()
    S = pd.DataFrame(scores, index=a.obs_names)
    # assign per leiden cluster by mean signature score
    lab = pd.Series(index=a.obs_names, dtype=object)
    for cl, idx in a.obs.groupby("leiden").groups.items():
        means = S.loc[idx].mean()
        lab.loc[idx] = means.idxmax()
    a.obs["cell_type"] = lab.values
    a.obs["cell_type_scbasecount"] = lab.values  # so build_summary prefers this
    if "_s" in a.obs: del a.obs["_s"]
    return a


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--max", type=int, default=14)
    args = ap.parse_args()
    if args.download:
        log.info("Fetching kidney samples…"); fetch(args.max)
    ANNOT.mkdir(parents=True, exist_ok=True)
    files = sorted(RAW.glob("*.h5ad"))
    log.info(f"Annotating {len(files)} kidney samples")
    for f in files:
        out = ANNOT / f.name
        if out.exists():
            log.info(f"  {f.stem}: already annotated"); continue
        try:
            a = annotate_one(f)
        except Exception as e:
            log.warning(f"  {f.stem}: FAILED ({type(e).__name__}: {e})"); continue
        if a is None:
            log.warning(f"  {f.stem}: too few cells"); continue
        vc = a.obs["cell_type"].value_counts()
        top = ", ".join(f"{k}:{v}" for k, v in vc.head(6).items())
        log.info(f"  {f.stem}: {a.n_obs} cells, {vc.size} types | {top}")
        a.write_h5ad(out)
    log.info("Done. Annotated files in annot/kidney/")


if __name__ == "__main__":
    main()
