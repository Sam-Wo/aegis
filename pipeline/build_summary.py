#!/usr/bin/env python3
"""
Aegis — Summary Builder
=======================
Distils single-cell healthy-tissue data into a compact per-(tissue, cell_type)
summary that powers the Aegis app.

For every surface protein it records, per cell type in every tissue:
  - mean expression (log-normalised)
  - % of cells expressing (raw count > 0)
And per cell type (gene-independent):
  - proliferation score (Tirosh S + G2M signature)
  - n_cells

This is a LIGHTWEIGHT pass — no PCA/UMAP/Leiden. It reads either raw
scBaseCount h5ads or already-annotated ones, and never writes big files.

Input:  a manifest CSV with columns: path,tissue,tier   (tier optional)
Output: data/aegis_summary.parquet   (long: tissue,cell_type,gene,mean_expr,pct_pos,n_cells)
        data/aegis_groups.parquet    (tissue,cell_type,n_cells,n_samples,prolif,tier,...)
        data/aegis_data.json         (compact, embedded-ready payload for the app)

Usage:
  python build_summary.py --manifest manifests/local.csv
  python build_summary.py --manifest manifests/all.csv --min-cells 20
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, json, logging, time
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("aegis")

HERE = Path(__file__).resolve().parent
AEGIS = HERE.parent
PROJECT = AEGIS.parent           # the parent ADC project (for surface list / data)
OUT = AEGIS / "data"
OUT.mkdir(parents=True, exist_ok=True)

# ── Surface protein list ──
def find_surface_csv():
    for name in ("master_surface_proteins_final.csv",
                 "master_surface_proteins_filtered.csv",
                 "master_surface_proteins.csv"):
        p = PROJECT / name
        if p.exists():
            return p
    # also allow a local copy inside aegis
    p = AEGIS / "data" / "surface_proteins.csv"
    return p if p.exists() else None

# ── Tissue safety tiers (dose-limiting-toxicity framework) ──
TISSUE_TIERS = {
    "bone_marrow": (1, "CRITICAL"), "heart": (1, "CRITICAL"), "liver": (1, "CRITICAL"),
    "kidney": (2, "SERIOUS"), "lung": (2, "SERIOUS"), "colon": (2, "SERIOUS"),
    "small_intestine": (2, "SERIOUS"), "esophagus": (2, "SERIOUS"), "eye": (2, "SERIOUS"),
    "stomach": (2, "SERIOUS"), "pancreas": (2, "SERIOUS"),
    "skin": (3, "MANAGEABLE"), "spleen": (3, "MANAGEABLE"), "thymus": (3, "MANAGEABLE"),
    "prostate": (4, "LOW"), "adipose": (4, "LOW"), "muscle": (4, "LOW"),
    "testis": (4, "LOW"), "vasculature": (4, "LOW"),
}
DEFAULT_TIER = (3, "OTHER")

# ── Cell-cycle signatures (Tirosh et al. 2016; standard scanpy set) ──
S_GENES = ("MCM5 PCNA TYMS FEN1 MCM2 MCM4 RRM1 UNG GINS2 MCM6 CDCA7 DTL PRIM1 "
    "HELLS RFC2 RPA2 NASP RAD51AP1 GMNN WDR76 SLBP CCNE2 UBR7 POLD3 MSH2 ATAD2 "
    "RAD51 RRM2 CDC45 CDC6 EXO1 TIPIN DSCC1 BLM CASP8AP2 USP1 CLSPN POLA1 CHAF1B "
    "BRIP1 E2F8").split()
G2M_GENES = ("HMGB2 CDK1 NUSAP1 UBE2C BIRC5 TPX2 TOP2A NDC80 CKS2 NUF2 CKS1B MKI67 "
    "TMPO CENPF TACC3 FAM64A SMC4 CCNB2 CKAP2L CKAP2 AURKB BUB1 KIF11 ANP32E TUBB4B "
    "GTSE1 KIF20B HJURP CDCA3 HN1 CDC20 TTK CDC25C KIF2C RANGAP1 NCAPD2 DLGAP5 "
    "CDCA2 CDCA8 ECT2 KIF23 HMMR AURKA PSRC1 ANLN LBR CKAP5 CENPE CTCF NEK2 G2E3 "
    "GAS2L3 CBX5 CENPA").split()


def load_surface_genes():
    p = find_surface_csv()
    if p is None:
        log.warning("No surface protein list found — using ALL genes.")
        return None, None
    df = pd.read_csv(p)
    genes = df["gene_name"].astype(str).tolist()
    # keep useful annotation columns for the app
    keep = [c for c in ("gene_name", "category", "confidence_tier",
                        "adc_clinical_status") if c in df.columns]
    ann = df[keep].drop_duplicates("gene_name").set_index("gene_name")
    log.info(f"Surface proteins: {len(genes)} ({p.name})")
    return genes, ann


def to_symbols(adata):
    """Ensure var_names are gene symbols (scBaseCount stores Ensembl in var_names)."""
    if "gene_symbols" in adata.var.columns:
        # var_names may already be symbols (annotated files) or Ensembl (raw)
        looks_ensembl = str(adata.var_names[0]).upper().startswith("ENSG")
        if looks_ensembl:
            adata.var_names = adata.var["gene_symbols"].astype(str).values
    adata.var_names_make_unique()
    return adata


def get_celltypes(adata):
    """Coalesced per-cell labels: native scBaseCount ontology first (captures
    parenchyma), then celltypist, then any cell_type, else 'unlabeled'."""
    order = ("cell_type_final", "cell_type_scbasecount", "cell_type_celltypist", "cell_type")
    have = [c for c in order if c in adata.obs.columns]
    if not have:
        return None
    out = pd.Series(["unlabeled"] * adata.n_obs, index=adata.obs.index, dtype=object)
    filled = pd.Series([False] * adata.n_obs, index=adata.obs.index)
    for c in have:
        vals = adata.obs[c].astype(str).str.strip()
        good = ~filled & vals.ne("") & ~vals.isin(["nan", "None", "unannotated", "unknown"])
        out[good] = vals[good]
        filled |= good
    return out


def get_lognorm_and_raw(adata):
    """Return (lognorm_X, raw_X) as sparse CSC, computing if needed."""
    # log-normalised expression
    if "log_normalized" in adata.layers:
        logX = adata.layers["log_normalized"]
    else:
        X = adata.X
        mx = (X[:200].toarray().max() if sp.issparse(X) else X[:200].max()) if adata.n_obs else 0
        if mx > 50:  # looks like raw counts
            adata.layers["_tmp_raw"] = X.copy()
            sc.pp.normalize_total(adata, target_sum=1e4)
            sc.pp.log1p(adata)
            logX = adata.X
        else:
            logX = X  # already normalised
    # raw counts (for detection / % positive)
    if "raw_counts" in adata.layers:
        rawX = adata.layers["raw_counts"]
    elif "_tmp_raw" in adata.layers:
        rawX = adata.layers["_tmp_raw"]
    else:
        rawX = logX  # log1p(0)==0, so >0 detection is identical to raw>0
    logX = sp.csc_matrix(logX)
    rawX = sp.csc_matrix(rawX)
    return logX, rawX


def light_qc(adata):
    adata.var["mt"] = adata.var_names.str.upper().str.startswith(("MT-", "MT."))
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None,
                               log1p=False, inplace=True)
    sc.pp.filter_cells(adata, min_genes=200)
    adata = adata[adata.obs["n_genes_by_counts"] < 8000, :].copy()
    adata = adata[adata.obs["pct_counts_mt"] < 25, :].copy()
    return adata


def proliferation_score(adata):
    """Per-cell proliferation = S_score + G2M_score on log-normalised data."""
    s = [g for g in S_GENES if g in adata.var_names]
    g2m = [g for g in G2M_GENES if g in adata.var_names]
    if len(s) < 5 or len(g2m) < 5:
        return np.zeros(adata.n_obs, dtype="float32"), np.zeros(adata.n_obs, dtype=bool)
    try:
        sc.tl.score_genes_cell_cycle(adata, s_genes=s, g2m_genes=g2m)
        prolif = (adata.obs["S_score"].values + adata.obs["G2M_score"].values).astype("float32")
        cycling = (adata.obs["phase"].values != "G1")
        return prolif, cycling
    except Exception as e:
        log.warning(f"  cell-cycle scoring failed: {e}")
        return np.zeros(adata.n_obs, dtype="float32"), np.zeros(adata.n_obs, dtype=bool)


class Accumulator:
    """Running per-(tissue, cell_type) sums so multiple samples merge correctly."""
    def __init__(self, genes):
        self.genes = genes
        self.gidx = {g: i for i, g in enumerate(genes)}
        self.n = len(genes)
        self.sum_expr = defaultdict(lambda: np.zeros(self.n, dtype="float64"))
        self.n_pos = defaultdict(lambda: np.zeros(self.n, dtype="float64"))
        self.n_cells = defaultdict(float)
        self.sum_prolif = defaultdict(float)
        self.n_cycling = defaultdict(float)
        self.samples = defaultdict(set)

    def add(self, tissue, ct, sample, logX_sub, rawX_sub, prolif, cycling, target_cols):
        """logX_sub/rawX_sub are (cells x len(target_cols)); scatter into width-N arrays."""
        key = (tissue, ct)
        self.sum_expr[key][target_cols] += np.asarray(logX_sub.sum(axis=0)).ravel()
        self.n_pos[key][target_cols] += np.asarray((rawX_sub > 0).sum(axis=0)).ravel()
        self.n_cells[key] += logX_sub.shape[0]
        self.sum_prolif[key] += float(prolif.sum())
        self.n_cycling[key] += float(cycling.sum())
        self.samples[key].add(sample)

    def to_frames(self, min_cells, tier_of, gene_order, kind_of=None):
        long_rows, group_rows = [], []
        for key in sorted(self.n_cells, key=lambda k: (k[0], k[1])):
            tissue, ct = key
            nc = self.n_cells[key]
            if nc < min_cells:
                continue
            kind = kind_of(tissue) if kind_of else "normal"
            tier, tlabel = (0, "TUMOR") if kind == "tumor" else tier_of(tissue)
            mean_expr = self.sum_expr[key] / nc
            pct_pos = self.n_pos[key] / nc
            group_rows.append(dict(
                tissue=tissue, cell_type=ct, tier=tier, tier_label=tlabel, kind=kind,
                n_cells=int(nc), n_samples=len(self.samples[key]),
                prolif=round(self.sum_prolif[key] / nc, 4),
                pct_cycling=round(self.n_cycling[key] / nc, 4),
            ))
            for g in gene_order:
                i = self.gidx[g]
                me, pp = mean_expr[i], pct_pos[i]
                if pp <= 0:      # skip genes undetected in this group (keeps output sparse)
                    continue
                long_rows.append(dict(
                    tissue=tissue, cell_type=ct, gene=g,
                    mean_expr=round(float(me), 4), pct_pos=round(float(pp), 4),
                    n_cells=int(nc),
                ))
        return pd.DataFrame(long_rows), pd.DataFrame(group_rows)


def build_json(long_df, group_df, gene_ann):
    """Compact payload keyed by gene for instant client-side lookup."""
    group_df = group_df.reset_index(drop=True)
    gkey = {(r.tissue, r.cell_type): i for i, r in group_df.iterrows()}
    groups = [dict(t=r.tissue, ct=r.cell_type, tier=int(r.tier), tl=r.tier_label,
                   kind=getattr(r, "kind", "normal"),
                   n=int(r.n_cells), ns=int(r.n_samples), pr=float(r.prolif),
                   cyc=float(r.pct_cycling))
              for r in group_df.itertuples(index=False)]
    expr = defaultdict(list)
    for r in long_df.itertuples(index=False):
        gi = gkey.get((r.tissue, r.cell_type))
        if gi is None:
            continue
        expr[r.gene].append([gi, r.mean_expr, r.pct_pos])
    genes = sorted(expr.keys())
    ann = {}
    if gene_ann is not None:
        for g in genes:
            if g in gene_ann.index:
                row = gene_ann.loc[g]
                ann[g] = {k: (None if pd.isna(row[k]) else str(row[k]))
                          for k in gene_ann.columns}
    norm = group_df[group_df.get("kind", "normal") != "tumor"] if "kind" in group_df else group_df
    tum = group_df[group_df.get("kind") == "tumor"] if "kind" in group_df else group_df.iloc[0:0]
    return dict(
        meta=dict(n_groups=len(groups), n_genes=len(genes),
                  total_cells=int(group_df.n_cells.sum()),
                  tissues=sorted(norm.tissue.unique().tolist()),
                  tumors=sorted(tum.tissue.unique().tolist())),
        groups=groups, genes=genes, expr=expr, ann=ann,
    )


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--manifest", required=True, help="CSV: path,tissue[,tier]")
    ap.add_argument("--min-cells", type=int, default=20,
                    help="Drop (tissue,cell_type) groups below this many cells")
    ap.add_argument("--all-genes", action="store_true",
                    help="Aggregate all genes instead of surface proteins only")
    ap.add_argument("--max-cells", type=int, default=0,
                    help="Subsample each sample to at most N cells (0 = all) — speeds the build")
    args = ap.parse_args()

    t0 = time.time()
    man = pd.read_csv(args.manifest)
    log.info(f"Manifest: {len(man)} samples")

    surface_genes, gene_ann = (None, None) if args.all_genes else load_surface_genes()

    acc = None
    gene_universe = None  # fixed gene order across samples (surface genes present anywhere)
    tissue_kind = {}      # tissue -> 'normal' | 'tumor'

    for i, row in man.iterrows():
        path = Path(row["path"])
        if not path.is_absolute():
            path = (PROJECT / path)
        tissue = str(row["tissue"])
        tissue_kind[tissue] = str(row.get("kind", "normal"))
        if not path.exists():
            log.warning(f"[{i+1}/{len(man)}] MISSING {path}")
            continue
        log.info(f"[{i+1}/{len(man)}] {tissue}: {path.name}")
        try:
            adata = sc.read_h5ad(path)
        except Exception as e:
            log.warning(f"  failed to read: {e}")
            continue
        adata = to_symbols(adata)
        if get_celltypes(adata) is None:
            log.warning("  no cell_type column — skipping")
            continue
        adata = light_qc(adata)
        if adata.n_obs < args.min_cells:
            log.warning(f"  only {adata.n_obs} cells after QC — skipping")
            continue
        if args.max_cells and adata.n_obs > args.max_cells:
            sc.pp.subsample(adata, n_obs=args.max_cells, random_state=0)

        prolif, cycling = proliferation_score(adata)
        logX, rawX = get_lognorm_and_raw(adata)

        # restrict to surface genes present
        if surface_genes is not None:
            present = [g for g in surface_genes if g in set(adata.var_names)]
        else:
            present = list(adata.var_names)
        if gene_universe is None:
            gene_universe = present if surface_genes is None else surface_genes
            acc = Accumulator(gene_universe)
        # sub-matrices restricted to the accumulator's gene order (present in this sample)
        var_pos = {g: j for j, g in enumerate(adata.var_names)}
        gpos = {g: k for k, g in enumerate(acc.genes)}
        keep_genes = [g for g in acc.genes if g in var_pos]
        src_cols = [var_pos[g] for g in keep_genes]        # cols in this adata
        target_cols = np.array([gpos[g] for g in keep_genes])  # cols in accumulator
        logX_s = logX[:, src_cols].tocsr()
        rawX_s = rawX[:, src_cols].tocsr()

        is_tumor = tissue_kind.get(tissue) == "tumor"
        if is_tumor:
            comp = adata.obs["compartment"].astype(str).values if "compartment" in adata.obs.columns else None
            mal = adata.obs["malignant"].values if "malignant" in adata.obs.columns else None
            if comp is not None and mal is not None:
                # malignant + adjacent-normal epithelium (same-sample batch-controlled ref); drop TME
                cts = np.where(mal, "malignant",
                       np.where(comp == "epithelial", "adjacent normal epithelium", "__skip__"))
            else:
                cf = get_celltypes(adata).values
                cts = np.where(cf == "malignant", "malignant", "__skip__")
        else:
            cts = get_celltypes(adata).values
        for ct in pd.unique(cts):
            if ct == "__skip__":
                continue
            m = np.where(cts == ct)[0]
            if len(m) == 0:
                continue
            acc.add(tissue, ct, path.stem,
                    logX_s[m], rawX_s[m], prolif[m], cycling[m], target_cols)
        # HER2-amplified malignant subset as its own group (for ERBB2 stratification)
        if is_tumor and "her2_amp" in adata.obs.columns and "malignant" in adata.obs.columns:
            ha = np.where(adata.obs["her2_amp"].values & adata.obs["malignant"].values)[0]
            if len(ha) > 0:
                acc.add(tissue, "malignant HER2-amp", path.stem,
                        logX_s[ha], rawX_s[ha], prolif[ha], cycling[ha], target_cols)
        log.info(f"  {adata.n_obs} cells, {len(pd.unique(cts))} cell types, "
                 f"{len(present)} surface genes present")
        del adata, logX, rawX, logX_s, rawX_s

    if acc is None:
        log.error("No data accumulated.")
        return

    def tier_of(t): return TISSUE_TIERS.get(t, DEFAULT_TIER)
    def kind_of(t): return tissue_kind.get(t, "normal")
    long_df, group_df = acc.to_frames(args.min_cells, tier_of, acc.genes, kind_of)
    log.info(f"Groups (>= {args.min_cells} cells): {len(group_df)}")
    log.info(f"Long rows (gene x group, detected): {len(long_df)}")

    long_df.to_parquet(OUT / "aegis_summary.parquet", index=False)
    group_df.to_parquet(OUT / "aegis_groups.parquet", index=False)
    payload = build_json(long_df, group_df, gene_ann)
    with open(OUT / "aegis_data.json", "w") as f:
        json.dump(payload, f, separators=(",", ":"))
    sz = (OUT / "aegis_data.json").stat().st_size / 1e6
    log.info(f"Wrote data/aegis_data.json ({sz:.1f} MB), "
             f"{payload['meta']['n_genes']} genes, {payload['meta']['n_groups']} groups, "
             f"{payload['meta']['total_cells']:,} cells")
    log.info(f"Done in {time.time()-t0:.0f}s")


if __name__ == "__main__":
    main()
