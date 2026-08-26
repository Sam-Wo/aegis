#!/usr/bin/env python3
"""
Aegis — HER2 (ERBB2) amplification stratification from infercnvpy CNV.
Reads each tumor annot file's obsm['X_cnv'] (genomic-window CNV), locates the
ERBB2 locus window (chr17q12) via the gene ordering infercnvpy used, and flags
malignant cells with focal amplification (CNV above the immune diploid ref).
Writes obs['erbb2_cnv'] and obs['her2_amp'] back into the file.

  python add_her2_amp.py --validate    # correlation check only
  python add_her2_amp.py               # annotate all tumor files in place
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, logging
from pathlib import Path
import numpy as np, pandas as pd, scanpy as sc, scipy.sparse as sp

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("aegis.her2")
AEGIS = Path(__file__).resolve().parent.parent
ANNOT_TUMOR = AEGIS / "annot_tumor"
POS = AEGIS / "refs" / "gene_positions.parquet"
WINDOW, STEP = 100, 10          # must match tumor_pipeline infercnv call
ERBB2 = ("ERBB2", "chr17")


def erbb2_cnv_per_cell(a, by_ens, by_sym):
    """Mean X_cnv over windows covering ERBB2, per cell. None if unmappable."""
    if "X_cnv" not in a.obsm or "cnv" not in a.uns:
        return None
    chr_pos = dict(a.uns["cnv"]["chr_pos"])
    if ERBB2[1] not in chr_pos:
        return None
    # genes infercnvpy used = var genes WITH a genomic position, sorted by (chrom,start)
    var = a.var
    chrom = var.get("chromosome"); start = var.get("start")
    if chrom is None or start is None:
        return None
    gp = pd.DataFrame({"gene": var.index, "chrom": chrom.values, "start": pd.to_numeric(start.values, errors="coerce")})
    gp = gp.dropna(subset=["chrom", "start"])
    gp = gp[gp["chrom"].astype(str).str.startswith("chr")]
    gp = gp.sort_values(["chrom", "start"])
    chrs = sorted(gp["chrom"].unique(), key=lambda c: list(chr_pos).index(c) if c in chr_pos else 999)
    g17 = gp[gp["chrom"] == ERBB2[1]].reset_index(drop=True)
    if ERBB2[0] not in set(g17["gene"]):
        return None
    r = int(g17.index[g17["gene"] == ERBB2[0]][0])        # rank of ERBB2 among chr17 genes
    n17 = len(g17)
    w17 = (chr_pos.get("chr18", a.obsm["X_cnv"].shape[1]) if "chr18" in chr_pos
           else a.obsm["X_cnv"].shape[1]) - chr_pos["chr17"]
    if w17 <= 0:
        return None
    i_center = int(round((r - WINDOW / 2) / STEP))
    i_center = max(0, min(w17 - 1, i_center))
    lo = max(0, i_center - 2); hi = min(w17, i_center + 3)
    cols = chr_pos["chr17"] + np.arange(lo, hi)
    X = a.obsm["X_cnv"]
    sub = X[:, cols].toarray() if sp.issparse(X) else np.asarray(X[:, cols])
    return sub.mean(axis=1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--validate", action="store_true")
    args = ap.parse_args()
    pos = pd.read_parquet(POS)
    by_ens = pos.dropna(subset=["ensembl_id"]).drop_duplicates("ensembl_id").set_index("ensembl_id")
    by_sym = pos.dropna(subset=["gene"]).drop_duplicates("gene").set_index("gene")

    files = sorted(ANNOT_TUMOR.glob("*/*.h5ad"))
    for f in files:
        a = sc.read_h5ad(f)
        ecnv = erbb2_cnv_per_cell(a, by_ens, by_sym)
        if ecnv is None:
            log.info(f"  {f.parent.name}/{f.stem}: ERBB2 locus unmappable — skip"); continue
        mal = a.obs["malignant"].values
        imm = a.obs["compartment"].values == "immune"
        # amplification threshold from immune (diploid) reference
        ref = ecnv[imm]
        thr = float(ref.mean() + 1.5 * ref.std(ddof=0)) if imm.sum() >= 30 else float(np.quantile(ecnv, 0.8))
        her2_amp = mal & (ecnv > thr)
        a.obs["erbb2_cnv"] = ecnv.astype("float32")
        a.obs["her2_amp"] = her2_amp
        if args.validate:
            # sanity: ERBB2-locus CNV should correlate with ERBB2 expression in malignant cells
            if "ERBB2" in a.var_names and mal.sum() > 50:
                log_ = a.layers.get("log_normalized", a.X)
                j = list(a.var_names).index("ERBB2")
                expr = (log_[:, j].toarray().ravel() if sp.issparse(log_) else np.asarray(log_[:, j]).ravel())
                m = mal
                cc = np.corrcoef(ecnv[m], expr[m])[0, 1] if m.sum() > 2 else np.nan
                log.info(f"  {f.parent.name}/{f.stem}: n_malig={int(mal.sum())} "
                         f"her2_amp={int(her2_amp.sum())} ({100*her2_amp.mean():.1f}% of all) "
                         f"corr(cnv,expr)={cc:.2f}")
        else:
            a.write_h5ad(f)
            log.info(f"  {f.parent.name}/{f.stem}: her2_amp {int(her2_amp.sum())}/{int(mal.sum())} malignant")
    log.info("validate done" if args.validate else "her2_amp written to all tumor files")


if __name__ == "__main__":
    main()
