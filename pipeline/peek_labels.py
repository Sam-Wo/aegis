#!/usr/bin/env python3
"""
Remotely peek at scBaseCount samples' native cell_type coverage WITHOUT
downloading whole files — HDF5 over gcsfs reads only the obs/cell_type arrays.
Lets us find natively-labeled samples (esp. for tissues like kidney that came
back unlabeled) before committing to multi-GB downloads.

Usage:
  python peek_labels.py --tissue kidney --max 30
  python peek_labels.py --tissue kidney --max 30 --write labeled_kidney.csv
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, logging
from pathlib import Path
import numpy as np, pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("aegis.peek")

AEGIS = Path(__file__).resolve().parent.parent
PROJECT = AEGIS.parent
META = PROJECT / "data" / "scbasecount" / "sample_metadata.parquet"

TISSUE_RX = {
    "kidney": r"\bkidney\b|renal|nephron",
    "heart": r"\bheart\b|myocard|cardiac",
    "lung": r"\blung\b|pulmonary|alveolar",
    "stomach": r"\bstomach\b|gastric",
}
BAD_TISSUE = (r"pbmc|periph|\bblood\b|bone marrow|organoid|cell line|cultured|ipsc|"
    r"lavage|sorted|enrich|mononuclear|-derived|\bislet|lymph node|fetal|embryo|"
    r"tumou?r|carcinoma|cancer|myxoma|adenoma|metasta|biopsy of|xenograft")
BAD_DISEASE = (r"cancer|carcinoma|tumou?r|neoplas|malign|adenoma|leukem|lymphoma|myeloma|"
    r"sarcoma|melanoma|metasta|cirrho|fibrosis|nephropathy|nephritis|glomerulo|"
    r"diabet|syndrome|\bdisease\b|infect|covid|sepsis|failure|injury|transplant|"
    r"rejection|scleros|dysplasia|carcinom")


def peek(fs, gcs_path):
    """Return (n_obs, labeled_fraction, top_types_list) reading only obs/cell_type."""
    import h5py
    with fs.open(gcs_path, "rb") as fh:
        with h5py.File(fh, "r") as f:
            if "obs" not in f or "cell_type" not in f["obs"]:
                return None
            node = f["obs"]["cell_type"]
            if isinstance(node, h5py.Group):  # categorical
                cats = [c.decode() if isinstance(c, bytes) else str(c)
                        for c in node["categories"][:]]
                codes = node["codes"][:]
            else:                              # plain string dataset
                raw = node[:]
                vals = [c.decode() if isinstance(c, bytes) else str(c) for c in raw]
                cats = sorted(set(vals))
                cmap = {c: i for i, c in enumerate(cats)}
                codes = np.array([cmap[v] for v in vals])
            n = len(codes)
            labels = np.array(["" if c < 0 else cats[c] for c in codes], dtype=object)
            good = np.array([str(x).strip() not in ("", "nan", "None", "unknown")
                             for x in labels])
            frac = good.mean() if n else 0.0
            vc = pd.Series(labels[good]).value_counts().head(4)
            return n, float(frac), list(vc.index)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tissue", required=True)
    ap.add_argument("--max", type=int, default=30)
    ap.add_argument("--lo", type=int, default=3000)
    ap.add_argument("--hi", type=int, default=90000)
    ap.add_argument("--write", type=str, default=None, help="CSV of labeled candidates")
    args = ap.parse_args()

    import gcsfs
    df = pd.read_parquet(META)
    if "organism" in df.columns:
        df = df[df["organism"].fillna("").str.contains("Homo sapiens|human", case=False, na=False)]
    df["tissue"] = df["tissue"].fillna(""); df["disease"] = df["disease"].fillna("")
    rx = TISSUE_RX.get(args.tissue, rf"\b{args.tissue}\b")
    m = df[df["tissue"].str.contains(rx, case=False, regex=True, na=False)].copy()
    m = m[~m["tissue"].str.contains(BAD_TISSUE, case=False, regex=True, na=False)]
    m = m[~m["disease"].str.contains(BAD_DISEASE, case=False, regex=True, na=False)]
    if "cell_line" in m.columns:
        cl = m["cell_line"].fillna("").astype(str).str.strip().str.lower()
        m = m[cl.isin({"", "nan", "none", "na", "unsure", "unknown", "not applicable", "other"})]
    m = m[(m["obs_count"] >= args.lo) & (m["obs_count"] <= args.hi)]
    m = m.sort_values("obs_count", ascending=False)
    key = "entrez_id" if "entrez_id" in m.columns else "srx_accession"
    m = m.drop_duplicates(subset=[key]).head(args.max)
    log.info(f"{args.tissue}: peeking {len(m)} candidates (of catalog match)")

    fs = gcsfs.GCSFileSystem(token="anon")
    rows = []
    for _, r in m.iterrows():
        try:
            res = peek(fs, r["file_path"])
        except Exception as e:
            log.warning(f"  {r['srx_accession']}: peek failed ({type(e).__name__})"); continue
        if res is None:
            log.info(f"  {r['srx_accession']:>14} {int(r['obs_count']):>7} cells  NO cell_type"); continue
        n, frac, top = res
        flag = "✓LABELED" if frac >= 0.5 else "·unlabeled"
        log.info(f"  {r['srx_accession']:>14} {n:>7} cells  {frac*100:>4.0f}% {flag}  "
                 f"{str(r['tissue'])[:28]:<28} | {', '.join(top)[:50]}")
        rows.append({"srx": r["srx_accession"], "gcs": r["file_path"], "tissue": args.tissue,
                     "cells": n, "labeled": round(frac, 3),
                     "disease": str(r["disease"])[:30], "top": "; ".join(top)})
    out = pd.DataFrame(rows)
    lab = out[out["labeled"] >= 0.5]
    log.info(f"\n{args.tissue}: {len(lab)}/{len(out)} labeled, "
             f"{int(lab['cells'].sum()):,} labeled cells total")
    if args.write and len(lab):
        p = AEGIS / "manifests" / args.write
        lab.to_csv(p, index=False)
        log.info(f"Wrote {p}")


if __name__ == "__main__":
    main()
