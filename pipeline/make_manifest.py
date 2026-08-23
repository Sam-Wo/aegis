#!/usr/bin/env python3
"""
Build the combined Aegis manifest: existing annotated normals + downloaded raw
samples that carry native scBaseCount cell-type labels (>= MIN_LABELED coverage).
Unlabeled samples are skipped (listed at the end) since we have no cross-tissue
annotator wired in yet.
"""
import warnings; warnings.filterwarnings("ignore")
import scanpy as sc
import pandas as pd
from pathlib import Path

AEGIS = Path(__file__).resolve().parent.parent
PROJECT = AEGIS.parent
RAW = AEGIS / "raw"
ANNOT = AEGIS / "annot"   # our own annotations (e.g. marker-based kidney)
MIN_LABELED = 0.5

rows, skipped = [], []

# 1) existing annotated normals (already have cell types)
for tissue, sub in [("liver", "normal/liver"), ("bone_marrow", "normal/bone_marrow")]:
    for f in sorted((PROJECT / "results" / sub).glob("*_annotated.h5ad")):
        rows.append({"path": f"results/{sub}/{f.name}", "tissue": tissue})

# 1b) our own annotated tissues (kidney via marker annotation, etc.)
if ANNOT.exists():
    for d in sorted(ANNOT.glob("*")):
        if d.is_dir():
            for f in sorted(d.glob("*.h5ad")):
                rows.append({"path": f"aegis/annot/{d.name}/{f.name}", "tissue": d.name})

# 2) downloaded raw samples with native labels
for d in sorted(RAW.glob("*")):
    if not d.is_dir():
        continue
    tissue = d.name
    for f in sorted(d.glob("*.h5ad")):
        a = sc.read_h5ad(f, backed="r")
        if "cell_type" not in a.obs.columns:
            skipped.append((tissue, f.stem, "no cell_type col")); continue
        ct = a.obs["cell_type"].astype(str).str.strip()
        frac = (~ct.isin(["", "nan", "None"])).mean()
        if frac >= MIN_LABELED:
            rows.append({"path": f"aegis/raw/{tissue}/{f.name}", "tissue": tissue})
        else:
            skipped.append((tissue, f.stem, f"{frac*100:.0f}% labeled"))

man = pd.DataFrame(rows)
out = AEGIS / "manifests" / "all.csv"
man.to_csv(out, index=False)
print(f"Wrote {out.relative_to(AEGIS)} — {len(man)} samples")
print(man["tissue"].value_counts().to_string())
if skipped:
    print("\nSKIPPED (unlabeled):")
    for t, s, why in skipped:
        print(f"  {t}/{s}: {why}")
