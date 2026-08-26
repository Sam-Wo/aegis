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
ANNOT = AEGIS / "annot"            # our own normal annotations (marker-based kidney)
ANNOT_TUMOR = AEGIS / "annot_tumor"  # CNV-refined tumor annotations
MIN_LABELED = 0.5

CANCER_LABEL = {"breast": "breast tumor", "luad": "lung adeno tumor", "crc": "colorectal tumor",
    "hcc": "HCC tumor", "pdac": "pancreatic tumor", "ovarian": "ovarian tumor",
    "gastric": "gastric tumor"}

rows, skipped = [], []

# 1) existing annotated normals (already have cell types)
for tissue, sub in [("liver", "normal/liver"), ("bone_marrow", "normal/bone_marrow")]:
    for f in sorted((PROJECT / "results" / sub).glob("*_annotated.h5ad")):
        rows.append({"path": f"results/{sub}/{f.name}", "tissue": tissue, "kind": "normal"})

# 1b) our own annotated normal tissues (kidney via marker annotation, etc.)
if ANNOT.exists():
    for d in sorted(ANNOT.glob("*")):
        if d.is_dir():
            for f in sorted(d.glob("*.h5ad")):
                rows.append({"path": f"aegis/annot/{d.name}/{f.name}", "tissue": d.name, "kind": "normal"})

# 1c) CNV-refined tumor annotations (malignant cells aggregated downstream)
if ANNOT_TUMOR.exists():
    for d in sorted(ANNOT_TUMOR.glob("*")):
        if d.is_dir():
            label = CANCER_LABEL.get(d.name, d.name + " tumor")
            for f in sorted(d.glob("*.h5ad")):
                rows.append({"path": f"aegis/annot_tumor/{d.name}/{f.name}", "tissue": label, "kind": "tumor"})

# 2) downloaded raw normal samples with native labels
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
            rows.append({"path": f"aegis/raw/{tissue}/{f.name}", "tissue": tissue, "kind": "normal"})
        else:
            skipped.append((tissue, f.stem, f"{frac*100:.0f}% labeled"))

man = pd.DataFrame(rows)
out = AEGIS / "manifests" / "all.csv"
man.to_csv(out, index=False)
print(f"Wrote {out.relative_to(AEGIS)} — {len(man)} samples "
      f"({(man['kind']=='tumor').sum()} tumor, {(man['kind']=='normal').sum()} normal)")
print(man.groupby(['kind','tissue']).size().to_string())
if skipped:
    print("\nSKIPPED (unlabeled):")
    for t, s, why in skipped:
        print(f"  {t}/{s}: {why}")
