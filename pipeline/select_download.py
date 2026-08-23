#!/usr/bin/env python3
"""
Aegis — Healthy-tissue sample selector / downloader
Selects healthy WHOLE-tissue scBaseCount samples for the target safety tissues
and (optionally) downloads them from GCS. Dry-run prints the plan first.

Usage:
  python select_download.py --dry-run
  python select_download.py --download --per-tissue 3
"""
import warnings; warnings.filterwarnings("ignore")
import argparse, logging, re
from pathlib import Path
import pandas as pd

logging.basicConfig(level=logging.INFO, format="%(asctime)s [%(levelname)s] %(message)s")
log = logging.getLogger("aegis.select")

AEGIS = Path(__file__).resolve().parent.parent
PROJECT = AEGIS.parent
META = PROJECT / "data" / "scbasecount" / "sample_metadata.parquet"
RAW = AEGIS / "raw"

# canonical tissue -> (regex on free-text tissue, safety tier)
TARGETS = {
    "heart":           (r"\bheart\b|myocard|cardiac|ventric|atri(um|al)", 1),
    "kidney":          (r"\bkidney\b|renal cortex|renal medulla|nephron|\brenal\b", 2),
    "lung":            (r"\blung\b|pulmonary|alveolar|lung parenchyma", 2),
    "colon":           (r"\bcolon\b|colonic|large intestine|sigmoid|rectum|rectal", 2),
    "small_intestine": (r"small intestine|\bileum\b|terminal ileum|jejunum|duoden", 2),
    "pancreas":        (r"\bpancrea", 2),
    "stomach":         (r"\bstomach\b|gastric (mucosa|tissue|corpus|antrum)|\bgastric\b", 2),
    "skin":            (r"\bskin\b|epiderm|dermis|dermal", 3),
}

# exclude non-whole-tissue / sorted / disease contexts in the free-text tissue label
BAD_TISSUE = (r"pbmc|periph|\bblood\b|bone marrow|organoid|cell line|cultured|ipsc|"
    r"lavage|sorted|enrich|mononuclear|phagocyte|dendritic|macrophage|myeloid|"
    r"lymphocyte|-derived|\bips\b|ipsc|lamina propria|\bislet|lymph node|fetal|embryo|"
    r"fibroblast|endothelial|arterial|epithelial cell|\bepithelium\b|"
    r"cd4\b|cd8\b|cd34|cd45|epcam|tumou?r|carcinoma|cancer|myxoma|adenoma|"
    r"metasta|biopsy of|xenograft|spheroid|barrett")
BAD_DISEASE = (r"cancer|carcinoma|tumou?r|neoplas|malign|adenoma|myxoma|polyp|leukem|leukaem|"
    r"lymphoma|myeloma|sarcoma|melanoma|metasta|cirrho|fibrosis|hepatitis|colitis|"
    r"crohn|diabet|syndrome|\bdisease\b|infect|covid|sepsis|failure|injury|"
    r"transplant|trisomy|dystrophy|sclero|barrett|metaplasia|dysplasia|"
    r"cardiomyopathy|dilated|stenosis|hypertension|hypoglycemia|"
    r"psoria|eczema|dermatitis|keratosis|atopic|\bscar\b|hypertrophic|\bwound\b|"
    r"vasculit|anca|neutrophil cytoplasm|immune-?relat|"
    r"\bild\b|ctd|amyloid|auto-?anti|auto-?immun|allerg|asthma|copd")
HEALTHY_HINT = r"normal|healthy|control|non-?disease|reference|adjacent normal|donor"
AMBIG = r"\bunsure\b"
NULL_TOKENS = {"", "nan", "none", "null", "na", "n/a", "unsure", "unknown",
               "not applicable", "not_applicable", "not collected", "other"}


def load():
    if not META.exists():
        raise SystemExit(f"Metadata not found: {META}")
    df = pd.read_parquet(META)
    if "organism" in df.columns:
        df = df[df["organism"].fillna("").str.contains("Homo sapiens|human", case=False, na=False)]
    df["tissue"] = df["tissue"].fillna("")
    df["disease"] = df["disease"].fillna("")
    return df


def select(df, per_tissue, lo=2000, hi=60000):
    plans = {}
    for tissue, (rx, tier) in TARGETS.items():
        m = df[df["tissue"].str.contains(rx, case=False, regex=True, na=False)].copy()
        m = m[~m["tissue"].str.contains(BAD_TISSUE, case=False, regex=True, na=False)]
        if "cell_line" in m.columns:
            cl = m["cell_line"].fillna("").astype(str).str.strip().str.lower()
            m = m[cl.isin(NULL_TOKENS)]
        dis = m["disease"].fillna("").astype(str).str.strip()
        not_bad = ~dis.str.contains(BAD_DISEASE, case=False, regex=True, na=False)
        not_ambig = ~dis.str.contains(AMBIG, case=False, regex=True, na=False)
        healthy = not_bad & not_ambig  # keep empty / explicit-healthy / benign-none
        m = m[healthy.astype(bool)]
        m = m[(m["obs_count"] >= lo) & (m["obs_count"] <= hi)]
        m = m.sort_values("obs_count", ascending=False)
        # diverse studies
        key = "entrez_id" if "entrez_id" in m.columns else "srx_accession"
        sel = m.drop_duplicates(subset=[key]).head(per_tissue)
        plans[tissue] = (tier, sel)
    return plans


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--download", action="store_true")
    ap.add_argument("--per-tissue", type=int, default=3)
    args = ap.parse_args()

    df = load()
    log.info(f"Human samples in catalog: {len(df):,}")
    plans = select(df, args.per_tissue)

    manifest_rows = []
    total_cells = tot = 0
    print(f"\n{'='*68}\nSELECTION PLAN ({args.per_tissue}/tissue)\n{'='*68}")
    for tissue, (tier, sel) in plans.items():
        print(f"\n── {tissue}  (Tier {tier})  [{len(sel)} samples, {sel['obs_count'].sum():,} cells]")
        if len(sel) == 0:
            print("   ⚠ no healthy whole-tissue sample matched")
            continue
        for _, r in sel.iterrows():
            print(f"   {r['srx_accession']:>14}  {r['obs_count']:>7,} cells  "
                  f"tissue='{str(r['tissue'])[:34]}'  disease='{str(r['disease'])[:22]}'")
            manifest_rows.append({"path": f"aegis/raw/{tissue}/{r['srx_accession']}.h5ad",
                                  "tissue": tissue, "srx": r["srx_accession"],
                                  "gcs": r["file_path"], "cells": r["obs_count"]})
            total_cells += int(r["obs_count"]); tot += 1
    print(f"\n{'='*68}\nTOTAL: {tot} samples, ~{total_cells:,} cells "
          f"(~{total_cells*12/1e6:.1f} GB raw est.)\n{'='*68}")

    man = pd.DataFrame(manifest_rows)
    (AEGIS / "manifests").mkdir(exist_ok=True)
    man[["path", "tissue"]].to_csv(AEGIS / "manifests" / "fetched.csv", index=False)
    man.to_csv(AEGIS / "manifests" / "fetched_full.csv", index=False)
    log.info(f"Wrote manifests/fetched.csv ({len(man)} samples)")

    if args.dry_run or not args.download:
        print("\n[DRY RUN] Re-run with --download to fetch.")
        return

    import gcsfs
    fs = gcsfs.GCSFileSystem(token="anon")
    ok = fail = 0
    for _, r in man.iterrows():
        out = AEGIS / "raw" / r["tissue"] / f"{r['srx']}.h5ad"
        if out.exists():
            log.info(f"  exists: {out.name}"); ok += 1; continue
        out.parent.mkdir(parents=True, exist_ok=True)
        try:
            fs.get(r["gcs"], str(out))
            log.info(f"  ↓ {r['tissue']}/{out.name} ({out.stat().st_size/1e6:.0f} MB)")
            ok += 1
        except Exception as e:
            log.warning(f"  FAILED {r['srx']}: {e}"); fail += 1
    log.info(f"Done: {ok} ok, {fail} failed")


if __name__ == "__main__":
    main()
