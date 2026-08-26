#!/usr/bin/env python3
"""
Build a cached gene genomic-position table (hg38) for infercnvpy CNV calling.
Keyed by both Ensembl ID and gene symbol. Saved to aegis/refs/gene_positions.parquet
"""
import warnings; warnings.filterwarnings("ignore")
import infercnvpy as cnv  # noqa: F401  (import check)
import pandas as pd
from pathlib import Path

AEGIS = Path(__file__).resolve().parent.parent
OUT = AEGIS / "refs"; OUT.mkdir(exist_ok=True)

from pybiomart import Server
print("Querying Ensembl (pybiomart)…")
server = Server(host="http://www.ensembl.org")
mart = server["ENSEMBL_MART_ENSEMBL"]["hsapiens_gene_ensembl"]
df = mart.query(attributes=[
    "ensembl_gene_id", "external_gene_name",
    "chromosome_name", "start_position", "end_position"])
df.columns = ["ensembl_id", "gene", "chromosome", "start", "end"]
# keep standard chromosomes only
keep = [str(c) for c in list(range(1, 23)) + ["X", "Y"]]
df = df[df["chromosome"].astype(str).isin(keep)].copy()
df["chromosome"] = "chr" + df["chromosome"].astype(str)
df["start"] = pd.to_numeric(df["start"], errors="coerce")
df["end"] = pd.to_numeric(df["end"], errors="coerce")
df = df.dropna(subset=["start", "end"])
df.to_parquet(OUT / "gene_positions.parquet", index=False)
print(f"Wrote {OUT/'gene_positions.parquet'} — {len(df):,} genes "
      f"({df['gene'].nunique():,} symbols, {df['ensembl_id'].nunique():,} ensembl ids)")
print(df.head(4).to_string(index=False))
