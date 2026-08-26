#!/bin/bash
# Re-annotate all tumor raw files with corrected gene positions (no re-download).
MAMBA=/c/Users/SamWo/miniforge3/condabin/mamba.bat
cd "$(dirname "$0")/../.."
for c in breast luad hcc pdac ovarian gastric crc; do
  echo "############## $c ##############"
  "$MAMBA" run -n adc-pipeline python aegis/pipeline/tumor_pipeline.py --cancer "$c" --n 14 || echo "  ($c errors, continuing)"
done
echo "REANNOTATE DONE"
