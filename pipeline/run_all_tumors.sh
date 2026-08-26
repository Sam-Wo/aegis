#!/bin/bash
# Run the CNV-refined tumor pipeline across the core ADC carcinomas.
set -e
MAMBA=/c/Users/SamWo/miniforge3/condabin/mamba.bat
cd "$(dirname "$0")/../.."
for c in breast luad hcc pdac ovarian gastric crc; do
  echo "############## $c ##############"
  "$MAMBA" run -n adc-pipeline python aegis/pipeline/tumor_pipeline.py --cancer "$c" --n 14 --download || echo "  ($c had errors, continuing)"
done
echo "ALL TUMORS DONE"
