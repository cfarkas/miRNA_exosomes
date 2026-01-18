#!/usr/bin/env bash
set -euo pipefail

# Convenience wrapper to reproduce the "without sample2" run.
#
# Usage:
#   bash scripts/run_without_sample2.sh data/miRNA-counts_example.csv results/without_sample2

INPUT_CSV="${1:-data/miRNA-counts_example.csv}"
OUTDIR="${2:-results/without_sample2}"

mkdir -p "${OUTDIR%/*}" || true

Rscript scripts/run_mirna_pipeline.R \
  --input "${INPUT_CSV}" \
  --outdir "${OUTDIR}" \
  --drop_samples Ctrl_2,Olig_2
