#!/usr/bin/env bash
set -euo pipefail

# Convenience wrapper to run the pipeline with all samples.
#
# Usage:
#   bash scripts/run_all_samples.sh data/miRNA-counts_example.csv results/all_samples

INPUT_CSV="${1:-data/miRNA-counts_example.csv}"
OUTDIR="${2:-results/all_samples}"

mkdir -p "${OUTDIR%/*}" || true

Rscript scripts/run_mirna_pipeline.R \
  --input "${INPUT_CSV}" \
  --outdir "${OUTDIR}"
