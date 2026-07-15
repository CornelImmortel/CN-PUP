#!/usr/bin/env bash
set -euo pipefail

nextflow run main.nf \
  -plugins nf-weblog \
  -profile conda \
  -w /scratch/cornelusp/work \
  -params-file docs/tzu_run/params.Pat_03_Pat_04_CTCs.monovar.exome.yml \
  -resume
