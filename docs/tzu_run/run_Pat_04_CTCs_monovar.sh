#!/usr/bin/env bash
set -euo pipefail

nextflow run main.nf \
  -plugins nf-weblog \
  -profile conda \
  -params-file docs/tzu_run/params.Pat_04_CTCs.monovar.yml \
  -resume
