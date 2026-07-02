#!/usr/bin/env bash
set -euo pipefail

nextflow run main.nf \
  -plugins nf-weblog \
  -profile conda \
  -params-file docs/tzu_run/params.Pat_03_CTCs.monovar.yml \
  -resume
