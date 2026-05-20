#!/usr/bin/env bash
set -euo pipefail

nextflow run main.nf \
  -plugins nf-weblog \
  -profile conda \
  -params-file docs/tzu_run/params.PAT-2026-03-00003.monovar.yml \
  -resume
