# Pipeline Overview

This document tracks the intended Nextflow workflow.

```text
patients.tsv
  -> CHECK_INPUTS
  -> RUN_MONOVAR
  -> SPLIT_MONOVAR
  -> PREPARE_GERMLINE_EXCLUSION
  -> FILTER_NORMALIZE_SUBTRACT
  -> VCF_TO_LONG_TABLE
  -> MERGE_LONG_TABLES
  -> PREPARE_VEP_INPUT
  -> RUN_VEP
  -> APPLY_POPULATION_COSMIC_FILTER
  -> BUILD_MATRICES
  -> WRITE_QC_REPORTS
```

The first implemented step is only `CHECK_INPUTS`.
