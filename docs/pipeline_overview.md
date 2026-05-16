# Pipeline Overview

This document tracks the intended Nextflow workflow.

```text
patients.tsv
  -> CHECK_INPUTS
  -> RUN_MONOVAR
  -> SPLIT_MONOVAR
  -> PREPARE_GERMLINE_EXCLUSION
  -> FILTER_NORMALIZE_SUBTRACT
  -> optional VEP population/COSMIC filtering
  -> MonoVar reports
```

Additional comparison branch:

```text
caller_vcfs.tsv
  -> STANDARDIZE_EXTERNAL_CALLERS
  -> comparable_calls/<caller>/*.comparable.vcf.gz
  -> VCF_TO_LONG_TABLE
  -> MERGE_LONG_TABLES
  -> FILL_MISSING_ANNOTATIONS
  -> BUILD_MATRICES
  -> WRITE_OVERLAP_QC
```

Optional internal SCcaller branch:

```text
cell_metadata + bulk_bam + hSNP VCF
  -> RUN_SCCALLER_FROM_BAM
  -> STANDARDIZE_INTERNAL_SCCALLER
  -> comparable_calls/sccaller/*.comparable.vcf.gz
```

The comparison branch reuses the validated downstream scripts vendored from
`mutation_matrix_pipeline` under `bin/mutation_matrix/`.
