# Patient 02 / 05 / 06 coverage-analysis configs (Su et al. 2019 / PRJNA448888)

Extends the same coverage-vs-tissue analysis already done for Patient 01 and
Patient 03 to three more patients from the same public cohort. No PTATO, no
caller-scorecard extension (only Patient 01/03 have that) -- this is coverage
only, using the DeepVariant tissue calls already present on
`/tzu-share-2/data/public/PRJNA448888/processed/sarek_variantcalling/`.

## Per-patient summary

| Patient | CTCs | Tissue | Germline mode | Germline source |
|---|---|---|---|---|
| 02 | 11 (CTC1-11) | primary + metastasis | `combined` | bulk-blood DeepVariant VCF (SRR8617527) + leukocyte BAM (SRR8617531), same pattern as Patient 01 |
| 05 | 12 (CTC1-12) | metastasis only | `joint_monovar_leukocyte` | leukocyte BAM only (SRR8617562) -- bulk-blood sample SRR8617564 has no precomputed caller VCF on this share, only a BAM, so this deliberately does not combine sources |
| 06 | 12 (CTC1-12) | primary only | `precomputed_vcf` | bulk-blood DeepVariant VCF only (SRR8617591) -- no leukocyte sample exists for this patient in the cohort at all |

## Run recipe (per patient, two Nextflow passes)

```bash
nextflow run main.nf -plugins nf-weblog -profile conda \
  -params-file configs/params.patient02_step1_monovar.yml -resume
nextflow run main.nf -plugins nf-weblog -profile conda \
  -params-file configs/params.patient02_step2_compare.yml -resume
```
Same for `patient05` / `patient06`. Step 2 must reuse step 1's `patients.tsv`
(it does, unchanged) so `-resume` finds the already-computed MonoVar split
VCFs instead of recalling them.

## What's verified vs. best-effort in these configs

**Verified against live schema/validation code in `main.nf`** (not guessed):
- `germline_mode` values and their required fields (`precomputed_vcf`,
  `joint_monovar_leukocyte`, `combined` -- read directly from the
  `CHECK_PATIENT_INPUTS` process's case statement, main.nf ~line 1713).
- `patients.tsv` / `cells.tsv` / `bams.txt` schemas -- copied from Patient
  01's and Patient 03's actual working configs.
- `tissue_vcfs.tsv` schema -- copied from
  `configs/tissue_vcfs.patient_03_wxs_meta_leuko.tsv`.
- All SRR accessions and DeepVariant VCF paths -- cross-checked against the
  actual `sarek_variantcalling`/`sarek/alignments` directory listings and
  the full PRJNA448888 SRA run-table workbook, not assumed.

**Best-effort, modeled on Patient 03's compare config but not dry-run
verified** (flagged individually in each `params.patientNN_step2_compare.yml`):
- `caller_vcfs.patientNN_compare.tsv`'s relative `results/...` paths for
  the MonoVar CTC splits. Patient 03's version of this file used absolute
  paths into a specific prior checkout (`CN-PUP2`); these use relative
  paths instead since the actual checkout location on the server wasn't
  known when writing this. If `run_comparison`/`run_snv_report` error on
  those rows, the fix is pointing them at wherever step 1's `outdir`
  actually resolves to on this checkout.
- Whether `run_comparison`/`run_snv_report` even need MonoVar rows in
  `caller_vcfs.tsv` at all, versus picking up the native MonoVar branch's
  own output automatically within the same run/DAG. Included them to
  match Patient 03's precedent; if they turn out to be redundant, that's
  harmless, not wrong.

Worth a quick sanity run on one patient (e.g. Patient 06, the simplest
germline setup) before launching all three, given the above two items
haven't been exercised end-to-end yet.
