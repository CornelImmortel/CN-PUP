# Pat_03_CTCs run config

Generated from `THESIS/analysis/tzu_runs/Pat_03_CTCs.txt` (a BFX ID list
supplied by the user): a 9-CTC subset of Patient 03 plus a WBC/germline
reference, run through the existing MonoVar WBC-subset branch of the
pipeline.

## Patient mapping

- Patient sheet `patient_id`: `PAT-2026-03-00003` (prostate_cancer).
- CTCs: BFX-2026-04-00022, 00024, 00026, 00027, 00028, 00029, 00030,
  00031, 00032.
- WBC/germline reference: BFX-2026-04-00033 — the same WBC used in the
  original full-cohort `PAT-2026-03-00003` run in this directory.

Patient identity and cell-type were resolved by cross-referencing
`THESIS/analysis/sample_tracking/raw_verbatim.txt` (lab tracking sheet)
with `docs/tzu_run/cnpup_patient_03_sample_manifest_with_bams.csv`. All
BAM paths for this run come from that manifest's `matched_by_id` rows,
i.e. they were confirmed against an actual server directory listing in a
prior session.

## Local environment limitation

`/tzu-share-2` is not mounted in this Windows/WSL environment, and no
conda/mamba is installed in WSL, so `CHECK_INPUTS` cannot pass and
`-profile conda` cannot build environments here. This config has been
structurally validated only: `nextflow run main.nf -params-file
docs/tzu_run/params.Pat_03_CTCs.monovar.yml` reaches `CHECK_INPUTS` and
fails at exactly the expected point (`ref_fasta not found under
/tzu-share-2`), confirming the patient sheet, cell metadata, and params
file are wired correctly. It has not been executed end-to-end.

## Run on the server

From the `CeciNestPasUnePipeline` repo root, with `/tzu-share-2` mounted
and the `conda` profile able to build environments:

```bash
bash docs/tzu_run/run_Pat_03_CTCs_monovar.sh
```

which runs:

```bash
nextflow run main.nf \
  -plugins nf-weblog \
  -profile conda \
  -params-file docs/tzu_run/params.Pat_03_CTCs.monovar.yml \
  -resume
```

Weblog is set to `http://localhost:8606` (the original full-cohort
Patient 03 run in this directory uses 8605, so both can run without
clashing). Output goes to `results_good_ctcs/PAT-2026-03-00003/...`, kept
separate from the original full-cohort run's `results/PAT-2026-03-00003/...`
via a dedicated `outdir`.

## Files

- `patients.Pat_03_CTCs.monovar.tsv` — patient sheet row.
- `Pat_03_CTCs.cells.tsv` — cell metadata (9 CTC + 1 WBC rows).
- `Pat_03_CTCs.monovar_bams.txt` — MonoVar CTC BAM list (CTCs only; WBC is
  handled separately via `leukocyte_bam` in the patient sheet).
- `params.Pat_03_CTCs.monovar.yml` — full params file (MonoVar, VEP,
  filtering, comparison, final/SNV reports, BAM QC, DelSIEVE stage 1 all
  enabled — mirrors the existing full-cohort Patient 03 config).
- `run_Pat_03_CTCs_monovar.sh` — launcher.
