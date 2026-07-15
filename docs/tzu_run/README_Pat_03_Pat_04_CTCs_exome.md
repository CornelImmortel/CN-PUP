# Pat_03_Pat_04_CTCs exome-restricted run config

Exome-restricted variant of `params.Pat_03_Pat_04_CTCs.monovar.yml` (the WGS-scope
combined run behind `results_good_ctcs/`) — same patients, same manifest, same WBC
identification, only the MonoVar candidate-site scope and BAM-QC target region change.

## Patient mapping (unchanged from the WGS combined run)

- `PAT-2026-03-00003` (prostate_cancer): 9 CTCs (BFX-2026-04-00022, 00024, 00026,
  00027, 00028, 00029, 00030, 00031, 00032) + WBC `BFX-2026-04-00051`.
- `PAT-2026-03-00004` (breast_cancer): 7 CTCs + WBC `BFX-2026-04-00050`.

Both WBC assignments were wet-lab-confirmed (see `STATUS.md`) and are already correct in
`patients.Pat_03_Pat_04_CTCs.monovar.tsv` and `Pat_0{3,4}_CTCs.cells.tsv` — reused as-is,
no changes needed for this exome variant.

## What actually changes vs. the WGS combined run

1. `outdir: results_exome_ctcs` (was `results_good_ctcs`) — dedicated output directory
   so this run doesn't overwrite/mix with the WGS-scope results.
2. `monovar_targets_bed` → GENCODE v50 merged exon BED (was `canonical_chroms.bed`,
   which is essentially whole-genome — see the comment in
   `params.Pat_03_Pat_04_CTCs.monovar.exome.yml` for the full explanation). Path:
   `/tzu-share-2/data/TZU/novovene_staging/wrappers/PTATO_pat03/exonic_coverage/gencode_v50_exons_merged.bed.gz`
   — 407,835 merged, coordinate-sorted, non-overlapping intervals, confirmed via that
   directory's own `README.md` as the file meant for exactly this use (not its sibling
   `gencode_v50_exons_raw.bed.gz`, which is unmerged/overlapping BED6, the wrong shape
   for a `-l` target).
3. `mosdepth_by` → the same exome BED (was empty/whole-genome), so BAM-QC coverage
   stats are scoped consistently with the restricted variant calls.

Everything else — filter thresholds, VEP config, `germline_mode: monovar_wbc_subset`,
DelSIEVE (off) — is identical to the WGS combined run.

## Run on the server

From the `CeciNestPasUnePipeline` repo root:

```bash
bash docs/tzu_run/run_Pat_03_Pat_04_CTCs_monovar_exome.sh
```

which runs:

```bash
nextflow run main.nf \
  -plugins nf-weblog \
  -profile conda \
  -w /scratch/cornelusp/work \
  -params-file docs/tzu_run/params.Pat_03_Pat_04_CTCs.monovar.exome.yml \
  -resume
```

Note the `-w /scratch/cornelusp/work` — Nextflow's work directory (intermediate
per-task files) is pointed at fast local scratch instead of the default (which would
land under the repo checkout on `/tzu-share-2`, network storage). Final results still
publish to `results_exome_ctcs/` on `/tzu-share-2` as normal — only the intermediate
work directory moves.

Weblog is set to `http://localhost:8605`, the shared node-wide listener, same as every
other run in this directory.

## Files

- `params.Pat_03_Pat_04_CTCs.monovar.exome.yml` — full params file.
- `run_Pat_03_Pat_04_CTCs_monovar_exome.sh` — launcher.
- `patients.Pat_03_Pat_04_CTCs.monovar.tsv`, `Pat_0{3,4}_CTCs.cells.tsv`,
  `Pat_0{3,4}_CTCs.monovar_bams.txt` — reused unchanged from the WGS combined run.
