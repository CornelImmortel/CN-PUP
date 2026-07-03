# Pat_04_CTCs run config

Generated from a user-supplied BFX ID list: a 7-CTC subset of Patient 04
plus a WBC/germline reference, run through the MonoVar WBC-subset branch
of the pipeline (same approach as `Pat_03_CTCs`).

## Patient mapping

- Patient sheet `patient_id`: `PAT-2026-03-00004` (breast_cancer).
- CTCs: BFX-2026-04-00037, 00038, 00039, 00040, 00041, 00044, 00045.
- WBC/germline reference: BFX-2026-04-00050.

## Verification

Unlike Patient 03, there was no dedicated `cnpup_patient_04_sample_manifest_with_bams.csv`
for this patient, so paths initially came only from the less-reliable
`cnpup_unmapped_bam_rows.csv`, which had a genuine ambiguity: both
BFX-2026-04-00043 and BFX-2026-04-00050 were listed under the same
folder index (`PAT-2026-03-00004_10`) — multiple samples can share one
numbered output folder, so the folder index alone doesn't disambiguate.

This was resolved two ways:
1. Cell-type identity (which BFX ID is the WBC vs. a CTC) was confirmed
   against the lab's specimen/gDNA tracking table, which explicitly
   flags BFX-2026-04-00050 as `WBC` and all others in this list as `CTC`.
2. The exact BAM path was confirmed with a direct `ls -la` on the server
   against `PAT-2026-03-00004_10/alignments/dna/`, showing
   `BFX-2026-04-00050.redux.bam` genuinely exists there (alongside
   BFX-2026-04-00043's files, which belong to a different, unrelated
   sample not used in this run).

All paths in this config are therefore verified, not inferred.

## Local environment limitation

Same as `Pat_03_CTCs`: `/tzu-share-2` is not mounted and no conda is
installed in this Windows/WSL environment, so this config has only been
structurally validated (`nextflow run main.nf -params-file
docs/tzu_run/params.Pat_04_CTCs.monovar.yml` reaches `CHECK_INPUTS` and
fails at exactly the expected `ref_fasta not found` point). It has not
been executed end-to-end here.

## Run on the server

From the `CeciNestPasUnePipeline` repo root, with `/tzu-share-2` mounted
and the `conda` profile able to build environments:

```bash
bash docs/tzu_run/run_Pat_04_CTCs_monovar.sh
```

which runs:

```bash
nextflow run main.nf \
  -plugins nf-weblog \
  -profile conda \
  -params-file docs/tzu_run/params.Pat_04_CTCs.monovar.yml \
  -resume
```

Weblog is set to `http://localhost:8605`, the shared node listener that
feeds the internal pipeline dashboard (same as `Pat_03_CTCs` — this is a
fixed port, not per-run). Output goes to
`results_good_ctcs/PAT-2026-03-00004/...`.

The DelSIEVE `applauncher`/`beast`/`treeannotator` params point directly
at `/tzu-share-2/users/students/cornelusp/beast/bin/` (confirmed to
contain all three binaries), carrying forward the fix discovered while
running `Pat_03_CTCs` — the original bare `applauncher`/`beast` command
names don't resolve because `DELSIEVE_DATA_COLLECTOR` and
`DELSIEVE_TREE_ANNOTATOR` have no `conda` directive and so run without
any environment activated.

## Files

- `patients.Pat_04_CTCs.monovar.tsv` — patient sheet row.
- `Pat_04_CTCs.cells.tsv` — cell metadata (7 CTC + 1 WBC rows).
- `Pat_04_CTCs.monovar_bams.txt` — MonoVar CTC BAM list (CTCs only; WBC
  is handled separately via `leukocyte_bam` in the patient sheet).
- `params.Pat_04_CTCs.monovar.yml` — full params file (MonoVar, VEP,
  filtering, comparison, final/SNV reports, BAM QC, DelSIEVE stage 1 all
  enabled — mirrors the `Pat_03_CTCs` config).
- `run_Pat_04_CTCs_monovar.sh` — launcher.
