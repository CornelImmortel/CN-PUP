# CeciNestPasUnePipeline

Nextflow pipeline under development for CTC single-cell variant calling and filtering.

The final goal is to automate:

1. MonoVar calling from patient-level BAM lists.
2. Germline/background exclusion using precomputed VCFs, bulk DeepVariant calls, leukocyte calls, or combined modes.
3. Quality filtering, SNV restriction, normalization, and subtraction.
4. VEP population AF / COSMIC filtering.
5. Long tables, matrices, Shiny app inputs, and QC reports.

## Step 0

By default the workflow validates the patient input sheet. MonoVar calling and per-cell MonoVar VCF splitting are implemented but opt-in with `--run_monovar true`.

Prepare a real config:

```bash
cp configs/patients.tsv.example configs/patients.tsv
# Optional: copy/edit the per-cell sheet if you do not want to use the example directly
cp configs/patient_01_cells.tsv.example configs/patient_01_cells.tsv
```

Then run:

```bash
nextflow run main.nf
```

Or with an explicit patient sheet:

```bash
nextflow run main.nf --patients configs/patients.tsv
```

Expected output:

```text
results/validation/<patient_id>.input_check.txt
```

## Germline modes planned

- `precomputed_vcf`
- `deepvariant_bulk_bam`
- `leukocyte_vcf`
- `joint_monovar_leukocyte`
- `combined`

Only validation is implemented for now. Real processes will be added step by step.

## Run MonoVar

After validation works on the server, launch MonoVar with:

```bash
nextflow run main.nf -profile conda --run_monovar true --monovar_script /tzu-share-2/users/students/cornelusp/monovar/monovar/src/monovar.py --monovar_threads 6
```

Outputs:

```text
results/<patient_id>/raw_calls/monovar/<patient_id>.monovar.vcf
results/<patient_id>/logs/<patient_id>.monovar.log
results/<patient_id>/split_calls/monovar/<cell_id>.monovar.split.vcf
results/<patient_id>/split_calls/monovar/<patient_id>.monovar.split_sample_map.tsv
results/<patient_id>/logs/<patient_id>.monovar.split.log
```

### Small region test

To test the MonoVar process without launching a full run, restrict `samtools mpileup` to one interval:

```bash
nextflow run main.nf \
  -profile conda \
  --run_monovar true \
  --monovar_script /tzu-share-2/users/students/cornelusp/monovar/monovar/src/monovar.py \
  --monovar_threads 2 \
  --monovar_region chr1:1-1000000
```

Use `-resume` after fixing failures:

```bash
nextflow run main.nf -profile conda --run_monovar true --monovar_region chr1:1-1000000 -resume
```


## Cell metadata

`configs/patients.tsv` stays one row per patient. The `cell_metadata` column points to a per-patient TSV with one row per cell:

```text
patient_id    cell_id    bam_path    cell_type
```

The splitter maps MonoVar sample columns to `cell_id` by matching SRR IDs in the BAM paths. If needed, it falls back to the row order in the cell metadata file.


## Step 4: background subtraction

After MonoVar splitting, the pipeline can produce comparable per-CTC MonoVar VCFs by filtering each split VCF and subtracting a background/exclusion VCF.

Implemented modes:

- `precomputed_vcf`: normalize a precomputed germline/bulk VCF from `germline_vcf`, then subtract it.
- `joint_monovar_leukocyte`: use the split MonoVar `leukocyte` sample from the same joint call as the exclusion set.

To test the leukocyte mode, use the leukocyte example sheet:

```bash
cp configs/patients.leukocyte.tsv.example configs/patients.tsv
nextflow run main.nf \
  -profile conda \
  -resume \
  --run_monovar true \
  --monovar_script /tzu-share-2/users/students/cornelusp/monovar/monovar/src/monovar.py \
  --monovar_threads 2 \
  --monovar_region chr1:1-1000000
```

Outputs:

```text
results/<patient_id>/germline_exclusion/<patient_id>.monovar_leukocyte.filtered.norm.vcf.gz
results/<patient_id>/processed_calls/monovar/<cell_id>.monovar.filtered.vcf
results/<patient_id>/normalized_calls/monovar/<cell_id>.monovar.filtered.norm.vcf.gz
results/<patient_id>/comparable_calls/monovar/<cell_id>.monovar.no_monovar_leukocyte.vcf.gz
```
