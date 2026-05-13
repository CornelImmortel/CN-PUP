# CeciNestPasUnePipeline

Nextflow pipeline under development for CTC single-cell variant calling and filtering.

The final goal is to automate:

1. MonoVar calling from patient-level BAM lists.
2. Germline/background exclusion using precomputed VCFs, bulk DeepVariant calls, leukocyte calls, or combined modes.
3. Quality filtering, SNV restriction, normalization, and subtraction.
4. VEP population AF / COSMIC filtering.
5. Long tables, matrices, Shiny app inputs, and QC reports.

## Step 0

By default the workflow validates the patient input sheet. MonoVar calling is now implemented but opt-in with `--run_monovar true`.

Prepare a real config:

```bash
cp configs/patients.tsv.example configs/patients.tsv
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
```
