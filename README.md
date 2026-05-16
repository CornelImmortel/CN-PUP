# CeciNestPasUnePipeline

Nextflow pipeline under development for CTC single-cell variant calling and filtering.

The final goal is to automate:

1. MonoVar calling from patient-level BAM lists.
2. Optional SCcaller calling from existing cell BAMs and matched bulk BAMs.
3. Germline/background exclusion using precomputed VCFs, bulk DeepVariant calls, leukocyte calls, or combined modes.
4. External caller ingestion for SCcaller, HaplotypeCaller, DeepVariant, and MonoVar VCFs.
5. Quality filtering, SNV restriction, normalization, and subtraction.
6. VEP population AF / COSMIC filtering.
7. Caller comparison: long tables, matrices, CTC-SCITE inputs, Shiny app inputs, and QC reports.

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

## Main operating modes

CN-PUP can now be used in three increasingly broad modes.

## Input Builder

CN-PUP includes a small static browser interface for building the required
input files:

```text
interface/index.html
```

Open it directly in a browser. It generates:

```text
patients.tsv
cells.tsv
monovar_bams.txt
caller_vcfs.tsv
cnpup_command.sh
```

The interface is intentionally static HTML/JavaScript. It does not upload data,
does not need a web server, and does not run Nextflow. It only helps users make
valid TSV inputs and a launch command. After downloading the files, place them
under `configs/` or adjust the generated command paths.

MonoVar only:

```bash
nextflow run main.nf \
  -profile conda \
  --run_monovar true \
  --monovar_script /path/to/MonoVar/src/monovar.py
```

Compare existing caller VCFs:

```bash
cp configs/caller_vcfs.tsv.example configs/caller_vcfs.tsv

nextflow run main.nf \
  -profile conda \
  --patients configs/patients.tsv \
  --run_external_callers true \
  --run_comparison true \
  --caller_vcfs configs/caller_vcfs.tsv
```

Run MonoVar, optionally run SCcaller, and compare with external Sarek callers:

```bash
nextflow run main.nf \
  -profile conda \
  --patients configs/patients.tsv \
  --run_monovar true \
  --run_sccaller true \
  --run_external_callers true \
  --run_comparison true \
  --monovar_script /path/to/MonoVar/src/monovar.py \
  --sccaller_script /path/to/sccaller_v2.0.0.py \
  --sccaller_hsnp_vcf /path/to/bulk.hsnp.biallelic.dbsnp.vcf \
  --caller_vcfs configs/caller_vcfs.tsv
```

The comparison path intentionally reuses the validated downstream code from
`mutation_matrix_pipeline`, vendored under `bin/mutation_matrix/`, instead of
reimplementing the filtering and matrix logic from scratch.

## External caller manifest

External caller VCFs are listed in `configs/caller_vcfs.tsv`:

```text
patient_id    cell_id    caller           vcf_path                         sccaller_mode
Patient_03    SRR8617653 haplotypecaller  /path/to/cell.haplotypecaller.vcf.gz
Patient_03    SRR8617653 deepvariant      /path/to/cell.deepvariant.vcf.gz
Patient_03    SRR8617653 sccaller         /path/to/cell.sccaller.vcf       external_bulk
```

Supported caller names:

- `monovar`
- `sccaller`
- `haplotypecaller`
- `deepvariant`

## Patient 03 WXS CTC-only MonoVar config

The repository includes a clean Patient 03 WXS CTC-only setup for rerunning
MonoVar without the leukocyte in the joint calling cohort:

```text
configs/patients.patient_03_wxs_ctc_only.tsv
configs/patient_03_wxs_ctc_only_cells.tsv
configs/patient_03_wxs_ctc_only_bams.txt
```

This uses:

```text
germline_mode = no_germline
```

That mode validates inputs and runs MonoVar/splitting without requiring an
immediate background VCF. Use it when you want clean CTC-only MonoVar split VCFs
first, then apply matched leukocyte/bulk subtraction later during comparison.

Example:

```bash
nextflow run main.nf \
  -profile conda \
  --patients configs/patients.patient_03_wxs_ctc_only.tsv \
  --run_monovar true \
  --monovar_script /path/to/monovar.py \
  --monovar_threads 6 \
  -resume
```

All caller branches are normalized into the same downstream contract:

```text
results/<patient_id>/comparable_calls/<caller>/<cell_id>.<caller>.comparable.vcf.gz
```

Those comparable VCFs are then converted into long tables, mutation matrices,
CTC-SCITE inputs, and overlap QC reports.

SCcaller supports the same comparison modes as `mutation_matrix_pipeline`:

- `so_true`: require SCcaller internal `SO=True`.
- `external_bulk`: ignore `SO`, then subtract the configured external germline/bulk VCF.
- `so_true_and_external_bulk`: require `SO=True`, then subtract external bulk.
- `no_bulk`: diagnostic mode; no internal `SO` requirement and no external subtraction.

## Optional SCcaller calling from BAMs

CN-PUP can call SCcaller from existing aligned BAMs. This first implementation
does not port the full SCcaller FASTQ-to-BAM preprocessing pipeline; it assumes
the cell BAMs and matched bulk BAM already exist.

Required inputs:

- `cell_metadata` in the patient sheet, with one row per cell BAM.
- `bulk_bam` in the patient sheet.
- `ref_fasta` in the patient sheet.
- `--sccaller_script /path/to/sccaller_v2.0.0.py`.
- `--sccaller_hsnp_vcf /path/to/bulk.hsnp.biallelic.dbsnp.vcf`.

The SCcaller process uses its own Conda environment:

```text
envs/sccaller.yml
```

It is separate from the MonoVar environment and includes Python 2.7, `pysam`
0.15.1, `numpy`, and `samtools`, matching the SCcaller core requirements.

Example:

```bash
nextflow run main.nf \
  -profile conda \
  --patients configs/patients.sccaller.tsv \
  --run_sccaller true \
  --run_comparison true \
  --sccaller_script /path/to/SCcaller-core/sccaller_v2.0.0.py \
  --sccaller_hsnp_vcf /path/to/bulk.hg38.hsnp.biallelic.dbsnp.vcf
```

Outputs:

```text
results/<patient_id>/raw_calls/sccaller/<cell_id>.sccaller.vcf
results/<patient_id>/processed_calls/sccaller/<cell_id>.somatic.snv.vcf
results/<patient_id>/comparable_calls/sccaller/<cell_id>.sccaller.comparable.vcf.gz
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


## Step 5: VEP population/COSMIC filtering

This optional stage annotates each comparable MonoVar VCF with VEP and keeps variants that either:

- have `MAX_AF <= --max_population_af` (default `0.001`), or
- have a COSMIC identifier in `Existing_variation` when `--keep_cosmic true`.

Run on the small test region:

```bash
nextflow run main.nf \
  -profile conda \
  -resume \
  --run_monovar true \
  --run_vep_filter true \
  --monovar_script /tzu-share-2/users/students/cornelusp/monovar/monovar/src/monovar.py \
  --monovar_threads 2 \
  --monovar_region chr1:1-1000000
```

Expected outputs:

```text
results/<patient_id>/vep/monovar/<cell_id>.monovar.vep.tsv
results/<patient_id>/vep/monovar/<cell_id>.monovar.population_cosmic.summary.tsv
results/<patient_id>/final_calls/monovar/<cell_id>.monovar.population_cosmic.vcf.gz
```


### VEP performance notes

The VEP stage follows the main nf-core/Sarek and Ensembl VEP recommendations that matter for this project:

- run VEP in offline/cache mode using `--vep_cache`, `--vep_species`, `--vep_cache_version`, and `--genome`;
- split each comparable VCF into VEP chunks with `--vep_sites_per_chunk`;
- run chunk-level VEP processes in parallel;
- use VEP internal parallelism within each chunk with `--vep_forks`;
- expose `--vep_buffer_size` so larger runs can trade memory for speed;
- optionally pass a FASTA with `--vep_fasta` if the cache does not auto-detect it.

This mirrors the scatter-gather idea used by nf-core's `vcf_annotate_ensemblvep_snpeff` subworkflow: chunks are annotated independently, then merged before applying the population/COSMIC filter.


## Step 6: tables and matrices

After final VEP population/COSMIC filtering, the pipeline builds app-ready outputs:

```text
results/<patient_id>/tables/<patient_id>.monovar.final.long.tsv
results/<patient_id>/tables/<patient_id>.monovar.final.summary.tsv
results/<patient_id>/matrices/<patient_id>.monovar.final.binary_matrix.tsv
results/<patient_id>/matrices/<patient_id>.monovar.final.altread_matrix.tsv
results/<patient_id>/matrices/<patient_id>.monovar.final.refread_matrix.tsv
```

The long table keeps one row per retained variant per CTC and includes `ctc_id`, `var_id`, `refread`, `altread`, genotype fields, and any gene/consequence annotation that is present in the final VCF.


## Step 7: QC and filter impact reports

The pipeline also writes compact audit reports:

```text
results/<patient_id>/reports/<patient_id>.monovar.filter_settings.tsv
results/<patient_id>/reports/<patient_id>.monovar.prefilter_depth_alt_qc.tsv
results/<patient_id>/reports/<patient_id>.monovar.filter_impact.tsv
```

These reports record the active thresholds, summarize how many calls were removed, and provide per-cell pre-filter distributions for total depth, alternative reads, and VAF across all MonoVar candidate sites and non-reference sites.

## Step 8: MultiQC HTML report

Following the nf-core/Sarek reporting pattern, the pipeline converts the CN-PUP QC tables into MultiQC custom content and writes a standalone HTML report:

```text
results/<patient_id>/reports/multiqc_custom/<patient_id>.*_mqc.json
results/<patient_id>/multiqc/<patient_id>.multiqc_report.html
results/<patient_id>/multiqc/<patient_id>.multiqc_report_data/
```

The report currently contains filter settings, total filter impact, per-cell final variant burden, pre-filter non-reference burden, median depth, median alt reads, median VAF, and per-cell retained/removed call counts. The plain TSV reports are still kept in `results/<patient_id>/reports/` for scripting and auditability.


The pipeline also writes Nextflow execution reports, similar to nf-core/Sarek:

```text
results/pipeline_info/execution_report.html
results/pipeline_info/execution_timeline.html
results/pipeline_info/execution_trace.txt
results/pipeline_info/pipeline_dag.dot
```

Future Sarek-like additions that require extra processes are alignment-level QC (`samtools stats`, `mosdepth` coverage distributions/contig coverage) and VCF-level QC (`bcftools stats`) for the final VCFs.


## Optional Sarek-like BAM and VCF QC

To add Sarek-like native MultiQC sections, enable BAM QC:

```bash
nextflow run main.nf \
  -profile conda \
  -resume \
  --run_monovar true \
  --run_vep_filter true \
  --run_bam_qc true \
  --monovar_script /path/to/MonoVar/src/monovar.py
```

This adds:

```text
results/<patient_id>/reports/samtools/<cell_id>/<cell_id>.samtools.stats.out
results/<patient_id>/reports/mosdepth/<cell_id>/<cell_id>.mosdepth.*
results/<patient_id>/reports/bcftools/<cell_id>.monovar.bcftools_stats.txt
```

`bcftools stats` runs on final VCFs whenever VEP filtering is enabled. `samtools stats` and `mosdepth` require `--run_bam_qc true` because full BAM coverage QC can be expensive. For exome/target panels, set `--mosdepth_by targets.bed` to get target-region coverage distributions.


## SNV biology report

When `--run_vep_filter true` is active, the pipeline also renders a patient-level SNV report from the final long table and matrices:

```text
results/<patient_id>/snv_report/<patient_id>.snv_report.html
results/<patient_id>/snv_report/<patient_id>.monovar.maf
```

The report reuses established R/Bioconductor visualization packages rather than hand-rolling every plot:

- `maftools` for cancer-style MAF summaries, oncoplots, Ti/Tv and VAF-style views;
- `UpSetR` for CTC overlap sets;
- `ComplexHeatmap` for clustered binary mutation matrices and Jaccard heatmaps;
- `circlize` for a Circos-like genome overview.

Set `--run_snv_report false` to skip this optional HTML report if the R reporting environment is not needed for a quick pipeline test.
