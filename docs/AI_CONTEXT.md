# AI Context Handoff

Last reviewed: 2026-07-01

2026-07-01 update: added `BUILD_CALLER_SCORECARD` (per-caller reliability
scorecard + statistical tests, gated by `params.run_caller_scorecard`) to
close the gap where `06_overlap_qc.py` only reports raw pairwise/all-way
Jaccard, not per-caller aggregate stats. This replaces a one-off analysis
that used to live in `Documents\CTC_logiciel\appv2\` (R Shiny workspace) --
that workspace was lost during an in-progress reorg into `Documents\THESIS\`
and is not recoverable; the equivalent analysis now lives in this pipeline
stage instead. See `THESIS\analysis\caller_comaprison\reconstructed_methods.md`
for the full history of what was lost/found during that reconstruction.

This document is a compact handoff for future AI agents working in this
repository. It summarizes the pipeline, the research objectives, the current
state of the codebase, and the practical traps that matter most.

## Project In One Sentence

CN-PUP (`CeciNestPasUnePipeline`) is a thesis-oriented Nextflow pipeline for
single-cell CTC DNA sequencing: it calls or imports SNVs, filters and normalizes
caller outputs, subtracts germline/background variants, compares callers and
samples, builds mutation matrices, and prepares downstream reports and
phylogenetic inputs, especially DelSIEVE.

## Scientific Objective

The biological goal is to evaluate whether circulating tumor cells capture
intra-patient tumor heterogeneity and how CTC SNVs relate to matched primary
tumor, metastasis, leukocyte, and/or bulk blood evidence.

The thesis questions driving the pipeline are:

- Do CTC SNVs overlap with matched primary or metastatic tissue?
- Can CTC sequencing alone reveal clonal heterogeneity?
- Can longitudinal or larger CTC sampling recover more subclonal structure?
- How sensitive are conclusions to caller choice, background subtraction, and
  recurrence filtering?
- Can filtered SNV read-count data feed phylogenetic tools such as DelSIEVE,
  CTC-SCITE, and eventually Phylonco/BEAST-style phylodynamics?

The main caution repeated throughout the docs is that a final comparable VCF is
a filtered hypothesis set, not ground truth. Interpretation should use caller
agreement, matched background subtraction, recurrence across CTCs, tissue
support, and QC.

## Repository Map

- `main.nf`: main Nextflow workflow. This is the source of process wiring.
- `nextflow.config`: default parameters, process executor/profile settings,
  Nextflow execution reports, weblog config.
- `configs/`: example and active patient sheets, caller manifests, params YAMLs,
  DelSIEVE XML template, tissue VCF manifests.
- `bin/`: helper scripts used by Nextflow.
- `bin/mutation_matrix/`: vendored downstream caller-comparison scripts from
  the earlier `mutation_matrix_pipeline`; these are central to long-table,
  matrix, overlap, CTC-SCITE, and flowchart generation.
- `assets/`: R Markdown report templates for SNV and final patient reports.
- `envs/`: Conda environments by tool family.
- `docs/`: thesis documentation, literature notes, example reports, current
  Patient 03 run setup, tree renders, flowcharts, and generated comparison
  artifacts.
- `interface/`: static HTML/JS input builder that creates patients/cells/caller
  TSVs and a launch command. It does not run Nextflow or upload data.
- `tools/DelSIEVE/`: vendored DelSIEVE source/distribution and examples.
- `results/`, `work/`, `.nextflow*`, `output/`: generated runtime artifacts.

## Key Documentation Read

- `README.md`: current user-facing usage and modes.
- `nextflow_pipeline_roadmap.txt`: original development plan and rationale.
- `docs/pipeline_overview.md`: minimal workflow summary.
- `docs/pipeline_thesis.qmd`: thesis-facing technical and biological
  explanation. This is the richest source of intent.
- `docs/cnpup_pipeline_flowchart.mmd`: compact workflow graph.
- `docs/Bibliography/*notes.txt`: literature notes on CTC/single-cell SNV
  analysis, recurrence filters, SIEVE/DelSIEVE candidate-site logic, and
  phylodynamic extensions.
- `docs/tzu_run/*`: concrete Tzu Patient `PAT-2026-03-00003` input mapping and
  run files.
- `docs/inspi_reports/BFX-2026-04-00007.orange.txt`: inspiration for a polished
  clinical-style report; not a CN-PUP output.
- `docs/example_reports/*`: static and R Markdown examples for final reports.
- `docs/tree_renders/*`: rendered Patient 03 DelSIEVE tree outputs and summary.

Binary/generated docs also present:

- PDF literature and ORANGE reports under `docs/Bibliography/` and
  `docs/inspi_reports/`.
- Many CNV JPEGs under `docs/CNVs/`.
- Rendered Quarto HTML and supporting libraries under `docs/pipeline_thesis*`.
- Excel workbooks in `docs/tzu_run/`; both contain 17-row sample manifests.

## Current Workflow Shape

The workflow starts from a patient sheet and branches depending on params:

1. `CHECK_INPUTS`
2. Optional caller generation:
   - `RUN_MONOVAR` then `SPLIT_MONOVAR`
   - `PREPARE_MONOVAR_WBC_SUBSET`, `RUN_MONOVAR_WBC_SUBSET`,
     `SPLIT_MONOVAR_WBC_SUBSET` for `germline_mode = monovar_wbc_subset`
   - `RUN_SCCALLER_FROM_BAM` then `STANDARDIZE_INTERNAL_SCCALLER`
   - `STANDARDIZE_EXTERNAL_CALLERS`
3. Background prep and subtraction:
   - `PREPARE_PRECOMPUTED_GERMLINE`
   - `PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION`
   - `FILTER_MONOVAR_AND_SUBTRACT_GERMLINE`
4. Optional VEP:
   - `MAKE_VEP_INPUT_CHUNKS`
   - `VEP_ANNOTATE_CHUNK`
   - `MERGE_VEP_AND_POPULATION_COSMIC_FILTER`
5. Tables, matrices, QC, reports:
   - `BUILD_MUTATION_MATRICES`
   - `BUILD_CALLER_COMPARISON`
   - `BUILD_CALLER_SCORECARD` (gated by `run_caller_scorecard`; per-caller
     reliability scorecard + statistical tests, consumes
     `BUILD_CALLER_COMPARISON.out` directly)
   - `BUILD_FINAL_REPORT_DATA`
   - `CN_PUP_FINAL_REPORT`
   - `SNV_HTML_REPORT`
   - `SUMMARIZE_FILTERS_QC`
   - `MAKE_MULTIQC_CUSTOM_CONTENT`
   - `MULTIQC_REPORT` or `MULTIQC_REPORT_WITH_BAM_QC`
6. Optional alignment/VCF QC:
   - `SAMTOOLS_STATS`
   - `MOSDEPTH_QC`
   - `BCFTOOLS_STATS`
7. Optional DelSIEVE:
   - `PREPARE_DELSIEVE_INPUTS`
   - `DELSIEVE_DATA_COLLECTOR`
   - `RUN_DELSIEVE_STAGE1`
   - `DELSIEVE_TREE_ANNOTATOR`
   - `DELSIEVE_TREE_PNG_STAGE1`
   - `DELSIEVE_VARIANT_CALLER_STAGE1`
   - `DELSIEVE_GENE_ANNOTATOR_STAGE1`
   - `DELSIEVE_GENE_TREE_PNG_STAGE1`

The downstream contract is:

```text
results/<patient_id>/comparable_calls/<caller>/<cell_id>.<caller>.comparable.vcf.gz
```

Once caller outputs are in this format, the comparison branch builds merged
long tables, matrices, overlap QC, CTC-SCITE inputs, variant flowcharts, and
report data.

## Data Model

Primary patient sheet columns:

```text
patient_id
ref_fasta
monovar_bam_list
cell_metadata
germline_mode
germline_vcf
bulk_bam
leukocyte_bam
leukocyte_vcf
```

Cell metadata columns:

```text
patient_id    cell_id    bam_path    cell_type
```

External caller manifest columns:

```text
patient_id    cell_id    caller    vcf_path    sccaller_mode
```

Supported caller names:

- `monovar`
- `sccaller`
- `haplotypecaller`
- `deepvariant`

Supported/recognized germline modes include:

- `precomputed_vcf`
- `deepvariant_bulk_bam` (validation exists; full internal bulk calling is
  still planned)
- `leukocyte_vcf`
- `joint_monovar_leukocyte`
- `monovar_wbc_subset`
- `combined`
- `no_germline`

**2026-07-03: `combined` mode fixed for the native MonoVar branch.** Previously
`combined` (bulk blood + leukocyte) only worked for
`STANDARDIZE_EXTERNAL_CALLERS`/`STANDARDIZE_INTERNAL_SCCALLER` (their bash
logic already merges both sources correctly via `bcftools view -G` + concat +
sort + `norm -d exact`). The native `RUN_MONOVAR`/`FILTER_MONOVAR_AND_SUBTRACT_GERMLINE`
branch only ever checked `germline_mode == 'precomputed_vcf'` XOR
`== 'joint_monovar_leukocyte'` -- never both, never `'combined'` -- so a
patient with `germline_mode: combined` and `run_monovar: true` got **no**
germline exclusion applied to MonoVar's own calls at all (this is why
`patients.patient_03_wxs_ctc_only_compare.tsv`, which does use `combined`,
sets `run_monovar: false` and imports MonoVar externally instead -- it never
hit this gap).

Fix: `precomputed_input_ch` and `leukocyte_exclusion_input_ch` now both also
fire for `germline_mode == 'combined'`; when both fire for the same patient,
a new process `MERGE_MONOVAR_GERMLINE_EXCLUSION` (added right after
`PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION` in `main.nf`) merges the bulk and
leukocyte exclusion VCFs into one (same `bcftools view -G` + concat + sort +
`norm -d exact` pattern already used by `STANDARDIZE_EXTERNAL_CALLERS`)
before `FILTER_MONOVAR_AND_SUBTRACT_GERMLINE` runs. Patients using
`precomputed_vcf`-only or `joint_monovar_leukocyte`-only still pass through
unchanged (implemented via `.join(..., remainder: true)` + `.branch{}` on
`PREPARE_PRECOMPUTED_GERMLINE.out`/`PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION.out`,
routing to `both`/`bulk_only`/`leuko_only`).

**This is new code, executed once (2026-07-06) and found to be correct at
the mechanics level, but with a sample-identity error in that first run --
see below.** Traced carefully against the existing
`STANDARDIZE_EXTERNAL_CALLERS` pattern and Nextflow channel semantics before
running. First run: `configs/patient01_ctc6_compare/` was switched to
`germline_mode: combined` for exactly this reason. Check
`results/patient01_ctc6_compare/germline_exclusion/*.combined_germline_exclusion.log`
and the resulting VCF's record count after running -- it should be roughly
(bulk count + leukocyte count − overlap), not wildly larger or smaller, as a
basic sanity check that the merge worked as intended.

**2026-07-06 correction: the leukocyte sample was wrong.** `SRR8617667` was
used as "patient01's leukocyte" based on it being the 7th sample column in
`vcf_bench/VCFs/patient01_ctcleuko_monovar.vcf` (labelled "leukocyte" in that
analysis's own report,
`vcf_bench/VCFs/comparisons/leukocyte_monovar_analysis/patient01_leukocyte_monovar_report.md`).
Checking SRA metadata directly (`Documents/Playground/SraRunTable.csv`)
shows `SRR8617667` is actually **Patient 02's CTC9** (`isolate=Patient 02`,
`tissue=Circulating tumor cells 9`) -- a pre-existing mislabeling in the
original `vcf_bench` analysis, inherited here without independently
cross-checking SRA metadata. Patient01's real leukocyte sample is
**`SRR8617610`** (`isolate=Patient 01`, `tissue=Single leokucyte [sic]`,
`Assay Type=WXS`). All three `configs/patient01_ctc6_compare/` files
(`bams.txt`, `cells.tsv`, `patients.tsv`) have been corrected to
`SRR8617610`. The 2026-07-06 Pass 1 run (MonoVar + germline exclusion) used
the wrong sample -- it subtracted Patient 02 tumor-cell variants from
Patient 01's CTC calls, not a real background/germline set -- and needs to
be rerun with the corrected config before Pass 2. Also check
`vcf_bench/VCFs/comparisons/leukocyte_monovar_analysis/` and any thesis text
built on it (including drafted methods text) for the same mislabeling before
citing its numbers.

Also worth double-checking: all patient01 samples used here are confirmed
WXS (not WGS) per the same SraRunTable -- `SRR8617606` (bulk),
`SRR8617611/608/609/612/613/526` (the 6 CTCs), and `SRR8617610` (leukocyte)
all show `Assay Type=WXS`. Note `PRJNA448888` also has WGS runs for some of
these same patient01 samples under different SRR accessions (e.g.
`SRR8617596`, `SRR8617605`) -- don't mix them in by accident when sourcing
future samples from this project.

That config needs a **two-pass run**: `params_leukocyte_prep.yml` first
(MonoVar-only, produces `split_calls/monovar/leuko.monovar.split.vcf`, which
`patients.tsv`'s `leukocyte_vcf` field points at), then `params.yml`
(`-resume` reuses the first pass's MonoVar result since both share the same
`monovar_bam_list`/`patients.tsv`). The external/SCcaller `combined`-mode
subtraction needs `leukocyte_vcf` to already be a real file, unlike the
native MonoVar branch which derives it live -- see the comment block at the
top of `configs/patient01_ctc6_compare/params.yml`.

## Important Parameters

Defaults live in `nextflow.config`; many are duplicated in params YAML files.

Filtering defaults:

```text
min_total_depth = 20
min_alt_reads = 4
min_vaf = 0.0
max_population_af = 0.001
keep_cosmic = true
```

Important run switches:

```text
run_monovar
run_sccaller
run_external_callers
run_comparison
run_caller_scorecard
run_vep_filter
run_snv_report
run_final_report
run_bam_qc
run_delsieve_prep
run_delsieve
run_delsieve_from_existing
```

Recent uncommitted changes added `monovar_targets_bed`, which is validated in
`CHECK_INPUTS` and passed to `samtools mpileup` as `-l` in both standard and
WBC-subset MonoVar processes. The Tzu params file points this at an HMF coding
coverage BED.

## Active Patient 03 Setups

### Public Patient 03 WXS CTC-only comparison

Key files:

```text
configs/patients.patient_03_wxs_ctc_only_compare.tsv
configs/caller_vcfs.patient_03_wxs_ctc_only_compare.tsv
configs/params.patient_03_wxs_ctc_only_compare.yml
configs/tissue_vcfs.patient_03_wxs_meta_leuko.tsv
```

Observed manifest counts:

- `configs/patient_03_wxs_ctc_only_cells.tsv`: 21 CTC rows.
- `configs/caller_vcfs.patient_03_wxs_ctc_only_compare.tsv`: 21 MonoVar,
  24 HaplotypeCaller, 24 DeepVariant rows.
- Tissue manifest includes metastasis and leukocyte VCFs for HaplotypeCaller
  and DeepVariant.

This setup currently uses `germline_mode = combined` with a DeepVariant
bulk/germline VCF and a MonoVar leukocyte split VCF.

**2026-07-01: SCcaller deferred for this patient (author's decision).**
`bulk_bam` in `configs/patients.patient_03_wxs_ctc_only_compare.tsv` was
updated to point at `SRR8617664.sorted.bam` (same SRR already used as the
DeepVariant germline/bulk VCF source) since it's harmless/more accurate
either way, but `run_sccaller` stays `false` in
`configs/params.patient_03_wxs_ctc_only_compare.yml`: unlike patient01
(which has a pre-built hSNP VCF for its bulk sample, see below),
patient_03's bulk sample (SRR8617664) has no equivalent hSNP VCF yet, and
building one needs a GATK HaplotypeCaller GVCF run + filtering (recipe in
`vcf_bench/SCcaller/sccaller_pipeline_1_v1.0.0.sh` sub-step 12) that wasn't
judged worth doing yet. `sccaller_script` is filled in (same tool as
patient01) for whenever this is revisited; `sccaller_hsnp_vcf` is
deliberately blank -- do not set `run_sccaller: true` here without setting
it to a real, verified hSNP VCF first. Patient_03's scorecard therefore
stays at 3 callers (MonoVar, HaplotypeCaller, DeepVariant) for now.

### Public Patient 01 6-CTC comparison (added 2026-07-01)

Key files:

```text
configs/patient01_ctc6_compare/patients.tsv
configs/patient01_ctc6_compare/caller_vcfs.tsv
configs/patient01_ctc6_compare/params.yml
configs/patient01_ctc6_compare/bams.txt
configs/patient01_ctc6_compare/cells.tsv
```

Reproduces, through CN-PUP natively, the same 6-CTC patient01 comparison
(`SRR8617611/608/609/612/613/526`) that the old `vcf_bench` ad-hoc scripts
originally produced (`vcf_bench\VCFs\comparisons\caller_choice_objective_report\`).
Intended as the first validation run for `BUILD_CALLER_SCORECARD`: if the
new pipeline-native scorecard roughly reproduces the old report's numbers
(e.g. MonoVar highest on `shared_all_n_ctcs`), that validates the ported
logic before trusting it on new data.

MonoVar and SCcaller run natively; HaplotypeCaller and DeepVariant are
imported. **The HaplotypeCaller/DeepVariant VCF paths and the germline VCF
path in these configs are inferred by analogy from the one confirmed
patient01 path pattern in `vcf_bench/variant_caller_comparison.qmd`
(`BULK_DV=.../sarek_variantcalling/variant_calling/deepvariant/patient_01/
<SRR>/<SRR>.deepvariant_snpEff.ann.vcf.gz`), not directly confirmed to
exist for all 6 cells -- verify over VPN before running.**

`sccaller_script` and `sccaller_hsnp_vcf`, unlike the HC/DV paths above,
**are** directly sourced from `vcf_bench/variant_caller_comparison.qmd`'s
`## SCcaller` section (the actual historical run for this exact patient),
not inferred: `sccaller_script =
/tzu-share-2/users/students/cornelusp/sccaller/SCcaller-core/sccaller_v2.0.0.py`,
`sccaller_hsnp_vcf =
/tzu-share-2/users/students/cornelusp/sccaller/ht/SRR8617606.hg38.hsnp.biallelic.dbsnp155common.vcf`
(a pre-built heterozygous-SNP VCF for the SRR8617606 bulk sample). Still
worth a quick existence check on the server since it lives outside any
git-tracked location and could have been cleaned up.

Kept minimal on purpose (`run_vep_filter`/`run_snv_report`/`run_final_report`/
`run_bam_qc`/`run_delsieve*` all `false`) so this run is fast and focused on
validating the scorecard stage, not a full patient report.

**2026-07-03: switched to `germline_mode: combined` (bulk + leukocyte).**
Patient01 has a leukocyte sample already jointly MonoVar-called with the 6
CTCs in `vcf_bench/VCFs/patient01_ctcleuko_monovar.vcf` -- missed in the
original config. Added it to `bams.txt`/`cells.tsv` (`cell_id: leuko`) and
switched `patients.tsv` off `precomputed_vcf` onto `combined`. This required
a `main.nf` fix (see the germline-mode section above) and makes this a
**two-pass run**: `configs/patient01_ctc6_compare/params_leukocyte_prep.yml`
first, then `params.yml`. Read the comment block at the top of `params.yml`
before running -- getting the pass order wrong will make the external-caller/
SCcaller side of the subtraction silently do less than intended (or the run
may fail outright if `leukocyte_vcf` doesn't resolve).

**2026-07-06 correction: `SRR8617610` is the correct leukocyte sample, not
`SRR8617667`.** `SRR8617667` (used above and in the first, now-invalid run)
is actually Patient 02's CTC9 per SRA metadata -- see the fuller note further
up this file. All three config files now use `SRR8617610`; the 2026-07-06
run needs to be redone from Pass 1.

### Tzu Patient PAT-2026-03-00003

Key files:

```text
docs/tzu_run/tzu_patient_03_samples.xlsx
docs/tzu_run/cnpup_patient_03_sample_manifest.csv
docs/tzu_run/cnpup_patient_03_sample_manifest_with_bams.csv
docs/tzu_run/cnpup_patient_03_sample_manifest_with_bams.xlsx
docs/tzu_run/PAT-2026-03-00003.cells.tsv
docs/tzu_run/PAT-2026-03-00003.monovar_bams.txt
docs/tzu_run/patients.PAT-2026-03-00003.monovar.tsv
docs/tzu_run/params.PAT-2026-03-00003.monovar.yml
docs/tzu_run/run_PAT-2026-03-00003_monovar.sh
```

Observed manifest counts:

- 15 CTC samples and 1 WBC sample.
- `patients.PAT-2026-03-00003.monovar.tsv` uses
  `germline_mode = monovar_wbc_subset`.
- The params file runs MonoVar, comparison, and final report; VEP, SNV report,
  and BAM QC are currently disabled.
- Weblog is enabled at `http://localhost:8605`.
- `monovar_targets_bed` restricts MonoVar/mpileup to the HMF coding coverage
  panel BED.

## Helper Script Roles

Top-level `bin/`:

- `split_monovar_vcf.py`: splits multi-sample MonoVar VCF into per-cell VCFs
  and maps VCF sample names to `cell_id` using SRR tokens or row order.
- `build_mutation_matrices.py`: older MonoVar final VCF to long table and
  matrices.
- `extract_existing_vep_from_vcf.py`: extracts VEP/CSQ-like annotations already
  embedded in VCFs.
- `filter_vcf_by_vep.py`: filters VCFs based on VEP max AF and COSMIC rescue.
- `merge_vep_tables.py`: concatenates chunk-level VEP TSVs.
- `summarize_filters_qc.py`: writes filter settings, depth/alt/VAF summaries,
  and filter-impact TSVs.
- `make_multiqc_custom_content.py`: converts CN-PUP QC tables into MultiQC
  custom content JSON.
- `prepare_delsieve_inputs.py`: selects candidate sites from long tables, piles
  original CTC BAMs at those sites, parses mpileup base strings, and writes
  DelSIEVE read-count/cell-name/mutation-map style inputs.
- `build_final_report_data.py`: produces compact tables for the final
  patient-level R Markdown report, including CTC QC, sharing spectra, top
  variants, tissue overlap, rarefaction, and warnings.

`bin/mutation_matrix/`:

- `02_process_raw_calls.sh`: standardizes raw caller VCFs into comparable
  filtered/normalized/subtracted callsets.
- `03_vcf_to_long_table.py`: VCF to per-sample/caller long table.
- `04_merge_long_tables.py`: merges long-table shards.
- `04_fill_missing_annotations.py`: fills missing annotation fields by `var_id`
  from rows where another caller has better annotation.
- `05_build_matrices.py`: builds binary, VAF, alt-read, and depth matrices.
- `06_overlap_qc.py`: pairwise and all-way overlap summaries.
- `07_variant_flowcharts.py`: per-cell/caller attrition counts and Mermaid
  flowcharts.
- `08_caller_scorecard.py` (added 2026-07-01): reads the comparison
  matrices and computes a per-caller reliability scorecard -- union/private/
  shared-across-all-CTCs variant counts, phylo-informative fraction,
  pairwise Jaccard/overlap-coefficient, cross-caller-support fraction,
  optional background(leukocyte/tissue)-overlap fraction, median depth/alt/
  VAF -- plus Kruskal-Wallis/Chi-square/Fisher-exact tests between callers
  and an auto-generated `scorecard_summary.md` conclusion. Generalized to N
  callers/M cells (infers both from the matrix column names, e.g.
  `CTC3__monovar`), no hardcoded caller list. CTC-vs-background sample
  classification uses the same `"ctc" in sample_id.lower()` heuristic as
  `build_final_report_data.py`. Deliberately has **no cancer-gene
  filtering** -- the scorecard measures caller reliability via cross-CTC
  sharing, not oncology-driver-gene enrichment; gene-level filtering would
  shrink N and could hide the sharing signal it's meant to measure. Needs
  `scipy` (added to `envs/python_reporting.yml`).

## Reports And Figures

`assets/snv_report.Rmd` uses `maftools`, `UpSetR`, `ComplexHeatmap`, `circlize`,
and `ggplot2` to summarize final MonoVar SNVs.

`assets/cnpup_final_report.Rmd` is the thesis-facing final patient report. It
organizes results around executive summary, QC, tissue representativity,
heterogeneity, longitudinal evolution, sampling saturation, per-CTC pages, and
methods/provenance.

The ORANGE report text/PDF in `docs/inspi_reports/` is only inspiration for
report polish and QC/driver-style presentation. It is not generated by CN-PUP.

## DelSIEVE Context

The DelSIEVE integration is not just a matrix export. DelSIEVE expects raw
read-count style inputs at selected candidate sites:

- total coverage per site/cell;
- counts for the three non-reference bases, ordered highest to lowest;
- optional background wildtype/acquisition-bias data;
- cell names;
- configured BEAST XML template.

The docs stress that SIEVE/DelSIEVE do not scan the whole genome/exome by
themselves during inference. Candidate-site definition is a separate critical
choice. The intended main strategy is to use recurrent CTC variants and/or
tissue-supported variants as a high-confidence backbone, with private variants
or DataFilter-style candidates reserved for sensitivity analyses.

Current rendered tree notes:

- Patient 03 DelSIEVE tree renders exist under `docs/tree_renders/`.
- Negative branch lengths in the MCC/NEXUS tree were set to zero for plotting.
- Tip labels were restored from the NEXUS `Taxlabels` block.
- Mutation labels use `gene_fsa` / `event_type_fsa` where available.
- Branch labels are compact summaries, and branches with fewer than 2 FSA
  variants are not labelled.

## Current Generated Outputs

`results/` currently contains mostly Nextflow pipeline-info reports from a
failed validation run. The latest `.nextflow.log` shows `CHECK_INPUTS
(patient_01)` failed because `/tzu-share-2/.../Homo_sapiens_assembly38.fasta`
does not exist in this local Windows/WSL environment. This is expected when
server paths are not mounted.

`docs/caller_comparison/` contains generated comparison artifacts, including:

- merged annotated long table (`.tsv.gz`);
- matrices for binary calls, alt reads, total depth, and VAF;
- CTC-SCITE input files;
- `reports/overlap_qc.tsv` with 2417 lines.

`output/playwright/` contains screenshots of example reports.

## Environment And Execution Notes

- The repository is on Windows at
  `C:\Users\CornelusPalsma\Documents\CeciNestPasUnePipeline`.
- Many configured data paths are Linux server paths under `/tzu-share-2/...`.
  They will not validate locally unless WSL/server mounts expose those paths.
- The latest logs use WSL-style work dirs under `/mnt/c/...`.
- Conda environments are split by tool family:
  - `monovar_py2.yml`: Python 2.7 MonoVar.
  - `sccaller.yml`: Python 2.7 SCcaller.
  - `bcftools.yml`: bcftools/samtools/htslib.
  - `vep.yml`: Ensembl VEP plus bcftools/samtools.
  - `alignment_qc.yml`: samtools/mosdepth.
  - `delsieve_prep.yml`: Python plus samtools/libdeflate.
  - `delsieve_tree_plot.yml`: FigTree.
  - `python_reporting.yml`: pandas/numpy/multiqc.
  - `snv_report.yml`: R reporting packages.

## Git State At Review Time

The worktree was dirty before this document was added. Existing modified files:

- `configs/params.example.yml`
- `docs/pipeline_thesis.qmd`
- `docs/tzu_run/params.PAT-2026-03-00003.monovar.yml`
- `main.nf`
- `nextflow.config`

There are many untracked generated/report files under `docs/`, `output/`,
`tools/`, and `.playwright-cli/`. Do not revert or delete them unless the user
explicitly asks.

The visible diff in code/config is mainly the addition of `monovar_targets_bed`
support and a large thesis-document extension.

## Practical Gotchas

- `rg` is not installed in this environment. Use PowerShell `Get-ChildItem`,
  `Select-String`, and `Get-Content`, or install/use another search tool if
  approved.
- Do not treat `patients.tsv` as the active analysis config by default; it
  points at server resources and lacks `cell_metadata`. Prefer params files for
  reproducible runs.
- `germline_mode = no_germline` is intended for producing clean CTC-only
  MonoVar split VCFs first, then subtracting backgrounds later in comparison.
- `combined` mode requires both a germline/bulk VCF and leukocyte evidence.
- For local dry runs, expect validation failures unless `/tzu-share-2` paths
  are accessible.
- **SCcaller on `Homo_sapiens_assembly38.fasta` needs `sccaller_head`/
  `sccaller_tail` set (e.g. `1`/`22`), or it silently fails.** This reference
  includes hundreds of HLA/ALT decoy contigs with `*`/`:` in their names
  (e.g. `HLA-DRB1*15:01:01:04`) that break `samtools mpileup`'s region
  parsing. Confirmed 2026-07-06 on patient01/CTC2: ran ~4.2 hours hitting
  these errors, never produced a VCF, and (before this was added) failed
  later with a confusing "cannot open file ... for reading" from the
  downstream `awk` step rather than a clear error at the source. Both params
  default to unset (scans every contig) -- set them per-patient based on
  that reference's actual chromosome ordering.
- Avoid editing generated Quarto/HTML support libraries under
  `docs/pipeline_thesis_files/`.
- `tools/DelSIEVE/` is vendored third-party code; changes there should be rare
  and deliberate.
- The thesis docs mix implemented behavior and planned extensions. Check
  `main.nf` before assuming a documented future feature exists.

## Best Next Engineering Tasks

1. Stabilize a small local/server smoke test using a params file, not default
   `configs/patients.tsv`.
2. Verify the recent `monovar_targets_bed` change on the intended WSL/server
   environment.
3. Make Patient `PAT-2026-03-00003` run end-to-end through MonoVar WBC-subset,
   comparison, and final report.
4. Promote primary/metastasis VCFs from flowchart-only context into first-class
   comparison samples/matrices if thesis analysis needs direct CTC-vs-tissue
   overlap.
5. Add explicit phylogeny candidate-site outputs:
   - candidate sites;
   - presence/absence matrix;
   - missingness matrix;
   - DelSIEVE/CTC-SCITE-ready exports;
   - reason/confidence tier per variant.
6. Keep recurrence/tissue-support confidence tiers separate from low-confidence
   private singleton CTC variants.

## Common Commands

Validation/default run:

```bash
nextflow run main.nf --patients configs/patients.tsv
```

Preferred Patient 03 comparison:

```bash
nextflow run main.nf \
  -profile conda \
  -params-file configs/params.patient_03_wxs_ctc_only_compare.yml \
  -resume
```

Tzu Patient MonoVar run:

```bash
nextflow run main.nf \
  -plugins nf-weblog \
  -profile conda \
  -params-file docs/tzu_run/params.PAT-2026-03-00003.monovar.yml \
  -resume
```

DelSIEVE from existing comparison long table:

```bash
nextflow run main.nf \
  -profile conda \
  -params-file configs/params.patient_03_wxs_ctc_only_delsieve_existing.yml \
  -resume
```

