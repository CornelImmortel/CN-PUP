# Repository Guidelines

## Project Structure & Module Organization

This repository contains a Nextflow workflow for CTC single-cell variant calling and filtering. The main workflow is `main.nf`; runtime defaults and profiles live in `nextflow.config`. Helper scripts are under `bin/`, including Python utilities and the `bin/mutation_matrix/` shell/Python pipeline. Reusable report templates are in `assets/`. Example and patient-specific inputs live in `configs/`; prefer adding `.example` files for shareable templates and keep local real inputs out of commits where possible. Documentation and diagrams are in `docs/`. The static input-builder UI is in `interface/` and can be opened directly via `interface/index.html`. Treat `.nextflow/`, `work/`, `results/`, and `.nextflow.log*` as generated run artifacts.

## Build, Test, and Development Commands

- `nextflow run main.nf`: validate the default `configs/patients.tsv` input sheet.
- `nextflow run main.nf --patients configs/patients.tsv`: run with an explicit patient sheet.
- `nextflow run main.nf -profile conda --run_monovar true --monovar_script /path/to/monovar.py`: run MonoVar-backed calling with Conda environments enabled.
- `nextflow run main.nf -profile conda --run_external_callers true --run_comparison true --caller_vcfs configs/caller_vcfs.tsv`: compare existing caller VCFs.
- `python bin/<script>.py --help`: inspect helper script arguments before changing or invoking them.

## Coding Style & Naming Conventions

Use existing style in each file. Keep Nextflow process and workflow names descriptive and uppercase where already established. Python helpers should use `snake_case` functions, explicit argument parsing, and clear TSV/VCF column names. Shell scripts should be POSIX/Bash compatible, fail early, and quote paths. Keep generated outputs under `results/` and intermediate data under `work/`, not beside source files.

## Testing Guidelines

Use the small fixtures in `test_data/` and `.example` config files for smoke tests. For workflow changes, run at least an input-validation pass and the smallest relevant mode, then inspect `results/pipeline_info/execution_trace.txt` and task logs on failure. For helper scripts, add or update tiny TSV/VCF fixtures and verify command-line behavior directly with `python bin/script.py ...`.

## Commit & Pull Request Guidelines

Recent commits use short, imperative subjects such as `Fix MonoVar target BED validation syntax` and `Recognize VEP gene aliases for DelSIEVE`. Follow that pattern: describe the behavior changed, not the work session. PRs should include the pipeline mode tested, exact command used, key output paths, and any data/config assumptions. Include screenshots only for `interface/` or report-rendering changes.

## Security & Configuration Tips

Do not commit patient-identifying data, full BAM/VCF outputs, local tool paths, or large Nextflow artifacts. Prefer `configs/*.example` and documented parameters in `nextflow.config` for reproducible configuration.
