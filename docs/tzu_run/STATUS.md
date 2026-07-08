# CN-PUP Patient 03 / Patient 04 CTC-subset runs — status summary

Last updated: 2026-07-03 (paused for the weekend here — user returning in ~2 days)

## Pick up here

Paste this whole file back to Claude and say what you want to do next.
Quick orientation on exactly where things stand:

1. **Patient 03**: pipeline succeeded through the core outputs (SNV
   report, final report, MultiQC, caller comparison) under a *resumed*
   run. DelSIEVE stage 1 is fixed but not re-verified (was explicitly
   left on hold). A **clean from-scratch rerun (no `-resume`)** was
   planned for full validation but **it's unclear whether it was ever
   started/finished** — check with the user first before assuming the
   analysis below reflects a clean run.
2. Results for Patient 03 were synced locally to
   `THESIS/analysis/tzu_runs/results_good_ctcs/PAT-2026-03-00003/` and
   analyzed in depth — see "Results analysis: Patient 03" below and the
   published report artifact.
3. **Open methodological question, not yet acted on**: whether the
   `min_total_depth=20`/`min_alt_reads=4` filter is too strict — see
   "Open question: filter thresholds" below. Two concrete follow-ups
   were floated but not started: (a) add the KCNJ5/CTC7 case study to
   the report artifact, (b) actually rerun with a looser DP threshold
   and diff the shared-variant results.
4. **Patient 04**: was running fresh on the server (`bash
   docs/tzu_run/run_Pat_04_CTCs_monovar.sh`) when this session paused.
   **Completion status unknown** — check `results_good_ctcs/PAT-2026-03-00004/`
   for how far it got (compare directory contents against the Patient 03
   layout listed in "Reference: where things live" to gauge progress —
   e.g. presence of `snv_report/`, `final_report/`, `delsieve_prep/`
   tells you roughly which stage it reached).
5. Nothing outstanding needs re-pushing — all config/pipeline fixes as
   of 2026-07-03 are committed and pushed to `origin/main`.

## Resuming after the computer was shut down

The actual pipeline runs are **unaffected** by the local machine's power
state — both are inside `screen` sessions on the remote server
(`tzu-node01`), detached with `Ctrl-A d`. `screen` sessions survive SSH
disconnects and keep running on the server regardless of what the local
laptop does. To check on them after reconnecting:
```bash
ssh cornelusp@tzu-node01
screen -ls              # lists active screen sessions
screen -r CNPUP          # reattach to see live progress / check if it finished
```

This Claude Code session itself is a local background job (its working
files live under `~/.claude/jobs/...` and `~/.claude/projects/...` on
this Windows machine) — it most likely does **not** keep actively
running while the computer is off, though the conversation history is
saved to disk and should still be there to resume once the machine is
back on. Safest assumption: nothing happens on the Claude Code side
while the computer is off; pick back up by either resuming this same
session from the Claude Code app/CLI, or — if that's not available for
whatever reason — starting a fresh session and pasting this whole file
in, then also pasting the output of the `screen -r` checks above so the
new session knows the current run state before doing anything else.

## Goal

Run the CN-PUP pipeline (MonoVar WBC-subset germline mode) for two
curated CTC subsets:

- **Patient 03** (`PAT-2026-03-00003`, prostate_cancer): 9 CTCs +
  1 WBC, from `THESIS/analysis/tzu_runs/Pat_03_CTCs.txt`.
- **Patient 04** (`PAT-2026-03-00004`, breast_cancer): 7 CTCs +
  1 WBC, from a user-supplied BFX ID list.

Work happens in two places:
- **This environment** (Windows/WSL, no `/tzu-share-2` mount, no conda):
  prepares and structurally validates configs, diagnoses failures from
  pasted server logs, edits/fixes `main.nf` and report assets, commits
  and pushes to `origin/main`.
- **The actual server** (`tzu-node01`, user-operated): pulls `main` and
  runs the pipeline for real via `bash docs/tzu_run/run_Pat_0X_CTCs_monovar.sh`.

## Patient 03 (`Pat_03_CTCs`)

**Config**: `docs/tzu_run/{patients,params,README}.Pat_03_CTCs*`,
`Pat_03_CTCs.cells.tsv`, `Pat_03_CTCs.monovar_bams.txt`,
`run_Pat_03_CTCs_monovar.sh`. Committed and pushed.

**WBC**: BFX-2026-04-00051 — corrected 2026-07-08 per direct confirmation
from the wet-lab team. A prior session had instead settled on
BFX-2026-04-00033 (matched against
`cnpup_patient_03_sample_manifest_with_bams.csv`'s `sample_type=WBC`
row), having rejected BFX-2026-04-00051 as coming from the
less-reliable "unmapped" source — that manifest-based assignment was
itself wrong. BFX-2026-04-00051's BAM path was re-verified against a
live `ls -la` on the server
(`PAT-2026-03-00003_10/alignments/dna/BFX-2026-04-00051.redux.bam`,
31.7GB — consistent with bulk WBC sequencing depth, not single-cell).
BFX-2026-04-00033 is dropped from the analysis entirely (it was never
one of the 9 "good" CTCs, only ever wired in as the WBC BAM).

**Run history / bugs found and fixed in `main.nf` and
`assets/snv_report.Rmd`** (all pushed to `origin/main`):

1. `SNV_HTML_REPORT`'s `flowcharts_dir` was resolved as a bare relative
   path (`normalizePath("variant_flowcharts")`) instead of an absolute
   path built from `Sys.getenv("PWD")` like every other input in the
   same params list. R's lazy argument evaluation + `rmarkdown::render()`
   changing the working directory meant this failed deterministically.
   Fixed to `normalizePath(file.path(Sys.getenv("PWD"), "${variant_flowcharts}"), mustWork=TRUE)`.
2. The germline-exclusion merge (`exclusion_join_ch`) used
   `join(remainder: true)` between `PREPARE_PRECOMPUTED_GERMLINE.out`
   (never fires for `germline_mode: monovar_wbc_subset`, so it's an
   entirely empty channel) and `PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION.out`.
   An entirely-empty side of a join collapses to a single `null` instead
   of one `null` per field, breaking the 9-parameter `.branch` closure
   (`Invalid method invocation doCall`). Rewrote as tag + `mix` +
   `groupTuple` to avoid the ambiguous join-remainder shape. This is a
   pre-existing bug that would also affect the original full-cohort
   Patient 03 config if it reaches that point.
3. The circos plot in `snv_report.Rmd` crashed
   (`circlize::circos.initialize` — "cell padding larger than sector
   width") on the smaller CTC subset, where some chromosomes have very
   few/small variant positions relative to others. Fixed by setting
   `cell.padding = c(0.02, 0, 0.02, 0)` in `circos.par()`, per
   circlize's own suggested fix.
4. `delsieve_applauncher`/`delsieve_beast`/`delsieve_treeannotator`
   were bare command names (`applauncher`, `beast`). `DELSIEVE_DATA_COLLECTOR`
   and `DELSIEVE_TREE_ANNOTATOR` have no `conda` directive, so they run
   with no environment activated and couldn't find these commands.
   Fixed by pointing all three at the confirmed BEAST2 install:
   `/tzu-share-2/users/students/cornelusp/beast/bin/{applauncher,beast}`.
   **Also affects the original full-cohort Patient 03 config**
   (`params.PAT-2026-03-00003.monovar.yml`), which was not touched.

**Other config fixes** (not bugs, just tuning):
- `outdir` changed from the shared `results` to a dedicated
  `results_good_ctcs`, so this subset run doesn't collide with the
  original full-cohort Patient 03 run's output.
- `weblog_url` corrected to `http://localhost:8605` — this is a fixed
  shared listener on the node that feeds an internal "All Users"
  pipeline dashboard, not a per-run port (an earlier guess of `8606` to
  avoid a perceived port clash was wrong and never showed up on the
  dashboard).

**Current state**: with all fixes applied, this run got all the way
through `CN_PUP_FINAL_REPORT`, `MULTIQC_REPORT_WITH_BAM_QC`,
`BUILD_CALLER_COMPARISON`, and `SNV_HTML_REPORT` successfully — the
core outputs (SNV report, final report, MultiQC, caller comparison
tables) exist under `results_good_ctcs/PAT-2026-03-00003/...`. It then
got stuck on the DelSIEVE branch (`applauncher: command not found`),
which is now fixed but not yet re-verified end to end.

**Explicitly on hold**: full DelSIEVE stage 1 (tree building) — user
said "let's leave DELSIEVE for now."

**In progress / decided just now**: user wants a **clean rerun without
`-resume`** for Patient 03, specifically to get full end-to-end
confidence that the fixed `main.nf` works correctly from scratch (not
just from `SNV_HTML_REPORT` onward, which is all that's actually been
exercised under the fixed code so far — everything upstream of it was
always served from cache within this session). Rationale discussed:
Nextflow's `-resume` is content-hash-based per task, so a differently-
composed CTC list from an earlier/unrelated run could **not** have
silently contaminated the current results even with `-resume` — but a
clean run is still valuable as full validation. Command to use (drop
`-resume`):
```bash
nextflow run main.nf \
  -plugins nf-weblog \
  -profile conda \
  -params-file docs/tzu_run/params.Pat_03_CTCs.monovar.yml
```
`run_Pat_03_CTCs_monovar.sh` itself still has `-resume` in it and has
not been changed — the user hasn't decided yet whether to make
no-resume the permanent default for that script.

## Results analysis: Patient 03

Results were synced from the server to
`THESIS/analysis/tzu_runs/results_good_ctcs/PAT-2026-03-00003/` and
analyzed thoroughly (following the reusable checklist in
`THESIS/analysis/tzu_runs/ANALYSIS_PROMPT.md`). Full write-up published
as an artifact: **filtering funnel → per-cell coverage/yield →
shared-vs-private variants → annotation profile → 6 notable HIGH-impact
findings → caveats**. The underlying HTML lives at
`THESIS/analysis/tzu_runs/PAT-2026-03-00003_analysis_report.html`
(re-publish with the Artifact tool to get a fresh link if the old one
has expired/changed — same file path redeploys to the same URL if the
artifact is republished in a new session, but a brand-new session has
no memory of the previous URL).

**Headline numbers** (all traceable to files under
`results_good_ctcs/PAT-2026-03-00003/`):
- 13,276 raw non-ref MonoVar calls (9 CTCs) → 422 after depth/alt-read/VAF
  filter (the big cut, −96.8%) → 359 after WBC germline subtraction →
  **335 final calls**, 240 distinct positions, 109 distinct genes.
- **Coverage is wildly uneven and explains almost all of the variant-count
  spread**: on-target (coding-panel) mean depth ranges from 0.28× (CTC3,
  zero final variants) to 7.08× (CTC7, 37 variants). The three
  worst-covered cells (CTC3/4/5) also had 3–4× higher alignment error
  rates than the well-covered ones — points to degraded/uneven WGA
  amplification in those specific cells, not just a sequencing-depth
  shortfall.
- **No variant is shared by all 9 CTCs.** 74% of variants are private to
  one cell; the most-shared variant (KCNJ5, chr11:128,916,399 G>A)
  reaches only 6 of 9 cells.
- 6 HIGH-impact stop-gained calls (ARID1A, MSH2, TGFBR2, FAT4, KMT2D,
  PABPC1). **PABPC1 is flagged as likely-benign**: its population AF is
  6.99% — 70× above the 0.1% filter threshold — and it only survived
  because of the `keep_cosmic=true` rescue; treat it with more
  skepticism than the other five.
- The gene panel (`HMF CoverageCodingPanel.38`) is a general pan-cancer
  panel, not prostate-specific — recurrent genes like TG/PABPC1/PREX2
  shouldn't be over-read as biologically significant without checking
  panel design.

## Open question: filter thresholds (raised by user, not yet acted on)

User asked whether the `min_total_depth=20`/`min_alt_reads=4` filter is
too strict and might be hiding additional shared (clonal) variants.
Investigated with a concrete example — the KCNJ5 shared variant's *raw*
(pre-filter) genotypes across all 9 cells:

| Cell | Raw GT | AD (ref,alt) | DP | Passed filter? |
|---|---|---|---|---|
| CTC3 | `./.` | — | — | No — literally zero reads at this position (true dropout) |
| CTC4 | `./.` | — | — | No — same, true dropout |
| CTC7 | 0/1 | 7,6 | 13 | **No — DP 13 < 20, but 46% VAF looks like a real call** |
| (other 6 cells) | 0/1 | — | 20–99 | Yes |

Conclusion reached: **don't just remove the filter** (would let in a lot
of low-VAF noise genome-wide), but the hard `DP≥20` cutoff is somewhat
blunt — CTC7's call looks real despite failing on raw depth alone.
Ideas discussed, none started yet:
1. Lower the DP floor a bit (e.g. 10–12) while adding a VAF sanity check
   (~≥20–25%) instead of relying on raw alt-read count alone.
2. Sensitivity analysis: rerun with a looser threshold, see how much the
   shared-variant histogram shifts, manually spot-check a sample of
   newly-recovered "shared" calls for artifacts.
3. **Probably the more principled fix**: this is exactly the kind of
   partial/uneven-coverage evidence problem DelSIEVE's probabilistic
   allelic-dropout model is meant to handle properly (vs. MonoVar's
   binary per-cell hard-filter) — worth prioritizing over hand-tuning
   the MonoVar filter, once DelSIEVE is un-paused.

## Patient 04 (`Pat_04_CTCs`)

**Config**: `docs/tzu_run/{patients,params,README}.Pat_04_CTCs*`,
`Pat_04_CTCs.cells.tsv`, `Pat_04_CTCs.monovar_bams.txt`,
`run_Pat_04_CTCs_monovar.sh`. Committed and pushed, carrying forward
all the fixes/tuning from Patient 03 (weblog port, dedicated outdir,
absolute BEAST2 paths) from the start.

**WBC**: BFX-2026-04-00050. Verification was more involved than
Patient 03 since there's no dedicated reliable manifest for Patient 04:
1. Cell-type identity (WBC vs. CTC) confirmed against the lab's
   specimen/gDNA tracking table, provided directly by the user.
2. The exact BAM path was independently confirmed via a live `ls -la`
   on the server against `PAT-2026-03-00004_10/alignments/dna/`, which
   showed `BFX-2026-04-00050.redux.bam` genuinely exists there. This
   resolved an apparent ambiguity in the less-reliable
   `cnpup_unmapped_bam_rows.csv` source, where BFX-2026-04-00043 and
   BFX-2026-04-00050 both listed the same folder index — turned out to
   be a non-issue, since multiple samples can share one numbered output
   folder; each BFX ID's `.redux.bam` filename is still unique.

**Current state**: was started fresh on the server (`bash
docs/tzu_run/run_Pat_04_CTCs_monovar.sh`, no prior cache for
`PAT-2026-03-00004`) and left running when the session paused for the
weekend. **Not yet checked how far it got or whether it hit any
errors** — this is the first thing to check when picking this back up.
As of the last directory listing (mid-run), only
`monovar_wbc_subset/config/` and partial `reports/mosdepth/`,
`reports/samtools/` existed — i.e. it had barely started. It almost
certainly needs a status check and possibly a `-resume` continuation,
not a fresh restart.

## Reference: where things live

- Server repo root: `/tzu-share-2/users/students/cornelusp/CN-PUP3`
  (this environment's local clone is
  `C:\Users\CornelusPalsma\Documents\CeciNestPasUnePipeline`, kept in
  sync via `origin/main` on GitHub — `CornelImmortel/CN-PUP`).
- Outputs: `results_good_ctcs/PAT-2026-03-0000{3,4}/...` on the server.
- Shared pipeline dashboard weblog: `http://localhost:8605` (fixed,
  shared across all runs on the node).
- BEAST2 install: `/tzu-share-2/users/students/cornelusp/beast/bin/`.
