#!/usr/bin/env bash
# Phase 2 of the explainability-layer pilot: replays CN-PUP's own RUN_MONOVAR
# command (main.nf's RUN_MONOVAR process, ~line 1916) directly against a tiny
# BED of hand-picked sites, instead of through Nextflow -- fast iteration for
# a 10-site diagnostic rerun rather than a full pipeline stage.
#
# Run from the CN-PUP repo root (paths below are relative to it, same as the
# real pipeline). Requires the monovar_py2 conda env (envs/monovar_py2.yml) --
# if you haven't got it under a name already, create it once:
#   conda env create -f envs/monovar_py2.yml -n monovar_py2_rerun
#
# Usage:
#   bin/rerun_monovar_on_sites.sh <03|04> <path/to/explainability_target.BED> [output_dir] [conda_env_name]
#
# Outputs (under output_dir, default: explainability_rerun/<patient_id>/):
#   rerun.<patient_id>.vcf            -- MonoVar's live joint-genotyping output at just these sites
#   rerun.<patient_id>.monovar.log    -- same log format RUN_MONOVAR itself produces
#   split_calls/<cell>.monovar.split.vcf  -- per-cell single-sample VCFs (via bin/split_monovar_vcf.py,
#                                             same script + same shape the real pipeline uses, so
#                                             build_monovar_variant_table.py's parse_vcf_records() can
#                                             read these completely unchanged)

set -euo pipefail

PATIENT_NUM="${1:?Usage: $0 <03|04> <bed_file> [output_dir] [conda_env_name]}"
BED_FILE="${2:?Usage: $0 <03|04> <bed_file> [output_dir] [conda_env_name]}"

case "$PATIENT_NUM" in
  03) PATIENT_ID="PAT-2026-03-00003" ;;
  04) PATIENT_ID="PAT-2026-03-00004" ;;
  *) echo "ERROR: patient must be 03 or 04, got: $PATIENT_NUM" >&2; exit 1 ;;
esac

OUTPUT_DIR="${3:-explainability_rerun/${PATIENT_ID}}"
CONDA_ENV="${4:-monovar_py2_rerun}"

REF_FASTA="/tzu-share-2/resources/genome/fasta_sarek/Homo_sapiens_assembly38.fasta"
MONOVAR_SCRIPT="/tzu-share-2/users/students/cornelusp/monovar/monovar/src/monovar.py"
BAM_LIST="docs/tzu_run/Pat_${PATIENT_NUM}_CTCs.monovar_bams.txt"
CELL_METADATA="docs/tzu_run/Pat_${PATIENT_NUM}_CTCs.cells.tsv"
THREADS=2

test -s "$REF_FASTA" || { echo "ERROR: ref_fasta not found: $REF_FASTA" >&2; exit 1; }
test -s "$MONOVAR_SCRIPT" || { echo "ERROR: monovar_script not found: $MONOVAR_SCRIPT" >&2; exit 1; }
test -s "$BAM_LIST" || { echo "ERROR: monovar_bam_list not found: $BAM_LIST (run from CN-PUP repo root)" >&2; exit 1; }
test -s "$CELL_METADATA" || { echo "ERROR: cell_metadata not found: $CELL_METADATA" >&2; exit 1; }
test -s "$BED_FILE" || { echo "ERROR: BED file not found: $BED_FILE" >&2; exit 1; }
command -v conda >/dev/null 2>&1 || { echo "ERROR: conda not found in PATH" >&2; exit 1; }

mkdir -p "$OUTPUT_DIR"
VCF_OUT="$OUTPUT_DIR/rerun.${PATIENT_ID}.vcf"
LOG_OUT="$OUTPUT_DIR/rerun.${PATIENT_ID}.monovar.log"

echo "[*] Patient $PATIENT_ID, $(grep -vc '^#' "$BED_FILE" 2>/dev/null || wc -l < "$BED_FILE") target sites, conda env: $CONDA_ENV"

source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate "$CONDA_ENV"

{
  python --version
  samtools --version | head -n 2
  echo "Starting live MonoVar rerun for $PATIENT_ID"
  echo "BAM list: $BAM_LIST"
  echo "Reference: $REF_FASTA"
  echo "MonoVar: $MONOVAR_SCRIPT"
  echo "Target BED: $BED_FILE"
  echo
} > "$LOG_OUT"

samtools mpileup \
  -BQ0 \
  -d10000 \
  -f "$REF_FASTA" \
  -q 40 \
  -l "$BED_FILE" \
  -b "$BAM_LIST" \
| python "$MONOVAR_SCRIPT" \
  -p 0.002 \
  -a 0.2 \
  -t 0.05 \
  -m "$THREADS" \
  -f "$REF_FASTA" \
  -b "$BAM_LIST" \
  -o "$VCF_OUT" \
>> "$LOG_OUT" 2>&1

echo "[+] Wrote $VCF_OUT"
conda deactivate

# Split into per-cell single-sample VCFs, same script + shape the real
# pipeline's SPLIT_MONOVAR process uses. split_monovar_vcf.py only needs
# stdlib (argparse/csv/re/pathlib) -- no need for the monovar_py2 env or any
# other conda env, plain python3 on PATH is enough.
python3 bin/split_monovar_vcf.py \
  --input "$VCF_OUT" \
  --cell-metadata "$CELL_METADATA" \
  --patient-id "$PATIENT_ID" \
  --outdir "$OUTPUT_DIR/split_calls" \
> "$OUTPUT_DIR/split_calls.log" 2>&1

echo "[+] Split per-cell VCFs written to $OUTPUT_DIR/split_calls/"
echo "[+] Done. Sanity-check row counts: $(ls "$OUTPUT_DIR"/split_calls/*.monovar.split.vcf 2>/dev/null | wc -l) per-cell files"
