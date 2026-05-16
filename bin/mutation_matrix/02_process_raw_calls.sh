#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$1"
REF_FASTA="$2"
SAMPLES_TSV="$3"

mkdir -p "${PROJECT_DIR}/processed_calls" "${PROJECT_DIR}/normalized_calls" "${PROJECT_DIR}/comparable_calls" "${PROJECT_DIR}/logs"

bulk_exclusion=""
if compgen -G "${PROJECT_DIR}/bulk_exclusion/*.norm.vcf.gz" > /dev/null; then
  bulk_exclusion="$(ls "${PROJECT_DIR}"/bulk_exclusion/*.norm.vcf.gz | head -n 1)"
fi

keep_filters="${KEEP_FILTERS:-PASS,.}"
min_dp="${MIN_TOTAL_DEPTH:-20}"
min_alt="${MIN_ALT_READS:-4}"
min_vaf="${MIN_VAF:-0}"
sccaller_somatic_mode="${SCCALLER_SOMATIC_MODE:-so_true}"

filter_expr_common="(FILTER=\"PASS\" || FILTER=\".\") && GT!=\"0/0\" && GT!=\"0|0\" && GT!=\"./.\" && GT!=\".|.\" && FORMAT/DP>=${min_dp} && FORMAT/AD[0:1]>=${min_alt}"

process_common() {
  local sample_id="$1" caller="$2" raw_vcf="$3"
  mkdir -p "${PROJECT_DIR}/processed_calls/${caller}" "${PROJECT_DIR}/normalized_calls/${caller}" "${PROJECT_DIR}/comparable_calls/${caller}"
  local filtered="${PROJECT_DIR}/processed_calls/${caller}/${sample_id}.${caller}.filtered.snvs.vcf.gz"
  local norm="${PROJECT_DIR}/normalized_calls/${caller}/${sample_id}.${caller}.filtered.snvs.norm.vcf.gz"
  local comp="${PROJECT_DIR}/comparable_calls/${caller}/${sample_id}.${caller}.comparable.vcf.gz"
  local norm_log="${PROJECT_DIR}/logs/${sample_id}.${caller}.norm.log"

  bcftools view -v snps -i "$filter_expr_common" "$raw_vcf" -Oz -o "$filtered"
  tabix -f -p vcf "$filtered"
  bcftools norm -f "$REF_FASTA" -m -any "$filtered" -Oz -o "$norm" 2> "$norm_log"
  tabix -f -p vcf "$norm"

  if [[ -n "$bulk_exclusion" ]]; then
    bcftools isec -C -w1 "$norm" "$bulk_exclusion" -Oz -o "$comp"
  else
    cp "$norm" "$comp"
  fi
  tabix -f -p vcf "$comp"
  echo "Comparable VCF: $comp"
}

process_sccaller() {
  local sample_id="$1" caller="$2" raw_vcf="$3" mode="${4:-$sccaller_somatic_mode}"
  mkdir -p "${PROJECT_DIR}/processed_calls/${caller}" "${PROJECT_DIR}/normalized_calls/${caller}" "${PROJECT_DIR}/comparable_calls/${caller}"

  local require_so="false"
  local subtract_bulk="false"
  case "$mode" in
    so_true)
      require_so="true"
      ;;
    external_bulk)
      subtract_bulk="true"
      ;;
    so_true_and_external_bulk)
      require_so="true"
      subtract_bulk="true"
      ;;
    no_bulk)
      ;;
    *)
      echo "Unsupported SCCALLER_SOMATIC_MODE: $mode" >&2
      echo "Allowed: so_true, external_bulk, so_true_and_external_bulk, no_bulk" >&2
      exit 1
      ;;
  esac

  if [[ "$subtract_bulk" == "true" && -z "$bulk_exclusion" ]]; then
    echo "SCcaller mode '$mode' requires config/bulk.tsv and a prepared bulk exclusion VCF." >&2
    exit 1
  fi

  local filtered_vcf="${PROJECT_DIR}/processed_calls/${caller}/${sample_id}.${caller}.filtered.${mode}.snvs.vcf"
  local filtered_gz="${filtered_vcf}.gz"
  local norm="${PROJECT_DIR}/normalized_calls/${caller}/${sample_id}.${caller}.filtered.${mode}.snvs.norm.vcf.gz"
  local comp="${PROJECT_DIR}/comparable_calls/${caller}/${sample_id}.${caller}.comparable.vcf.gz"
  local norm_log="${PROJECT_DIR}/logs/${sample_id}.${caller}.${mode}.norm.log"

  awk -v min_dp="$min_dp" -v min_alt="$min_alt" -v require_so="$require_so" 'BEGIN{FS=OFS="\t"}
    /^##contig=<ID=HLA-/ {next}
    /^#/ {print; next}
    {
      split($9, fmt, ":"); split($10, val, ":"); delete f; delete ad
      for (i=1; i<=length(fmt); i++) f[fmt[i]]=val[i]
      split(f["AD"], ad, ",")
      dp=0; for (i in ad) dp+=ad[i]
      passdot=($7=="." || $7=="PASS")
      nonref=(f["GT"]!="0/0" && f["GT"]!="0|0" && f["GT"]!="./." && f["GT"]!=".|.")
      snv=(length($4)==1 && length($5)==1 && $4 ~ /^[ACGT]$/ && $5 ~ /^[ACGT]$/)
      so_ok=(require_so!="true" || f["SO"]=="True")
      if (passdot && nonref && so_ok && snv && dp>=min_dp && ad[2]>=min_alt) print
    }' "$raw_vcf" > "$filtered_vcf"

  bgzip -f "$filtered_vcf"
  tabix -f -p vcf "$filtered_gz"
  bcftools norm -f "$REF_FASTA" -m -any "$filtered_gz" -Oz -o "$norm" 2> "$norm_log"
  tabix -f -p vcf "$norm"

  if [[ "$subtract_bulk" == "true" ]]; then
    bcftools isec -C -w1 "$norm" "$bulk_exclusion" -Oz -o "$comp"
  else
    cp "$norm" "$comp"
  fi
  tabix -f -p vcf "$comp"
  echo "Comparable VCF: $comp"
  echo "SCcaller mode: $mode"
}
tail -n +2 "$SAMPLES_TSV" | while IFS=$'\t' read -r sample_id caller raw_vcf sccaller_mode; do
  [[ -z "${sample_id:-}" ]] && continue
  [[ "${sample_id:0:1}" == "#" ]] && continue
  caller_lc="$(echo "$caller" | tr '[:upper:]' '[:lower:]')"
  case "$caller_lc" in
    monovar|deepvariant|haplotypecaller) process_common "$sample_id" "$caller_lc" "$raw_vcf" ;;
    sccaller) process_sccaller "$sample_id" "$caller_lc" "$raw_vcf" "${sccaller_mode:-}" ;;
    *) echo "Unsupported caller: $caller" >&2; exit 1 ;;
  esac
done
