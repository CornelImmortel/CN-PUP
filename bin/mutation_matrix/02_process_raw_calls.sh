#!/usr/bin/env bash
set -euo pipefail

PROJECT_DIR="$1"
REF_FASTA="$2"
SAMPLES_TSV="$3"
PROJECT_CODE_DIR="${4:-}"

mkdir -p "${PROJECT_DIR}/processed_calls" "${PROJECT_DIR}/normalized_calls" "${PROJECT_DIR}/comparable_calls" "${PROJECT_DIR}/logs"

bulk_exclusion=""
if compgen -G "${PROJECT_DIR}/bulk_exclusion/*.external_exclusion.norm.vcf.gz" > /dev/null; then
  bulk_exclusion="$(ls "${PROJECT_DIR}"/bulk_exclusion/*.external_exclusion.norm.vcf.gz | head -n 1)"
elif compgen -G "${PROJECT_DIR}/bulk_exclusion/*.norm.vcf.gz" > /dev/null; then
  bulk_exclusion="$(ls "${PROJECT_DIR}"/bulk_exclusion/*.norm.vcf.gz | head -n 1)"
fi

keep_filters="${KEEP_FILTERS:-PASS,.}"
min_dp="${MIN_TOTAL_DEPTH:-20}"
min_alt="${MIN_ALT_READS:-4}"
min_vaf="${MIN_VAF:-0}"
run_existing_vep_filter="${RUN_EXISTING_VEP_FILTER:-false}"
max_population_af="${MAX_POPULATION_AF:-0.001}"
keep_cosmic="${KEEP_COSMIC:-true}"
genome="${GENOME:-GRCh38}"
vep_species="${VEP_SPECIES:-homo_sapiens}"
vep_cache="${VEP_CACHE:-$HOME/.vep}"
vep_cache_version="${VEP_CACHE_VERSION:-}"
vep_fasta="${VEP_FASTA:-}"
vep_forks="${VEP_FORKS:-2}"
vep_buffer_size="${VEP_BUFFER_SIZE:-5000}"
sccaller_somatic_mode="${SCCALLER_SOMATIC_MODE:-so_true}"

filter_expr_common="(FILTER=\"PASS\" || FILTER=\".\") && GT!=\"0/0\" && GT!=\"0|0\" && GT!=\"./.\" && GT!=\".|.\" && FORMAT/DP>=${min_dp} && FORMAT/AD[0:1]>=${min_alt}"

prepare_vcf_header() {
  local sample_id="$1" caller="$2" raw_vcf="$3"
  local header="${PROJECT_DIR}/logs/${sample_id}.${caller}.raw_header.txt"
  local fixed="${PROJECT_DIR}/processed_calls/${caller}/${sample_id}.${caller}.with_contigs.vcf"

  bcftools view -h "$raw_vcf" > "$header"
  if grep -q '^##contig=' "$header"; then
    printf '%s\n' "$raw_vcf"
    return
  fi

  if [[ ! -s "${REF_FASTA}.fai" ]]; then
    samtools faidx "$REF_FASTA"
  fi
  bcftools reheader -f "${REF_FASTA}.fai" -o "$fixed" "$raw_vcf"
  printf '%s\n' "$fixed"
}

process_common() {
  local sample_id="$1" caller="$2" raw_vcf="$3"
  mkdir -p "${PROJECT_DIR}/processed_calls/${caller}" "${PROJECT_DIR}/normalized_calls/${caller}" "${PROJECT_DIR}/comparable_calls/${caller}"
  local input_vcf
  input_vcf="$(prepare_vcf_header "$sample_id" "$caller" "$raw_vcf")"
  local filtered="${PROJECT_DIR}/processed_calls/${caller}/${sample_id}.${caller}.filtered.snvs.vcf.gz"
  local norm="${PROJECT_DIR}/normalized_calls/${caller}/${sample_id}.${caller}.filtered.snvs.norm.vcf.gz"
  local comp="${PROJECT_DIR}/comparable_calls/${caller}/${sample_id}.${caller}.comparable.vcf.gz"
  local norm_log="${PROJECT_DIR}/logs/${sample_id}.${caller}.norm.log"

  bcftools view -v snps -i "$filter_expr_common" "$input_vcf" -Oz -o "$filtered"
  tabix -f -p vcf "$filtered"
  bcftools norm -f "$REF_FASTA" -m -any "$filtered" -Oz -o "$norm" 2> "$norm_log"
  tabix -f -p vcf "$norm"

  if [[ -n "$bulk_exclusion" ]]; then
    bcftools isec -C -w1 "$norm" "$bulk_exclusion" -Oz -o "$comp"
  else
    cp "$norm" "$comp"
  fi
  tabix -f -p vcf "$comp"

  if [[ "$run_existing_vep_filter" == "true" && -n "$PROJECT_CODE_DIR" ]]; then
    local vep_dir="${PROJECT_DIR}/vep_existing/${caller}"
    local final="${PROJECT_DIR}/comparable_calls/${caller}/${sample_id}.${caller}.population_cosmic.comparable.vcf.gz"
    local vep_tsv="${vep_dir}/${sample_id}.${caller}.existing_vep.tsv"
    local annotations="${vep_dir}/${sample_id}.${caller}.population_cosmic.annotations.tsv"
    local summary="${vep_dir}/${sample_id}.${caller}.population_cosmic.summary.tsv"
    local filter_log="${PROJECT_DIR}/logs/${sample_id}.${caller}.existing_vep_filter.log"
    mkdir -p "$vep_dir"
    if bcftools view -h "$input_vcf" > "${PROJECT_DIR}/logs/${sample_id}.${caller}.prepared_header.txt" && grep -q '##INFO=<ID=CSQ' "${PROJECT_DIR}/logs/${sample_id}.${caller}.prepared_header.txt"; then
      python "$PROJECT_CODE_DIR/bin/extract_existing_vep_from_vcf.py" \
        --vcf "$input_vcf" \
        --output "$vep_tsv" \
        > "$filter_log" 2>&1
    else
      echo "No CSQ annotation found in $raw_vcf; running VEP on comparable VCF." > "$filter_log"
      local vep_input="${vep_dir}/${sample_id}.${caller}.vep_input.tsv"
      bcftools query \
        -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\n' \
        "$comp" \
      | awk 'BEGIN { OFS="\t" } { print $1, $2, $3, $4 "/" $5, "+", $1 ":" $2 "_" $4 "/" $5 }' \
      | sort -k1,1V -k2,2n > "$vep_input"
      if [[ ! -s "$vep_input" ]]; then
        printf 'Uploaded_variation\tLocation\tAllele\tGene\tSYMBOL\tFeature\tConsequence\tIMPACT\tHGVSc\tHGVSp\tExisting_variation\tMAX_AF\tMAX_AF_POPS\tAF\tgnomADe_AF\tgnomADg_AF\n' > "$vep_tsv"
      else
        local fasta_args=()
        local cache_version_args=()
        if [[ -n "$vep_fasta" ]]; then
          fasta_args=(--fasta "$vep_fasta")
        fi
        if [[ -n "$vep_cache_version" ]]; then
          cache_version_args=(--cache_version "$vep_cache_version")
        fi
        vep \
          --input_file "$vep_input" \
          --output_file "$vep_tsv" \
          --tab \
          --force_overwrite \
          --offline \
          --cache \
          --assembly "$genome" \
          --species "$vep_species" \
          --dir_cache "$vep_cache" \
          "${cache_version_args[@]}" \
          "${fasta_args[@]}" \
          --fork "$vep_forks" \
          --buffer_size "$vep_buffer_size" \
          --symbol \
          --canonical \
          --numbers \
          --hgvs \
          --af \
          --af_1kg \
          --af_gnomade \
          --af_gnomadg \
          --max_af \
          --no_stats \
          --fields Uploaded_variation,Location,Allele,Gene,SYMBOL,Feature,Consequence,IMPACT,HGVSc,HGVSp,Existing_variation,MAX_AF,MAX_AF_POPS,AF,gnomADe_AF,gnomADg_AF \
          >> "$filter_log" 2>&1
      fi
    fi

    if [[ -s "$vep_tsv" ]]; then
      python "$PROJECT_CODE_DIR/bin/filter_vcf_by_vep.py" \
        --vcf "$comp" \
        --vep "$vep_tsv" \
        --output-vcf "${final%.gz}" \
        --output-annotations "$annotations" \
        --summary "$summary" \
        --max-af "$max_population_af" \
        $([[ "$keep_cosmic" == "true" ]] && printf '%s' '--keep-cosmic') \
        >> "$filter_log" 2>&1
      bgzip -f -c "${final%.gz}" > "$final"
      tabix -f -p vcf "$final"
      mv -f "$final" "$comp"
      mv -f "${final}.tbi" "${comp}.tbi"
    fi
  fi
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
