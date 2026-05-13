#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * CeciNestPasUnePipeline
 * Input validation, optional MonoVar calling, per-cell splitting, and first
 * background-excluded comparable MonoVar VCFs.
 */

params.patients = params.patients ?: "configs/patients.tsv"
params.outdir = params.outdir ?: "results"
params.monovar_script = params.monovar_script ?: ""
params.monovar_threads = params.monovar_threads ?: 2
params.monovar_region = params.monovar_region ?: ""
params.run_monovar = params.run_monovar ?: false
params.run_vep_filter = params.run_vep_filter ?: false
params.run_bam_qc = params.run_bam_qc ?: false
params.mosdepth_threads = params.mosdepth_threads ?: 2
params.mosdepth_by = params.mosdepth_by ?: ""
params.vep_species = params.vep_species ?: "homo_sapiens"
params.vep_cache_version = params.vep_cache_version ?: ""
params.vep_fasta = params.vep_fasta ?: ""
params.vep_forks = params.vep_forks ?: 2
params.vep_buffer_size = params.vep_buffer_size ?: 5000
params.vep_sites_per_chunk = params.vep_sites_per_chunk ?: 5000

def resolveInputPath(value) {
    if (value == null) {
        return ''
    }
    def v = value.toString()
    if (v == '' || v == 'NA') {
        return v
    }
    return file(v).toString()
}

def isLeukocyteCell(cellId) {
    def x = cellId.toString().toLowerCase()
    return x == 'leukocyte' || x == 'leuko' || x == 'wbc' || x.contains('leuk')
}


def cellMetadataRows(patientId, cellMetadata) {
    def rows = []
    file(cellMetadata).readLines().drop(1).each { line ->
        if (line.trim()) {
            def cols = line.split('\t', -1)
            if (cols.size() >= 4 && cols[0] == patientId) {
                rows << tuple(patientId, cols[1], resolveInputPath(cols[2]), cols[3])
            }
        }
    }
    return rows
}

workflow {
    patients_ch = Channel
        .fromPath(params.patients, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.patient_id,
                resolveInputPath(row.ref_fasta),
                resolveInputPath(row.monovar_bam_list),
                resolveInputPath(row.cell_metadata),
                row.germline_mode,
                resolveInputPath(row.germline_vcf),
                resolveInputPath(row.bulk_bam),
                resolveInputPath(row.leukocyte_bam),
                resolveInputPath(row.leukocyte_vcf)
            )
        }

    checked_ch = CHECK_INPUTS(patients_ch)
    checked_ch.view { "Validated patient: ${it[0]} -> ${it[9]}" }

    if (params.run_bam_qc) {
        bam_qc_input_ch = checked_ch.flatMap { patient_id, ref_fasta, monovar_bam_list, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_bam, leukocyte_vcf, input_check ->
            cellMetadataRows(patient_id, cell_metadata)
        }
        SAMTOOLS_STATS(bam_qc_input_ch)
        MOSDEPTH_QC(bam_qc_input_ch)
    }

    if (params.run_monovar) {
        filter_info_ch = checked_ch.map { patient_id, ref_fasta, monovar_bam_list, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_bam, leukocyte_vcf, input_check ->
            tuple(patient_id, ref_fasta, germline_mode, germline_vcf)
        }

        precomputed_input_ch = checked_ch.filter { patient_id, ref_fasta, monovar_bam_list, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_bam, leukocyte_vcf, input_check -> germline_mode == 'precomputed_vcf' }
        PREPARE_PRECOMPUTED_GERMLINE(precomputed_input_ch)
        PREPARE_PRECOMPUTED_GERMLINE.out.view { "Precomputed germline exclusion VCF: ${it[1]}" }

        RUN_MONOVAR(checked_ch)
        RUN_MONOVAR.out.view { "MonoVar VCF: ${it[1]}" }

        SPLIT_MONOVAR(RUN_MONOVAR.out)
        SPLIT_MONOVAR.out.view { "Split MonoVar VCFs: ${it[1]}" }

        split_monovar_vcfs_ch = SPLIT_MONOVAR.out.flatMap { patient_id, split_vcfs, sample_map, split_log ->
            def files = split_vcfs instanceof List ? split_vcfs : [split_vcfs]
            files.collect { split_vcf ->
                def cell_id = split_vcf.getBaseName().replaceFirst(/\.monovar\.split$/, '')
                tuple(patient_id, cell_id, split_vcf)
            }
        }

        split_with_info_ch = split_monovar_vcfs_ch.combine(filter_info_ch, by: 0)

        leukocyte_exclusion_input_ch = split_with_info_ch.filter { patient_id, cell_id, split_vcf, ref_fasta, germline_mode, germline_vcf ->
            germline_mode == 'joint_monovar_leukocyte' && isLeukocyteCell(cell_id)
        }
        PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION(leukocyte_exclusion_input_ch)
        PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION.out.view { "MonoVar leukocyte exclusion VCF: ${it[1]}" }

        exclusion_ch = PREPARE_PRECOMPUTED_GERMLINE.out.mix(PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION.out)

        target_split_ch = split_with_info_ch.filter { patient_id, cell_id, split_vcf, ref_fasta, germline_mode, germline_vcf ->
            !isLeukocyteCell(cell_id)
        }

        monovar_background_ch = target_split_ch.combine(exclusion_ch, by: 0)
        FILTER_MONOVAR_AND_SUBTRACT_GERMLINE(monovar_background_ch)
        FILTER_MONOVAR_AND_SUBTRACT_GERMLINE.out.view { "Comparable MonoVar VCF: ${it[3]}" }

        if (params.run_vep_filter) {
            MAKE_VEP_INPUT_CHUNKS(FILTER_MONOVAR_AND_SUBTRACT_GERMLINE.out)

            vep_chunk_ch = MAKE_VEP_INPUT_CHUNKS.out.flatMap { patient_id, cell_id, comparable_vcf, chunks ->
                def files = chunks instanceof List ? chunks : [chunks]
                files.collect { chunk -> tuple(patient_id, cell_id, chunk) }
            }
            comparable_for_vep_ch = MAKE_VEP_INPUT_CHUNKS.out.map { patient_id, cell_id, comparable_vcf, chunks ->
                tuple(patient_id, cell_id, comparable_vcf)
            }

            VEP_ANNOTATE_CHUNK(vep_chunk_ch)
            vep_annotated_chunks_ch = VEP_ANNOTATE_CHUNK.out.groupTuple(by: [0, 1])
            vep_merge_input_ch = vep_annotated_chunks_ch.join(comparable_for_vep_ch, by: [0, 1])

            MERGE_VEP_AND_POPULATION_COSMIC_FILTER(vep_merge_input_ch)
            MERGE_VEP_AND_POPULATION_COSMIC_FILTER.out.view { "Population/COSMIC filtered VCF: ${it[3]}" }

            final_vcf_cell_ch = MERGE_VEP_AND_POPULATION_COSMIC_FILTER.out
                .map { patient_id, cell_id, vep_tsv, final_vcf, summary_tsv, log_file -> tuple(patient_id, cell_id, final_vcf) }
            BCFTOOLS_STATS(final_vcf_cell_ch)

            final_vcfs_by_patient_ch = final_vcf_cell_ch
                .map { patient_id, cell_id, final_vcf -> tuple(patient_id, final_vcf) }
                .groupTuple(by: 0)
            BUILD_MUTATION_MATRICES(final_vcfs_by_patient_ch)
            BUILD_MUTATION_MATRICES.out.view { "Mutation matrix long table: ${it[1]}" }

            split_vcfs_by_patient_ch = SPLIT_MONOVAR.out.map { patient_id, split_vcfs, sample_map, split_log -> tuple(patient_id, split_vcfs) }
            qc_inputs_ch = split_vcfs_by_patient_ch.combine(BUILD_MUTATION_MATRICES.out, by: 0).combine(filter_info_ch, by: 0)
            SUMMARIZE_FILTERS_QC(qc_inputs_ch)
            SUMMARIZE_FILTERS_QC.out.view { "Filter/QC summary: ${it[1]}" }

            MAKE_MULTIQC_CUSTOM_CONTENT(SUMMARIZE_FILTERS_QC.out)

            bcftools_stats_by_patient_ch = BCFTOOLS_STATS.out
                .map { patient_id, cell_id, stats_file -> tuple(patient_id, stats_file) }
                .groupTuple(by: 0)

            if (params.run_bam_qc) {
                samtools_stats_by_patient_ch = SAMTOOLS_STATS.out
                    .map { patient_id, cell_id, stats_file -> tuple(patient_id, stats_file) }
                    .groupTuple(by: 0)
                mosdepth_by_patient_ch = MOSDEPTH_QC.out
                    .map { patient_id, cell_id, mosdepth_files -> tuple(patient_id, mosdepth_files) }
                    .groupTuple(by: 0)
                multiqc_full_input_ch = MAKE_MULTIQC_CUSTOM_CONTENT.out
                    .combine(bcftools_stats_by_patient_ch, by: 0)
                    .combine(samtools_stats_by_patient_ch, by: 0)
                    .combine(mosdepth_by_patient_ch, by: 0)
                MULTIQC_REPORT_WITH_BAM_QC(multiqc_full_input_ch)
                MULTIQC_REPORT_WITH_BAM_QC.out.view { "MultiQC report: ${it[1]}" }
            } else {
                multiqc_input_ch = MAKE_MULTIQC_CUSTOM_CONTENT.out
                    .combine(bcftools_stats_by_patient_ch, by: 0)
                MULTIQC_REPORT(multiqc_input_ch)
                MULTIQC_REPORT.out.view { "MultiQC report: ${it[1]}" }
            }
        } else {
            log.info "VEP population/COSMIC filtering skipped. Add --run_vep_filter true when ready."
        }
    } else {
        log.info "MonoVar calling skipped. Re-run with --run_monovar true when ready."
    }
}

process CHECK_INPUTS {
    tag "$patient_id"
    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    tuple val(patient_id), val(ref_fasta), val(monovar_bam_list), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_bam), val(leukocyte_vcf)

    output:
    tuple val(patient_id), val(ref_fasta), val(monovar_bam_list), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_bam), val(leukocyte_vcf), path("${patient_id}.input_check.txt")

    script:
    """
    set -euo pipefail

    report="${patient_id}.input_check.txt"
    : > "\$report"

    echo "patient_id\t${patient_id}" >> "\$report"
    echo "ref_fasta\t${ref_fasta}" >> "\$report"
    echo "monovar_bam_list\t${monovar_bam_list}" >> "\$report"
    echo "cell_metadata\t${cell_metadata}" >> "\$report"
    echo "germline_mode\t${germline_mode}" >> "\$report"
    echo "germline_vcf\t${germline_vcf}" >> "\$report"
    echo "bulk_bam\t${bulk_bam}" >> "\$report"
    echo "leukocyte_bam\t${leukocyte_bam}" >> "\$report"
    echo "leukocyte_vcf\t${leukocyte_vcf}" >> "\$report"
    echo "monovar_script\t${params.monovar_script}" >> "\$report"
    echo >> "\$report"

    test -n "${patient_id}" || { echo "ERROR: missing patient_id" >&2; exit 1; }
    test -s "${ref_fasta}" || { echo "ERROR: ref_fasta not found: ${ref_fasta}" >&2; exit 1; }
    test -s "${monovar_bam_list}" || { echo "ERROR: monovar_bam_list not found: ${monovar_bam_list}" >&2; exit 1; }
    test -s "${cell_metadata}" || { echo "ERROR: cell_metadata not found: ${cell_metadata}" >&2; exit 1; }

    awk -F '\t' -v patient_id="${patient_id}" '
      NR == 1 { next }
      NF < 4 || \$1 == "" || \$2 == "" || \$3 == "" || \$4 == "" { print "ERROR: invalid cell_metadata row " NR > "/dev/stderr"; exit 1 }
      \$1 != patient_id { print "ERROR: cell_metadata patient_id mismatch on row " NR ": " \$1 > "/dev/stderr"; exit 1 }
      { count++ }
      END { if (count < 1) { print "ERROR: no cells found in cell_metadata" > "/dev/stderr"; exit 1 } }
    ' "${cell_metadata}"

    case "${germline_mode}" in
      precomputed_vcf)
        test -s "${germline_vcf}" || { echo "ERROR: germline_vcf not found: ${germline_vcf}" >&2; exit 1; }
        ;;
      deepvariant_bulk_bam)
        test -s "${bulk_bam}" || { echo "ERROR: bulk_bam not found: ${bulk_bam}" >&2; exit 1; }
        ;;
      leukocyte_vcf)
        test -s "${leukocyte_vcf}" || { echo "ERROR: leukocyte_vcf not found: ${leukocyte_vcf}" >&2; exit 1; }
        ;;
      joint_monovar_leukocyte)
        test -s "${leukocyte_bam}" || { echo "ERROR: leukocyte_bam not found: ${leukocyte_bam}" >&2; exit 1; }
        ;;
      combined)
        test -s "${germline_vcf}" || { echo "ERROR: germline_vcf not found for combined mode: ${germline_vcf}" >&2; exit 1; }
        test -s "${leukocyte_bam}" || test -s "${leukocyte_vcf}" || { echo "ERROR: combined mode requires leukocyte_bam or leukocyte_vcf" >&2; exit 1; }
        ;;
      *)
        echo "ERROR: unsupported germline_mode '${germline_mode}'" >&2
        echo "Allowed: precomputed_vcf, deepvariant_bulk_bam, leukocyte_vcf, joint_monovar_leukocyte, combined" >&2
        exit 1
        ;;
    esac

    if [[ "${params.run_monovar}" == "true" ]]; then
      test -n "${params.monovar_script}" || { echo "ERROR: --monovar_script is required when --run_monovar true" >&2; exit 1; }
      test -s "${params.monovar_script}" || { echo "ERROR: monovar_script not found: ${params.monovar_script}" >&2; exit 1; }
    fi

    echo "status\tOK" >> "\$report"
    """
}

process RUN_MONOVAR {
    tag "$patient_id"
    cpus params.monovar_threads
    conda "envs/monovar_py2.yml"
    publishDir { "${params.outdir}/${patient_id}/raw_calls/monovar" }, mode: 'copy', pattern: "*.vcf"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.log"

    input:
    tuple val(patient_id), val(ref_fasta), val(monovar_bam_list), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_bam), val(leukocyte_vcf), path(input_check)

    output:
    tuple val(patient_id), path("${patient_id}.monovar.vcf"), val(cell_metadata), val(ref_fasta), path("${patient_id}.monovar.log")

    script:
    """
    set -euo pipefail

    test -s "${params.monovar_script}" || { echo "ERROR: monovar_script not found: ${params.monovar_script}" >&2; exit 1; }
    test -s "${ref_fasta}" || { echo "ERROR: ref_fasta not found: ${ref_fasta}" >&2; exit 1; }
    test -s "${monovar_bam_list}" || { echo "ERROR: monovar_bam_list not found: ${monovar_bam_list}" >&2; exit 1; }
    command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }
    python --version > "${patient_id}.monovar.log" 2>&1 || true
    samtools --version | head -n 2 >> "${patient_id}.monovar.log" 2>&1 || true
    echo "Starting MonoVar for ${patient_id}" >> "${patient_id}.monovar.log"
    echo "BAM list: ${monovar_bam_list}" >> "${patient_id}.monovar.log"
    echo "Reference: ${ref_fasta}" >> "${patient_id}.monovar.log"
    echo "MonoVar: ${params.monovar_script}" >> "${patient_id}.monovar.log"
    echo "Region: ${params.monovar_region}" >> "${patient_id}.monovar.log"
    echo >> "${patient_id}.monovar.log"

    region_args=()
    if [[ -n "${params.monovar_region}" ]]; then
      region_args=(-r "${params.monovar_region}")
    fi

    samtools mpileup \
      -BQ0 \
      -d10000 \
      -f "${ref_fasta}" \
      -q 40 \
      -b "${monovar_bam_list}" \
      "\${region_args[@]}" \
    | python "${params.monovar_script}" \
      -p 0.002 \
      -a 0.2 \
      -t 0.05 \
      -m ${params.monovar_threads} \
      -f "${ref_fasta}" \
      -b "${monovar_bam_list}" \
      -o "${patient_id}.monovar.vcf" \
    >> "${patient_id}.monovar.log" 2>&1
    """
}

process SPLIT_MONOVAR {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/split_calls/monovar" }, mode: 'copy', pattern: "*.monovar.split.vcf"
    publishDir { "${params.outdir}/${patient_id}/split_calls/monovar" }, mode: 'copy', pattern: "*.split_sample_map.tsv"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.split.log"

    input:
    tuple val(patient_id), path(monovar_vcf), val(cell_metadata), val(ref_fasta), path(monovar_log)

    output:
    tuple val(patient_id), path("*.monovar.split.vcf"), path("${patient_id}.monovar.split_sample_map.tsv"), path("${patient_id}.monovar.split.log")

    script:
    """
    set -euo pipefail

    python "${projectDir}/bin/split_monovar_vcf.py" \
      --input "${monovar_vcf}" \
      --cell-metadata "${cell_metadata}" \
      --patient-id "${patient_id}" \
      --outdir . \
    > "${patient_id}.monovar.split.log" 2>&1
    """
}

process PREPARE_PRECOMPUTED_GERMLINE {
    tag "$patient_id"
    conda "envs/bcftools.yml"
    publishDir { "${params.outdir}/${patient_id}/germline_exclusion" }, mode: 'copy', pattern: "*.vcf.gz*"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.germline_exclusion.log"

    input:
    tuple val(patient_id), val(ref_fasta), val(monovar_bam_list), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_bam), val(leukocyte_vcf), path(input_check)

    output:
    tuple val(patient_id), path("${patient_id}.precomputed_germline.pass.snvs.norm.vcf.gz"), path("${patient_id}.precomputed_germline.pass.snvs.norm.vcf.gz.tbi"), val(ref_fasta), val("precomputed_germline")

    script:
    """
    set -euo pipefail

    bcftools view \
      -v snps \
      -f PASS \
      "${germline_vcf}" \
      -Ou \
    | bcftools norm \
      -f "${ref_fasta}" \
      -m -any \
      -Oz \
      -o "${patient_id}.precomputed_germline.pass.snvs.norm.vcf.gz" \
      2> "${patient_id}.precomputed_germline.germline_exclusion.log"

    tabix -f -p vcf "${patient_id}.precomputed_germline.pass.snvs.norm.vcf.gz"
    """
}

process PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION {
    tag { "${patient_id}:${cell_id}" }
    conda "envs/bcftools.yml"
    publishDir { "${params.outdir}/${patient_id}/germline_exclusion" }, mode: 'copy', pattern: "*.vcf.gz*"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.leukocyte_exclusion.log"

    input:
    tuple val(patient_id), val(cell_id), path(split_vcf), val(ref_fasta), val(germline_mode), val(germline_vcf)

    output:
    tuple val(patient_id), path("${patient_id}.monovar_leukocyte.filtered.norm.vcf.gz"), path("${patient_id}.monovar_leukocyte.filtered.norm.vcf.gz.tbi"), val(ref_fasta), val("monovar_leukocyte")

    script:
    """
    set -euo pipefail

    filtered="${patient_id}.monovar_leukocyte.filtered.vcf"
    filtered_gz="${patient_id}.monovar_leukocyte.filtered.vcf.gz"
    norm="${patient_id}.monovar_leukocyte.filtered.norm.vcf.gz"
    log="${patient_id}.monovar_leukocyte.leukocyte_exclusion.log"

    awk -v min_dp="${params.min_total_depth}" -v min_alt="${params.min_alt_reads}" -v min_vaf="${params.min_vaf}" '
      BEGIN { FS=OFS="\t" }
      /^#/ { print; next }
      {
        split(\$9, fmt, ":")
        split(\$10, val, ":")
        delete f
        for (i = 1; i <= length(fmt); i++) f[fmt[i]] = val[i]
        split(f["AD"], ad, ",")
        gt = f["GT"]
        dp = f["DP"] + 0
        ref = ad[1] + 0
        alt = ad[2] + 0
        vaf = (ref + alt) > 0 ? alt / (ref + alt) : 0
        if ((gt != "0/0") && (gt != "0|0") && (gt != "./.") && (gt != ".|.") && dp >= min_dp && alt >= min_alt && vaf >= min_vaf) print
      }
    ' "${split_vcf}" > "\$filtered"

    echo "Input leukocyte split VCF: ${split_vcf}" > "\$log"
    echo "min_total_depth=${params.min_total_depth}" >> "\$log"
    echo "min_alt_reads=${params.min_alt_reads}" >> "\$log"
    echo "min_vaf=${params.min_vaf}" >> "\$log"
    echo -n "Leukocyte exclusion variants before normalization: " >> "\$log"
    grep -vc '^#' "\$filtered" >> "\$log"

    bgzip -f -c "\$filtered" > "\$filtered_gz"
    tabix -f -p vcf "\$filtered_gz"

    bcftools view -v snps "\$filtered_gz" -Ou \
    | bcftools norm \
      -f "${ref_fasta}" \
      -m -any \
      -Oz \
      -o "\$norm" \
      2>> "\$log"
    tabix -f -p vcf "\$norm"

    echo -n "Leukocyte exclusion variants after normalization: " >> "\$log"
    bcftools view -H "\$norm" | wc -l >> "\$log"
    """
}

process FILTER_MONOVAR_AND_SUBTRACT_GERMLINE {
    tag { "${patient_id}:${cell_id}:${exclusion_source}" }
    conda "envs/bcftools.yml"
    publishDir { "${params.outdir}/${patient_id}/processed_calls/monovar" }, mode: 'copy', pattern: "*.filtered.vcf"
    publishDir { "${params.outdir}/${patient_id}/normalized_calls/monovar" }, mode: 'copy', pattern: "*.norm.vcf.gz*"
    publishDir { "${params.outdir}/${patient_id}/comparable_calls/monovar" }, mode: 'copy', pattern: "*.no_*.vcf.gz*"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.monovar_filter.log"

    input:
    tuple val(patient_id), val(cell_id), path(split_vcf), val(ref_fasta), val(germline_mode), val(germline_vcf), path(exclusion_vcf), path(exclusion_tbi), val(exclusion_ref_fasta), val(exclusion_source)

    output:
    tuple val(patient_id), val(cell_id), path("*.filtered.norm.vcf.gz"), path("*.no_*.vcf.gz"), path("*.monovar_filter.log")

    script:
    """
    set -euo pipefail

    filtered="${cell_id}.monovar.filtered.vcf"
    filtered_gz="${cell_id}.monovar.filtered.vcf.gz"
    norm="${cell_id}.monovar.filtered.norm.vcf.gz"
    comparable="${cell_id}.monovar.no_${exclusion_source}.vcf.gz"
    log="${cell_id}.monovar_filter.log"

    awk -v min_dp="${params.min_total_depth}" -v min_alt="${params.min_alt_reads}" -v min_vaf="${params.min_vaf}" '
      BEGIN { FS=OFS="\t" }
      /^#/ { print; next }
      {
        split(\$9, fmt, ":")
        split(\$10, val, ":")
        delete f
        for (i = 1; i <= length(fmt); i++) f[fmt[i]] = val[i]
        split(f["AD"], ad, ",")
        gt = f["GT"]
        dp = f["DP"] + 0
        ref = ad[1] + 0
        alt = ad[2] + 0
        vaf = (ref + alt) > 0 ? alt / (ref + alt) : 0
        if ((gt != "0/0") && (gt != "0|0") && (gt != "./.") && (gt != ".|.") && dp >= min_dp && alt >= min_alt && vaf >= min_vaf) print
      }
    ' "${split_vcf}" > "\$filtered"

    echo "Input split VCF: ${split_vcf}" > "\$log"
    echo "Exclusion source: ${exclusion_source}" >> "\$log"
    echo "Exclusion VCF: ${exclusion_vcf}" >> "\$log"
    echo "min_total_depth=${params.min_total_depth}" >> "\$log"
    echo "min_alt_reads=${params.min_alt_reads}" >> "\$log"
    echo "min_vaf=${params.min_vaf}" >> "\$log"
    echo -n "Filtered variants: " >> "\$log"
    grep -vc '^#' "\$filtered" >> "\$log"

    bgzip -f -c "\$filtered" > "\$filtered_gz"
    tabix -f -p vcf "\$filtered_gz"

    bcftools view -v snps "\$filtered_gz" -Ou \
    | bcftools norm \
      -f "${ref_fasta}" \
      -m -any \
      -Oz \
      -o "\$norm" \
      2>> "\$log"
    tabix -f -p vcf "\$norm"

    bcftools isec \
      -C \
      -w1 \
      "\$norm" \
      "${exclusion_vcf}" \
      -Oz \
      -o "\$comparable" \
      2>> "\$log"
    tabix -f -p vcf "\$comparable"

    echo -n "Comparable variants after background subtraction: " >> "\$log"
    bcftools view -H "\$comparable" | wc -l >> "\$log"
    """
}

process MAKE_VEP_INPUT_CHUNKS {
    tag { "${patient_id}:${cell_id}" }
    conda "envs/bcftools.yml"
    publishDir { "${params.outdir}/${patient_id}/vep/monovar/chunks" }, mode: 'copy', pattern: "vep_chunks/*.tsv"

    input:
    tuple val(patient_id), val(cell_id), path(norm_vcf), path(comparable_vcf), path(filter_log)

    output:
    tuple val(patient_id), val(cell_id), path("${cell_id}.comparable.input.vcf.gz"), path("vep_chunks/*.tsv")

    script:
    """
    set -euo pipefail

    mkdir -p vep_chunks
    cp "${comparable_vcf}" "${cell_id}.comparable.input.vcf.gz"

    bcftools query \
      -f '%CHROM\t%POS\t%POS\t%REF\t%ALT\n' \
      "${comparable_vcf}" \
    | awk 'BEGIN { OFS="\t" } { print \$1, \$2, \$3, \$4 "/" \$5, "+", \$1 ":" \$2 "_" \$4 "/" \$5 }' \
    | sort -k1,1V -k2,2n > "${cell_id}.vep_input.all.tsv"

    if [[ ! -s "${cell_id}.vep_input.all.tsv" ]]; then
      : > "vep_chunks/${cell_id}.chunk_000000.tsv"
    else
      split \
        -l "${params.vep_sites_per_chunk}" \
        -d \
        -a 6 \
        --additional-suffix=.tsv \
        "${cell_id}.vep_input.all.tsv" \
        "vep_chunks/${cell_id}.chunk_"
    fi
    """
}

process VEP_ANNOTATE_CHUNK {
    tag { "${patient_id}:${cell_id}:${chunk.simpleName}" }
    cpus params.vep_forks
    conda "envs/vep.yml"
    publishDir { "${params.outdir}/${patient_id}/vep/monovar/chunk_annotations" }, mode: 'copy', pattern: "*.vep.tsv"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.vep_chunk.log"

    input:
    tuple val(patient_id), val(cell_id), path(chunk)

    output:
    tuple val(patient_id), val(cell_id), path("*.vep.tsv")

    script:
    """
    set -euo pipefail

    out="${chunk.simpleName}.vep.tsv"
    log="${chunk.simpleName}.vep_chunk.log"

    if [[ ! -s "${chunk}" ]]; then
      printf '#Uploaded_variation\tLocation\tAllele\tGene\tSYMBOL\tFeature\tConsequence\tIMPACT\tHGVSc\tHGVSp\tExisting_variation\tMAX_AF\tMAX_AF_POPS\tAF\tgnomADe_AF\tgnomADg_AF\n' > "\$out"
      echo "Empty VEP chunk: ${chunk}" > "\$log"
      exit 0
    fi

    fasta_args=()
    if [[ -n "${params.vep_fasta}" ]]; then
      fasta_args=(--fasta "${params.vep_fasta}")
    fi

    cache_version_args=()
    if [[ -n "${params.vep_cache_version}" ]]; then
      cache_version_args=(--cache_version "${params.vep_cache_version}")
    fi

    vep \
      --input_file "${chunk}" \
      --output_file "\$out" \
      --tab \
      --force_overwrite \
      --offline \
      --cache \
      --assembly "${params.genome}" \
      --species "${params.vep_species}" \
      --dir_cache "${params.vep_cache}" \
      "\${cache_version_args[@]}" \
      "\${fasta_args[@]}" \
      --fork "${params.vep_forks}" \
      --buffer_size "${params.vep_buffer_size}" \
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
      > "\$log" 2>&1
    """
}

process MERGE_VEP_AND_POPULATION_COSMIC_FILTER {
    tag { "${patient_id}:${cell_id}" }
    conda "envs/vep.yml"
    publishDir { "${params.outdir}/${patient_id}/vep/monovar" }, mode: 'copy', pattern: "*.vep.tsv"
    publishDir { "${params.outdir}/${patient_id}/vep/monovar" }, mode: 'copy', pattern: "*.population_cosmic.*.tsv"
    publishDir { "${params.outdir}/${patient_id}/final_calls/monovar" }, mode: 'copy', pattern: "*.population_cosmic.vcf.gz*"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.vep_filter.log"

    input:
    tuple val(patient_id), val(cell_id), path(vep_chunks), path(comparable_vcf)

    output:
    tuple val(patient_id), val(cell_id), path("*.vep.tsv"), path("*.population_cosmic.vcf.gz"), path("*.population_cosmic.summary.tsv"), path("*.vep_filter.log")

    script:
    """
    set -euo pipefail

    prefix="${cell_id}.monovar"
    merged_vep="\${prefix}.vep.tsv"
    kept_vcf="\${prefix}.population_cosmic.vcf"
    kept_gz="\${prefix}.population_cosmic.vcf.gz"
    summary="\${prefix}.population_cosmic.summary.tsv"
    annotations="\${prefix}.population_cosmic.annotations.tsv"
    log="\${prefix}.vep_filter.log"

    echo "Input comparable VCF: ${comparable_vcf}" > "\$log"
    echo "VEP chunks: ${vep_chunks}" >> "\$log"
    echo "max_population_af=${params.max_population_af}" >> "\$log"
    echo "keep_cosmic=${params.keep_cosmic}" >> "\$log"
    echo "vep_sites_per_chunk=${params.vep_sites_per_chunk}" >> "\$log"

    python "${projectDir}/bin/merge_vep_tables.py" \
      --output "\$merged_vep" \
      ${vep_chunks}

    python "${projectDir}/bin/filter_vcf_by_vep.py" \
      --vcf "${comparable_vcf}" \
      --vep "\$merged_vep" \
      --output-vcf "\$kept_vcf" \
      --output-annotations "\$annotations" \
      --summary "\$summary" \
      --max-af "${params.max_population_af}" \
      ${params.keep_cosmic ? '--keep-cosmic' : ''} \
      >> "\$log" 2>&1

    bgzip -f -c "\$kept_vcf" > "\$kept_gz"
    tabix -f -p vcf "\$kept_gz"
    """
}

process BUILD_MUTATION_MATRICES {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/tables" }, mode: 'copy', pattern: "*.tsv"
    publishDir { "${params.outdir}/${patient_id}/matrices" }, mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(patient_id), path(final_vcfs)

    output:
    tuple val(patient_id), path("*.long.tsv"), path("*.binary_matrix.tsv"), path("*.altread_matrix.tsv"), path("*.refread_matrix.tsv"), path("*.summary.tsv")

    script:
    """
    set -euo pipefail

    python "${projectDir}/bin/build_mutation_matrices.py" \
      --patient-id "${patient_id}" \
      --caller monovar \
      --out-prefix "${patient_id}.monovar.final" \
      ${final_vcfs}
    """
}

process SUMMARIZE_FILTERS_QC {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/reports" }, mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(patient_id), path(split_vcfs), path(long_table), path(binary_matrix), path(altread_matrix), path(refread_matrix), path(final_summary), val(ref_fasta), val(germline_mode), val(germline_vcf)

    output:
    tuple val(patient_id), path("*.filter_settings.tsv"), path("*.prefilter_depth_alt_qc.tsv"), path("*.filter_impact.tsv")

    script:
    """
    set -euo pipefail

    python "${projectDir}/bin/summarize_filters_qc.py" \
      --patient-id "${patient_id}" \
      --caller monovar \
      --out-prefix "${patient_id}.monovar" \
      --final-summary "${final_summary}" \
      --min-total-depth "${params.min_total_depth}" \
      --min-alt-reads "${params.min_alt_reads}" \
      --min-vaf "${params.min_vaf}" \
      --max-population-af "${params.max_population_af}" \
      --keep-cosmic "${params.keep_cosmic}" \
      --germline-mode "${germline_mode}" \
      --vep-sites-per-chunk "${params.vep_sites_per_chunk}" \
      ${split_vcfs}
    """
}


process SAMTOOLS_STATS {
    tag { "${patient_id}:${cell_id}" }
    conda "envs/alignment_qc.yml"
    publishDir { "${params.outdir}/${patient_id}/reports/samtools/${cell_id}" }, mode: 'copy', pattern: "*.samtools.stats.out"

    input:
    tuple val(patient_id), val(cell_id), val(bam), val(cell_type)

    output:
    tuple val(patient_id), val(cell_id), path("*.samtools.stats.out")

    script:
    """
    set -euo pipefail

    samtools stats "${bam}" > "${cell_id}.samtools.stats.out"
    """
}

process MOSDEPTH_QC {
    tag { "${patient_id}:${cell_id}" }
    cpus params.mosdepth_threads
    conda "envs/alignment_qc.yml"
    publishDir { "${params.outdir}/${patient_id}/reports/mosdepth/${cell_id}" }, mode: 'copy'

    input:
    tuple val(patient_id), val(cell_id), val(bam), val(cell_type)

    output:
    tuple val(patient_id), val(cell_id), path("${cell_id}*")

    script:
    """
    set -euo pipefail

    by_args=()
    if [[ -n "${params.mosdepth_by}" ]]; then
      by_args=(--by "${params.mosdepth_by}")
    fi

    mosdepth \
      --threads ${task.cpus} \
      "\${by_args[@]}" \
      "${cell_id}" \
      "${bam}"
    """
}

process BCFTOOLS_STATS {
    tag { "${patient_id}:${cell_id}" }
    conda "envs/bcftools.yml"
    publishDir { "${params.outdir}/${patient_id}/reports/bcftools" }, mode: 'copy', pattern: "*.bcftools_stats.txt"

    input:
    tuple val(patient_id), val(cell_id), path(vcf)

    output:
    tuple val(patient_id), val(cell_id), path("*.bcftools_stats.txt")

    script:
    """
    set -euo pipefail

    bcftools stats "${vcf}" > "${cell_id}.monovar.bcftools_stats.txt"
    """
}

process MAKE_MULTIQC_CUSTOM_CONTENT {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/reports/multiqc_custom" }, mode: 'copy', pattern: "*_mqc.json"

    input:
    tuple val(patient_id), path(filter_settings), path(prefilter_qc), path(filter_impact)

    output:
    tuple val(patient_id), path("*_mqc.json")

    script:
    """
    set -euo pipefail

    python "${projectDir}/bin/make_multiqc_custom_content.py" \
      --patient-id "${patient_id}" \
      --settings "${filter_settings}" \
      --prefilter-qc "${prefilter_qc}" \
      --filter-impact "${filter_impact}" \
      --outdir .
    """
}

process MULTIQC_REPORT {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/multiqc" }, mode: 'copy'

    input:
    tuple val(patient_id), path(multiqc_custom_files), path(bcftools_stats_files)

    output:
    tuple val(patient_id), path("${patient_id}.multiqc_report.html"), path("${patient_id}.multiqc_report_data")

    script:
    """
    set -euo pipefail

    mkdir -p multiqc_inputs
    cp ${multiqc_custom_files} multiqc_inputs/
    cp ${bcftools_stats_files} multiqc_inputs/
    multiqc multiqc_inputs --force --outdir . --filename "${patient_id}.multiqc_report.html"
    """
}

process MULTIQC_REPORT_WITH_BAM_QC {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/multiqc" }, mode: 'copy'

    input:
    tuple val(patient_id), path(multiqc_custom_files), path(bcftools_stats_files), path(samtools_stats_files), path(mosdepth_files)

    output:
    tuple val(patient_id), path("${patient_id}.multiqc_report.html"), path("${patient_id}.multiqc_report_data")

    script:
    """
    set -euo pipefail

    mkdir -p multiqc_inputs
    cp ${multiqc_custom_files} multiqc_inputs/
    cp ${bcftools_stats_files} multiqc_inputs/
    cp ${samtools_stats_files} multiqc_inputs/
    cp ${mosdepth_files} multiqc_inputs/
    multiqc multiqc_inputs --force --outdir . --filename "${patient_id}.multiqc_report.html"
    """
}
