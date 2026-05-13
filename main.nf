#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * CeciNestPasUnePipeline
 * Step 1: input validation plus optional MonoVar calling and splitting.
 *
 * Default run only validates inputs. Add --run_monovar true to launch MonoVar and split per-cell VCFs.
 */

params.patients = params.patients ?: "configs/patients.tsv"
params.outdir = params.outdir ?: "results"
params.monovar_script = params.monovar_script ?: ""
params.monovar_threads = params.monovar_threads ?: 2
params.monovar_region = params.monovar_region ?: ""
params.run_monovar = params.run_monovar ?: false

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

    if (params.run_monovar) {
        RUN_MONOVAR(checked_ch)
        RUN_MONOVAR.out.view { "MonoVar VCF: ${it[1]}" }
        SPLIT_MONOVAR(RUN_MONOVAR.out)
        SPLIT_MONOVAR.out.view { "Split MonoVar VCFs: ${it[1]}" }
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
    tuple val(patient_id), path("${patient_id}.monovar.vcf"), val(cell_metadata), path("${patient_id}.monovar.log")

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
    tuple val(patient_id), path(monovar_vcf), val(cell_metadata), path(monovar_log)

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
