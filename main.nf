#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

/*
 * CeciNestPasUnePipeline
 * Step 0: minimal skeleton and input validation.
 *
 * This first workflow intentionally does not run MonoVar yet.
 * It reads a patient sheet, checks that required files exist, and creates
 * a small per-patient validation report. Once this works, RUN_MONOVAR will
 * be added as the first real process.
 */

params.patients = params.patients ?: "configs/patients.tsv"
params.outdir = params.outdir ?: "results"
params.monovar_script = params.monovar_script ?: ""

workflow {
    patients_ch = Channel
        .fromPath(params.patients, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            tuple(
                row.patient_id,
                row.ref_fasta,
                row.monovar_bam_list,
                row.germline_mode,
                row.germline_vcf,
                row.bulk_bam,
                row.leukocyte_bam,
                row.leukocyte_vcf
            )
        }

    CHECK_INPUTS(patients_ch)

    CHECK_INPUTS.out.view { "Validated patient: ${it[0]} -> ${it[1]}" }
}

process CHECK_INPUTS {
    tag "$patient_id"
    publishDir "${params.outdir}/validation", mode: 'copy'

    input:
    tuple val(patient_id), val(ref_fasta), val(monovar_bam_list), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_bam), val(leukocyte_vcf)

    output:
    tuple val(patient_id), path("${patient_id}.input_check.txt")

    script:
    """
    set -euo pipefail

    report="${patient_id}.input_check.txt"
    : > "\$report"

    echo "patient_id\t${patient_id}" >> "\$report"
    echo "ref_fasta\t${ref_fasta}" >> "\$report"
    echo "monovar_bam_list\t${monovar_bam_list}" >> "\$report"
    echo "germline_mode\t${germline_mode}" >> "\$report"
    echo "germline_vcf\t${germline_vcf}" >> "\$report"
    echo "bulk_bam\t${bulk_bam}" >> "\$report"
    echo "leukocyte_bam\t${leukocyte_bam}" >> "\$report"
    echo "leukocyte_vcf\t${leukocyte_vcf}" >> "\$report"
    echo >> "\$report"

    test -n "${patient_id}" || { echo "ERROR: missing patient_id" >&2; exit 1; }
    test -s "${ref_fasta}" || { echo "ERROR: ref_fasta not found: ${ref_fasta}" >&2; exit 1; }
    test -s "${monovar_bam_list}" || { echo "ERROR: monovar_bam_list not found: ${monovar_bam_list}" >&2; exit 1; }

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

    echo "status\tOK" >> "\$report"
    """
}
