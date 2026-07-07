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
params.monovar_targets_bed = params.monovar_targets_bed ?: ""
params.monovar_background_ctc_count = params.monovar_background_ctc_count ?: 5
params.monovar_background_seed = params.monovar_background_seed ?: 13
params.run_monovar = params.run_monovar ?: false
params.run_sccaller = params.run_sccaller ?: false
params.run_external_callers = params.run_external_callers ?: false
params.run_comparison = params.run_comparison ?: false
params.run_caller_scorecard = params.run_caller_scorecard ?: false
params.run_delsieve_prep = params.run_delsieve_prep ?: false
params.run_delsieve = params.run_delsieve ?: false
params.run_delsieve_from_existing = params.run_delsieve_from_existing ?: false
params.run_delsieve_tree_annotator = params.run_delsieve_tree_annotator == null ? true : params.run_delsieve_tree_annotator
params.caller_vcfs = params.caller_vcfs ?: "configs/caller_vcfs.tsv"
params.sccaller_script = params.sccaller_script ?: ""
params.sccaller_cpu = params.sccaller_cpu ?: 8
params.sccaller_bias = params.sccaller_bias ?: 0.75
params.sccaller_engine = params.sccaller_engine ?: "samtools"
params.sccaller_mode = params.sccaller_mode ?: "external_bulk"
params.sccaller_hsnp_vcf = params.sccaller_hsnp_vcf ?: ""
params.sccaller_head = params.sccaller_head ?: ""
params.sccaller_tail = params.sccaller_tail ?: ""
params.run_vep_filter = params.run_vep_filter ?: false
params.run_snv_report = params.run_snv_report == null ? true : params.run_snv_report
params.run_final_report = params.run_final_report == null ? true : params.run_final_report
params.final_report_top_variants = params.final_report_top_variants ?: 200
params.final_report_rarefaction_iterations = params.final_report_rarefaction_iterations ?: 200
params.run_bam_qc = params.run_bam_qc ?: false
params.mosdepth_threads = params.mosdepth_threads ?: 2
params.mosdepth_by = params.mosdepth_by ?: ""
params.vep_species = params.vep_species ?: "homo_sapiens"
params.vep_cache_version = params.vep_cache_version ?: ""
params.vep_fasta = params.vep_fasta ?: ""
params.vep_forks = params.vep_forks ?: 2
params.vep_buffer_size = params.vep_buffer_size ?: 5000
params.vep_sites_per_chunk = params.vep_sites_per_chunk ?: 5000
params.delsieve_min_ctc_support = params.delsieve_min_ctc_support ?: 2
params.delsieve_callers = params.delsieve_callers ?: "all"
params.delsieve_min_base_quality = params.delsieve_min_base_quality ?: 0
params.delsieve_min_mapping_quality = params.delsieve_min_mapping_quality ?: 0
params.delsieve_existing_long_table = params.delsieve_existing_long_table ?: ""
params.delsieve_applauncher = params.delsieve_applauncher ?: "applauncher"
params.delsieve_beast = params.delsieve_beast ?: "beast"
params.delsieve_treeannotator = params.delsieve_treeannotator ?: "applauncher"
params.delsieve_version_file = params.delsieve_version_file ?: ""
params.delsieve_template_stage1 = params.delsieve_template_stage1 ?: ""
params.delsieve_threads = params.delsieve_threads ?: 4
params.delsieve_datatype = params.delsieve_datatype ?: 1
params.delsieve_missing_threshold = params.delsieve_missing_threshold ?: 1
params.delsieve_sample_cells = params.delsieve_sample_cells ?: 0
params.delsieve_sample_loci = params.delsieve_sample_loci ?: 0
params.delsieve_ignore_sex = params.delsieve_ignore_sex ?: false
params.delsieve_beast_extra_args = params.delsieve_beast_extra_args ?: ""
params.delsieve_datacollector_extra_args = params.delsieve_datacollector_extra_args ?: ""
params.delsieve_tree_burnin = params.delsieve_tree_burnin ?: 10
params.delsieve_tree_heights = params.delsieve_tree_heights ?: "median"
params.run_delsieve_variant_caller = params.run_delsieve_variant_caller == null ? true : params.run_delsieve_variant_caller
params.delsieve_variant_burnin = params.delsieve_variant_burnin ?: params.delsieve_tree_burnin
params.delsieve_variant_estimates = params.delsieve_variant_estimates ?: "median"
params.delsieve_variant_extra_args = params.delsieve_variant_extra_args ?: ""
params.run_delsieve_tree_png = params.run_delsieve_tree_png == null ? true : params.run_delsieve_tree_png
params.delsieve_figtree = params.delsieve_figtree ?: "figtree"
params.delsieve_figtree_extra_args = params.delsieve_figtree_extra_args ?: ""
params.run_delsieve_gene_annotator = params.run_delsieve_gene_annotator == null ? true : params.run_delsieve_gene_annotator
params.delsieve_mutation_map = params.delsieve_mutation_map ?: ""
params.delsieve_geneannotator_subst = params.delsieve_geneannotator_subst ?: 1
params.delsieve_mutation_map_annovar = params.delsieve_mutation_map_annovar ?: false
params.delsieve_annovar_separator = params.delsieve_annovar_separator ?: 0
params.delsieve_annovar_filter_col = params.delsieve_annovar_filter_col == null ? -1 : params.delsieve_annovar_filter_col
params.delsieve_annovar_filter_keywords = params.delsieve_annovar_filter_keywords ?: ""
params.delsieve_gene_filter = params.delsieve_gene_filter ?: ""
params.delsieve_gene_filter_sep = params.delsieve_gene_filter_sep ?: 0
params.delsieve_gene_filter_col = params.delsieve_gene_filter_col ?: 0
params.delsieve_geneannotator_extra_args = params.delsieve_geneannotator_extra_args ?: ""

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

def isMonovarSubsetBackgroundMode(germlineMode) {
    return germlineMode == 'monovar_wbc_subset'
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

    checked_info_ch = checked_ch.map { patient_id, ref_fasta, monovar_bam_list, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_bam, leukocyte_vcf, input_check ->
        tuple(patient_id, ref_fasta, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_vcf)
    }

    comparison_inputs = Channel.empty()

    delsieve_source_ch = Channel.empty()

    if (params.run_delsieve_from_existing) {
        def existing_long_table = resolveInputPath(params.delsieve_existing_long_table)
        if (!existing_long_table || existing_long_table == 'NA') {
            error "delsieve_existing_long_table is required when run_delsieve_from_existing is true"
        }
        delsieve_existing_ch = checked_info_ch.map { patient_id, ref_fasta, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_vcf ->
            tuple(patient_id, file(existing_long_table), ref_fasta, cell_metadata)
        }
        delsieve_source_ch = delsieve_source_ch.mix(delsieve_existing_ch)
    }

    if (params.run_external_callers) {
        external_rows_ch = Channel
            .fromPath(params.caller_vcfs, checkIfExists: true)
            .splitCsv(header: true, sep: '\t')
            .filter { row -> row.patient_id && !row.patient_id.toString().startsWith('#') }
            .map { row ->
                tuple(
                    row.patient_id.toString(),
                    [
                        row.cell_id.toString(),
                        row.caller.toString(),
                        resolveInputPath(row.vcf_path),
                        row.sccaller_mode == null ? '' : row.sccaller_mode.toString()
                    ]
                )
            }
            .groupTuple(by: 0)

        STANDARDIZE_EXTERNAL_CALLERS(external_rows_ch.join(checked_info_ch, by: 0))
        comparison_inputs = comparison_inputs.mix(STANDARDIZE_EXTERNAL_CALLERS.out)
        STANDARDIZE_EXTERNAL_CALLERS.out.view { "External comparable VCFs for ${it[0]}" }
    }

    if (params.run_bam_qc) {
        bam_qc_input_ch = checked_ch.flatMap { patient_id, ref_fasta, monovar_bam_list, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_bam, leukocyte_vcf, input_check ->
            cellMetadataRows(patient_id, cell_metadata)
        }
        SAMTOOLS_STATS(bam_qc_input_ch)
        MOSDEPTH_QC(bam_qc_input_ch)

        mosdepth_by_patient_ch = MOSDEPTH_QC.out
            .map { patient_id, cell_id, mosdepth_files -> tuple(patient_id, mosdepth_files) }
            .groupTuple(by: 0)
            .map { patient_id, mosdepth_file_groups -> tuple(patient_id, mosdepth_file_groups.flatten()) }
        SUMMARIZE_BREADTH_QC(mosdepth_by_patient_ch)
        SUMMARIZE_BREADTH_QC.out.view { "Breadth-of-coverage summary: ${it[1]}" }
    }

    if (params.run_monovar) {
        filter_info_ch = checked_ch.map { patient_id, ref_fasta, monovar_bam_list, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_bam, leukocyte_vcf, input_check ->
            tuple(patient_id, ref_fasta, germline_mode, germline_vcf)
        }

        precomputed_input_ch = checked_ch.filter { patient_id, ref_fasta, monovar_bam_list, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_bam, leukocyte_vcf, input_check -> germline_mode == 'precomputed_vcf' || germline_mode == 'combined' }
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

        if (params.run_monovar) {
            monovar_subset_bg_input_ch = checked_ch.filter { patient_id, ref_fasta, monovar_bam_list, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_bam, leukocyte_vcf, input_check ->
                isMonovarSubsetBackgroundMode(germline_mode)
            }
            PREPARE_MONOVAR_WBC_SUBSET(monovar_subset_bg_input_ch)
            PREPARE_MONOVAR_WBC_SUBSET.out.view { "MonoVar WBC subset input for ${it[0]}: ${it[2]}" }

            RUN_MONOVAR_WBC_SUBSET(PREPARE_MONOVAR_WBC_SUBSET.out)
            RUN_MONOVAR_WBC_SUBSET.out.view { "MonoVar WBC subset VCF: ${it[1]}" }

            SPLIT_MONOVAR_WBC_SUBSET(RUN_MONOVAR_WBC_SUBSET.out)
            SPLIT_MONOVAR_WBC_SUBSET.out.view { "Split MonoVar WBC subset VCFs: ${it[1]}" }
        }

        leukocyte_exclusion_input_ch = split_with_info_ch.filter { patient_id, cell_id, split_vcf, ref_fasta, germline_mode, germline_vcf ->
            (germline_mode == 'joint_monovar_leukocyte' || germline_mode == 'combined') && isLeukocyteCell(cell_id)
        }
        subset_leukocyte_exclusion_input_ch = SPLIT_MONOVAR_WBC_SUBSET.out
            .flatMap { patient_id, split_vcfs, sample_map, split_log, ref_fasta ->
                def files = split_vcfs instanceof List ? split_vcfs : [split_vcfs]
                files.collect { split_vcf ->
                    def cell_id = split_vcf.getBaseName().replaceFirst(/\.monovar\.split$/, '')
                    tuple(patient_id, cell_id, split_vcf, ref_fasta, 'monovar_wbc_subset', 'NA')
                }
            }
            .filter { patient_id, cell_id, split_vcf, ref_fasta, germline_mode, germline_vcf ->
                isLeukocyteCell(cell_id)
            }
        leukocyte_exclusion_all_ch = leukocyte_exclusion_input_ch.mix(subset_leukocyte_exclusion_input_ch)
        PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION(leukocyte_exclusion_all_ch)
        PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION.out.view { "MonoVar leukocyte exclusion VCF: ${it[1]}" }

        // 'combined' mode makes both PREPARE_PRECOMPUTED_GERMLINE and
        // PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION fire for the same patient; a plain
        // mix()+combine(by:0) would pair each cell with EACH exclusion source
        // separately (two independent single-source subtractions), not a real
        // bulk+leukocyte combined exclusion. Join and merge those cases instead;
        // 'precomputed_vcf'-only and 'joint_monovar_leukocyte'-only patients pass
        // through unchanged (only one side of the join is populated for them).
        // join(remainder:true) collapses an entirely-empty upstream channel (e.g.
        // PREPARE_PRECOMPUTED_GERMLINE never firing for a germline_mode that isn't
        // 'precomputed_vcf'/'combined') into a single null instead of one null per
        // missing field, which breaks positional destructuring. Tag + mix +
        // groupTuple instead so the branch closure sees a predictable shape
        // regardless of which side(s) actually emitted for a given patient.
        exclusion_join_ch = PREPARE_PRECOMPUTED_GERMLINE.out
            .map { patient_id, vcf, tbi, ref_fasta, source -> tuple(patient_id, 'bulk', vcf, tbi, ref_fasta, source) }
            .mix(PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION.out
                .map { patient_id, vcf, tbi, ref_fasta, source -> tuple(patient_id, 'leuko', vcf, tbi, ref_fasta, source) })
            .groupTuple(by: 0)
            .branch { patient_id, kinds, vcfs, tbis, refs, sources ->
                def bulk_idx = kinds.indexOf('bulk')
                def leuko_idx = kinds.indexOf('leuko')
                both: bulk_idx >= 0 && leuko_idx >= 0
                    return tuple(patient_id, vcfs[bulk_idx], tbis[bulk_idx], vcfs[leuko_idx], tbis[leuko_idx], refs[bulk_idx])
                bulk_only: bulk_idx >= 0
                    return tuple(patient_id, vcfs[bulk_idx], tbis[bulk_idx], refs[bulk_idx], sources[bulk_idx])
                leuko_only: true
                    return tuple(patient_id, vcfs[leuko_idx], tbis[leuko_idx], refs[leuko_idx], sources[leuko_idx])
            }

        MERGE_MONOVAR_GERMLINE_EXCLUSION(exclusion_join_ch.both)
        MERGE_MONOVAR_GERMLINE_EXCLUSION.out.view { "Combined bulk+leukocyte germline exclusion VCF: ${it[1]}" }

        exclusion_ch = exclusion_join_ch.bulk_only
            .mix(exclusion_join_ch.leuko_only)
            .mix(MERGE_MONOVAR_GERMLINE_EXCLUSION.out)

        target_split_ch = split_with_info_ch.filter { patient_id, cell_id, split_vcf, ref_fasta, germline_mode, germline_vcf ->
            !isLeukocyteCell(cell_id)
        }

        monovar_background_ch = target_split_ch.combine(exclusion_ch, by: 0)
        FILTER_MONOVAR_AND_SUBTRACT_GERMLINE(monovar_background_ch)
        FILTER_MONOVAR_AND_SUBTRACT_GERMLINE.out.view { "Comparable MonoVar VCF: ${it[3]}" }
        monovar_split_by_patient_ch = SPLIT_MONOVAR.out
            .map { patient_id, split_vcfs, sample_map, split_log -> tuple(patient_id, split_vcfs) }
        monovar_norm_by_patient_ch = FILTER_MONOVAR_AND_SUBTRACT_GERMLINE.out
            .map { patient_id, cell_id, norm_vcf, comparable_vcf, log_file -> tuple(patient_id, norm_vcf) }
            .groupTuple(by: 0)
        monovar_comparable_by_patient_ch = FILTER_MONOVAR_AND_SUBTRACT_GERMLINE.out
            .map { patient_id, cell_id, norm_vcf, comparable_vcf, log_file -> tuple(patient_id, comparable_vcf) }
            .groupTuple(by: 0)
        monovar_exclusion_by_patient_ch = PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION.out
            .map { patient_id, exclusion_vcf, exclusion_tbi, ref_fasta, exclusion_source -> tuple(patient_id, exclusion_vcf) }
            .groupTuple(by: 0)
        monovar_flowchart_input_ch = monovar_split_by_patient_ch
            .combine(monovar_norm_by_patient_ch, by: 0)
            .combine(monovar_comparable_by_patient_ch, by: 0)
            .combine(monovar_exclusion_by_patient_ch, by: 0)
        MONOVAR_VARIANT_FLOWCHARTS(monovar_flowchart_input_ch)
        MONOVAR_VARIANT_FLOWCHARTS.out.view { "MonoVar variant flowcharts: ${it[1]}" }

        monovar_comparison_ch = FILTER_MONOVAR_AND_SUBTRACT_GERMLINE.out
            .map { patient_id, cell_id, norm_vcf, comparable_vcf, log_file -> tuple(patient_id, comparable_vcf) }
            .groupTuple(by: 0)
        comparison_inputs = comparison_inputs.mix(monovar_comparison_ch)

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
            final_vep_by_patient_ch = MERGE_VEP_AND_POPULATION_COSMIC_FILTER.out
                .map { patient_id, cell_id, vep_tsv, final_vcf, summary_tsv, log_file -> tuple(patient_id, vep_tsv) }
                .groupTuple(by: 0)
            final_matrix_input_ch = final_vcfs_by_patient_ch.combine(final_vep_by_patient_ch, by: 0)
            BUILD_MUTATION_MATRICES(final_matrix_input_ch)
            BUILD_MUTATION_MATRICES.out.view { "Mutation matrix long table: ${it[1]}" }

            if (params.run_snv_report) {
                snv_report_breadth_ch = params.run_bam_qc
                    ? SUMMARIZE_BREADTH_QC.out.map { patient_id, breadth_file -> tuple(patient_id, breadth_file.toString()) }
                    : BUILD_MUTATION_MATRICES.out.map { patient_id, a, b, c, d, e -> tuple(patient_id, "") }
                snv_report_input_ch = BUILD_MUTATION_MATRICES.out
                    .combine(MONOVAR_VARIANT_FLOWCHARTS.out, by: 0)
                    .combine(snv_report_breadth_ch, by: 0)
                SNV_HTML_REPORT(snv_report_input_ch)
                SNV_HTML_REPORT.out.view { "SNV HTML report: ${it[1]}" }
            }

            split_vcfs_by_patient_ch = SPLIT_MONOVAR.out.map { patient_id, split_vcfs, sample_map, split_log -> tuple(patient_id, split_vcfs) }
            qc_inputs_ch = split_vcfs_by_patient_ch.combine(BUILD_MUTATION_MATRICES.out, by: 0).combine(filter_info_ch, by: 0)
            SUMMARIZE_FILTERS_QC(qc_inputs_ch)
            SUMMARIZE_FILTERS_QC.out.view { "Filter/QC summary: ${it[1]}" }

            bcftools_stats_by_patient_ch = BCFTOOLS_STATS.out
                .map { patient_id, cell_id, stats_file -> tuple(patient_id, stats_file) }
                .groupTuple(by: 0)

            if (params.run_bam_qc) {
                samtools_stats_by_patient_ch = SAMTOOLS_STATS.out
                    .map { patient_id, cell_id, stats_file -> tuple(patient_id, stats_file) }
                    .groupTuple(by: 0)
                // mosdepth_by_patient_ch and SUMMARIZE_BREADTH_QC are already
                // built once, earlier, in the top-level `if (params.run_bam_qc)`
                // block right after MOSDEPTH_QC -- reused here rather than
                // recomputed to avoid a duplicate process invocation.
                breadth_for_multiqc_ch = SUMMARIZE_BREADTH_QC.out
                    .map { patient_id, breadth_file -> tuple(patient_id, breadth_file.toString()) }
                multiqc_custom_input_ch = SUMMARIZE_FILTERS_QC.out.combine(breadth_for_multiqc_ch, by: 0)
                MAKE_MULTIQC_CUSTOM_CONTENT(multiqc_custom_input_ch)

                multiqc_full_input_ch = MAKE_MULTIQC_CUSTOM_CONTENT.out
                    .combine(bcftools_stats_by_patient_ch, by: 0)
                    .combine(samtools_stats_by_patient_ch, by: 0)
                    .combine(mosdepth_by_patient_ch, by: 0)
                MULTIQC_REPORT_WITH_BAM_QC(multiqc_full_input_ch)
                MULTIQC_REPORT_WITH_BAM_QC.out.view { "MultiQC report: ${it[1]}" }
            } else {
                multiqc_custom_input_ch = SUMMARIZE_FILTERS_QC.out
                    .map { patient_id, settings, prefilter_qc, impact -> tuple(patient_id, settings, prefilter_qc, impact, "") }
                MAKE_MULTIQC_CUSTOM_CONTENT(multiqc_custom_input_ch)

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

    if (params.run_sccaller) {
        sccaller_input_ch = checked_info_ch.flatMap { patient_id, ref_fasta, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_vcf ->
            cellMetadataRows(patient_id, cell_metadata).findAll { row -> !isLeukocyteCell(row[1]) }.collect { row ->
                tuple(row[0], row[1], row[2], ref_fasta, bulk_bam, resolveInputPath(params.sccaller_hsnp_vcf))
            }
        }
        RUN_SCCALLER_FROM_BAM(sccaller_input_ch)
        sccaller_by_patient_ch = RUN_SCCALLER_FROM_BAM.out
            .map { patient_id, cell_id, raw_vcf -> tuple(patient_id, [cell_id, raw_vcf]) }
            .groupTuple(by: 0)
        STANDARDIZE_INTERNAL_SCCALLER(sccaller_by_patient_ch.join(checked_info_ch, by: 0))
        comparison_inputs = comparison_inputs.mix(STANDARDIZE_INTERNAL_SCCALLER.out)
        STANDARDIZE_INTERNAL_SCCALLER.out.view { "Internal SCcaller comparable VCFs for ${it[0]}" }
    }

    if (params.run_comparison) {
        comparison_by_patient_ch = comparison_inputs
            .flatMap { patient_id, vcfs ->
                def files = vcfs instanceof List ? vcfs : [vcfs]
                files.collect { vcf -> tuple(patient_id, vcf) }
            }
            .groupTuple(by: 0)
        BUILD_CALLER_COMPARISON(comparison_by_patient_ch)
        BUILD_CALLER_COMPARISON.out.view { "Caller comparison matrices: ${it[2]}" }

        if (params.run_caller_scorecard) {
            BUILD_CALLER_SCORECARD(BUILD_CALLER_COMPARISON.out)
            BUILD_CALLER_SCORECARD.out.view { "Caller scorecard: ${it[1]}" }
        }

        if (params.run_final_report) {
            if (params.run_vep_filter) {
                final_report_base_ch = BUILD_MUTATION_MATRICES.out
                    .join(checked_info_ch, by: 0)
                if (params.run_bam_qc) {
                    final_report_samtools_ch = SAMTOOLS_STATS.out
                        .map { patient_id, cell_id, stats_file -> tuple(patient_id, stats_file) }
                        .groupTuple(by: 0)
                    final_report_mosdepth_ch = MOSDEPTH_QC.out
                        .map { patient_id, cell_id, mosdepth_files -> tuple(patient_id, mosdepth_files) }
                        .groupTuple(by: 0)
                        .map { patient_id, mosdepth_file_groups -> tuple(patient_id, mosdepth_file_groups.flatten()) }
                    final_report_with_qc_ch = final_report_base_ch
                        .combine(final_report_samtools_ch, by: 0)
                        .combine(final_report_mosdepth_ch, by: 0)
                    BUILD_FINAL_REPORT_DATA_FROM_FINAL_CALLS_WITH_ALIGNMENT_QC(final_report_with_qc_ch)
                    final_report_data_ch = BUILD_FINAL_REPORT_DATA_FROM_FINAL_CALLS_WITH_ALIGNMENT_QC.out
                } else {
                    BUILD_FINAL_REPORT_DATA_FROM_FINAL_CALLS(final_report_base_ch)
                    final_report_data_ch = BUILD_FINAL_REPORT_DATA_FROM_FINAL_CALLS.out
                }
            } else {
                final_report_base_ch = BUILD_CALLER_COMPARISON.out
                    .join(checked_info_ch, by: 0)
                if (params.run_bam_qc) {
                    final_report_samtools_ch = SAMTOOLS_STATS.out
                        .map { patient_id, cell_id, stats_file -> tuple(patient_id, stats_file) }
                        .groupTuple(by: 0)
                    final_report_mosdepth_ch = MOSDEPTH_QC.out
                        .map { patient_id, cell_id, mosdepth_files -> tuple(patient_id, mosdepth_files) }
                        .groupTuple(by: 0)
                        .map { patient_id, mosdepth_file_groups -> tuple(patient_id, mosdepth_file_groups.flatten()) }
                    final_report_with_qc_ch = final_report_base_ch
                        .combine(final_report_samtools_ch, by: 0)
                        .combine(final_report_mosdepth_ch, by: 0)
                    BUILD_FINAL_REPORT_DATA_WITH_ALIGNMENT_QC(final_report_with_qc_ch)
                    final_report_data_ch = BUILD_FINAL_REPORT_DATA_WITH_ALIGNMENT_QC.out
                } else {
                    BUILD_FINAL_REPORT_DATA(final_report_base_ch)
                    final_report_data_ch = BUILD_FINAL_REPORT_DATA.out
                }
            }
            CN_PUP_FINAL_REPORT(final_report_data_ch)
            CN_PUP_FINAL_REPORT.out.view { "CN-PUP final report: ${it[1]}" }
        }

        if (params.run_delsieve_prep || params.run_delsieve) {
            delsieve_patient_info_ch = checked_info_ch.map { patient_id, ref_fasta, cell_metadata, germline_mode, germline_vcf, bulk_bam, leukocyte_vcf ->
                tuple(patient_id, ref_fasta, cell_metadata)
            }
            delsieve_from_comparison_ch = BUILD_CALLER_COMPARISON.out
                .join(delsieve_patient_info_ch, by: 0)
                .map { patient_id, long_table, matrices, ctc_scite_inputs, overlap_qc, ref_fasta, cell_metadata ->
                    tuple(patient_id, long_table, ref_fasta, cell_metadata)
                }
            delsieve_source_ch = delsieve_source_ch.mix(delsieve_from_comparison_ch)
        }
    }

    if (params.run_delsieve_prep || params.run_delsieve || params.run_delsieve_from_existing) {
        PREPARE_DELSIEVE_INPUTS(delsieve_source_ch)
        PREPARE_DELSIEVE_INPUTS.out.view { "DelSIEVE prep inputs: ${it[1]}" }

        if (params.run_delsieve) {
            DELSIEVE_DATA_COLLECTOR(PREPARE_DELSIEVE_INPUTS.out)
            DELSIEVE_DATA_COLLECTOR.out.view { "DelSIEVE stage 1 XML: ${it[2]}" }

            RUN_DELSIEVE_STAGE1(DELSIEVE_DATA_COLLECTOR.out)
            RUN_DELSIEVE_STAGE1.out.view { "DelSIEVE stage 1 run: ${it[1]}" }

            if (params.run_delsieve_tree_annotator) {
                DELSIEVE_TREE_ANNOTATOR(RUN_DELSIEVE_STAGE1.out)
                DELSIEVE_TREE_ANNOTATOR.out.view { "DelSIEVE annotated tree: ${it[1]}" }

                if (params.run_delsieve_tree_png) {
                    DELSIEVE_TREE_PNG_STAGE1(DELSIEVE_TREE_ANNOTATOR.out)
                    DELSIEVE_TREE_PNG_STAGE1.out.view { "DelSIEVE stage 1 tree PNG: ${it[1]}" }
                }

                if (params.run_delsieve_variant_caller) {
                    DELSIEVE_VARIANT_CALLER_STAGE1(DELSIEVE_TREE_ANNOTATOR.out)
                    DELSIEVE_VARIANT_CALLER_STAGE1.out.view { "DelSIEVE stage 1 variants: ${it[1]}" }

                    if (params.run_delsieve_gene_annotator) {
                        DELSIEVE_GENE_ANNOTATOR_STAGE1(DELSIEVE_VARIANT_CALLER_STAGE1.out)
                        DELSIEVE_GENE_ANNOTATOR_STAGE1.out.view { "DelSIEVE stage 1 mutation-annotated tree: ${it[1]}" }

                        if (params.run_delsieve_tree_png) {
                            DELSIEVE_GENE_TREE_PNG_STAGE1(DELSIEVE_GENE_ANNOTATOR_STAGE1.out)
                            DELSIEVE_GENE_TREE_PNG_STAGE1.out.view { "DelSIEVE stage 1 mutation-annotated PNG: ${it[1]}" }
                        }
                    }
                }
            }
        }
    }
}

process STANDARDIZE_EXTERNAL_CALLERS {
    tag "$patient_id"
    conda { params.run_vep_filter ? "envs/vep.yml" : "envs/bcftools.yml" }
    publishDir { "${params.outdir}/${patient_id}/processed_calls" }, mode: 'copy', pattern: "comparison_project/processed_calls/**"
    publishDir { "${params.outdir}/${patient_id}/normalized_calls" }, mode: 'copy', pattern: "comparison_project/normalized_calls/**"
    publishDir { "${params.outdir}/${patient_id}/comparable_calls" }, mode: 'copy', pattern: "comparison_project/comparable_calls/**"
    publishDir { "${params.outdir}/${patient_id}/vep/external_existing" }, mode: 'copy', pattern: "comparison_project/vep_existing/**"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "comparison_project/logs/*"
    publishDir { "${params.outdir}/${patient_id}/reports/variant_flowcharts" }, mode: 'copy', pattern: "comparison_project/reports/variant_flowcharts/**"

    input:
    tuple val(patient_id), val(rows), val(ref_fasta), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_vcf)

    output:
    tuple val(patient_id), path("comparison_project/comparable_calls/*/*.comparable.vcf.gz")

    script:
    def tissue_vcfs_arg = params.tissue_vcfs ? "--tissue-vcfs ${resolveInputPath(params.tissue_vcfs)}" : ""
    def samples_tsv = rows.collect { r ->
        def mode = r[3] ?: params.sccaller_mode
        "${r[0]}\t${r[1].toString().toLowerCase()}\t${r[2]}\t${mode}"
    }.join('\n')
    """
    set -euo pipefail

    mkdir -p comparison_project/config comparison_project/bulk_exclusion comparison_project/logs
    cat > comparison_project/config/samples.tsv <<'EOF'
sample_id	caller	raw_vcf	sccaller_mode
${samples_tsv}
EOF

    exclusion_inputs=()
    if [[ ( "${germline_mode}" == "precomputed_vcf" || "${germline_mode}" == "combined" ) && -n "${germline_vcf}" && "${germline_vcf}" != "NA" ]]; then
      if [[ ! -s "${ref_fasta}.fai" ]]; then
        samtools faidx "${ref_fasta}"
      fi
      bcftools reheader -f "${ref_fasta}.fai" \\
        -o "comparison_project/bulk_exclusion/${patient_id}.bulk_blood.with_contigs.vcf" \\
        "${germline_vcf}"
      bcftools view -v snps -f PASS,. "comparison_project/bulk_exclusion/${patient_id}.bulk_blood.with_contigs.vcf" -Ou \\
      | bcftools norm -f "${ref_fasta}" -m -any -Oz \\
        -o "comparison_project/bulk_exclusion/${patient_id}.bulk_blood.norm.vcf.gz" \\
        2> "comparison_project/logs/${patient_id}.bulk_blood.norm.log"
      tabix -f -p vcf "comparison_project/bulk_exclusion/${patient_id}.bulk_blood.norm.vcf.gz"
      exclusion_inputs+=("comparison_project/bulk_exclusion/${patient_id}.bulk_blood.norm.vcf.gz")
    fi

    if [[ ( "${germline_mode}" == "leukocyte_vcf" || "${germline_mode}" == "combined" ) && -n "${leukocyte_vcf}" && "${leukocyte_vcf}" != "NA" ]]; then
      if [[ ! -s "${ref_fasta}.fai" ]]; then
        samtools faidx "${ref_fasta}"
      fi
      bcftools reheader -f "${ref_fasta}.fai" \\
        -o "comparison_project/bulk_exclusion/${patient_id}.leukocyte.with_contigs.vcf" \\
        "${leukocyte_vcf}"
      bcftools view -v snps -f PASS,. "comparison_project/bulk_exclusion/${patient_id}.leukocyte.with_contigs.vcf" -Ou \\
      | bcftools norm -f "${ref_fasta}" -m -any -Oz \\
        -o "comparison_project/bulk_exclusion/${patient_id}.leukocyte.norm.vcf.gz" \\
        2> "comparison_project/logs/${patient_id}.leukocyte.norm.log"
      tabix -f -p vcf "comparison_project/bulk_exclusion/${patient_id}.leukocyte.norm.vcf.gz"
      exclusion_inputs+=("comparison_project/bulk_exclusion/${patient_id}.leukocyte.norm.vcf.gz")
    fi

    if (( \${#exclusion_inputs[@]} == 1 )); then
      cp "\${exclusion_inputs[0]}" "comparison_project/bulk_exclusion/${patient_id}.external_exclusion.norm.vcf.gz"
      cp "\${exclusion_inputs[0]}.tbi" "comparison_project/bulk_exclusion/${patient_id}.external_exclusion.norm.vcf.gz.tbi"
    elif (( \${#exclusion_inputs[@]} > 1 )); then
      site_only_inputs=()
      for exclusion_vcf in "\${exclusion_inputs[@]}"; do
        site_only="\${exclusion_vcf%.vcf.gz}.sites_only.vcf.gz"
        bcftools view -G "\$exclusion_vcf" -Oz -o "\$site_only"
        tabix -f -p vcf "\$site_only"
        site_only_inputs+=("\$site_only")
      done
      bcftools concat -a "\${site_only_inputs[@]}" -Ou \\
      | bcftools sort -Ou \\
      | bcftools norm -d exact -Oz \\
        -o "comparison_project/bulk_exclusion/${patient_id}.external_exclusion.norm.vcf.gz" \\
        2> "comparison_project/logs/${patient_id}.external_exclusion.merge.log"
      tabix -f -p vcf "comparison_project/bulk_exclusion/${patient_id}.external_exclusion.norm.vcf.gz"
    fi

    export MIN_TOTAL_DEPTH="${params.min_total_depth}"
    export MIN_ALT_READS="${params.min_alt_reads}"
    export MIN_VAF="${params.min_vaf}"
    export SCCALLER_SOMATIC_MODE="${params.sccaller_mode}"
    export RUN_EXISTING_VEP_FILTER="${params.run_vep_filter}"
    export MAX_POPULATION_AF="${params.max_population_af}"
    export KEEP_COSMIC="${params.keep_cosmic}"
    export GENOME="${params.genome}"
    export VEP_SPECIES="${params.vep_species}"
    export VEP_CACHE="${params.vep_cache}"
    export VEP_CACHE_VERSION="${params.vep_cache_version}"
    export VEP_FASTA="${params.vep_fasta}"
    export VEP_FORKS="${params.vep_forks}"
    export VEP_BUFFER_SIZE="${params.vep_buffer_size}"

    bash "${projectDir}/bin/mutation_matrix/02_process_raw_calls.sh" \\
      "\$PWD/comparison_project" \\
      "${ref_fasta}" \\
      "\$PWD/comparison_project/config/samples.tsv" \\
      "${projectDir}"

    python "${projectDir}/bin/mutation_matrix/07_variant_flowcharts.py" \\
      --project "\$PWD/comparison_project" \\
      --samples-tsv "\$PWD/comparison_project/config/samples.tsv" \\
      --output-dir "\$PWD/comparison_project/reports/variant_flowcharts" ${tissue_vcfs_arg}
    """
}

process RUN_SCCALLER_FROM_BAM {
    tag { "${patient_id}:${cell_id}" }
    cpus params.sccaller_cpu
    conda "envs/sccaller.yml"
    publishDir { "${params.outdir}/${patient_id}/raw_calls/sccaller" }, mode: 'copy', pattern: "*.sccaller.vcf"
    publishDir { "${params.outdir}/${patient_id}/processed_calls/sccaller" }, mode: 'copy', pattern: "*.somatic.snv.vcf"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.sccaller.log"

    input:
    tuple val(patient_id), val(cell_id), val(cell_bam), val(ref_fasta), val(bulk_bam), val(hsnp_vcf)

    output:
    tuple val(patient_id), val(cell_id), path("${cell_id}.sccaller.vcf")

    script:
    def sccaller_extra_args = ""
    if (params.sccaller_head) {
        sccaller_extra_args += " --head ${params.sccaller_head}"
    }
    if (params.sccaller_tail) {
        sccaller_extra_args += " --tail ${params.sccaller_tail}"
    }
    """
    set -euo pipefail

    if [[ -z "${params.sccaller_script}" ]]; then
      echo "--sccaller_script is required when --run_sccaller true" >&2
      exit 1
    fi
    if [[ -z "${hsnp_vcf}" || "${hsnp_vcf}" == "NA" ]]; then
      echo "--sccaller_hsnp_vcf is required when --run_sccaller true" >&2
      exit 1
    fi
    if [[ -z "${bulk_bam}" || "${bulk_bam}" == "NA" ]]; then
      echo "bulk_bam in the patient sheet is required when --run_sccaller true" >&2
      exit 1
    fi

    mkdir -p sccaller_work
    python "${params.sccaller_script}" \\
      --bam "${cell_bam}" \\
      --fasta "${ref_fasta}" \\
      --bulk "${bulk_bam}" \\
      --output "${cell_id}.sccaller.vcf" \\
      --snp_type hsnp \\
      --snp_in "${hsnp_vcf}" \\
      --cpu_num "${params.sccaller_cpu}" \\
      --bias "${params.sccaller_bias}" \\
      --wkdir sccaller_work \\
      --engine "${params.sccaller_engine}"${sccaller_extra_args} \\
      > "${cell_id}.sccaller.log" 2>&1

    if [[ ! -s "${cell_id}.sccaller.vcf" ]]; then
      echo "ERROR: ${params.sccaller_script} did not produce ${cell_id}.sccaller.vcf -- see ${cell_id}.sccaller.log. A common cause on this reference is SCcaller choking on HLA/ALT decoy contigs; try setting --sccaller_head/--sccaller_tail to restrict it to the primary chromosomes (e.g. 1/22 for Homo_sapiens_assembly38.fasta)." >&2
      exit 1
    fi

    awk 'BEGIN{FS=OFS="\\t"} /^#/ {print; next} /0\\/1/ && /True/ && \$7=="." && length(\$5)==1 {split(\$10,a,":"); split(a[2],ad,","); if (ad[1]+ad[2]>=20) print}' \\
      "${cell_id}.sccaller.vcf" > "${cell_id}.somatic.snv.vcf"
    """
}

process STANDARDIZE_INTERNAL_SCCALLER {
    tag "$patient_id"
    conda { params.run_vep_filter ? "envs/vep.yml" : "envs/bcftools.yml" }
    publishDir { "${params.outdir}/${patient_id}/normalized_calls/sccaller" }, mode: 'copy', pattern: "comparison_project/normalized_calls/sccaller/*"
    publishDir { "${params.outdir}/${patient_id}/comparable_calls/sccaller" }, mode: 'copy', pattern: "comparison_project/comparable_calls/sccaller/*"
    publishDir { "${params.outdir}/${patient_id}/vep/external_existing" }, mode: 'copy', pattern: "comparison_project/vep_existing/**"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "comparison_project/logs/*"
    publishDir { "${params.outdir}/${patient_id}/reports/variant_flowcharts" }, mode: 'copy', pattern: "comparison_project/reports/variant_flowcharts/**"

    input:
    tuple val(patient_id), val(rows), val(ref_fasta), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_vcf)

    output:
    tuple val(patient_id), path("comparison_project/comparable_calls/sccaller/*.comparable.vcf.gz")

    script:
    def tissue_vcfs_arg = params.tissue_vcfs ? "--tissue-vcfs ${resolveInputPath(params.tissue_vcfs)}" : ""
    def samples_tsv = rows.collect { r ->
        "${r[0]}\tsccaller\t${r[1]}\t${params.sccaller_mode}"
    }.join('\n')
    """
    set -euo pipefail

    mkdir -p comparison_project/config comparison_project/bulk_exclusion comparison_project/logs
    cat > comparison_project/config/samples.tsv <<'EOF'
sample_id	caller	raw_vcf	sccaller_mode
${samples_tsv}
EOF

    exclusion_inputs=()
    if [[ ( "${germline_mode}" == "precomputed_vcf" || "${germline_mode}" == "combined" ) && -n "${germline_vcf}" && "${germline_vcf}" != "NA" ]]; then
      if [[ ! -s "${ref_fasta}.fai" ]]; then
        samtools faidx "${ref_fasta}"
      fi
      bcftools reheader -f "${ref_fasta}.fai" \\
        -o "comparison_project/bulk_exclusion/${patient_id}.bulk_blood.with_contigs.vcf" \\
        "${germline_vcf}"
      bcftools view -v snps -f PASS,. "comparison_project/bulk_exclusion/${patient_id}.bulk_blood.with_contigs.vcf" -Ou \\
      | bcftools norm -f "${ref_fasta}" -m -any -Oz \\
        -o "comparison_project/bulk_exclusion/${patient_id}.bulk_blood.norm.vcf.gz" \\
        2> "comparison_project/logs/${patient_id}.bulk_blood.norm.log"
      tabix -f -p vcf "comparison_project/bulk_exclusion/${patient_id}.bulk_blood.norm.vcf.gz"
      exclusion_inputs+=("comparison_project/bulk_exclusion/${patient_id}.bulk_blood.norm.vcf.gz")
    fi

    if [[ ( "${germline_mode}" == "leukocyte_vcf" || "${germline_mode}" == "combined" ) && -n "${leukocyte_vcf}" && "${leukocyte_vcf}" != "NA" ]]; then
      if [[ ! -s "${ref_fasta}.fai" ]]; then
        samtools faidx "${ref_fasta}"
      fi
      bcftools reheader -f "${ref_fasta}.fai" \\
        -o "comparison_project/bulk_exclusion/${patient_id}.leukocyte.with_contigs.vcf" \\
        "${leukocyte_vcf}"
      bcftools view -v snps -f PASS,. "comparison_project/bulk_exclusion/${patient_id}.leukocyte.with_contigs.vcf" -Ou \\
      | bcftools norm -f "${ref_fasta}" -m -any -Oz \\
        -o "comparison_project/bulk_exclusion/${patient_id}.leukocyte.norm.vcf.gz" \\
        2> "comparison_project/logs/${patient_id}.leukocyte.norm.log"
      tabix -f -p vcf "comparison_project/bulk_exclusion/${patient_id}.leukocyte.norm.vcf.gz"
      exclusion_inputs+=("comparison_project/bulk_exclusion/${patient_id}.leukocyte.norm.vcf.gz")
    fi

    if (( \${#exclusion_inputs[@]} == 1 )); then
      cp "\${exclusion_inputs[0]}" "comparison_project/bulk_exclusion/${patient_id}.external_exclusion.norm.vcf.gz"
      cp "\${exclusion_inputs[0]}.tbi" "comparison_project/bulk_exclusion/${patient_id}.external_exclusion.norm.vcf.gz.tbi"
    elif (( \${#exclusion_inputs[@]} > 1 )); then
      site_only_inputs=()
      for exclusion_vcf in "\${exclusion_inputs[@]}"; do
        site_only="\${exclusion_vcf%.vcf.gz}.sites_only.vcf.gz"
        bcftools view -G "\$exclusion_vcf" -Oz -o "\$site_only"
        tabix -f -p vcf "\$site_only"
        site_only_inputs+=("\$site_only")
      done
      bcftools concat -a "\${site_only_inputs[@]}" -Ou \\
      | bcftools sort -Ou \\
      | bcftools norm -d exact -Oz \\
        -o "comparison_project/bulk_exclusion/${patient_id}.external_exclusion.norm.vcf.gz" \\
        2> "comparison_project/logs/${patient_id}.external_exclusion.merge.log"
      tabix -f -p vcf "comparison_project/bulk_exclusion/${patient_id}.external_exclusion.norm.vcf.gz"
    fi

    export MIN_TOTAL_DEPTH="${params.min_total_depth}"
    export MIN_ALT_READS="${params.min_alt_reads}"
    export MIN_VAF="${params.min_vaf}"
    export SCCALLER_SOMATIC_MODE="${params.sccaller_mode}"
    export RUN_EXISTING_VEP_FILTER="${params.run_vep_filter}"
    export MAX_POPULATION_AF="${params.max_population_af}"
    export KEEP_COSMIC="${params.keep_cosmic}"
    export GENOME="${params.genome}"
    export VEP_SPECIES="${params.vep_species}"
    export VEP_CACHE="${params.vep_cache}"
    export VEP_CACHE_VERSION="${params.vep_cache_version}"
    export VEP_FASTA="${params.vep_fasta}"
    export VEP_FORKS="${params.vep_forks}"
    export VEP_BUFFER_SIZE="${params.vep_buffer_size}"

    bash "${projectDir}/bin/mutation_matrix/02_process_raw_calls.sh" \\
      "\$PWD/comparison_project" \\
      "${ref_fasta}" \\
      "\$PWD/comparison_project/config/samples.tsv" \\
      "${projectDir}"

    python "${projectDir}/bin/mutation_matrix/07_variant_flowcharts.py" \\
      --project "\$PWD/comparison_project" \\
      --samples-tsv "\$PWD/comparison_project/config/samples.tsv" \\
      --output-dir "\$PWD/comparison_project/reports/variant_flowcharts" ${tissue_vcfs_arg}
    """
}

process BUILD_CALLER_COMPARISON {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/caller_comparison" }, mode: 'copy'

    input:
    tuple val(patient_id), path(comparable_vcfs)

    output:
    tuple val(patient_id), path("cache/merged/all_mutations_long.annotated.tsv.gz"), path("matrices/*.tsv"), path("ctc_scite_inputs/*.tsv"), path("reports/overlap_qc.tsv")

    script:
    """
    set -euo pipefail

    mkdir -p comparable_calls cache/per_sample cache/merged matrices ctc_scite_inputs reports
    for vcf in ${comparable_vcfs}; do
      base=\$(basename "\$vcf")
      caller=""
      staged_base=""
      if [[ "\$base" =~ ^(.+)\\.([^.]+)\\.comparable\\.vcf\\.gz\$ ]]; then
        caller="\${BASH_REMATCH[2]}"
        staged_base="\$base"
      elif [[ "\$base" =~ ^(.+)\\.monovar\\.no_.+\\.vcf\\.gz\$ ]]; then
        caller="monovar"
        staged_base="\${BASH_REMATCH[1]}.monovar.comparable.vcf.gz"
      elif [[ "\$base" =~ ^(.+)\\.monovar\\.population_cosmic\\.vcf\\.gz\$ ]]; then
        caller="monovar"
        staged_base="\${BASH_REMATCH[1]}.monovar.comparable.vcf.gz"
      else
        echo "ERROR: cannot infer sample/caller from comparable VCF name: \$base" >&2
        exit 1
      fi
      mkdir -p "comparable_calls/\$caller"
      cp "\$vcf" "comparable_calls/\$caller/\$staged_base"
      if [[ -f "\${vcf}.tbi" ]]; then
        cp "\${vcf}.tbi" "comparable_calls/\$caller/\$staged_base.tbi"
      fi
    done

    python "${projectDir}/bin/mutation_matrix/03_vcf_to_long_table.py" \\
      --project "\$PWD" \\
      --input-root "\$PWD/comparable_calls" \\
      --output-dir "\$PWD/cache/per_sample"

    python "${projectDir}/bin/mutation_matrix/04_merge_long_tables.py" \\
      --input-dir "\$PWD/cache/per_sample" \\
      --output "\$PWD/cache/merged/all_mutations_long.tsv.gz"

    python "${projectDir}/bin/mutation_matrix/04_fill_missing_annotations.py" \\
      --input "\$PWD/cache/merged/all_mutations_long.tsv.gz" \\
      --output "\$PWD/cache/merged/all_mutations_long.annotated.tsv.gz"

    python "${projectDir}/bin/mutation_matrix/05_build_matrices.py" \\
      --input "\$PWD/cache/merged/all_mutations_long.annotated.tsv.gz" \\
      --matrix-dir "\$PWD/matrices" \\
      --ctc-scite-dir "\$PWD/ctc_scite_inputs"

    python "${projectDir}/bin/mutation_matrix/06_overlap_qc.py" \\
      --input "\$PWD/cache/merged/all_mutations_long.annotated.tsv.gz" \\
      --output "\$PWD/reports/overlap_qc.tsv"
    """
}

process BUILD_CALLER_SCORECARD {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/caller_comparison/reports" }, mode: 'copy'

    input:
    tuple val(patient_id), path(long_table), path(matrices), path(ctc_scite_inputs), path(overlap_qc)

    output:
    tuple val(patient_id), path("caller_scorecard.tsv"), path("statistical_tests.tsv"), path("scorecard_summary.md")

    script:
    """
    set -euo pipefail

    mkdir -p matrix_dir
    cp ${matrices} matrix_dir/

    python "${projectDir}/bin/mutation_matrix/08_caller_scorecard.py" \\
      --matrix-dir "\$PWD/matrix_dir" \\
      --output-dir "\$PWD"
    """
}

process BUILD_FINAL_REPORT_DATA {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/final_report/data" }, mode: 'copy', pattern: "final_report_data/*"

    input:
    tuple val(patient_id), path(long_table), path(matrices), path(ctc_scite_inputs), path(overlap_qc), val(ref_fasta), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_vcf)

    output:
    tuple val(patient_id), path("final_report_data")

    script:
    """
    set -euo pipefail

    mkdir -p final_report_data final_report_matrices
    cp ${matrices} final_report_matrices/

    python "${projectDir}/bin/build_final_report_data.py" \
      --patient-id "${patient_id}" \
      --long-table "${long_table}" \
      --matrices-dir "\$PWD/final_report_matrices" \
      --cell-metadata "${cell_metadata}" \
      --output-dir "\$PWD/final_report_data" \
      --top-variant-limit "${params.final_report_top_variants}" \
      --rarefaction-iterations "${params.final_report_rarefaction_iterations}"
    """
}

process BUILD_FINAL_REPORT_DATA_WITH_ALIGNMENT_QC {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/final_report/data" }, mode: 'copy', pattern: "final_report_data/*"

    input:
    tuple val(patient_id), path(long_table), path(matrices), path(ctc_scite_inputs), path(overlap_qc), val(ref_fasta), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_vcf), path(samtools_stats_files), path(mosdepth_files)

    output:
    tuple val(patient_id), path("final_report_data")

    script:
    """
    set -euo pipefail

    mkdir -p final_report_data final_report_matrices final_alignment_qc
    cp ${matrices} final_report_matrices/
    cp ${samtools_stats_files} final_alignment_qc/
    cp ${mosdepth_files} final_alignment_qc/

    python "${projectDir}/bin/build_final_report_data.py" \
      --patient-id "${patient_id}" \
      --long-table "${long_table}" \
      --matrices-dir "\$PWD/final_report_matrices" \
      --cell-metadata "${cell_metadata}" \
      --alignment-qc-dir "\$PWD/final_alignment_qc" \
      --output-dir "\$PWD/final_report_data" \
      --top-variant-limit "${params.final_report_top_variants}" \
      --rarefaction-iterations "${params.final_report_rarefaction_iterations}"
    """
}

process BUILD_FINAL_REPORT_DATA_FROM_FINAL_CALLS {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/final_report/data" }, mode: 'copy', pattern: "final_report_data/*"

    input:
    tuple val(patient_id), path(long_table), path(binary_matrix), path(altread_matrix), path(refread_matrix), path(summary_table), val(ref_fasta), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_vcf)

    output:
    tuple val(patient_id), path("final_report_data")

    script:
    """
    set -euo pipefail

    mkdir -p final_report_data final_report_matrices
    cp "${binary_matrix}" final_report_matrices/mutation_binary.tsv

    python "${projectDir}/bin/build_final_report_data.py" \
      --patient-id "${patient_id}" \
      --long-table "${long_table}" \
      --matrices-dir "\$PWD/final_report_matrices" \
      --cell-metadata "${cell_metadata}" \
      --output-dir "\$PWD/final_report_data" \
      --top-variant-limit "${params.final_report_top_variants}" \
      --rarefaction-iterations "${params.final_report_rarefaction_iterations}"
    """
}

process BUILD_FINAL_REPORT_DATA_FROM_FINAL_CALLS_WITH_ALIGNMENT_QC {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/final_report/data" }, mode: 'copy', pattern: "final_report_data/*"

    input:
    tuple val(patient_id), path(long_table), path(binary_matrix), path(altread_matrix), path(refread_matrix), path(summary_table), val(ref_fasta), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_vcf), path(samtools_stats_files), path(mosdepth_files)

    output:
    tuple val(patient_id), path("final_report_data")

    script:
    """
    set -euo pipefail

    mkdir -p final_report_data final_report_matrices final_alignment_qc
    cp "${binary_matrix}" final_report_matrices/mutation_binary.tsv
    cp ${samtools_stats_files} final_alignment_qc/
    cp ${mosdepth_files} final_alignment_qc/

    python "${projectDir}/bin/build_final_report_data.py" \
      --patient-id "${patient_id}" \
      --long-table "${long_table}" \
      --matrices-dir "\$PWD/final_report_matrices" \
      --cell-metadata "${cell_metadata}" \
      --alignment-qc-dir "\$PWD/final_alignment_qc" \
      --output-dir "\$PWD/final_report_data" \
      --top-variant-limit "${params.final_report_top_variants}" \
      --rarefaction-iterations "${params.final_report_rarefaction_iterations}"
    """
}

process CN_PUP_FINAL_REPORT {
    tag "$patient_id"
    conda "envs/snv_report.yml"
    publishDir { "${params.outdir}/${patient_id}/final_report" }, mode: 'copy'

    input:
    tuple val(patient_id), path(report_data)

    output:
    tuple val(patient_id), path("${patient_id}.cnpup_final_report.html")

    script:
    """
    set -euo pipefail

    mkdir -p report_data
    cp -r ${report_data}/* report_data/

    Rscript -e 'rmarkdown::render(
      "${projectDir}/assets/cnpup_final_report.Rmd",
      output_file="${patient_id}.cnpup_final_report.html",
      output_dir=".",
      knit_root_dir=getwd(),
      params=list(
        patient_id="${patient_id}",
        data_dir=file.path(Sys.getenv("PWD"), "report_data"),
        out_prefix="${patient_id}.cnpup_final"
      )
    )'
    """
}

process PREPARE_DELSIEVE_INPUTS {
    tag "$patient_id"
    conda "envs/delsieve_prep.yml"
    publishDir { "${params.outdir}/${patient_id}/delsieve_prep" }, mode: 'copy', pattern: "delsieve_prep/*"

    input:
    tuple val(patient_id), path(long_table), val(ref_fasta), val(cell_metadata)

    output:
    tuple val(patient_id), path("delsieve_prep")

    script:
    """
    set -euo pipefail

    # DelSIEVE DataCollector consumes DataFilter-style sections, not a flat numeric matrix.
    mkdir -p delsieve_prep

    python "${projectDir}/bin/prepare_delsieve_inputs.py" select-candidates \\
      --long-table "${long_table}" \\
      --cell-metadata "${cell_metadata}" \\
      --output-dir "\$PWD/delsieve_prep" \\
      --min-ctc-support "${params.delsieve_min_ctc_support}" \\
      --min-total-depth "${params.min_total_depth}" \\
      --min-alt-reads "${params.min_alt_reads}" \\
      --min-vaf "${params.min_vaf}" \\
      --callers "${params.delsieve_callers}"

    if [[ -s "\$PWD/delsieve_prep/candidate_sites.bed" ]]; then
      samtools mpileup \\
        -aa -A -B \\
        -Q "${params.delsieve_min_base_quality}" \\
        -q "${params.delsieve_min_mapping_quality}" \\
        -f "${ref_fasta}" \\
        -l "\$PWD/delsieve_prep/candidate_sites.bed" \\
        -b "\$PWD/delsieve_prep/bams.txt" \\
        > "\$PWD/delsieve_prep/candidate_sites.mpileup"

      python "${projectDir}/bin/prepare_delsieve_inputs.py" parse-mpileup \\
        --mpileup "\$PWD/delsieve_prep/candidate_sites.mpileup" \\
        --output-dir "\$PWD/delsieve_prep"
    else
      touch "\$PWD/delsieve_prep/candidate_sites.mpileup"
      sample_count=\$(awk '{print NF; exit}' "\$PWD/delsieve_prep/cell_names")
      {
        printf "=numSamples=\\n\$sample_count\\n"
        printf "=numCandidateMutatedSites=\\n0\\n"
        printf "=numBackgroundSites=\\n0\\n"
        printf "=mutations=\\n=background=\\n"
        printf "0,0\\n0,0\\n0,0\\n0,0\\n0,0\\n"
      } > "\$PWD/delsieve_prep/read_counts.full_support_coverage.tsv"
      printf "var_id\tchrom\tpos\tref\tcandidate_alt\tcell_id\tA\tC\tG\tT\talt1\talt2\talt3\tref_count\tcoverage\\n" > "\$PWD/delsieve_prep/read_counts.human_readable.tsv"
    fi
    """
}

process DELSIEVE_DATA_COLLECTOR {
    tag "$patient_id"
    publishDir { "${params.outdir}/${patient_id}/delsieve_stage1" }, mode: 'copy', pattern: "delsieve_stage1/*"

    input:
    tuple val(patient_id), path(delsieve_prep)

    output:
    tuple val(patient_id), path(delsieve_prep), path("delsieve_stage1/${patient_id}.stage1.xml"), path("delsieve_stage1/${patient_id}.datacollector.log")

    script:
    def sample_args = (params.delsieve_sample_cells.toInteger() > 0 && params.delsieve_sample_loci.toInteger() > 0) ? "-sample ${params.delsieve_sample_cells} ${params.delsieve_sample_loci}" : ""
    def ignore_sex_arg = params.delsieve_ignore_sex ? "-ignoreSex" : ""
    def version_file_args = params.delsieve_version_file ? "-version_file ${resolveInputPath(params.delsieve_version_file)}" : ""
    """
    set -euo pipefail

    unset DISPLAY
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.awt.headless=true"

    mkdir -p delsieve_stage1

    if [[ -z "${params.delsieve_template_stage1}" ]]; then
      echo "delsieve_template_stage1 is required when run_delsieve is true" >&2
      exit 1
    fi
    if [[ ! -s "${params.delsieve_template_stage1}" ]]; then
      echo "DelSIEVE stage 1 template not found: ${params.delsieve_template_stage1}" >&2
      exit 1
    fi
    if [[ \$(wc -l < "${delsieve_prep}/read_counts.full_support_coverage.tsv") -eq 0 ]]; then
      echo "No DelSIEVE candidate sites were found. Check ${delsieve_prep}/candidate_summary.tsv." >&2
      exit 1
    fi
    cell_count=\$(awk '{print NF; exit}' "${delsieve_prep}/cell_names")
    if [[ "${params.delsieve_datatype}" != "1" ]]; then
      echo "CN-PUP currently writes DelSIEVE datatype 1 input; got delsieve_datatype=${params.delsieve_datatype}." >&2
      exit 1
    fi
    for marker in '=numSamples=' '=numCandidateMutatedSites=' '=numBackgroundSites=' '=mutations=' '=background='; do
      if ! grep -qx "\$marker" "${delsieve_prep}/read_counts.full_support_coverage.tsv"; then
        echo "DelSIEVE data file is missing required DataFilter marker: \$marker" >&2
        exit 1
      fi
    done
    declared_cells=\$(awk 'prev == "=numSamples=" { print; exit } { prev=\$0 }' "${delsieve_prep}/read_counts.full_support_coverage.tsv")
    if [[ "\$declared_cells" -ne "\$cell_count" ]]; then
      echo "DelSIEVE data declares \$declared_cells cells, but ${delsieve_prep}/cell_names contains \$cell_count cells." >&2
      exit 1
    fi
    declared_sites=\$(awk 'prev == "=numCandidateMutatedSites=" { print; exit } { prev=\$0 }' "${delsieve_prep}/read_counts.full_support_coverage.tsv")
    if [[ "\$declared_sites" -eq 0 ]]; then
      echo "No DelSIEVE candidate sites were found. Check ${delsieve_prep}/candidate_summary.tsv." >&2
      exit 1
    fi
    usable_sites=\$(awk -v ignore_sex="${params.delsieve_ignore_sex}" '
      \$0 == "=mutations=" { in_mut=1; next }
      \$0 == "=background=" { in_mut=0 }
      in_mut {
        chrom = toupper(\$1)
        sub(/^CHR/, "", chrom)
        if (chrom == "M" || chrom == "MT") next
        if (ignore_sex == "true" && (chrom == "X" || chrom == "Y")) next
        n++
      }
      END { print n + 0 }
    ' "${delsieve_prep}/read_counts.full_support_coverage.tsv")
    if [[ "\$usable_sites" -eq 0 ]]; then
      echo "No DelSIEVE candidate sites remain after DelSIEVE chromosome filtering. Check mitochondrial/sex chromosome settings." >&2
      exit 1
    fi
    bad_mutation_row=\$(awk -v expected="\$((cell_count + 4))" '
      \$0 == "=mutations=" { in_mut=1; next }
      \$0 == "=background=" { in_mut=0 }
      in_mut && NF != expected { print NR ":" NF; exit }
    ' "${delsieve_prep}/read_counts.full_support_coverage.tsv")
    if [[ -n "\$bad_mutation_row" ]]; then
      echo "DelSIEVE mutation row has the wrong number of columns at row \$bad_mutation_row; expected \$((cell_count + 4))." >&2
      echo "Rows must be: chrom pos ref alt followed by one 'alt1,alt2,alt3;count1,count2,count3,coverage' field per cell." >&2
      exit 1
    fi
    background_rows=\$(awk '
      \$0 == "=background=" { in_bg=1; next }
      in_bg && NF > 0 { n++ }
      END { print n + 0 }
    ' "${delsieve_prep}/read_counts.full_support_coverage.tsv")
    if [[ "\$background_rows" -ne 5 ]]; then
      echo "DelSIEVE datatype 1 requires five background rows; found \$background_rows." >&2
      exit 1
    fi

    set +e
    "${params.delsieve_applauncher}" ${version_file_args} DataCollectorLauncher \\
      -cell "${delsieve_prep}/cell_names" \\
      -data "${delsieve_prep}/read_counts.full_support_coverage.tsv" \\
      -datatype "${params.delsieve_datatype}" \\
      -template "${params.delsieve_template_stage1}" \\
      -out "\$PWD/delsieve_stage1/${patient_id}.stage1.xml" \\
      -miss "${params.delsieve_missing_threshold}" \\
      ${sample_args} \\
      ${ignore_sex_arg} \\
      ${params.delsieve_datacollector_extra_args} \\
      > "\$PWD/delsieve_stage1/${patient_id}.datacollector.log" 2>&1
    datacollector_status=\$?
    set -e

    if [[ "\$datacollector_status" -ne 0 ]]; then
      echo "DataCollectorLauncher exited with status \$datacollector_status" >&2
      echo "DataCollector log:" >&2
      cat "\$PWD/delsieve_stage1/${patient_id}.datacollector.log" >&2
      echo "Files produced under delsieve_stage1:" >&2
      find "\$PWD/delsieve_stage1" -maxdepth 2 -type f -print >&2
      exit "\$datacollector_status"
    fi

    if [[ ! -s "\$PWD/delsieve_stage1/${patient_id}.stage1.xml" ]]; then
      echo "DataCollectorLauncher finished but did not create delsieve_stage1/${patient_id}.stage1.xml" >&2
      echo "DataCollector log:" >&2
      cat "\$PWD/delsieve_stage1/${patient_id}.datacollector.log" >&2
      echo "Files produced under delsieve_stage1:" >&2
      find "\$PWD/delsieve_stage1" -maxdepth 2 -type f -print >&2
      exit 1
    fi
    """
}

process RUN_DELSIEVE_STAGE1 {
    tag "$patient_id"
    cpus params.delsieve_threads
    publishDir { "${params.outdir}/${patient_id}/delsieve_stage1_run" }, mode: 'copy', pattern: "delsieve_stage1_run/*"

    input:
    tuple val(patient_id), path(delsieve_prep), path(stage1_xml), path(datacollector_log)

    output:
    tuple val(patient_id), path("delsieve_stage1_run"), path(stage1_xml), path(delsieve_prep)

    script:
    """
    set -euo pipefail

    unset DISPLAY
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.awt.headless=true"

    mkdir -p delsieve_stage1_run

    "${params.delsieve_beast}" \\
      -overwrite \\
      -threads "${params.delsieve_threads}" \\
      -prefix "\$PWD/delsieve_stage1_run/" \\
      ${params.delsieve_beast_extra_args} \\
      "${stage1_xml}" \\
      > "\$PWD/delsieve_stage1_run/${patient_id}.beast.log" 2>&1

    find "\$PWD" -maxdepth 1 -type f \\( -name "*.trees" -o -name "*.tree" -o -name "*.log" \\) \\
      ! -name ".command.*" \\
      -exec cp -n {} "\$PWD/delsieve_stage1_run/" \\;

    cp "${stage1_xml}" "\$PWD/delsieve_stage1_run/${patient_id}.stage1.xml"
    cp "${datacollector_log}" "\$PWD/delsieve_stage1_run/${patient_id}.datacollector.log"
    """
}

process DELSIEVE_TREE_ANNOTATOR {
    tag "$patient_id"
    publishDir { "${params.outdir}/${patient_id}/delsieve_stage1_tree" }, mode: 'copy', pattern: "delsieve_stage1_tree/*"

    input:
    tuple val(patient_id), path(stage1_run), path(stage1_xml), path(delsieve_prep)

    output:
    tuple val(patient_id), path("delsieve_stage1_tree/*"), path(stage1_run), path(delsieve_prep)

    script:
    def version_file_args = params.delsieve_version_file ? "-version_file ${resolveInputPath(params.delsieve_version_file)}" : ""
    """
    set -euo pipefail

    unset DISPLAY
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.awt.headless=true"

    mkdir -p delsieve_stage1_tree

    tree_file=\$(find -L "${stage1_run}" -maxdepth 1 -type f \\( -name "*.trees" -o -name "*.tree" \\) | head -n 1)
    if [[ -z "\$tree_file" ]]; then
      echo "No BEAST tree file found in ${stage1_run}" >&2
      exit 1
    fi

    "${params.delsieve_treeannotator}" ${version_file_args} ScsTreeAnnotatorLauncher \\
      -burnin "${params.delsieve_tree_burnin}" \\
      -heights "${params.delsieve_tree_heights}" \\
      "\$tree_file" \\
      "\$PWD/delsieve_stage1_tree/${patient_id}.stage1.mcc.tree" \\
      > "\$PWD/delsieve_stage1_tree/${patient_id}.treeannotator.log" 2>&1

    cp "\$tree_file" "\$PWD/delsieve_stage1_tree/\$(basename "\$tree_file")"
    cp "${stage1_xml}" "\$PWD/delsieve_stage1_tree/${patient_id}.stage1.xml"
    cp "${delsieve_prep}/candidate_sites.tsv" "\$PWD/delsieve_stage1_tree/candidate_sites.tsv"
    cp "${delsieve_prep}/cell_names.tsv" "\$PWD/delsieve_stage1_tree/cell_names.tsv"
    """
}

process DELSIEVE_TREE_PNG_STAGE1 {
    tag "$patient_id"
    conda "envs/delsieve_tree_plot.yml"
    publishDir { "${params.outdir}/${patient_id}/delsieve_stage1_tree" }, mode: 'copy', pattern: "delsieve_stage1_tree_png/*"

    input:
    tuple val(patient_id), path(stage1_tree_files), path(stage1_run), path(delsieve_prep)

    output:
    tuple val(patient_id), path("delsieve_stage1_tree_png/*")

    script:
    """
    set -euo pipefail

    unset DISPLAY
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.awt.headless=true"

    mkdir -p delsieve_stage1_tree_png

    tree_file=\$(find -L . -maxdepth 4 -type f -name "${patient_id}.stage1.mcc.tree" | head -n 1)
    if [[ -z "\$tree_file" ]]; then
      echo "No DelSIEVE MCC tree found for PNG rendering" >&2
      find -L . -maxdepth 4 -type f -print >&2
      exit 1
    fi

    plot_tree="\$PWD/delsieve_stage1_tree_png/${patient_id}.stage1.mcc.plot_clipped.tree"
    python "${projectDir}/bin/clip_tree_branch_lengths.py" \
      --input "\$tree_file" \
      --output "\$plot_tree" \
      --summary "\$PWD/delsieve_stage1_tree_png/${patient_id}.stage1.mcc.branch_length_clipping.tsv"

    "${params.delsieve_figtree}" -graphic PNG ${params.delsieve_figtree_extra_args} \\
      "\$plot_tree" \\
      "\$PWD/delsieve_stage1_tree_png/${patient_id}.stage1.mcc.figtree.png" \\
      > "\$PWD/delsieve_stage1_tree_png/${patient_id}.stage1.mcc.figtree.log" 2>&1

    cp "\$tree_file" "\$PWD/delsieve_stage1_tree_png/${patient_id}.stage1.mcc.tree"
    """
}

process DELSIEVE_VARIANT_CALLER_STAGE1 {
    tag "$patient_id"
    cpus params.delsieve_threads
    publishDir { "${params.outdir}/${patient_id}/delsieve_stage1_variants" }, mode: 'copy', pattern: "delsieve_stage1_variants/*"

    input:
    tuple val(patient_id), path(stage1_tree_files), path(stage1_run), path(delsieve_prep)

    output:
    tuple val(patient_id), path("delsieve_stage1_variants")

    script:
    def version_file_args = params.delsieve_version_file ? "-version_file ${resolveInputPath(params.delsieve_version_file)}" : ""
    """
    set -euo pipefail

    unset DISPLAY
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.awt.headless=true"

    mkdir -p delsieve_stage1_variants

    mcc_tree=\$(find -L . -maxdepth 4 -type f -name "${patient_id}.stage1.mcc.tree" | head -n 1)
    if [[ -z "\$mcc_tree" ]]; then
      echo "No DelSIEVE MCC tree found for ${patient_id}" >&2
      find -L . -maxdepth 4 -type f -print >&2
      exit 1
    fi

    mcmc_log=\$(find -L "${stage1_run}" -maxdepth 1 -type f -name "*.log" ! -name "*.beast.log" ! -name "*.datacollector.log" | head -n 1)
    if [[ -z "\$mcmc_log" ]]; then
      echo "No DelSIEVE MCMC log found in ${stage1_run}" >&2
      find -L "${stage1_run}" -maxdepth 1 -type f -print >&2
      exit 1
    fi

    stage1_xml_file=\$(find -L . -maxdepth 4 -type f -name "${patient_id}.stage1.xml" | head -n 1)
    if [[ -z "\$stage1_xml_file" ]]; then
      echo "No DelSIEVE stage 1 XML found for ${patient_id}" >&2
      find -L . -maxdepth 4 -type f -print >&2
      exit 1
    fi

    "${params.delsieve_applauncher}" ${version_file_args} VariantCallerLauncher \\
      -threads "${params.delsieve_threads}" \\
      -prefix "\$PWD/delsieve_stage1_variants/" \\
      -burnin "${params.delsieve_variant_burnin}" \\
      -estimates "${params.delsieve_variant_estimates}" \\
      -mcmclog "\$mcmc_log" \\
      -config "\$stage1_xml_file" \\
      -tree "\$mcc_tree" \\
      -details \\
      ${params.delsieve_variant_extra_args} \\
      > "\$PWD/delsieve_stage1_variants/${patient_id}.variantcaller.log" 2>&1

    cp "\$mcc_tree" "\$PWD/delsieve_stage1_variants/${patient_id}.stage1.mcc.tree"
    cp "\$mcmc_log" "\$PWD/delsieve_stage1_variants/\$(basename "\$mcmc_log")"
    cp "\$stage1_xml_file" "\$PWD/delsieve_stage1_variants/${patient_id}.stage1.xml"
    cp "${delsieve_prep}/candidate_sites.tsv" "\$PWD/delsieve_stage1_variants/candidate_sites.tsv"
    cp "${delsieve_prep}/cell_names.tsv" "\$PWD/delsieve_stage1_variants/cell_names.tsv"
    awk 'BEGIN { FS=OFS="\t" }
      NR == 1 {
        for (i = 1; i <= NF; i++) idx[\$i] = i
        next
      }
      {
        chrom = \$(idx["chrom"])
        pos = \$(idx["pos"])
        var_id = \$(idx["var_id"])
        gene = (idx["gene_symbol"] ? \$(idx["gene_symbol"]) : "")
        gsub(/[;,][[:space:]]*/, ",", gene)
        if (gene == "" || gene == "." || tolower(gene) == "unknown" || tolower(gene) == "unannotated" || tolower(gene) == "na" || tolower(gene) == "n/a") {
          gene = var_id
        }
        print chrom, pos, pos, gene
      }
    ' "${delsieve_prep}/candidate_sites.tsv" > "\$PWD/delsieve_stage1_variants/mutation_map.tsv"

    {
      echo "patient_id	stage	mcc_tree	mcmc_log	config_xml	variantcaller_log"
      echo "${patient_id}	stage1	${patient_id}.stage1.mcc.tree	\$(basename "\$mcmc_log")	${patient_id}.stage1.xml	${patient_id}.variantcaller.log"
      echo
      echo "Produced files:"
      find "\$PWD/delsieve_stage1_variants" -maxdepth 1 -type f -printf "%f\\n" | sort
    } > "\$PWD/delsieve_stage1_variants/${patient_id}.stage1_variantcaller_manifest.tsv"
    """
}

process DELSIEVE_GENE_ANNOTATOR_STAGE1 {
    tag "$patient_id"
    publishDir { "${params.outdir}/${patient_id}/delsieve_stage1_gene_tree" }, mode: 'copy', pattern: "delsieve_stage1_gene_tree/*"

    input:
    tuple val(patient_id), path(delsieve_stage1_variants)

    output:
    tuple val(patient_id), path("delsieve_stage1_gene_tree/*")

    script:
    def version_file_args = params.delsieve_version_file ? "-version_file ${resolveInputPath(params.delsieve_version_file)}" : ""
    def supplied_mutation_map = params.delsieve_mutation_map ? resolveInputPath(params.delsieve_mutation_map) : ""
    def annovar_args = params.delsieve_mutation_map_annovar ? "-anv ${params.delsieve_annovar_separator} -anvfiltercol ${params.delsieve_annovar_filter_col}" : ""
    def annovar_filter_args = params.delsieve_annovar_filter_keywords ? "-anvfilterkws ${params.delsieve_annovar_filter_keywords}" : ""
    def gene_filter_args = params.delsieve_gene_filter ? "-filter ${resolveInputPath(params.delsieve_gene_filter)} -sep ${params.delsieve_gene_filter_sep} -col ${params.delsieve_gene_filter_col}" : ""
    """
    set -euo pipefail

    unset DISPLAY
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.awt.headless=true"

    mkdir -p delsieve_stage1_gene_tree

    mutation_map="${supplied_mutation_map}"
    if [[ -z "\$mutation_map" ]]; then
      mutation_map="${delsieve_stage1_variants}/mutation_map.tsv"
    fi
    test -s "\$mutation_map" || { echo "DelSIEVE mutation map not found or empty: \$mutation_map" >&2; exit 1; }

    details_input="\$PWD/delsieve_stage1_geneannotator_details"
    mkdir -p "\$details_input"
    for detail_file in "${delsieve_stage1_variants}"/*; do
      if [[ -f "\$detail_file" ]]; then
        case "\$(basename "\$detail_file")" in
          *.loci_info|*.ternary)
            awk 'NF > 0' "\$detail_file" > "\$details_input/\$(basename "\$detail_file")"
            ;;
          *)
            cp "\$detail_file" "\$details_input/\$(basename "\$detail_file")"
            ;;
        esac
      fi
    done

    loci_info=\$(find "\$details_input" -maxdepth 1 -type f -name "*.loci_info" | head -n 1)
    test -s "\$loci_info" || { echo "DelSIEVE GeneAnnotator loci_info detail file is missing or empty" >&2; exit 1; }
    awk 'NF != 4 { print "Malformed DelSIEVE loci_info row " NR ": " \$0 > "/dev/stderr"; exit 1 }' "\$loci_info"

    "${params.delsieve_applauncher}" ${version_file_args} GeneAnnotatorLauncher \\
      -details "\$details_input" \\
      -subst "${params.delsieve_geneannotator_subst}" \\
      -map "\$mutation_map" \\
      ${annovar_args} \\
      ${annovar_filter_args} \\
      ${gene_filter_args} \\
      ${params.delsieve_geneannotator_extra_args} \\
      -out "\$PWD/delsieve_stage1_gene_tree/${patient_id}.stage1.mutation_annotated.tree" \\
      > "\$PWD/delsieve_stage1_gene_tree/${patient_id}.geneannotator.log" 2>&1

    expected_tree="\$PWD/delsieve_stage1_gene_tree/${patient_id}.stage1.mutation_annotated.tree"
    if [[ ! -s "\$expected_tree" ]]; then
      produced_tree=\$(find "\$PWD" -maxdepth 2 -type f \\( -name "*.tree" -o -name "*.trees" \\) ! -name "${patient_id}.stage1.mcc.tree" | head -n 1)
      if [[ -n "\$produced_tree" ]]; then
        cp "\$produced_tree" "\$expected_tree"
      else
        cp "${delsieve_stage1_variants}/${patient_id}.stage1.mcc.tree" "\$expected_tree"
        {
          echo "WARNING: GeneAnnotatorLauncher completed but did not create an annotated tree."
          echo "CN-PUP copied the MCC tree to ${patient_id}.stage1.mutation_annotated.tree so downstream PNG rendering can continue."
          echo "Inspect ${patient_id}.geneannotator.log to confirm whether GeneAnnotator found branch mutation labels."
        } > "\$PWD/delsieve_stage1_gene_tree/${patient_id}.geneannotator.warning.txt"
      fi
    fi

    cp "${delsieve_stage1_variants}/${patient_id}.stage1.mcc.tree" "\$PWD/delsieve_stage1_gene_tree/${patient_id}.stage1.mcc.tree"
    cp "${delsieve_stage1_variants}/candidate_sites.tsv" "\$PWD/delsieve_stage1_gene_tree/candidate_sites.tsv"
    cp "${delsieve_stage1_variants}/cell_names.tsv" "\$PWD/delsieve_stage1_gene_tree/cell_names.tsv"
    cp "\$mutation_map" "\$PWD/delsieve_stage1_gene_tree/mutation_map.tsv"
    cp "\$details_input"/*.loci_info "\$PWD/delsieve_stage1_gene_tree/" 2>/dev/null || true
    cp "\$details_input"/*.ternary "\$PWD/delsieve_stage1_gene_tree/" 2>/dev/null || true
    """
}

process DELSIEVE_GENE_TREE_PNG_STAGE1 {
    tag "$patient_id"
    conda "envs/delsieve_tree_plot.yml"
    publishDir { "${params.outdir}/${patient_id}/delsieve_stage1_gene_tree" }, mode: 'copy', pattern: "delsieve_stage1_gene_tree_png/*"

    input:
    tuple val(patient_id), path(gene_tree_files)

    output:
    tuple val(patient_id), path("delsieve_stage1_gene_tree_png/*")

    script:
    """
    set -euo pipefail

    unset DISPLAY
    export JAVA_TOOL_OPTIONS="\${JAVA_TOOL_OPTIONS:-} -Djava.awt.headless=true"

    mkdir -p delsieve_stage1_gene_tree_png

    tree_file=\$(find -L . -maxdepth 4 -type f -name "${patient_id}.stage1.mutation_annotated.tree" | head -n 1)
    if [[ -z "\$tree_file" ]]; then
      echo "No DelSIEVE mutation-annotated tree found for PNG rendering" >&2
      find -L . -maxdepth 4 -type f -print >&2
      exit 1
    fi

    plot_tree="\$PWD/delsieve_stage1_gene_tree_png/${patient_id}.stage1.mutation_annotated.plot_clipped.tree"
    python "${projectDir}/bin/clip_tree_branch_lengths.py" \
      --input "\$tree_file" \
      --output "\$plot_tree" \
      --summary "\$PWD/delsieve_stage1_gene_tree_png/${patient_id}.stage1.mutation_annotated.branch_length_clipping.tsv"

    "${params.delsieve_figtree}" -graphic PNG ${params.delsieve_figtree_extra_args} \\
      "\$plot_tree" \\
      "\$PWD/delsieve_stage1_gene_tree_png/${patient_id}.stage1.mutation_annotated.figtree.png" \\
      > "\$PWD/delsieve_stage1_gene_tree_png/${patient_id}.stage1.mutation_annotated.figtree.log" 2>&1

    cp "\$tree_file" "\$PWD/delsieve_stage1_gene_tree_png/${patient_id}.stage1.mutation_annotated.tree"
    """
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
    echo "monovar_region\t${params.monovar_region}" >> "\$report"
    echo "monovar_targets_bed\t${params.monovar_targets_bed}" >> "\$report"
    echo >> "\$report"

    test -n "${patient_id}" || { echo "ERROR: missing patient_id" >&2; exit 1; }
    test -s "${ref_fasta}" || { echo "ERROR: ref_fasta not found: ${ref_fasta}" >&2; exit 1; }
    if [[ "${params.run_monovar}" == "true" ]]; then
      test -s "${monovar_bam_list}" || { echo "ERROR: monovar_bam_list not found: ${monovar_bam_list}" >&2; exit 1; }
    fi

    if [[ "${params.run_monovar}" == "true" || "${params.run_sccaller}" == "true" || "${params.run_bam_qc}" == "true" ]]; then
      test -s "${cell_metadata}" || { echo "ERROR: cell_metadata not found: ${cell_metadata}" >&2; exit 1; }
      awk -F '\t' -v patient_id="${patient_id}" '
        NR == 1 { next }
        NF < 4 || \$1 == "" || \$2 == "" || \$3 == "" || \$4 == "" { print "ERROR: invalid cell_metadata row " NR > "/dev/stderr"; exit 1 }
        \$1 != patient_id { print "ERROR: cell_metadata patient_id mismatch on row " NR ": " \$1 > "/dev/stderr"; exit 1 }
        { count++ }
        END { if (count < 1) { print "ERROR: no cells found in cell_metadata" > "/dev/stderr"; exit 1 } }
      ' "${cell_metadata}"
    fi

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
      monovar_wbc_subset)
        test -s "${leukocyte_bam}" || { echo "ERROR: leukocyte_bam not found: ${leukocyte_bam}" >&2; exit 1; }
        ;;
      combined)
        test -s "${germline_vcf}" || { echo "ERROR: germline_vcf not found for combined mode: ${germline_vcf}" >&2; exit 1; }
        test -s "${leukocyte_bam}" || test -s "${leukocyte_vcf}" || { echo "ERROR: combined mode requires leukocyte_bam or leukocyte_vcf" >&2; exit 1; }
        ;;
      no_germline)
        echo "No germline/background exclusion configured; MonoVar split-only or external comparison-only mode." >> "\$report"
        ;;
      *)
        echo "ERROR: unsupported germline_mode '${germline_mode}'" >&2
        echo "Allowed: precomputed_vcf, deepvariant_bulk_bam, leukocyte_vcf, joint_monovar_leukocyte, monovar_wbc_subset, combined, no_germline" >&2
        exit 1
        ;;
    esac

    if [[ "${params.run_monovar}" == "true" ]]; then
      test -n "${params.monovar_script}" || { echo "ERROR: --monovar_script is required when --run_monovar true" >&2; exit 1; }
      test -s "${params.monovar_script}" || { echo "ERROR: monovar_script not found: ${params.monovar_script}" >&2; exit 1; }
      if [[ -n "${params.monovar_targets_bed}" ]]; then
        test -s "${params.monovar_targets_bed}" || { echo "ERROR: monovar_targets_bed not found: ${params.monovar_targets_bed}" >&2; exit 1; }
        zcat -f "${params.monovar_targets_bed}" | awk 'BEGIN{ok=0} /^[[:space:]]*#/ {next} NF == 0 {next} NF < 3 { print "ERROR: malformed monovar_targets_bed row " NR ": " \$0 > "/dev/stderr"; exit 1 } {ok=1} END { if (!ok) { print "ERROR: monovar_targets_bed has no BED intervals" > "/dev/stderr"; exit 1 } }'
      fi
    fi

    if [[ "${params.run_sccaller}" == "true" ]]; then
      test -n "${params.sccaller_script}" || { echo "ERROR: --sccaller_script is required when --run_sccaller true" >&2; exit 1; }
      test -s "${params.sccaller_script}" || { echo "ERROR: sccaller_script not found: ${params.sccaller_script}" >&2; exit 1; }
      test -n "${params.sccaller_hsnp_vcf}" || { echo "ERROR: --sccaller_hsnp_vcf is required when --run_sccaller true" >&2; exit 1; }
      test -s "${params.sccaller_hsnp_vcf}" || { echo "ERROR: sccaller_hsnp_vcf not found: ${params.sccaller_hsnp_vcf}" >&2; exit 1; }
      test -s "${bulk_bam}" || { echo "ERROR: bulk_bam is required when --run_sccaller true" >&2; exit 1; }
    fi

    echo "status\tOK" >> "\$report"
    """
}

process PREPARE_MONOVAR_WBC_SUBSET {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/monovar_wbc_subset/config" }, mode: 'copy'

    input:
    tuple val(patient_id), val(ref_fasta), val(monovar_bam_list), val(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_bam), val(leukocyte_vcf), path(input_check)

    output:
    tuple val(patient_id), val(ref_fasta), path("${patient_id}.monovar_wbc_subset_bams.txt"), path("${patient_id}.monovar_wbc_subset_cells.tsv"), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_bam), val(leukocyte_vcf), path("${patient_id}.monovar_wbc_subset_selection.tsv")

    script:
    """
    set -euo pipefail

    python - <<'PY'
import csv
import random
from pathlib import Path

patient_id = "${patient_id}"
cell_metadata = Path("${cell_metadata}")
leukocyte_bam = "${leukocyte_bam}"
ctc_count = int("${params.monovar_background_ctc_count}")
seed = "${params.monovar_background_seed}"

with cell_metadata.open(newline="") as handle:
    rows = [row for row in csv.DictReader(handle, delimiter="\\t") if row.get("patient_id") == patient_id]

def is_wbc(row):
    cell_id = row.get("cell_id", "").lower()
    cell_type = row.get("cell_type", "").lower()
    return cell_id in {"wbc", "leukocyte", "leuko"} or "leuk" in cell_id or cell_type in {"wbc", "leukocyte", "leuko"} or "leuk" in cell_type

ctcs = [row for row in rows if not is_wbc(row)]
wbcs = [row for row in rows if is_wbc(row)]
if leukocyte_bam and leukocyte_bam != "NA":
    matching = [row for row in wbcs if row.get("bam_path") == leukocyte_bam]
    if matching:
        wbcs = matching
if not wbcs:
    raise SystemExit("No WBC/leukocyte row found in cell_metadata for monovar_wbc_subset mode")
if not ctcs:
    raise SystemExit("No CTC rows found in cell_metadata for monovar_wbc_subset mode")

rng = random.Random(f"{seed}:{patient_id}")
selected_ctcs = sorted(rng.sample(ctcs, min(ctc_count, len(ctcs))), key=lambda row: rows.index(row))
selected = selected_ctcs + [wbcs[0]]

with open(f"{patient_id}.monovar_wbc_subset_bams.txt", "w", newline="") as out:
    for row in selected:
        out.write(row["bam_path"] + "\\n")

fieldnames = ["patient_id", "cell_id", "bam_path", "cell_type"]
with open(f"{patient_id}.monovar_wbc_subset_cells.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\\t", fieldnames=fieldnames, lineterminator="\\n")
    writer.writeheader()
    for row in selected:
        writer.writerow({key: row.get(key, "") for key in fieldnames})

with open(f"{patient_id}.monovar_wbc_subset_selection.tsv", "w", newline="") as out:
    writer = csv.DictWriter(out, delimiter="\\t", fieldnames=fieldnames + ["selection_role"], lineterminator="\\n")
    writer.writeheader()
    for row in selected_ctcs:
        payload = {key: row.get(key, "") for key in fieldnames}
        payload["selection_role"] = "random_ctc_background_context"
        writer.writerow(payload)
    payload = {key: wbcs[0].get(key, "") for key in fieldnames}
    payload["selection_role"] = "wbc_background"
    writer.writerow(payload)
PY
    """
}

process RUN_MONOVAR_WBC_SUBSET {
    tag "$patient_id"
    cpus params.monovar_threads
    conda "envs/monovar_py2.yml"
    publishDir { "${params.outdir}/${patient_id}/raw_calls/monovar_wbc_subset" }, mode: 'copy', pattern: "*.vcf"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.log"

    input:
    tuple val(patient_id), val(ref_fasta), path(monovar_bam_list), path(cell_metadata), val(germline_mode), val(germline_vcf), val(bulk_bam), val(leukocyte_bam), val(leukocyte_vcf), path(selection_tsv)

    output:
    tuple val(patient_id), path("${patient_id}.monovar_wbc_subset.vcf"), path(cell_metadata), val(ref_fasta), path("${patient_id}.monovar_wbc_subset.log")

    script:
    """
    set -euo pipefail

    test -s "${params.monovar_script}" || { echo "ERROR: monovar_script not found: ${params.monovar_script}" >&2; exit 1; }
    test -s "${ref_fasta}" || { echo "ERROR: ref_fasta not found: ${ref_fasta}" >&2; exit 1; }
    command -v samtools >/dev/null 2>&1 || { echo "ERROR: samtools not found in PATH" >&2; exit 1; }
    python --version > "${patient_id}.monovar_wbc_subset.log" 2>&1 || true
    samtools --version | head -n 2 >> "${patient_id}.monovar_wbc_subset.log" 2>&1 || true
    echo "Starting MonoVar WBC subset background run for ${patient_id}" >> "${patient_id}.monovar_wbc_subset.log"
    echo "BAM list: ${monovar_bam_list}" >> "${patient_id}.monovar_wbc_subset.log"
    echo "Selection TSV: ${selection_tsv}" >> "${patient_id}.monovar_wbc_subset.log"
    echo "Reference: ${ref_fasta}" >> "${patient_id}.monovar_wbc_subset.log"
    echo "MonoVar: ${params.monovar_script}" >> "${patient_id}.monovar_wbc_subset.log"
    echo "Region: ${params.monovar_region}" >> "${patient_id}.monovar_wbc_subset.log"
    echo "Targets BED: ${params.monovar_targets_bed}" >> "${patient_id}.monovar_wbc_subset.log"
    echo >> "${patient_id}.monovar_wbc_subset.log"

    region_args=()
    if [[ -n "${params.monovar_region}" ]]; then
      region_args=(-r "${params.monovar_region}")
    fi
    target_args=()
    if [[ -n "${params.monovar_targets_bed}" ]]; then
      test -s "${params.monovar_targets_bed}" || { echo "ERROR: monovar_targets_bed not found: ${params.monovar_targets_bed}" >&2; exit 1; }
      zcat -f "${params.monovar_targets_bed}" > "\$PWD/monovar_targets.bed"
      test -s "\$PWD/monovar_targets.bed" || { echo "ERROR: monovar_targets_bed produced no intervals: ${params.monovar_targets_bed}" >&2; exit 1; }
      target_args=(-l "\$PWD/monovar_targets.bed")
    fi

    samtools mpileup \
      -BQ0 \
      -d10000 \
      -f "${ref_fasta}" \
      -q 40 \
      "\${target_args[@]}" \
      -b "${monovar_bam_list}" \
      "\${region_args[@]}" \
    | python "${params.monovar_script}" \
      -p 0.002 \
      -a 0.2 \
      -t 0.05 \
      -m ${params.monovar_threads} \
      -f "${ref_fasta}" \
      -b "${monovar_bam_list}" \
      -o "${patient_id}.monovar_wbc_subset.vcf" \
    >> "${patient_id}.monovar_wbc_subset.log" 2>&1
    """
}

process SPLIT_MONOVAR_WBC_SUBSET {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/split_calls/monovar_wbc_subset" }, mode: 'copy', pattern: "*.monovar.split.vcf"
    publishDir { "${params.outdir}/${patient_id}/split_calls/monovar_wbc_subset" }, mode: 'copy', pattern: "*.split_sample_map.tsv"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.split.log"

    input:
    tuple val(patient_id), path(monovar_vcf), path(cell_metadata), val(ref_fasta), path(monovar_log)

    output:
    tuple val(patient_id), path("*.monovar.split.vcf"), path("${patient_id}.monovar_wbc_subset.split_sample_map.tsv"), path("${patient_id}.monovar_wbc_subset.split.log"), val(ref_fasta)

    script:
    """
    set -euo pipefail

    python "${projectDir}/bin/split_monovar_vcf.py" \
      --input "${monovar_vcf}" \
      --cell-metadata "${cell_metadata}" \
      --patient-id "${patient_id}" \
      --outdir . \
    > "${patient_id}.monovar_wbc_subset.split.log" 2>&1
    mv "${patient_id}.monovar.split_sample_map.tsv" "${patient_id}.monovar_wbc_subset.split_sample_map.tsv"
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
    echo "Targets BED: ${params.monovar_targets_bed}" >> "${patient_id}.monovar.log"
    echo >> "${patient_id}.monovar.log"

    region_args=()
    if [[ -n "${params.monovar_region}" ]]; then
      region_args=(-r "${params.monovar_region}")
    fi
    target_args=()
    if [[ -n "${params.monovar_targets_bed}" ]]; then
      test -s "${params.monovar_targets_bed}" || { echo "ERROR: monovar_targets_bed not found: ${params.monovar_targets_bed}" >&2; exit 1; }
      zcat -f "${params.monovar_targets_bed}" > "\$PWD/monovar_targets.bed"
      test -s "\$PWD/monovar_targets.bed" || { echo "ERROR: monovar_targets_bed produced no intervals: ${params.monovar_targets_bed}" >&2; exit 1; }
      target_args=(-l "\$PWD/monovar_targets.bed")
    fi

    samtools mpileup \
      -BQ0 \
      -d10000 \
      -f "${ref_fasta}" \
      -q 40 \
      "\${target_args[@]}" \
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

    if [[ ! -s "${ref_fasta}.fai" ]]; then
      samtools faidx "${ref_fasta}"
    fi
    bcftools reheader -f "${ref_fasta}.fai" \
      -o "${patient_id}.precomputed_germline.with_contigs.vcf" \
      "${germline_vcf}"

    bcftools view \
      -v snps \
      -f PASS \
      "${patient_id}.precomputed_germline.with_contigs.vcf" \
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

    # NOTE: this no longer drops low-depth/low-alt leukocyte calls -- it tags
    # them LOWCONF instead and keeps them. bcftools isec (used downstream in
    # FILTER_MONOVAR_AND_SUBTRACT_GERMLINE) matches on CHROM/POS/REF/ALT only
    # and does not read the FILTER column, so every non-ref leukocyte call
    # (PASS or LOWCONF) now counts as excludable germline evidence -- not
    # only the ones that independently clear the same depth bar as somatic
    # calling. This is intentional: it closes the germline-leakage gap where
    # a real germline SNP in shallow WBC coverage escaped exclusion entirely.
    awk -v min_dp="${params.min_total_depth}" -v min_alt="${params.min_alt_reads}" -v min_vaf="${params.min_vaf}" '
      BEGIN { FS=OFS="\t" }
      /^#CHROM/ { print "##FILTER=<ID=LOWCONF,Description=\"Non-reference call below the depth/alt-read/VAF confidence threshold; kept for visibility, not dropped\">" }
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
        if ((gt != "0/0") && (gt != "0|0") && (gt != "./.") && (gt != ".|.")) {
          \$7 = (dp >= min_dp && alt >= min_alt && vaf >= min_vaf) ? "PASS" : "LOWCONF"
          print
        }
      }
    ' "${split_vcf}" > "\$filtered"

    echo "Input leukocyte split VCF: ${split_vcf}" > "\$log"
    echo "min_total_depth=${params.min_total_depth}" >> "\$log"
    echo "min_alt_reads=${params.min_alt_reads}" >> "\$log"
    echo "min_vaf=${params.min_vaf}" >> "\$log"
    echo -n "Leukocyte exclusion variants before normalization: " >> "\$log"
    grep -vc '^#' "\$filtered" >> "\$log"
    echo -n "Leukocyte exclusion variants tagged LOWCONF (still used for exclusion): " >> "\$log"
    awk -F'\t' '\$7=="LOWCONF"' "\$filtered" | wc -l >> "\$log"

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

process MERGE_MONOVAR_GERMLINE_EXCLUSION {
    tag "$patient_id"
    conda "envs/bcftools.yml"
    publishDir { "${params.outdir}/${patient_id}/germline_exclusion" }, mode: 'copy', pattern: "*.vcf.gz*"
    publishDir { "${params.outdir}/${patient_id}/logs" }, mode: 'copy', pattern: "*.combined_germline_exclusion.log"

    input:
    tuple val(patient_id), path(bulk_vcf), path(bulk_tbi), path(leuko_vcf), path(leuko_tbi), val(ref_fasta)

    output:
    tuple val(patient_id), path("${patient_id}.combined_bulk_leukocyte.vcf.gz"), path("${patient_id}.combined_bulk_leukocyte.vcf.gz.tbi"), val(ref_fasta), val("combined_bulk_leukocyte")

    script:
    """
    set -euo pipefail

    log="${patient_id}.combined_germline_exclusion.log"
    echo "Bulk exclusion source: ${bulk_vcf}" > "\$log"
    echo "Leukocyte exclusion source: ${leuko_vcf}" >> "\$log"

    # bulk_vcf and leuko_vcf come from different single-sample callsets (bulk
    # DeepVariant vs. leukocyte MonoVar split), so they don't share sample
    # columns, and both are already normalized upstream (PREPARE_PRECOMPUTED_GERMLINE
    # / PREPARE_MONOVAR_LEUKOCYTE_EXCLUSION each run bcftools norm -f ref -m -any).
    # Same pattern as STANDARDIZE_EXTERNAL_CALLERS' multi-source exclusion merge:
    # drop genotypes to a pure sites union, concat, sort, dedupe exact duplicates.
    bcftools view -G "${bulk_vcf}" -Oz -o bulk.sites.vcf.gz
    tabix -f -p vcf bulk.sites.vcf.gz
    bcftools view -G "${leuko_vcf}" -Oz -o leuko.sites.vcf.gz
    tabix -f -p vcf leuko.sites.vcf.gz

    bcftools concat -a bulk.sites.vcf.gz leuko.sites.vcf.gz -Ou \
    | bcftools sort -Ou \
    | bcftools norm -d exact -Oz -o "${patient_id}.combined_bulk_leukocyte.vcf.gz" \
      2>> "\$log"
    tabix -f -p vcf "${patient_id}.combined_bulk_leukocyte.vcf.gz"

    echo -n "Combined bulk+leukocyte exclusion variants: " >> "\$log"
    bcftools view -H "${patient_id}.combined_bulk_leukocyte.vcf.gz" | wc -l >> "\$log"
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

    # NOTE: this no longer drops calls below the depth/alt/VAF bar -- it
    # tags them LOWCONF and keeps them, so they stay visible in the long
    # table/matrices/reports instead of silently disappearing. Downstream
    # bcftools norm/isec and the VEP merge all pass FILTER through
    # unchanged, so PASS/LOWCONF survives to the final per-cell VCF and the
    # long table's FILTER column.
    awk -v min_dp="${params.min_total_depth}" -v min_alt="${params.min_alt_reads}" -v min_vaf="${params.min_vaf}" '
      BEGIN { FS=OFS="\t" }
      /^#CHROM/ { print "##FILTER=<ID=LOWCONF,Description=\"Non-reference call below the depth/alt-read/VAF confidence threshold; kept for visibility, not dropped\">" }
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
        if ((gt != "0/0") && (gt != "0|0") && (gt != "./.") && (gt != ".|.")) {
          \$7 = (dp >= min_dp && alt >= min_alt && vaf >= min_vaf) ? "PASS" : "LOWCONF"
          print
        }
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
    echo -n "Of which tagged LOWCONF (kept, low depth/alt): " >> "\$log"
    awk -F'\t' '\$7=="LOWCONF"' "\$filtered" | wc -l >> "\$log"

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

process MONOVAR_VARIANT_FLOWCHARTS {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/reports/variant_flowcharts" }, mode: 'copy'

    input:
    tuple val(patient_id), path(split_vcfs), path(norm_vcfs), path(comparable_vcfs), path(exclusion_vcfs)

    output:
    tuple val(patient_id), path("variant_flowcharts")

    script:
    """
    set -euo pipefail

    mkdir -p flow_project/config flow_project/bulk_exclusion flow_project/normalized_calls/monovar flow_project/comparable_calls/monovar variant_flowcharts

    samples_tsv="\$PWD/flow_project/config/samples.tsv"
    printf "sample_id\tcaller\traw_vcf\tsccaller_mode\n" > "\$samples_tsv"

    for raw_vcf in ${split_vcfs}; do
      cell_id=\$(basename "\$raw_vcf" .monovar.split.vcf)
      case "\$(echo "\$cell_id" | tr '[:upper:]' '[:lower:]')" in
        wbc|leukocyte|leuko|*leuk*) continue ;;
      esac
      printf "%s\tmonovar\t%s\t\n" "\$cell_id" "\$raw_vcf" >> "\$samples_tsv"
    done

    for norm_vcf in ${norm_vcfs}; do
      cell_id=\$(basename "\$norm_vcf" .monovar.filtered.norm.vcf.gz)
      cp "\$norm_vcf" "flow_project/normalized_calls/monovar/\${cell_id}.monovar.filtered.snvs.norm.vcf.gz"
      if [[ -f "\${norm_vcf}.tbi" ]]; then
        cp "\${norm_vcf}.tbi" "flow_project/normalized_calls/monovar/\${cell_id}.monovar.filtered.snvs.norm.vcf.gz.tbi"
      fi
    done

    for comparable_vcf in ${comparable_vcfs}; do
      cell_id=\$(basename "\$comparable_vcf" | sed -E 's/\\.monovar\\.no_.+\\.vcf\\.gz\$//')
      cp "\$comparable_vcf" "flow_project/comparable_calls/monovar/\${cell_id}.monovar.comparable.vcf.gz"
      if [[ -f "\${comparable_vcf}.tbi" ]]; then
        cp "\${comparable_vcf}.tbi" "flow_project/comparable_calls/monovar/\${cell_id}.monovar.comparable.vcf.gz.tbi"
      fi
    done

    first_exclusion=""
    for exclusion_vcf in ${exclusion_vcfs}; do
      first_exclusion="\$exclusion_vcf"
      break
    done
    if [[ -n "\$first_exclusion" ]]; then
      cp "\$first_exclusion" "flow_project/bulk_exclusion/${patient_id}.leukocyte.norm.vcf.gz"
      if [[ -f "\${first_exclusion}.tbi" ]]; then
        cp "\${first_exclusion}.tbi" "flow_project/bulk_exclusion/${patient_id}.leukocyte.norm.vcf.gz.tbi"
      fi
    fi

    python "${projectDir}/bin/mutation_matrix/07_variant_flowcharts.py" \
      --project "\$PWD/flow_project" \
      --samples-tsv "\$samples_tsv" \
      --output-dir "\$PWD/variant_flowcharts"
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

    out="${chunk.name}.vep.tsv"
    log="${chunk.name}.vep_chunk.log"

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
    tuple val(patient_id), path(final_vcfs), path(vep_tsvs)

    output:
    tuple val(patient_id), path("*.long.tsv"), path("*.binary_matrix.tsv"), path("*.altread_matrix.tsv"), path("*.refread_matrix.tsv"), path("*.summary.tsv")

    script:
    def vep_args = vep_tsvs.collect { "--vep-tsv ${it}" }.join(' ')
    """
    set -euo pipefail

    python "${projectDir}/bin/build_mutation_matrices.py" \
      --patient-id "${patient_id}" \
      --caller monovar \
      --out-prefix "${patient_id}.monovar.final" \
      ${vep_args} \
      ${final_vcfs}
    """
}


process SNV_HTML_REPORT {
    tag "$patient_id"
    conda "envs/snv_report.yml"
    publishDir { "${params.outdir}/${patient_id}/snv_report" }, mode: 'copy'

    input:
    tuple val(patient_id), path(long_table), path(binary_matrix), path(altread_matrix), path(refread_matrix), path(summary_table), path(variant_flowcharts), val(breadth_summary_path)

    output:
    tuple val(patient_id), path("${patient_id}.snv_report.html"), path("${patient_id}.monovar.maf")

    script:
    """
    set -euo pipefail

    mkdir -p report_inputs
    cp "${long_table}" report_inputs/long.tsv
    cp "${binary_matrix}" report_inputs/binary_matrix.tsv
    cp "${altread_matrix}" report_inputs/altread_matrix.tsv
    cp "${refread_matrix}" report_inputs/refread_matrix.tsv
    cp "${summary_table}" report_inputs/summary.tsv

    Rscript -e 'rmarkdown::render(
      "${projectDir}/assets/snv_report.Rmd",
      output_file="${patient_id}.snv_report.html",
      output_dir=".",
      knit_root_dir=getwd(),
      params=list(
        patient_id="${patient_id}",
        long_table=file.path(Sys.getenv("PWD"), "report_inputs", "long.tsv"),
        binary_matrix=file.path(Sys.getenv("PWD"), "report_inputs", "binary_matrix.tsv"),
        altread_matrix=file.path(Sys.getenv("PWD"), "report_inputs", "altread_matrix.tsv"),
        refread_matrix=file.path(Sys.getenv("PWD"), "report_inputs", "refread_matrix.tsv"),
        summary_table=file.path(Sys.getenv("PWD"), "report_inputs", "summary.tsv"),
        flowcharts_dir=normalizePath(file.path(Sys.getenv("PWD"), "${variant_flowcharts}"), mustWork=TRUE),
        breadth_summary="${breadth_summary_path}",
        out_prefix="${patient_id}.monovar"
      )
    )'
    """
}

process SUMMARIZE_BREADTH_QC {
    tag "$patient_id"
    conda "envs/python_reporting.yml"
    publishDir { "${params.outdir}/${patient_id}/reports" }, mode: 'copy', pattern: "*.tsv"

    input:
    tuple val(patient_id), path(mosdepth_files)

    output:
    tuple val(patient_id), path("*.breadth_summary.tsv")

    script:
    """
    set -euo pipefail

    python "${projectDir}/bin/summarize_breadth_qc.py" \
      --patient-id "${patient_id}" \
      --out-prefix "${patient_id}.monovar" \
      --thresholds "${params.breadth_thresholds}" \
      ${mosdepth_files}
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
      --long-table "${long_table}" \
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
    tuple val(patient_id), path(filter_settings), path(prefilter_qc), path(filter_impact), val(breadth_summary_path)

    output:
    tuple val(patient_id), path("*_mqc.json")

    script:
    def breadth_arg = breadth_summary_path ? "--breadth-summary \"${breadth_summary_path}\"" : ""
    """
    set -euo pipefail

    python "${projectDir}/bin/make_multiqc_custom_content.py" \
      --patient-id "${patient_id}" \
      --settings "${filter_settings}" \
      --prefilter-qc "${prefilter_qc}" \
      --filter-impact "${filter_impact}" \
      ${breadth_arg} \
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
