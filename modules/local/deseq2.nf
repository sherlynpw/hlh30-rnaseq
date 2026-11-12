process DESEQ2 {
    tag { group == 'all' ? 'deseq2' : "deseq2:${group}" }
    label 'process_medium'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/deseq2${group == 'all' ? '' : '/' + group}", mode: params.publish_dir_mode

    input:
    path count_matrix
    path sample_meta
    path gene_db
    path ruv_factors   // file or 'NO_RUV' sentinel
    tuple val(group), val(reference)   // group name ('all' = ungrouped) + its reference level

    output:
    path "results/*.csv",                  emit: results
    path "contrast_results/*.csv",         emit: contrast_results, optional: true
    path "dds.rds",                        emit: dds_rds
    path "vst_counts.csv",                 emit: vst_counts
    path "vst_counts_raw.csv",             emit: vst_counts_raw
    path "normalized_counts.csv",          emit: normalized_counts
    path "size_factors.csv",               emit: size_factors

    script:
    def batch_arg     = params.deseq2.batch ? "--batch ${params.deseq2.batch}" : ''
    def db_arg        = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    def ruv_arg       = ruv_factors.name != 'NO_RUV' ? "--ruv_factors ${ruv_factors}" : ''
    def correct_arg   = params.clustering?.correct_vst ? '--apply_batch_correction' : ''
    // Multi-factor mode: --design_formula overrides condition/pairwise logic.
    // When set, also pass through --reference_levels for per-factor refs.
    // Both must be quoted because the formula contains '+' and spaces.
    def design_arg    = params.deseq2.design_formula ? "--design_formula '${params.deseq2.design_formula}'" : ''
    def reflevels_arg = params.deseq2.reference_levels ? "--reference_levels '${params.deseq2.reference_levels}'" : ''
    // Compound contrasts: serialize the user's map to JSON the R script reads.
    def contrasts_map = params.deseq2.contrasts ?: [:]
    def contrasts_json = contrasts_map ? groovy.json.JsonOutput.toJson(contrasts_map) : ''
    def contrasts_setup = contrasts_json ? "cat > deseq2_contrasts.json <<'_EOF_'\n${contrasts_json}\n_EOF_" : ''
    def contrasts_arg = contrasts_json ? "--contrasts_json deseq2_contrasts.json" : ''
    // condition_from / reference are ignored when design_formula is set, but we
    // still pass them so optparse doesn't fail (they have safe defaults).
    """
    ${contrasts_setup}
    deseq2.R \\
        --counts ${count_matrix} \\
        --metadata ${sample_meta} \\
        --condition_from ${params.deseq2.condition_from} \\
        --reference ${reference} \\
        --padj ${params.deseq2.padj_threshold} \\
        --lfc ${params.deseq2.lfc_threshold} \\
        --shrinkage ${params.deseq2.shrinkage} \\
        ${batch_arg} \\
        ${ruv_arg} \\
        ${correct_arg} \\
        ${design_arg} \\
        ${reflevels_arg} \\
        ${contrasts_arg} \\
        ${db_arg} \\
        --outdir .
    """
}
