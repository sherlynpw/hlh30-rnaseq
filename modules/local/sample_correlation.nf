process SAMPLE_CORRELATION {
    tag "sample_correlation"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/exploration/sample_correlation", mode: params.publish_dir_mode

    input:
    path vst_counts
    path sample_meta

    output:
    path "sample_correlation.csv",            emit: matrix
    path "sample_correlation.pdf",            emit: plot
    path "sample_correlation_summary.tsv",    emit: summary, optional: true
    path "*",                                 emit: results

    script:
    // top_n = 0 -> use all genes (recommended for cross-group similarity questions)
    def top_n    = params.exploration?.correlation_top_n   != null ? params.exploration.correlation_top_n : 0
    def sel      = params.exploration?.correlation_selection ?: 'expressed'
    def method   = params.exploration?.correlation_method  ?: 'spearman'
    def group_by = params.exploration?.correlation_group_by ?: 'condition'
    """
    sample_correlation.R \\
        --vst ${vst_counts} \\
        --metadata ${sample_meta} \\
        --top_n ${top_n} \\
        --selection ${sel} \\
        --method ${method} \\
        --group_by ${group_by} \\
        --outdir .
    """
}
