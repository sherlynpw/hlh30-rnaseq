process RUVSEQ {
    tag "ruvseq"
    label 'process_medium'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/ruvseq", mode: params.publish_dir_mode

    input:
    path count_matrix
    path sample_meta
    path anchor_groups   // file or 'NO_ANCHOR_GROUPS' sentinel

    output:
    path "w_factors.tsv",        emit: w_factors
    path "ruv_diagnostics.tsv",  emit: diagnostics
    path "ruv_diagnostics.pdf",  emit: plots

    script:
    def anchor_arg = anchor_groups.name != 'NO_ANCHOR_GROUPS' ? "--anchor_groups ${anchor_groups}" : ''
    """
    ruvseq.R \\
        --counts ${count_matrix} \\
        --metadata ${sample_meta} \\
        --method ${params.ruvseq.method} \\
        --k ${params.ruvseq.k} \\
        --condition_from ${params.deseq2.condition_from} \\
        ${anchor_arg} \\
        --outdir .
    """
}
