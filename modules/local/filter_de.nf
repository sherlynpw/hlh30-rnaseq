process FILTER_DE {
    tag "filter_de"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/filtered", mode: params.publish_dir_mode

    input:
    path de_results

    output:
    path "filtered/",  emit: filtered

    script:
    """
    filter_de.R \\
        --input_dir . \\
        --padj ${params.deseq2.padj_threshold} \\
        --lfc ${params.deseq2.lfc_threshold} \\
        --outdir filtered
    """
}
