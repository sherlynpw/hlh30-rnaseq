process WORMCAT {
    tag "wormcat_${source_label}"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/enrichment/wormcat/${source_label}", mode: params.publish_dir_mode

    input:
    path input_data
    val source_label
    val filter_by
    val threshold

    output:
    path "**/*.csv", emit: tables, optional: true
    path "**/*.svg", emit: plots, optional: true
    path "**",       emit: results, optional: true

    when:
    params.enrichment.wormcat.enabled

    script:
    def filter_arg = filter_by ? "--filter_by ${filter_by} --threshold ${threshold}" : ''
    """
    wormcat.R \\
        --input . \\
        --cpus ${task.cpus} \\
        ${filter_arg} \\
        --outdir .
    """
}
