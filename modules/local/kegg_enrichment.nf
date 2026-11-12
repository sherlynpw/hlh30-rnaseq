process KEGG_ENRICHMENT {
    tag "kegg_${source_label}"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/enrichment/kegg/${source_label}", mode: params.publish_dir_mode

    input:
    path input_data
    path gene_db
    path universe      // background gene list (enrichment universe)
    val source_label
    val filter_by
    val threshold

    output:
    path "**/*.csv",  emit: tables, optional: true
    path "**/*.pdf",  emit: plots, optional: true
    path "**",        emit: results, optional: true

    when:
    params.enrichment.kegg.enabled

    script:
    def db_arg = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    def filter_arg = filter_by ? "--filter_by ${filter_by} --threshold ${threshold}" : ''
    def lfc_arg = params.enrichment.kegg.lfc ? "--lfc ${params.enrichment.kegg.lfc}" : ''
    """
    kegg_enrichment.R \\
        --input . \\
        --kegg_org ${params.organism.kegg_org} \\
        --universe ${universe} \\
        --cpus ${task.cpus} \\
        ${filter_arg} \\
        ${lfc_arg} \\
        ${db_arg} \\
        --outdir .
    """
}
