process GO_ENRICHMENT {
    tag "go_${source_label}"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/enrichment/go/${source_label}", mode: params.publish_dir_mode

    input:
    path input_data
    path gene_db
    path universe      // background gene list (enrichment universe)
    val source_label
    val filter_by      // 'padj', 'pvalue', or '' (for gene list .txt files)
    val threshold

    output:
    path "**/*.csv",  emit: tables,  optional: true
    path "**/*.pdf",  emit: plots,   optional: true
    // Emit TOP-LEVEL entries only (the per-gene-list subdirs go_enrichment.R
    // produces). Recursive ** would expand every CSV/PDF as a separate flat
    // path, colliding on duplicate basenames (multiple gene lists each have
    // go_BP.csv) when staged into downstream processes like REVIGO.
    path "*",         emit: results, optional: true

    when:
    params.enrichment.go.enabled

    script:
    def db_arg = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    def filter_arg = filter_by ? "--filter_by ${filter_by} --threshold ${threshold}" : ''
    def lfc_arg = params.enrichment.go.lfc ? "--lfc ${params.enrichment.go.lfc}" : ''
    """
    go_enrichment.R \\
        --input . \\
        --org_db ${params.organism.org_db} \\
        --universe ${universe} \\
        --cpus ${task.cpus} \\
        ${filter_arg} \\
        ${lfc_arg} \\
        ${db_arg} \\
        --outdir .
    """
}
