process GSEA {
    tag "gsea"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/enrichment/gsea", mode: params.publish_dir_mode

    input:
    path de_results
    path gene_db
    path manifest_json
    path collection_files
    path ortholog_map

    output:
    path "**/*.csv",  emit: tables, optional: true
    path "**/*.pdf",  emit: plots, optional: true
    path "**",        emit: results, optional: true

    when:
    params.enrichment.gsea.enabled

    script:
    def db_arg = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    def builtin_list = (params.enrichment.gsea.builtin ?: []).join(',')
    def builtin_arg = builtin_list ? "--builtin ${builtin_list}" : ''
    def manifest_arg = manifest_json.name != 'NO_COLLECTIONS_MANIFEST' ? "--collections_manifest ${manifest_json}" : ''
    def min_size_arg = params.enrichment.gsea.min_size ? "--min_size ${params.enrichment.gsea.min_size}" : ''
    def max_size_arg = params.enrichment.gsea.max_size ? "--max_size ${params.enrichment.gsea.max_size}" : ''
    def ortholog_map_arg = ortholog_map.name != 'NO_ORTHOLOG_MAP' ? "--ortholog_map ${ortholog_map}" : ''
    def ortholog_from_arg = params.enrichment.gsea.ortholog_from ? "--ortholog_from ${params.enrichment.gsea.ortholog_from}" : ''
    def ortholog_to_arg = params.enrichment.gsea.ortholog_to ? "--ortholog_to ${params.enrichment.gsea.ortholog_to}" : ''
    """
    gsea.R \\
        --de_results_dir . \\
        --collections_dir . \\
        --ranking ${params.enrichment.gsea.ranking} \\
        --org_db ${params.organism.org_db} \\
        --kegg_org ${params.organism.kegg_org ?: 'cel'} \\
        ${builtin_arg} \\
        ${manifest_arg} \\
        ${min_size_arg} \\
        ${max_size_arg} \\
        ${ortholog_map_arg} \\
        ${ortholog_from_arg} \\
        ${ortholog_to_arg} \\
        ${db_arg} \\
        --outdir .
    """
}
