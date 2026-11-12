process REVIGO {
    tag "revigo_${source_label}"
    label 'process_medium'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/enrichment/revigo/${source_label}", mode: params.publish_dir_mode

    input:
    path go_results
    val source_label

    output:
    path "**/revigo_*_reduced.csv", emit: tables,  optional: true
    path "**/revigo_*.pdf",         emit: plots,   optional: true
    // Top-level entries (per-gene-list subdirs); mirrors GO_ENRICHMENT.
    path "*",                       emit: results, optional: true

    when:
    params.enrichment?.go?.enabled && (params.enrichment?.revigo?.enabled != false)

    script:
    def threshold  = params.enrichment?.revigo?.threshold ?: 0.7
    def ontologies = params.enrichment?.revigo?.ontologies ?
        params.enrichment.revigo.ontologies.join(',') : 'BP,MF,CC'
    def min_terms  = params.enrichment?.revigo?.min_terms ?: 3
    """
    revigo_reduce.R \\
        --input . \\
        --org_db ${params.organism.org_db} \\
        --threshold ${threshold} \\
        --ontologies ${ontologies} \\
        --min_terms ${min_terms} \\
        --cpus ${task.cpus} \\
        --outdir .
    """
}
