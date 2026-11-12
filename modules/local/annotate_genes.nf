process ANNOTATE_GENES {
    tag "annotate"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/annotation", mode: params.publish_dir_mode

    input:
    path filtered
    path gene_db

    output:
    path "annotated_genes.csv",  emit: results

    script:
    def db_arg = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    """
    annotate_genes.R \\
        --genes ${filtered} \\
        ${db_arg} \\
        --outdir .
    """
}
