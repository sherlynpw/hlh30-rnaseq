process FUZZY_CLUSTERING {
    tag "fuzzy_clustering"
    label 'process_medium'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/clustering/fuzzy", mode: params.publish_dir_mode

    input:
    path vst_counts
    path filtered
    path sample_meta
    path gene_db

    output:
    path "cluster_assignments.csv",  emit: assignments
    path "core_genes.csv",           emit: core_genes
    path "plots/",                   emit: plots
    path "*",                        emit: results

    script:
    def db_arg = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    def order_arg = params.clustering.condition_order ? "--condition_order ${params.clustering.condition_order}" : ''
    """
    fuzzy_clustering.R \\
        --vst ${vst_counts} \\
        --genes ${filtered} \\
        --metadata ${sample_meta} \\
        --num_clusters ${params.clustering.num_clusters} \\
        --sort_by ${params.clustering.sort_by} \\
        --condition_from ${params.deseq2.condition_from} \\
        ${order_arg} \\
        ${db_arg} \\
        --outdir .
    """
}
