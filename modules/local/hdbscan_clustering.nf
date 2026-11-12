process HDBSCAN_CLUSTERING {
    tag "hdbscan"
    label 'process_medium'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/clustering/hdbscan", mode: params.publish_dir_mode

    input:
    path vst_counts
    path filtered
    path sample_meta
    path gene_db

    output:
    path "cluster_assignments.csv",  emit: assignments
    path "plots/",                   emit: plots
    path "*",                        emit: results

    script:
    def db_arg = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    def order_arg = params.clustering.condition_order ? "--condition_order ${params.clustering.condition_order}" : ''
    def min_pts_arg = params.clustering.min_pts ? "--min_pts ${params.clustering.min_pts}" : ''
    """
    hdbscan_clustering.R \\
        --vst ${vst_counts} \\
        --genes ${filtered} \\
        --metadata ${sample_meta} \\
        --condition_from ${params.deseq2.condition_from} \\
        ${min_pts_arg} \\
        ${order_arg} \\
        ${db_arg} \\
        --outdir .
    """
}
