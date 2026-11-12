process KMEANS_CLUSTERING {
    tag "kmeans"
    label 'process_medium'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/clustering/kmeans", mode: params.publish_dir_mode

    input:
    path vst_counts
    path filtered
    path sample_meta
    path gene_db

    output:
    path "cluster_assignments.csv",  emit: assignments
    path "cluster_summary.csv",      emit: summary
    path "plots/",                   emit: plots
    path "*",                        emit: results

    script:
    def db_arg = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    def order_arg = params.clustering.condition_order ? "--condition_order ${params.clustering.condition_order}" : ''
    def num_clusters = params.clustering.kmeans_num_clusters ?: params.clustering.num_clusters
    """
    kmeans_clustering.R \\
        --vst ${vst_counts} \\
        --genes ${filtered} \\
        --metadata ${sample_meta} \\
        --num_clusters ${num_clusters} \\
        --cluster_range ${params.clustering.cluster_range} \\
        --condition_from ${params.deseq2.condition_from} \\
        ${order_arg} \\
        ${db_arg} \\
        --outdir .
    """
}
