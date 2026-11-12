process HIERARCHICAL_CLUSTERING {
    tag "hierarchical"
    label 'process_medium'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/clustering/hierarchical", mode: params.publish_dir_mode

    input:
    path vst_counts
    path filtered
    path sample_meta
    path gene_db

    output:
    path "cluster_assignments.csv",  emit: assignments
    path "cluster_summary.csv",      emit: summary
    path "hclust.rds",               emit: model
    path "plots/",                   emit: plots
    path "*",                        emit: results

    script:
    def db_arg = gene_db.name != 'NO_GENE_DB' ? "--gene_db ${gene_db}" : ''
    def order_arg = params.clustering.condition_order ? "--condition_order ${params.clustering.condition_order}" : ''
    def num_clusters = params.clustering.hierarchical_num_clusters ?: params.clustering.num_clusters
    def linkage = params.clustering.hierarchical_linkage ?: 'ward.D2'
    def distance = params.clustering.hierarchical_distance ?: 'euclidean'
    def max_plot = params.clustering.hierarchical_max_plot_genes ?: 2000
    def auto_method = params.clustering.hierarchical_auto_method ?: 'dynamic_tree_cut'
    def deep_split = params.clustering.hierarchical_deep_split ?: 2
    def min_size = params.clustering.hierarchical_min_cluster_size ?: 30
    """
    hierarchical_clustering.R \\
        --vst ${vst_counts} \\
        --genes ${filtered} \\
        --metadata ${sample_meta} \\
        --num_clusters ${num_clusters} \\
        --cluster_range ${params.clustering.cluster_range} \\
        --condition_from ${params.deseq2.condition_from} \\
        --linkage ${linkage} \\
        --distance ${distance} \\
        --auto_method ${auto_method} \\
        --deep_split ${deep_split} \\
        --min_cluster_size ${min_size} \\
        --max_plot_genes ${max_plot} \\
        ${order_arg} \\
        ${db_arg} \\
        --outdir .
    """
}
