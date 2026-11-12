process UMAP {
    tag "umap"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/exploration/umap", mode: params.publish_dir_mode

    input:
    path vst_counts
    path vst_counts_raw
    path sample_meta

    output:
    path "umap_coordinates.csv",  emit: coordinates
    path "umap_*.pdf",            emit: plots
    path "*",                     emit: results

    script:
    def color_by = params.exploration.pca_color_by ?
        "--color_by ${params.exploration.pca_color_by.join(',')}" : ''
    def raw_arg = params.clustering?.correct_vst ? "--vst_raw ${vst_counts_raw}" : ''
    """
    umap.R \\
        --vst ${vst_counts} \\
        ${raw_arg} \\
        --metadata ${sample_meta} \\
        ${color_by} \\
        --outdir .
    """
}
