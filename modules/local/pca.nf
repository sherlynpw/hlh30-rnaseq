process PCA {
    tag "pca"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/exploration/pca", mode: params.publish_dir_mode

    input:
    path vst_counts
    path vst_counts_raw
    path sample_meta

    output:
    path "pca_coordinates.csv",  emit: coordinates
    path "pca_*.pdf",            emit: plots
    path "*",                    emit: results

    script:
    def color_by = params.exploration.pca_color_by ?
        "--color_by ${params.exploration.pca_color_by.join(',')}" : ''
    // Only request before/after when correction was actually applied
    def raw_arg = params.clustering?.correct_vst ? "--vst_raw ${vst_counts_raw}" : ''
    """
    pca.R \\
        --vst ${vst_counts} \\
        ${raw_arg} \\
        --metadata ${sample_meta} \\
        ${color_by} \\
        --outdir .
    """
}
