process VST_ALL {
    tag "vst_all"
    label 'process_medium'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/vst_all", mode: params.publish_dir_mode

    input:
    path count_matrix
    path sample_meta

    output:
    path "vst_counts.csv",     emit: vst_counts
    path "vst_counts_raw.csv", emit: vst_counts_raw

    script:
    """
    vst_all.R \\
        --counts ${count_matrix} \\
        --metadata ${sample_meta} \\
        --outdir .
    """
}
