process CV_DENSITY {
    tag "cv_density"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/exploration/cv_density", mode: params.publish_dir_mode

    input:
    path vst_counts

    output:
    path "cv_density_plot.pdf",  emit: plots
    path "*",                    emit: results

    script:
    """
    cv_density.R \\
        --counts ${vst_counts} \\
        --outdir .
    """
}
