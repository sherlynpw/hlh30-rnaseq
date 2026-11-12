process MERGE_COUNTS {
    tag "merge_counts"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/counts", mode: params.publish_dir_mode

    input:
    path count_files  // collected featureCounts output files
    path sample_meta  // sample metadata YAML

    output:
    path "gene_counts.tsv",   emit: count_matrix
    path "sample_meta.tsv",   emit: sample_meta

    script:
    """
    merge_counts.R \\
        --counts_dir . \\
        --metadata ${sample_meta} \\
        --outdir .
    """
}
