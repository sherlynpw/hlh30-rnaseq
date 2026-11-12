process BACKGROUND_GENES {
    tag "background_genes"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    publishDir "${params.outdir}/${params.id}/enrichment", mode: params.publish_dir_mode

    input:
    path de_results   // DESeq2 result CSVs (all groups pooled)

    output:
    path "background.txt", emit: genes

    script:
    """
    background_genes.R --input_dir . --outdir .
    """
}
