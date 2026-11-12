process FEATURECOUNTS {
    tag "$meta.id"
    label 'process_low'

    container 'rayandrew/celegans-r-analysis:0.11.0'

    input:
    tuple val(meta), path(bam)
    path gtf

    output:
    tuple val(meta), path("${meta.id}.featureCounts.txt"), emit: counts
    tuple val(meta), path("${meta.id}.featureCounts.txt.summary"), emit: summary
    val(meta), emit: meta

    script:
    def strandedness = meta.strandedness ?: 'auto'
    def strand_val = strandedness == 'forward' ? '1' :
                     strandedness == 'reverse' ? '2' : '0'
    def paired_val = meta.input_type == 'fastq_standard' ? 'TRUE' : 'FALSE'

    """
    featurecounts.R ${bam} ${gtf} ${meta.id} ${strand_val} ${paired_val} ${task.cpus}
    """
}
