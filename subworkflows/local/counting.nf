include { FEATURECOUNTS } from '../../modules/local/featurecounts'
include { MERGE_COUNTS  } from '../../modules/local/merge_counts'

workflow COUNTING {

    take:
    ch_bams       // [meta, bam] - all BAMs (Plasmidsaurus + nf-core/rnaseq output)
    ch_gtf        // path to Ensembl 115 GTF (channel)
    samplesheet   // path to original YAML samplesheet

    main:
    // Run featureCounts on each BAM with consistent GTF
    FEATURECOUNTS(ch_bams, ch_gtf)

    // Collect all per-sample counts and merge into matrix
    MERGE_COUNTS(
        FEATURECOUNTS.out.counts.map { meta, counts -> counts }.collect(),
        file(samplesheet)
    )

    emit:
    count_matrix = MERGE_COUNTS.out.count_matrix   // gene_counts.tsv
    sample_meta  = MERGE_COUNTS.out.sample_meta     // sample_meta.tsv
    per_sample   = FEATURECOUNTS.out.counts          // [meta, counts] per sample
    summaries    = FEATURECOUNTS.out.summary          // [meta, summary] per sample
}
