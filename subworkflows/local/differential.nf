include { DESEQ2    } from '../../modules/local/deseq2'
include { FILTER_DE } from '../../modules/local/filter_de'
include { RUVSEQ    } from '../../modules/local/ruvseq'

workflow DIFFERENTIAL {

    take:
    ch_counts   // count matrix
    ch_meta     // sample metadata
    ch_gene_db  // gene_db.sqlite

    main:
    // Optional RUVSeq - emits W factors when enabled, sentinel otherwise
    if (params.ruvseq?.enabled) {
        ch_anchor = params.ruvseq.anchor_groups
            ? Channel.value(file(params.ruvseq.anchor_groups))
            : Channel.value(file("${projectDir}/assets/NO_ANCHOR_GROUPS"))
        RUVSEQ(ch_counts, ch_meta, ch_anchor)
        ch_ruv = RUVSEQ.out.w_factors
    } else {
        ch_ruv = Channel.value(file("${projectDir}/assets/NO_RUV"))
    }

    // One DESeq2 run over all samples, tagged 'all'. The reference may be absent in multi-factor configs (design_formula); it
    // is ignored downstream in that mode, so a placeholder is fine.
    ch_groups = Channel.of(tuple('all', (params.deseq2.reference ?: 'NA').toString()))

    // Path inputs are value channels, reused across every group
    DESEQ2(ch_counts, ch_meta, ch_gene_db, ch_ruv, ch_groups)
    FILTER_DE(DESEQ2.out.results)

    emit:
    de_results        = DESEQ2.out.results.flatten().collect()  // all groups pooled
    contrast_results  = DESEQ2.out.contrast_results   // optional compound contrast CSVs
    dds_rds           = DESEQ2.out.dds_rds            // fitted DESeqDataSet(s)
    vst_counts        = DESEQ2.out.vst_counts         // per-group when split
    vst_counts_raw    = DESEQ2.out.vst_counts_raw
    filtered          = FILTER_DE.out.filtered        // filter_de gene lists (padj + lfc)
    ruv               = ch_ruv                        // W factors path or sentinel
}
