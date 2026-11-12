include { FUZZY_CLUSTERING        } from '../../modules/local/fuzzy_clustering'
include { HDBSCAN_CLUSTERING      } from '../../modules/local/hdbscan_clustering'
include { KMEANS_CLUSTERING       } from '../../modules/local/kmeans_clustering'
include { HIERARCHICAL_CLUSTERING } from '../../modules/local/hierarchical_clustering'

workflow CLUSTERING {

    take:
    ch_vst         // VST-normalized counts
    ch_de_results  // DESeq2 CSVs (not used here)
    ch_filtered    // filtered DE gene lists
    ch_meta        // sample metadata
    ch_gene_db     // gene_db.sqlite

    main:
    FUZZY_CLUSTERING(ch_vst, ch_filtered, ch_meta, ch_gene_db)
    HDBSCAN_CLUSTERING(ch_vst, ch_filtered, ch_meta, ch_gene_db)
    KMEANS_CLUSTERING(ch_vst, ch_filtered, ch_meta, ch_gene_db)
    HIERARCHICAL_CLUSTERING(ch_vst, ch_filtered, ch_meta, ch_gene_db)

    ch_fuzzy   = FUZZY_CLUSTERING.out.results
    ch_hdbscan = HDBSCAN_CLUSTERING.out.results
    ch_kmeans  = KMEANS_CLUSTERING.out.results
    ch_hier    = HIERARCHICAL_CLUSTERING.out.results

    emit:
    fuzzy        = ch_fuzzy
    hdbscan      = ch_hdbscan
    kmeans       = ch_kmeans
    hierarchical = ch_hier
}
