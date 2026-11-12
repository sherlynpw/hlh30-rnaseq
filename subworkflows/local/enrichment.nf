include { GO_ENRICHMENT as GO_DESEQ         } from '../../modules/local/go_enrichment'
include { GO_ENRICHMENT as GO_FILTERED       } from '../../modules/local/go_enrichment'
include { GO_ENRICHMENT as GO_FUZZY          } from '../../modules/local/go_enrichment'
include { GO_ENRICHMENT as GO_HDBSCAN        } from '../../modules/local/go_enrichment'
include { GO_ENRICHMENT as GO_KMEANS         } from '../../modules/local/go_enrichment'
include { GO_ENRICHMENT as GO_HIERARCHICAL   } from '../../modules/local/go_enrichment'
include { GO_ENRICHMENT as GO_COMPARE        } from '../../modules/local/go_enrichment'

include { REVIGO as REVIGO_DESEQ             } from '../../modules/local/revigo'
include { REVIGO as REVIGO_FILTERED          } from '../../modules/local/revigo'
include { REVIGO as REVIGO_FUZZY             } from '../../modules/local/revigo'
include { REVIGO as REVIGO_HDBSCAN           } from '../../modules/local/revigo'
include { REVIGO as REVIGO_KMEANS            } from '../../modules/local/revigo'
include { REVIGO as REVIGO_HIERARCHICAL      } from '../../modules/local/revigo'
include { REVIGO as REVIGO_COMPARE           } from '../../modules/local/revigo'

include { KEGG_ENRICHMENT as KEGG_DESEQ     } from '../../modules/local/kegg_enrichment'
include { KEGG_ENRICHMENT as KEGG_FILTERED   } from '../../modules/local/kegg_enrichment'
include { KEGG_ENRICHMENT as KEGG_FUZZY      } from '../../modules/local/kegg_enrichment'
include { KEGG_ENRICHMENT as KEGG_HDBSCAN    } from '../../modules/local/kegg_enrichment'
include { KEGG_ENRICHMENT as KEGG_KMEANS     } from '../../modules/local/kegg_enrichment'
include { KEGG_ENRICHMENT as KEGG_HIERARCHICAL } from '../../modules/local/kegg_enrichment'
include { KEGG_ENRICHMENT as KEGG_COMPARE    } from '../../modules/local/kegg_enrichment'

include { WORMCAT as WORMCAT_DESEQ           } from '../../modules/local/wormcat'
include { WORMCAT as WORMCAT_FILTERED        } from '../../modules/local/wormcat'
include { WORMCAT as WORMCAT_FUZZY           } from '../../modules/local/wormcat'
include { WORMCAT as WORMCAT_HDBSCAN         } from '../../modules/local/wormcat'
include { WORMCAT as WORMCAT_KMEANS          } from '../../modules/local/wormcat'
include { WORMCAT as WORMCAT_HIERARCHICAL    } from '../../modules/local/wormcat'
include { WORMCAT as WORMCAT_COMPARE         } from '../../modules/local/wormcat'

include { GSEA            } from '../../modules/local/gsea'
include { SSGSEA          } from '../../modules/local/ssgsea'
include { ANNOTATE_GENES  } from '../../modules/local/annotate_genes'
include { COMPARE_DE      } from '../../modules/local/compare_de'
include { BACKGROUND_GENES } from '../../modules/local/background_genes'

workflow ENRICHMENT {

    take:
    ch_de_results    // DESeq2 result CSVs (full, unfiltered)
    ch_filtered      // filter_de output (gene list .txt files)
    ch_fuzzy_genes   // fuzzy clustering output (gene list .txt files)
    ch_hdbscan_genes // HDBSCAN output (gene list .txt files)
    ch_kmeans_genes  // k-means output (gene list .txt files)
    ch_hierarchical_genes // hierarchical output (gene list .txt files)
    ch_vst_counts    // VST counts CSV (RUV/batch-corrected when available) — for ssGSEA
    ch_sample_meta   // sample metadata TSV — for ssGSEA
    ch_contrast_results // optional DESeq2 compound-contrast CSVs — for ssGSEA gene-driver plots
    ch_gene_db       // gene_db.sqlite

    main:

    ch_ortholog_map = Channel.value(file(params.enrichment.gsea.ortholog_map ?: 'NO_ORTHOLOG_MAP'))

    // GSEA collections - arbitrary list of external GMTs
    def gsea_collections = params.enrichment.gsea.collections ?: []
    def collection_files = gsea_collections.collect { file(it.file) }
    def manifest_entries = gsea_collections.collect { entry ->
        [
            name        : entry.name,
            file        : file(entry.file).name,
            format      : entry.format ?: 'gmt',
            map_ortholog: entry.map_ortholog ?: false
        ]
    }
    def manifest_file
    if (manifest_entries) {
        def tmp_dir = file("${projectDir}/.nextflow_tmp")
        tmp_dir.mkdirs()
        manifest_file = file("${tmp_dir}/gsea_collections_manifest.${params.id}.json")
        manifest_file.text = groovy.json.JsonOutput.toJson(manifest_entries)
    } else {
        manifest_file = file("${projectDir}/assets/NO_COLLECTIONS_MANIFEST")
    }
    ch_gsea_manifest = Channel.value(manifest_file)
    ch_gsea_collection_files = Channel.value(collection_files ?: [file("${projectDir}/assets/NO_COLLECTION_FILES")])

    // Compare DE (produces gene list .txt files) - optional
    def compare_enabled = params.enrichment.compare_de != false
    if (compare_enabled) {
        COMPARE_DE(ch_de_results)
    }

    // Enrichment universe = genes tested in this analysis, not the whole
    // annotated genome (omitting it inflates broad, spurious GO/KEGG terms).
    BACKGROUND_GENES(ch_de_results)
    ch_background = BACKGROUND_GENES.out.genes

    // --- GO enrichment ---
    // From DESeq2 CSVs: filter by padj
    GO_DESEQ(ch_de_results, ch_gene_db, ch_background, 'deseq2',
             params.enrichment.go.padj_filter ?: 'padj',
             params.enrichment.go.padj)
    // From gene list sources: no filtering (already filtered)
    GO_FILTERED(ch_filtered, ch_gene_db, ch_background, 'filter_de', '', 0)
    GO_FUZZY(ch_fuzzy_genes, ch_gene_db, ch_background, 'fuzzy_clustering', '', 0)
    GO_HDBSCAN(ch_hdbscan_genes, ch_gene_db, ch_background, 'hdbscan', '', 0)
    GO_KMEANS(ch_kmeans_genes, ch_gene_db, ch_background, 'kmeans', '', 0)
    GO_HIERARCHICAL(ch_hierarchical_genes, ch_gene_db, ch_background, 'hierarchical', '', 0)
    if (compare_enabled) {
        GO_COMPARE(COMPARE_DE.out.results, ch_gene_db, ch_background, 'compare_de', '', 0)
    }

    // --- REVIGO reduction over every GO enrichment output ---
    // Applied per source so cluster contexts stay separated. Internally,
    // revigo_reduce.R further splits each input by (ontology, gene_list).
    REVIGO_DESEQ(GO_DESEQ.out.results,         'deseq2')
    REVIGO_FILTERED(GO_FILTERED.out.results,   'filter_de')
    REVIGO_FUZZY(GO_FUZZY.out.results,         'fuzzy_clustering')
    REVIGO_HDBSCAN(GO_HDBSCAN.out.results,     'hdbscan')
    REVIGO_KMEANS(GO_KMEANS.out.results,       'kmeans')
    REVIGO_HIERARCHICAL(GO_HIERARCHICAL.out.results, 'hierarchical')
    if (compare_enabled) {
        REVIGO_COMPARE(GO_COMPARE.out.results, 'compare_de')
    }

    // --- KEGG enrichment ---
    KEGG_DESEQ(ch_de_results, ch_gene_db, ch_background, 'deseq2',
               params.enrichment.kegg.padj_filter ?: 'padj',
               params.enrichment.kegg.padj)
    KEGG_FILTERED(ch_filtered, ch_gene_db, ch_background, 'filter_de', '', 0)
    KEGG_FUZZY(ch_fuzzy_genes, ch_gene_db, ch_background, 'fuzzy_clustering', '', 0)
    KEGG_HDBSCAN(ch_hdbscan_genes, ch_gene_db, ch_background, 'hdbscan', '', 0)
    KEGG_KMEANS(ch_kmeans_genes, ch_gene_db, ch_background, 'kmeans', '', 0)
    KEGG_HIERARCHICAL(ch_hierarchical_genes, ch_gene_db, ch_background, 'hierarchical', '', 0)
    if (compare_enabled) {
        KEGG_COMPARE(COMPARE_DE.out.results, ch_gene_db, ch_background, 'compare_de', '', 0)
    }

    // --- WormCat ---
    // WormCat uses pvalue (more lenient, matches old pipeline)
    WORMCAT_DESEQ(ch_de_results, 'deseq2',
                  params.enrichment.wormcat.padj_filter ?: 'pvalue',
                  params.enrichment.wormcat.padj)
    WORMCAT_FILTERED(ch_filtered, 'filter_de', '', 0)
    WORMCAT_FUZZY(ch_fuzzy_genes, 'fuzzy_clustering', '', 0)
    WORMCAT_HDBSCAN(ch_hdbscan_genes, 'hdbscan', '', 0)
    WORMCAT_KMEANS(ch_kmeans_genes, 'kmeans', '', 0)
    WORMCAT_HIERARCHICAL(ch_hierarchical_genes, 'hierarchical', '', 0)
    if (compare_enabled) {
        WORMCAT_COMPARE(COMPARE_DE.out.results, 'compare_de', '', 0)
    }

    // --- GSEA on full ranked DESeq2 results (no filtering) ---
    GSEA(ch_de_results, ch_gene_db, ch_gsea_manifest, ch_gsea_collection_files, ch_ortholog_map)

    // --- ssGSEA per-sample pathway scores on VST counts ---
    // Gene-driver plots need both default and compound-contrast DE CSVs together;
    // mix() concatenates the two channels, collect() folds into a single staged dir.
    ch_de_for_ssgsea = ch_de_results.mix(ch_contrast_results).collect()
    SSGSEA(ch_vst_counts, ch_sample_meta, ch_gene_db,
           ch_gsea_manifest, ch_gsea_collection_files, ch_ortholog_map,
           ch_de_for_ssgsea)

    // --- Annotate all gene lists ---
    ANNOTATE_GENES(ch_filtered, ch_gene_db)

    emit:
    go_deseq       = GO_DESEQ.out.results
    go_filtered    = GO_FILTERED.out.results
    go_fuzzy       = GO_FUZZY.out.results
    go_hdbscan     = GO_HDBSCAN.out.results
    go_kmeans      = GO_KMEANS.out.results
    go_hierarchical = GO_HIERARCHICAL.out.results
    go_compare     = compare_enabled ? GO_COMPARE.out.results       : Channel.empty()
    revigo_deseq     = REVIGO_DESEQ.out.results
    revigo_filtered  = REVIGO_FILTERED.out.results
    revigo_fuzzy     = REVIGO_FUZZY.out.results
    revigo_hdbscan   = REVIGO_HDBSCAN.out.results
    revigo_kmeans    = REVIGO_KMEANS.out.results
    revigo_hierarchical = REVIGO_HIERARCHICAL.out.results
    revigo_compare   = compare_enabled ? REVIGO_COMPARE.out.results  : Channel.empty()
    kegg_deseq     = KEGG_DESEQ.out.results
    kegg_filtered  = KEGG_FILTERED.out.results
    kegg_fuzzy     = KEGG_FUZZY.out.results
    kegg_hdbscan   = KEGG_HDBSCAN.out.results
    kegg_kmeans    = KEGG_KMEANS.out.results
    kegg_hierarchical = KEGG_HIERARCHICAL.out.results
    kegg_compare   = compare_enabled ? KEGG_COMPARE.out.results     : Channel.empty()
    gsea           = GSEA.out.results
    ssgsea         = SSGSEA.out.results
    annotated      = ANNOTATE_GENES.out.results
    comparison     = compare_enabled ? COMPARE_DE.out.results       : Channel.empty()
}
