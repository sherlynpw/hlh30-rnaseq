#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

/*
 * C. elegans RNA-seq Pipeline
 *
 * Input paths:
 *   1. Plasmidsaurus UMI-dedup BAMs → featureCounts (re-count with Ensembl 115)
 *   2. nf-core/rnaseq output BAMs   → featureCounts (same GTF)
 *      (nf-core/rnaseq is run separately via run.sh wrapper)
 *
 * All paths use the same GTF for counting, ensuring comparable results.
 * If no GTF is provided, it is auto-downloaded from Ensembl.
 */

include { INPUT_CHECK  } from './subworkflows/local/input_check'
include { COUNTING     } from './subworkflows/local/counting'
include { DIFFERENTIAL } from './subworkflows/local/differential'
include { EXPLORATION  } from './subworkflows/local/exploration'
include { CLUSTERING   } from './subworkflows/local/clustering'
include { ENRICHMENT   } from './subworkflows/local/enrichment'
include { VST_ALL      } from './modules/local/vst_all'

// Optional: path to nf-core/rnaseq output directory (from run.sh step 1)
params.rnaseq_outdir = null

workflow {

    // Validate inputs
    if (!params.id)           { error "Please provide --id or set 'id' in params YAML (e.g. id: hlh30-all-ruv)" }
    if (!params.organism.gtf) { error "Please provide organism.gtf - run ./scripts/upload_refs.sh first" }

    ch_gtf = Channel.value(file(params.organism.gtf))

    // Gene database (optional but recommended - avoids BioMart API calls)
    ch_gene_db = params.organism.gene_db
        ? Channel.value(file(params.organism.gene_db))
        : Channel.value(file("${projectDir}/assets/NO_GENE_DB"))

    if (!params.input) { error "Please provide --input (samplesheet YAML)" }

    // Parse samplesheet and split by input_type
    INPUT_CHECK(params.input)

    // Collect BAMs from Plasmidsaurus (direct from samplesheet)
    ch_plasmidsaurus_bams = INPUT_CHECK.out.plasmidsaurus_bams  // [meta, bam]

    // Collect BAMs from nf-core/rnaseq output (if FASTQ samples were aligned).
    // rnaseq_outdir is expected to be a shared pool that may contain BAMs from
    // unrelated runs, so we filter the glob against sample_ids declared in
    // this run's samplesheet. Without this filter, COUNTING would pick up
    // stranger samples whose names happen to live in the pool.
    ch_nfcore_bams = Channel.empty()
    if (params.rnaseq_outdir) {
        // Keyed by sample_id so .join() with the pool BAMs intersects to just
        // the BAMs declared in this run's samplesheet.
        ch_expected_std = INPUT_CHECK.out.fastq_standard.map { meta, _fq -> [meta.id, 'fastq_standard'] }
        ch_expected_tp  = INPUT_CHECK.out.fastq_3prime.map  { meta, _fq -> [meta.id, 'fastq_3prime']  }

        // Standard RNA-seq BAMs
        ch_standard_bams = Channel
            .fromPath("${params.rnaseq_outdir}/standard/star_salmon/*.markdup.sorted.bam", checkIfExists: false)
            .map { bam ->
                def sample_id = bam.baseName.replaceAll(/\.markdup\.sorted$/, '')
                [sample_id, bam]
            }
            .join(ch_expected_std)
            .map { sample_id, bam, input_type -> [[id: sample_id, input_type: input_type], bam] }

        // 3' RNA-seq BAMs
        ch_3prime_bams = Channel
            .fromPath("${params.rnaseq_outdir}/3prime/star_salmon/*.markdup.sorted.bam", checkIfExists: false)
            .map { bam ->
                def sample_id = bam.baseName.replaceAll(/\.markdup\.sorted$/, '')
                [sample_id, bam]
            }
            .join(ch_expected_tp)
            .map { sample_id, bam, input_type -> [[id: sample_id, input_type: input_type], bam] }

        ch_nfcore_bams = ch_standard_bams.mix(ch_3prime_bams)
    }

    // Merge all BAMs into one channel for unified counting
    ch_all_bams = ch_plasmidsaurus_bams.mix(ch_nfcore_bams)

    // Count with consistent GTF - same for all inputs
    COUNTING(ch_all_bams, ch_gtf, params.input)

    // DESeq2 is only needed for clustering or enrichment. A pure
    // PCA/exploration run skips it and uses VST_ALL for the combined VST.
    def need_differential = !params.skip_clustering || !params.skip_enrichment
    if (need_differential) {
        DIFFERENTIAL(COUNTING.out.count_matrix, COUNTING.out.sample_meta, ch_gene_db)
    }

    // VST for exploration/visualisation. Use the all-sample VST_ALL when
    // DESeq2 is skipped; otherwise the DESeq2 VST.
    if (!need_differential) {
        VST_ALL(COUNTING.out.count_matrix, COUNTING.out.sample_meta)
        ch_vst_main     = VST_ALL.out.vst_counts
        ch_vst_raw_main = VST_ALL.out.vst_counts_raw
    } else {
        ch_vst_main     = DIFFERENTIAL.out.vst_counts
        ch_vst_raw_main = DIFFERENTIAL.out.vst_counts_raw
    }

    // Exploration - PCA, UMAP, CV density (run in parallel)
    if (!params.skip_exploration) {
        EXPLORATION(ch_vst_main, ch_vst_raw_main, COUNTING.out.sample_meta)
    }

    // Clustering - VST-based (fuzzy/HDBSCAN/k-means/hierarchical).
    if (!params.skip_clustering) {
        CLUSTERING(ch_vst_main, DIFFERENTIAL.out.de_results, DIFFERENTIAL.out.filtered, COUNTING.out.sample_meta, ch_gene_db)
    }

    // Enrichment depends on clustering outputs, so it's gated on both flags.
    if (!params.skip_enrichment && !params.skip_clustering) {
        ENRICHMENT(
            DIFFERENTIAL.out.de_results,     // full DESeq2 CSVs (unfiltered) → GO/KEGG/WormCat filter internally
            DIFFERENTIAL.out.filtered,       // filter_de gene lists (padj + lfc)
            CLUSTERING.out.fuzzy,            // fuzzy cluster gene lists
            CLUSTERING.out.hdbscan,          // HDBSCAN cluster gene lists
            CLUSTERING.out.kmeans,           // k-means cluster gene lists
            CLUSTERING.out.hierarchical,     // hierarchical cluster gene lists
            ch_vst_main,                     // VST → ssGSEA
            COUNTING.out.sample_meta,        // sample metadata → ssGSEA limma design
            DIFFERENTIAL.out.contrast_results, // compound-contrast DE CSVs → ssGSEA gene-driver plots
            ch_gene_db
        )
    }
}

workflow.onComplete {
    log.info ""
    log.info "Pipeline completed: ${workflow.success ? 'SUCCESS' : 'FAILED'}"
    log.info "Duration:  ${workflow.duration}"
    log.info "Results:   ${params.outdir}/${params.id}"
    log.info ""
}
