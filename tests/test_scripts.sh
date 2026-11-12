#!/usr/bin/env bash
set -euo pipefail

# Test each R script locally with synthetic data
# Usage: ./tests/test_scripts.sh
#
# Requires R with all pipeline packages installed (or run inside the container):
#   docker run --rm -v $(pwd):/work -w /work rayandrew/celegans-r-analysis:0.4.0 ./tests/test_scripts.sh

cd "$(dirname "$0")/.."

TESTDIR="tests/data"
RESULTS="tests/results"
PASS=0
FAIL=0
SKIP=0

rm -rf "$RESULTS"
mkdir -p "$RESULTS"

run_test() {
    local name="$1"
    shift
    echo -n "Testing ${name}... "
    if "$@" > "${RESULTS}/${name}.log" 2>&1; then
        echo "PASS"
        PASS=$((PASS + 1))
    else
        echo "FAIL (see ${RESULTS}/${name}.log)"
        tail -5 "${RESULTS}/${name}.log" 2>/dev/null
        FAIL=$((FAIL + 1))
    fi
}

skip_test() {
    echo "SKIP: $1 ($2)"
    SKIP=$((SKIP + 1))
}

# Step 0: Generate test data
echo "=== Generating test data ==="
Rscript tests/generate_test_data.R

# Step 0b: Build gene_db if go-kegg cache is available
GENE_DB=""
GOKEGG_CACHE="${GOKEGG_CACHE:-.cache/go-kegg}"
if [[ -d "$GOKEGG_CACHE" ]]; then
    echo "Building gene_db.sqlite from go-kegg cache..."
    run_test "build_gene_db" \
        Rscript bin/build_gene_db.R \
            --cache "$GOKEGG_CACHE" \
            --output "${RESULTS}/gene_db.sqlite"
    if [[ -f "${RESULTS}/gene_db.sqlite" ]]; then
        GENE_DB="${RESULTS}/gene_db.sqlite"
        echo "  gene_db available for testing"
    fi
else
    echo "  go-kegg cache not found, skipping gene_db tests"
fi

# Step 1: merge_counts
run_test "merge_counts" \
    Rscript bin/merge_counts.R \
        --counts_dir "$TESTDIR" \
        --metadata "${TESTDIR}/samplesheet_test.yaml" \
        --outdir "${RESULTS}/counts"

# Step 2: deseq2
if [[ -f "${RESULTS}/counts/gene_counts.tsv" ]]; then
    run_test "deseq2" \
        Rscript bin/deseq2.R \
            --counts "${RESULTS}/counts/gene_counts.tsv" \
            --metadata "${RESULTS}/counts/sample_meta.tsv" \
            --condition_from "genotype_timepoint" \
            --reference "WT_Day1" \
            --padj 0.05 \
            --lfc 1.0 \
            --shrinkage "apeglm" \
            --outdir "${RESULTS}/deseq2"
else
    skip_test "deseq2" "merge_counts failed"
fi

# Step 3: filter_de
if [[ -d "${RESULTS}/deseq2/results" ]]; then
    run_test "filter_de" \
        Rscript bin/filter_de.R \
            --input_dir "${RESULTS}/deseq2/results" \
            --padj 0.05 \
            --lfc 1.0 \
            --outdir "${RESULTS}/filtered"
else
    skip_test "filter_de" "deseq2 failed"
fi

# Step 4: PCA
if [[ -f "${RESULTS}/deseq2/vst_counts.csv" ]]; then
    run_test "pca" \
        Rscript bin/pca.R \
            --vst "${RESULTS}/deseq2/vst_counts.csv" \
            --metadata "${RESULTS}/counts/sample_meta.tsv" \
            --outdir "${RESULTS}/pca"
else
    skip_test "pca" "deseq2 failed"
fi

# Step 5: UMAP
if [[ -f "${RESULTS}/deseq2/vst_counts.csv" ]]; then
    run_test "umap" \
        Rscript bin/umap.R \
            --vst "${RESULTS}/deseq2/vst_counts.csv" \
            --metadata "${RESULTS}/counts/sample_meta.tsv" \
            --outdir "${RESULTS}/umap"
else
    skip_test "umap" "deseq2 failed"
fi

# Step 6: CV density
if [[ -f "${RESULTS}/deseq2/vst_counts.csv" ]]; then
    run_test "cv_density" \
        Rscript bin/cv_density.R \
            --counts "${RESULTS}/deseq2/vst_counts.csv" \
            --metadata "${RESULTS}/counts/sample_meta.tsv" \
            --outdir "${RESULTS}/cv_density"
else
    skip_test "cv_density" "deseq2 failed"
fi

# Step 7: Fuzzy clustering
if [[ -d "${RESULTS}/filtered" && -f "${RESULTS}/deseq2/vst_counts.csv" ]]; then
    FUZZY_ARGS="--vst ${RESULTS}/deseq2/vst_counts.csv --genes ${RESULTS}/filtered --num_clusters 3 --outdir ${RESULTS}/fuzzy"
    [[ -f "${RESULTS}/counts/sample_meta.tsv" ]] && FUZZY_ARGS="$FUZZY_ARGS --metadata ${RESULTS}/counts/sample_meta.tsv"
    [[ -n "$GENE_DB" ]] && FUZZY_ARGS="$FUZZY_ARGS --gene_db $GENE_DB"
    run_test "fuzzy_clustering" Rscript bin/fuzzy_clustering.R $FUZZY_ARGS
else
    skip_test "fuzzy_clustering" "filter_de or deseq2 failed"
fi

# Step 8: HDBSCAN
if [[ -d "${RESULTS}/filtered" && -f "${RESULTS}/deseq2/vst_counts.csv" ]]; then
    HDBSCAN_ARGS="--vst ${RESULTS}/deseq2/vst_counts.csv --genes ${RESULTS}/filtered --outdir ${RESULTS}/hdbscan"
    [[ -f "${RESULTS}/counts/sample_meta.tsv" ]] && HDBSCAN_ARGS="$HDBSCAN_ARGS --metadata ${RESULTS}/counts/sample_meta.tsv"
    [[ -n "$GENE_DB" ]] && HDBSCAN_ARGS="$HDBSCAN_ARGS --gene_db $GENE_DB"
    run_test "hdbscan" Rscript bin/hdbscan_clustering.R $HDBSCAN_ARGS
else
    skip_test "hdbscan" "filter_de or deseq2 failed"
fi

# Step 8b: KMeans clustering
if [[ -d "${RESULTS}/filtered" && -f "${RESULTS}/deseq2/vst_counts.csv" ]]; then
    KMEANS_ARGS="--vst ${RESULTS}/deseq2/vst_counts.csv --genes ${RESULTS}/filtered --outdir ${RESULTS}/kmeans --num_clusters auto --cluster_range 2,6"
    [[ -f "${RESULTS}/counts/sample_meta.tsv" ]] && KMEANS_ARGS="$KMEANS_ARGS --metadata ${RESULTS}/counts/sample_meta.tsv"
    [[ -n "$GENE_DB" ]] && KMEANS_ARGS="$KMEANS_ARGS --gene_db $GENE_DB"
    run_test "kmeans" Rscript bin/kmeans_clustering.R $KMEANS_ARGS
    if [[ -d "${RESULTS}/kmeans" ]]; then
        for kmeans_output in \
            cluster_assignments.csv \
            cluster_summary.csv \
            k_selection_summary.csv \
            plots/elbow_selection.pdf
        do
            if [[ ! -f "${RESULTS}/kmeans/${kmeans_output}" ]]; then
                echo "Missing expected kmeans output: ${RESULTS}/kmeans/${kmeans_output}" >> "${RESULTS}/kmeans.log"
                echo "FAIL (missing ${kmeans_output})"
                FAIL=$((FAIL + 1))
            fi
        done
        grep -q "cluster" "${RESULTS}/kmeans/cluster_assignments.csv" || { echo "cluster column missing" >> "${RESULTS}/kmeans.log"; FAIL=$((FAIL + 1)); }
        grep -q "tot_withinss" "${RESULTS}/kmeans/k_selection_summary.csv" || { echo "tot_withinss column missing" >> "${RESULTS}/kmeans.log"; FAIL=$((FAIL + 1)); }
        find "${RESULTS}/kmeans" -maxdepth 1 -name 'cluster_*_genes.txt' | grep -q . || { echo "no kmeans gene lists found" >> "${RESULTS}/kmeans.log"; FAIL=$((FAIL + 1)); }
    fi
else
    skip_test "kmeans" "filter_de or deseq2 failed"
fi

# Step 9: Compare DE
if [[ -d "${RESULTS}/deseq2/results" ]]; then
    run_test "compare_de" \
        Rscript bin/compare_de.R \
            --input_dir "${RESULTS}/deseq2/results" \
            --outdir "${RESULTS}/compare_de"
else
    skip_test "compare_de" "deseq2 failed"
fi

# Step 10: GO enrichment (from DESeq2 results - filters by padj internally)
if [[ -d "${RESULTS}/deseq2/results" ]]; then
    GO_ARGS="--input ${RESULTS}/deseq2/results --filter_by padj --threshold 0.05 --org_db org.Ce.eg.db --outdir ${RESULTS}/go"
    [[ -n "$GENE_DB" ]] && GO_ARGS="$GO_ARGS --gene_db $GENE_DB"
    run_test "go_enrichment" Rscript bin/go_enrichment.R $GO_ARGS
else
    skip_test "go_enrichment" "deseq2 failed"
fi

# Step 10b: GO enrichment from gene list .txt files (filter_de)
if [[ -d "${RESULTS}/filtered" ]]; then
    GO_ARGS2="--input ${RESULTS}/filtered --org_db org.Ce.eg.db --outdir ${RESULTS}/go_filtered"
    [[ -n "$GENE_DB" ]] && GO_ARGS2="$GO_ARGS2 --gene_db $GENE_DB"
    run_test "go_enrichment_filtered" Rscript bin/go_enrichment.R $GO_ARGS2
else
    skip_test "go_enrichment_filtered" "filter_de failed"
fi

# Step 10c: GO enrichment from kmeans cluster gene lists
if [[ -d "${RESULTS}/kmeans" ]]; then
    GO_KMEANS_ARGS="--input ${RESULTS}/kmeans --org_db org.Ce.eg.db --outdir ${RESULTS}/go_kmeans"
    [[ -n "$GENE_DB" ]] && GO_KMEANS_ARGS="$GO_KMEANS_ARGS --gene_db $GENE_DB"
    run_test "go_enrichment_kmeans" Rscript bin/go_enrichment.R $GO_KMEANS_ARGS
else
    skip_test "go_enrichment_kmeans" "kmeans failed"
fi

# Step 11: KEGG enrichment (from DESeq2 results)
if [[ -d "${RESULTS}/deseq2/results" ]]; then
    KEGG_ARGS="--input ${RESULTS}/deseq2/results --filter_by padj --threshold 0.05 --kegg_org cel --outdir ${RESULTS}/kegg"
    [[ -n "$GENE_DB" ]] && KEGG_ARGS="$KEGG_ARGS --gene_db $GENE_DB"
    run_test "kegg_enrichment" Rscript bin/kegg_enrichment.R $KEGG_ARGS
else
    skip_test "kegg_enrichment" "deseq2 failed"
fi

# Step 12: GSEA
if [[ -d "${RESULTS}/deseq2/results" ]]; then
    GSEA_ARGS="--de_results_dir ${RESULTS}/deseq2/results --ranking stat --org_db org.Ce.eg.db --collections hallmark,curated,custom --hallmark_gmt ${TESTDIR}/hallmark_test.gmt --curated_gmt ${TESTDIR}/curated_test.gmt --custom_gene_set_file ${TESTDIR}/neuron_genes.txt --custom_gene_set_name neuron_genes --custom_gene_set_format list --ortholog_map ${TESTDIR}/ortholog_map.tsv --ortholog_from wormbase_gene --ortholog_to human_symbol --outdir ${RESULTS}/gsea"
    [[ -n "$GENE_DB" ]] && GSEA_ARGS="$GSEA_ARGS --gene_db $GENE_DB"
    run_test "gsea" Rscript bin/gsea.R $GSEA_ARGS
else
    skip_test "gsea" "deseq2 failed"
fi

# Step 13: Annotate genes
if [[ -d "${RESULTS}/filtered" ]]; then
    ANNO_ARGS="--genes ${RESULTS}/filtered --outdir ${RESULTS}/annotation"
    [[ -n "$GENE_DB" ]] && ANNO_ARGS="$ANNO_ARGS --gene_db $GENE_DB"
    run_test "annotate_genes" Rscript bin/annotate_genes.R $ANNO_ARGS
else
    skip_test "annotate_genes" "filter_de failed"
fi

# Step 14: featurecounts.R (needs real BAM + GTF)
TEST_BAM="data/test_sample_dedup-mapped-reads.bam"
TEST_GTF="/refs/genome.gtf"  # mounted via docker -v
if [[ -f "$TEST_GTF" ]] || [[ -f "/tmp/Caenorhabditis_elegans.WBcel235.115.gtf" ]]; then
    [[ ! -f "$TEST_GTF" ]] && TEST_GTF="/tmp/Caenorhabditis_elegans.WBcel235.115.gtf"
    if [[ -f "$TEST_BAM" ]] || [[ -L "$TEST_BAM" ]]; then
        run_test "featurecounts" \
            Rscript bin/featurecounts.R \
                "$TEST_BAM" "$TEST_GTF" "test_featurecounts" 0 FALSE 4
    else
        skip_test "featurecounts" "BAM not found: $TEST_BAM"
    fi
else
    skip_test "featurecounts" "GTF not found (mount with -v /tmp/...gtf:/refs/genome.gtf)"
fi

# Step 15: WormCat from gene list files
if [[ -d "${RESULTS}/filtered" ]] && ls "${RESULTS}/filtered"/*genes*.txt >/dev/null 2>&1; then
    run_test "wormcat" \
        Rscript bin/wormcat.R \
            --input "${RESULTS}/filtered" \
            --outdir "${RESULTS}/wormcat"
else
    skip_test "wormcat" "no gene lists with real WormBase IDs"
fi

# Step 15b: WormCat from DESeq2 results (pvalue filter)
if [[ -d "${RESULTS}/deseq2/results" ]]; then
    run_test "wormcat_deseq" \
        Rscript bin/wormcat.R \
            --input "${RESULTS}/deseq2/results" \
            --filter_by pvalue \
            --threshold 0.05 \
            --outdir "${RESULTS}/wormcat_deseq"
else
    skip_test "wormcat_deseq" "deseq2 failed"
fi

echo ""
echo "=== Results ==="
echo "PASS: ${PASS}"
echo "FAIL: ${FAIL}"
echo "SKIP: ${SKIP}"
echo ""

if [[ $FAIL -gt 0 ]]; then
    echo "Failed tests - check logs in ${RESULTS}/"
    exit 1
else
    echo "All executed tests passed!"
fi
