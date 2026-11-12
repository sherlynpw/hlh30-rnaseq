#!/usr/bin/env Rscript

# Generate small synthetic test data for local pipeline testing

set.seed(42)
outdir <- "tests/data"
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# 4 samples, 2 conditions, 2 replicates
samples <- c("WT_Day1_1", "WT_Day1_2", "mut_Day1_1", "mut_Day1_2")

# 1000 genes - enough for vst(), with 100 DE
n_genes <- 1000
gene_ids <- paste0("WBGene", sprintf("%08d", 1:n_genes))

# Base expression
base_expr <- rpois(n_genes, lambda = 100)

# 100 DE genes: 50 up, 50 down in mut
de_up   <- 1:50
de_down <- 51:100

for (i in seq_along(samples)) {
  counts <- base_expr + rpois(n_genes, lambda = 10)

  if (grepl("mut", samples[i])) {
    counts[de_up]   <- counts[de_up] * 4
    counts[de_down] <- as.integer(counts[de_down] / 4)
  }

  counts <- pmax(counts + round(rnorm(n_genes, 0, 5)), 0)

  # Write featureCounts format (comment line + TSV)
  outfile <- file.path(outdir, paste0(samples[i], ".featureCounts.txt"))
  header_line <- "# Program:featureCounts v2.0.6; Fake test data"

  df <- data.frame(
    Geneid = gene_ids,
    Chr    = "I",
    Start  = seq(1, by = 1000, length.out = n_genes),
    End    = seq(999, by = 1000, length.out = n_genes),
    Strand = "+",
    Length = 1000,
    Count  = as.integer(counts)
  )

  con <- file(outfile, "w")
  writeLines(header_line, con)
  close(con)
  suppressWarnings(
    write.table(df, outfile, sep = "\t", row.names = FALSE, quote = FALSE, append = TRUE)
  )
}

# Write samplesheet
samplesheet <- list(samples = lapply(seq_along(samples), function(i) {
  list(
    sample_id  = samples[i],
    input_type = "plasmidsaurus_bam",
    genotype   = sub("_D.*", "", samples[i]),
    timepoint  = "Day1",
    source     = "test",
    replicate  = as.integer(sub(".*_", "", samples[i]))
  )
}))
yaml::write_yaml(samplesheet, file.path(outdir, "samplesheet_test.yaml"))

# Write ortholog map and example gene sets for GSEA testing
ortholog_map <- data.frame(
  wormbase_gene = gene_ids[1:200],
  human_symbol = paste0("HUMAN", seq_len(200)),
  stringsAsFactors = FALSE
)
write.table(
  ortholog_map,
  file.path(outdir, "ortholog_map.tsv"),
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

writeLines(
  c(
    paste(c("HALLMARK_TEST_UP", "na", paste0("HUMAN", 1:25)), collapse = "\t"),
    paste(c("HALLMARK_TEST_DOWN", "na", paste0("HUMAN", 51:75)), collapse = "\t")
  ),
  file.path(outdir, "hallmark_test.gmt")
)

writeLines(
  c(
    paste(c("CURATED_TEST_MIXED", "na", paste0("HUMAN", 10:35)), collapse = "\t"),
    paste(c("CURATED_TEST_DOWN", "na", paste0("HUMAN", 76:100)), collapse = "\t")
  ),
  file.path(outdir, "curated_test.gmt")
)

writeLines(gene_ids[1:30], file.path(outdir, "neuron_genes.txt"))

# Write params
params <- list(
  input = "tests/data/samplesheet_test.yaml",
  organism = list(
    genome = "WBcel235", gtf = NULL, ensembl_release = 115,
    kegg_org = "cel", org_db = "org.Ce.eg.db",
    biomart_dataset = "celegans_gene_ensembl", gene_db = NULL
  ),
  deseq2 = list(
    padj_threshold = 0.05, lfc_threshold = 1.0, shrinkage = "apeglm",
    condition_from = "genotype_timepoint", reference = "WT_Day1", batch = NULL
  ),
  groups = list(all = list(genotype = c("WT", "mut"))),
  clustering = list(num_clusters = "3", cluster_range = "2,6", sort_by = "pvalue", correct_vst = FALSE),
  enrichment = list(
    go      = list(enabled = TRUE,  padj = 0.05, lfc = 1.0),
    kegg    = list(enabled = TRUE,  padj = 0.05, lfc = 1.0),
    wormcat = list(enabled = FALSE, padj = 0.1),
    gsea    = list(
      enabled = TRUE,
      ranking = "stat",
      collections = list("hallmark", "curated", "custom"),
      hallmark_gmt = "tests/data/hallmark_test.gmt",
      curated_gmt = "tests/data/curated_test.gmt",
      custom_gene_set_file = "tests/data/neuron_genes.txt",
      custom_gene_set_name = "neuron_genes",
      custom_gene_set_format = "list",
      custom_gene_set_map_ortholog = FALSE,
      ortholog_map = "tests/data/ortholog_map.tsv",
      ortholog_from = "wormbase_gene",
      ortholog_to = "human_symbol"
    )
  ),
  exploration = list(pca_color_by = list()),
  outdir = "tests/results"
)
yaml::write_yaml(params, file.path(outdir, "params_test.yaml"))

message("Test data generated in ", outdir)
message("  4 samples, ", n_genes, " genes, ", length(c(de_up, de_down)), " DE genes")
