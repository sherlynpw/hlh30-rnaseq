#!/usr/bin/env Rscript
# Supplementary tables: how many genes are detected at each stage.
#
# The pipeline never produces this. Its only filter is rowSums(counts) >= 10
# applied ONCE across all samples (bin/deseq2.R:207), which defines the DE/PCA
# gene universe. Split by stage that rule becomes replicate-count dependent
# (Day 3 has 12 libraries, L1 has 6), so the per-stage figure quoted here is
# the per-SAMPLE count averaged within a stage, which is immune to unequal n.
#
# Depth is not a confounder at these library sizes (r ~ 0.06 between depth and
# genes detected; ~21 genes per extra 1M reads), so the raw >=10 threshold is
# safe to use per sample without normalising.
#
#   Rscript gene-detection-tables.R

suppressPackageStartupMessages({ library(optparse) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--counts", type = "character",
              default = "../kmeans-umap/cache/counts/gene_counts.tsv"),
  make_option("--meta",   type = "character",
              default = "../kmeans-umap/cache/hlh30-all-ruv/sample_meta.tsv"),
  make_option("--gene-db", type = "character", dest = "gene_db",
              default = "../../../gene_db.sqlite"),
  make_option("--outdir", type = "character", default = "out")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

cnt  <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
meta <- read.delim(opt$meta, stringsAsFactors = FALSE)
meta <- meta[match(colnames(cnt), meta$sample_id), ]
stopifnot(!any(is.na(meta$sample_id)))

STAGE_LAB <- c(L1 = "L1", L4 = "L4", Day3 = "Day 3", Day6 = "Day 6")
GENO_LAB  <- c(N2 = "WT", WT = "WT", `hlh-30` = "hlh-30(tm1978)", tm1978 = "hlh-30(tm1978)")
stages <- intersect(c("L1","L4","Day3","Day6"), unique(meta$stage))

depth <- colSums(cnt)
cpm   <- t(t(cnt) / depth) * 1e6

# ---- per-sample table --------------------------------------------------------
per_sample <- data.frame(
  sample_id   = meta$sample_id,
  genotype    = unname(GENO_LAB[meta$genotype]),
  stage       = unname(STAGE_LAB[meta$stage]),
  replicate   = meta$replicate,
  total_reads = as.integer(depth),
  genes_ge10_counts = as.integer(colSums(cnt >= 10)),
  genes_cpm_ge1     = as.integer(colSums(cpm >= 1)),
  stringsAsFactors = FALSE
)
per_sample <- per_sample[order(match(per_sample$stage, STAGE_LAB),
                               per_sample$genotype, per_sample$replicate), ]
f1 <- file.path(opt$outdir, "supp_gene_detection_per_sample.tsv")
write.table(per_sample, f1, sep = "\t", quote = FALSE, row.names = FALSE)

# ---- per-stage summary, overall and split by genotype ------------------------
summarise <- function(idx, stage, geno) {
  sub <- cnt[, idx, drop = FALSE]; subc <- cpm[, idx, drop = FALSE]
  data.frame(
    stage = stage, genotype = geno, n_samples = length(idx),
    mean_total_reads = round(mean(depth[idx])),
    # headline figure: per-sample >=10, averaged within the group
    mean_genes_ge10_per_sample = round(mean(colSums(sub >= 10)), 1),
    sd_genes_ge10_per_sample   = round(sd(colSums(sub >= 10)), 1),
    # set-based alternatives, for reference
    genes_ge10_summed_over_group = sum(rowSums(sub) >= 10),
    genes_cpm_ge1_in_all_reps    = sum(rowSums(subc >= 1) == length(idx)),
    genes_detected_any_read      = sum(rowSums(sub > 0) >= 1),
    stringsAsFactors = FALSE)
}
rows <- list()
for (s in stages) {
  k <- which(meta$stage == s)
  rows[[length(rows)+1]] <- summarise(k, STAGE_LAB[[s]], "both genotypes")
  for (g in unique(meta$genotype[k])) {
    kk <- which(meta$stage == s & meta$genotype == g)
    rows[[length(rows)+1]] <- summarise(kk, STAGE_LAB[[s]], unname(GENO_LAB[g]))
  }
}
summary_df <- do.call(rbind, rows)
f2 <- file.path(opt$outdir, "supp_gene_detection_per_stage.tsv")
write.table(summary_df, f2, sep = "\t", quote = FALSE, row.names = FALSE)

message("wrote ", f1, "\nwrote ", f2)
message("depth vs genes>=10 across all libraries: r = ",
        round(cor(depth, colSums(cnt >= 10)), 3))
