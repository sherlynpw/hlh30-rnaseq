#!/usr/bin/env Rscript
# Is hlh-30(tm1978) different from WT at each stage, for the fig S3 gene sets?
# Two tests, both with the biological replicate as the unit:
#   camera  limma gene-set test on the VST matrix, one contrast per stage; it
#           allows for the correlation between the genes of a set (the sperm
#           genes move together, so a test that treats genes as independent
#           gives p-values that are far too small)
#   score   per sample, the mean over the set of each gene's z-score across
#           samples; Welch t-test per stage
# BH over the 8 tests (2 sets x 4 stages) for each method.
#
#   Rscript germ-stats.R --out out-stats/germ_stats.tsv

suppressPackageStartupMessages({ library(optparse); library(limma); library(RSQLite) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--vst",     type = "character", default = "cache/vst.csv"),
  make_option("--meta",    type = "character", default = "cache/meta.tsv"),
  make_option("--sets",    type = "character", default = "germ=cache/germ_genes.txt,sperm=cache/sperm_genes.txt"),
  make_option("--gene-db", type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--out",     type = "character", default = "out-stats/germ_stats.tsv")
)))
dir.create(dirname(opt$out), recursive = TRUE, showWarnings = FALSE)
STAGES <- c("L1", "L4", "Day3", "Day6")

con <- dbConnect(SQLite(), opt$gene_db)
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name, sequence_name FROM annotations")
dbDisconnect(con)
sym2wb <- c(setNames(ann$wormbase_gene, ann$gene_name), setNames(ann$wormbase_gene, ann$sequence_name))

meta <- read.delim(opt$meta, stringsAsFactors = FALSE)
vst  <- as.matrix(read.csv(opt$vst, row.names = 1, check.names = FALSE))[, meta$sample_id]
vst  <- vst[rowMeans(vst) > min(vst) + 1, ]   # drop genes at the VST floor

kv <- strsplit(strsplit(opt$sets, ",")[[1]], "=")
sets <- setNames(lapply(kv, function(x) intersect(na.omit(sym2wb[readLines(x[2])]), rownames(vst))),
                 sapply(kv, `[`, 1))

grp <- factor(paste(meta$genotype, meta$stage, sep = "_"))
X <- model.matrix(~ 0 + grp); colnames(X) <- levels(grp)
z <- t(scale(t(vst)))

res <- do.call(rbind, lapply(STAGES, function(st) {
  cm <- makeContrasts(contrasts = sprintf("tm1978_%s - WT_%s", st, st), levels = X)
  cam <- camera(vst, sets, X, contrast = cm[, 1], inter.gene.cor = NA)
  do.call(rbind, lapply(names(sets), function(s) {
    sc <- colMeans(z[sets[[s]], , drop = FALSE])
    a <- sc[meta$stage == st & meta$genotype == "WT"]; b <- sc[meta$stage == st & meta$genotype == "tm1978"]
    data.frame(set = s, stage = st, n_genes = length(sets[[s]]), n_wt = length(a), n_hlh30 = length(b),
               score_diff = mean(b) - mean(a), score_p = t.test(b, a)$p.value,
               camera_direction = cam[s, "Direction"], camera_gene_cor = cam[s, "Correlation"],
               camera_p = cam[s, "PValue"])
  }))
}))
res$score_q  <- p.adjust(res$score_p, "BH")
res$camera_q <- p.adjust(res$camera_p, "BH")
write.table(res, opt$out, sep = "\t", quote = FALSE, row.names = FALSE)
print(within(res, { score_diff <- round(score_diff, 2); camera_gene_cor <- round(camera_gene_cor, 2) }), digits = 2)
