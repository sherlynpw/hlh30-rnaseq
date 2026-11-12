#!/usr/bin/env Rscript
# The fig S3 gene set: every gene in a germ-cell GO term that is enriched
# (adj. p < 0.05) in any compare-DE gene list of the three adult transitions
# (L4 -> D3, L4 -> D6, D3 -> D6). Built from the GO tables so the set follows
# the current analysis; it used to be a hand-kept file and went stale.
# --kw and --exclude pick the terms, so the same script makes the oocyte and
# germline set (S3B) and the sperm set (S3C, the pseudopodium term of the
# major sperm protein genes, whose name has no germ-cell word).
#
#   Rscript germ-genes.R --go cache/go --out cache/germ_genes.txt

suppressPackageStartupMessages({ library(optparse); library(readr); library(dplyr); library(stringr) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--go",  type = "character", default = "cache/go", help = "compare-DE GO dirs, one per gene list"),
  make_option("--out", type = "character", default = "cache/germ_genes.txt"),
  make_option("--padj", type = "double", default = 0.05),
  make_option("--kw", type = "character",
              default = "germ|gamete|oogenesis|spermato|oocyte|\\bsperm\\b|meiosis|meiotic|P granule|P-granule|reproduc|gonad|fertiliz|sex determination|sex differ",
              help = "regex a term name must match"),
  make_option("--exclude", type = "character", default = "",
              help = "regex a term name must not match (empty = none)")
)))
files <- Sys.glob(file.path(opt$go, "*", "go_*.csv"))
files <- files[str_detect(basename(dirname(files)),
                          "hlh-30_(Day3|Day6)_vs_hlh-30_(L4|Day3)__X__N2_(Day3|Day6)_vs_N2_(L4|Day3)")]
stopifnot(length(files) > 0)
terms <- bind_rows(lapply(files, function(f) {
  read_csv(f, show_col_types = FALSE) %>%
    filter(p.adjust < opt$padj, str_detect(Description, regex(opt$kw, ignore_case = TRUE)),
           !nzchar(opt$exclude) | !str_detect(Description, regex(opt$exclude, ignore_case = TRUE))) %>%
    transmute(list = sub(".*_genes_", "", basename(dirname(f))), ontology = sub("go_(.*)\\.csv", "\\1", basename(f)),
              term = Description, p.adjust, geneID)
}))
genes <- sort(unique(unlist(str_split(terms$geneID, "/"))))
writeLines(genes, opt$out)
write_tsv(terms, sub("\\.txt$", "_terms.tsv", opt$out))
message(length(genes), " genes from ", nrow(terms), " term rows in ", length(files), " GO tables -> ", opt$out)
