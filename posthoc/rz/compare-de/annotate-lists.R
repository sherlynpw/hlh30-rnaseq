#!/usr/bin/env Rscript
# Turn compare_de's bare WBGene lists into complete annotated DEG tables:
# gene names, both contrasts' LFC/padj, category, sorted by padj ascending.
# Split into up/down so collaborators get the full post-DESeq2 list, not just
# the genes that happened to fall into a GO term.

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr); library(RSQLite) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--compared", type = "character", help = "*_compared_genes.csv from compare_de.R"),
  make_option("--gene-db",  type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--outdir",   type = "character", default = "out")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

d <- read_csv(opt$compared, show_col_types = FALSE)
lfc <- grep("^log2FoldChange\\.", names(d), value = TRUE)
pad <- grep("^padj\\.", names(d), value = TRUE)
labs <- sub("^log2FoldChange\\.", "", lfc)
message("contrasts: ", paste(labs, collapse = " | "))

con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name, sequence_name, concise_description FROM annotations")

d <- d %>% left_join(ann, by = c("gene_id" = "wormbase_gene"))

write_set <- function(df, defining_padj, defining_lfc, name) {
  if (nrow(df) == 0) { message("  (empty) ", name); return(invisible()) }
  out <- df %>%
    mutate(direction = ifelse(.data[[defining_lfc]] > 0, "up", "down")) %>%
    arrange(.data[[defining_padj]])
  for (dir_ in c("up", "down")) {
    sub <- out %>% filter(direction == dir_) %>%
      select(gene_id, gene_name, sequence_name, category, direction,
             all_of(c(lfc, pad)), lfc_diff, concise_description)
    f <- file.path(opt$outdir, sprintf("%s_%s.tsv", name, dir_))
    write_tsv(sub, f)
    message(sprintf("  %-46s %5d genes", basename(f), nrow(sub)))
  }
}

for (i in seq_along(labs)) {
  cat_name <- paste0(labs[i], "_exclusive")
  write_set(d %>% filter(category == cat_name), pad[i], lfc[i], cat_name)
}
# shared / differential: sort by the more significant of the two padj
d2 <- d %>% mutate(padj_min = pmin(.data[[pad[1]]], .data[[pad[2]]], na.rm = TRUE))
for (cat_name in c("shared", "differential")) {
  sub <- d2 %>% filter(category == cat_name)
  if (nrow(sub)) write_set(sub, "padj_min", lfc[1], cat_name)
}

# one combined table, everything that is not 'unchanged'
allc <- d2 %>% filter(category != "not_de") %>% arrange(padj_min) %>%
  select(gene_id, gene_name, sequence_name, category,
         all_of(c(lfc, pad)), lfc_diff, concise_description)
write_tsv(allc, file.path(opt$outdir, "all_categories.tsv"))
message(sprintf("  %-46s %5d genes", "all_categories.tsv", nrow(allc)))
