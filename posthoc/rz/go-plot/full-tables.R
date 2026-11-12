#!/usr/bin/env Rscript
# One spreadsheet per GO paper panel with every significant term of its gene
# list and ontology. The terms the panel shows are bold, and the "shown in
# panel" column keeps that mark when the sheet is saved as CSV. The terms a
# panel showed before the 2026-10-12 forward-strand rerun that are no longer
# significant (--lost) follow in grey, with their current values from the
# no-cutoff GO rerun in --lost-dir (`make lost`) and their old adjusted P.
#
#   Rscript full-tables.R --config targets-paper-figures.tsv --outdir out-full

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(stringr); library(openxlsx2)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "go-plot-common.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--bucket", type = "character", default = "s3://celegans-pipeline-rayandrew"),
  make_option("--cache",  type = "character", default = "cache"),
  make_option("--lost",   type = "character", default = NULL,
              help = "TSV of terms shown before the rerun that are no longer significant"),
  make_option("--lost-dir", dest = "lost_dir", type = "character", default = "out-lost",
              help = "GO tables with no cutoff for the --lost gene lists (rerun-lost.sh)"),
  make_option("--old-prefix", dest = "old_prefix", type = "character", default = "_unstranded-2026-10-12",
              help = "results/ prefix that holds the run the old panels came from")
)))
stopifnot(!is.null(opt$config), !is.null(opt$outdir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

targets <- read_tsv(opt$config, show_col_types = FALSE)
lost <- if (is.null(opt$lost)) tibble(fig = character()) else read_tsv(opt$lost, comment = "#", show_col_types = FALSE)
OLD <- "no longer significant (shown before 2026-10-12)"
for (i in seq_len(nrow(targets))) {
  row <- targets[i, ]
  fig <- str_extract(row$outfile, "(?<=\\[)Fig [^\\]]+")
  csv <- ensure_csv(row$run, row$path, row$ontology, opt$bucket, opt$cache)
  go  <- read_csv(csv, show_col_types = FALSE)
  # the pipeline writes only terms with adjusted p < 0.05
  stopifnot(all(go$p.adjust < 0.05))
  # the label the panel prints: the full name when it fits in 2 lines
  shown <- select_terms(go, row$terms, 20, fig) %>%
    mutate(display = gsub("\n", " ", paper_label_for(Description, display, row$placed_in))) %>%
    select(Description, display)

  tab <- go %>%
    left_join(shown, by = "Description") %>%
    arrange(p.adjust) %>%
    transmute(`GO ID` = ID, Term = Description,
              `Shown in panel` = ifelse(is.na(display), "", "yes"),
              `Panel label` = ifelse(is.na(display) | display == Description, "", display),
              Ontology = row$ontology, Genes = Count, GeneRatio, BgRatio,
              `Fold enrichment` = round(FoldEnrichment, 2),
              `P value` = signif(pvalue, 3), `Adjusted P` = signif(p.adjust, 3),
              `Gene names` = str_replace_all(geneID, "/", ", "))

  # current values from the no-cutoff rerun of the same gene list; the old
  # adjusted P from the run the author's panel was made from
  gone <- lost %>% filter(fig == !!fig)
  old <- bind_rows(lapply(split(gone, seq_len(nrow(gone))), function(g) {
    stopifnot(g$run == row$run, g$path == row$path, !g$term %in% go$Description)
    rerun <- read_csv(file.path(opt$lost_dir, g$run, basename(g$path), sprintf("go_%s.csv", g$ontology)),
                      show_col_types = FALSE)
    # the rerun must reproduce the pipeline's own table
    chk <- inner_join(go, rerun, by = "ID")
    stopifnot(nrow(chk) == nrow(go), isTRUE(all.equal(chk$p.adjust.x, chk$p.adjust.y)))
    f <- ensure_csv(file.path(opt$old_prefix, g$old_run), g$old_path, g$ontology, opt$bucket, opt$cache)
    was <- read_csv(f, show_col_types = FALSE) %>% filter(Description == g$term)
    hit <- rerun %>% filter(Description == g$term)
    stopifnot(nrow(hit) == 1, nrow(was) == 1)
    transmute(hit, `GO ID` = ID, Term = Description, `Shown in panel` = OLD, `Panel label` = "",
              Ontology = g$ontology, Genes = Count, GeneRatio, BgRatio,
              `Fold enrichment` = round(FoldEnrichment, 2),
              `P value` = signif(pvalue, 3), `Adjusted P` = signif(p.adjust, 3),
              `Gene names` = str_replace_all(geneID, "/", ", "),
              `Adjusted P before 2026-10-12` = signif(was$p.adjust, 3))
  }))
  tab <- bind_rows(tab, old)

  n <- ncol(tab)
  bold <- which(tab$`Shown in panel` == "yes") + 4
  grey <- which(tab$`Shown in panel` == OLD) + 4
  wb <- wb_workbook()$add_worksheet("terms")$
    add_data(x = paste0(sprintf("%s: all %d significant GO %s terms of the panel's gene list (BH-adjusted P < 0.05). Bold = shown in the panel.",
                                fig, nrow(tab) - nrow(old), row$ontology),
                        if (nrow(old)) " Grey = shown in the panel before the 2026-10-12 forward-strand rerun and no longer significant; values from the same gene list with no cutoff, last column = adjusted P before the rerun." else ""),
             start_row = 1)$
    add_data(x = sprintf("Gene list: %s / %s", row$run, row$path), start_row = 2)$
    add_data(x = tab, start_row = 4)$
    add_font(dims = wb_dims(rows = 4, cols = seq_len(n)), bold = TRUE)$
    set_col_widths(cols = seq_len(n), widths = c(12, 45, 9, 22, 8, 6, 9, 10, 9, 10, 10, 60, 12)[seq_len(n)])$
    freeze_pane(first_active_row = 5)
  for (r in bold) wb$add_font(dims = wb_dims(rows = r, cols = seq_len(n)), bold = TRUE)
  for (r in grey) wb$add_font(dims = wb_dims(rows = r, cols = seq_len(n)), color = wb_color(hex = "FF808080"))
  out <- file.path(opt$outdir, sprintf("%s-full.xlsx", str_replace(tolower(fig), " ", "-")))
  wb$save(out)
  message(sprintf("[wrote] %s  (%d terms, %d shown, %d no longer significant)", out, nrow(tab) - nrow(old), length(bold), nrow(old)))
}
