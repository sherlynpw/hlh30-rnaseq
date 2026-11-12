#!/usr/bin/env Rscript
# What fraction of each timepoint's DEGs have experimental HLH-30 ChIP-seq
# binding near their TSS (ChIP-Atlas), versus the expressed-gene background.
#
# Uses the ChIP-Atlas HLH-30 peaks, keyed on DESeq2 DEG lists.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2); library(RSQLite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--chipatlas", type = "character", help = "ChIP-Atlas hlh-30.<dist>.tsv"),
  make_option("--dist",      type = "character", default = "5", help = "kb window, for labels"),
  make_option("--deg-dir",   type = "character", dest = "deg_dir",
              help = "dir with <contrast>_genes_{de,up,down}.txt"),
  make_option("--de-dir",    type = "character", dest = "de_dir",
              help = "dir with DESeq2 <contrast>.csv (defines the tested background)"),
  make_option("--contrasts", type = "character",
              help = "comma-separated contrast names, e.g. tm1978_L1_vs_WT_L1,..."),
  make_option("--gene-db",   type = "character", dest = "gene_db"),
  make_option("--score-min", type = "double", default = 50, dest = "score_min",
              help = "min hlh-30|Average (-log10 q) to count as bound"),
  make_option("--font-scale", type = "double", default = 1, dest = "font_scale",
              help = "Multiply every font size; render at the ORIGINAL canvas."),
  make_option("--outdir",    type = "character")
)))
stopifnot(!is.null(opt$chipatlas), !is.null(opt$deg_dir), !is.null(opt$de_dir),
          !is.null(opt$contrasts), !is.null(opt$outdir))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

contrasts <- trimws(strsplit(opt$contrasts, ",")[[1]])

ca <- read_tsv(opt$chipatlas, show_col_types = FALSE)
avg <- grep("^hlh-30\\|Average$", names(ca), value = TRUE)[1]
ca <- ca %>% rename(symbol = Target_genes, score = !!sym(avg)) %>% filter(!is.na(score))
bound_syms <- ca$symbol[ca$score > opt$score_min]
message("ChIP-Atlas genes: ", nrow(ca), "; bound at score > ", opt$score_min, ": ", length(bound_syms))

con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name, sequence_name FROM annotations")
ann$lookup <- ifelse(!is.na(ann$gene_name) & nzchar(ann$gene_name),
                     ann$gene_name, ann$sequence_name)
wb2sym <- setNames(ann$lookup, ann$wormbase_gene)

stage_of <- function(c) {
  s <- sub("^tm1978_", "", sub("_vs_WT_.*$", "", c))
  c(L1 = "L1", L4 = "Day 0", Day3 = "Day 3", Day6 = "Day 6")[s]
}

rows <- list()
for (ct in contrasts) {
  de <- read_csv(file.path(opt$de_dir, paste0(ct, ".csv")), show_col_types = FALSE)
  bg_wb  <- de$gene_id[!is.na(de$padj)]
  bg_sym <- unique(na.omit(wb2sym[bg_wb])); bg_sym <- bg_sym[nzchar(bg_sym)]
  bg_bound   <- intersect(bg_sym, bound_syms)
  # Precomputed: inside tibble() a later expression would resolve `bg_bound`
  # to the column just defined (a scalar), not this vector.
  n_bg       <- length(bg_sym)
  n_bg_bound <- length(bg_bound)
  pct_bg     <- 100 * n_bg_bound / n_bg

  for (dir_ in c("de", "down", "up")) {
    f <- file.path(opt$deg_dir, sprintf("%s_genes_%s.txt", ct, dir_))
    if (!file.exists(f)) next
    wb <- readLines(f); wb <- wb[nzchar(wb)]
    sym <- unique(na.omit(wb2sym[wb])); sym <- intersect(sym[nzchar(sym)], bg_sym)
    hit <- intersect(sym, bound_syms)

    # DEGs bound vs background bound, Fisher on the 2x2
    m <- matrix(c(length(hit), length(sym) - length(hit),
                  length(bg_bound) - length(hit),
                  length(bg_sym) - length(sym) - (length(bg_bound) - length(hit))), nrow = 2)
    ft <- tryCatch(fisher.test(m), error = function(e) list(p.value = NA, estimate = NA))

    rows[[length(rows) + 1]] <- tibble(
      contrast = ct, stage = stage_of(ct), direction = dir_,
      n_deg = length(sym), n_bound = length(hit),
      pct_bound = 100 * length(hit) / max(length(sym), 1),
      bg_n = n_bg, bg_bound = n_bg_bound, bg_pct = pct_bg,
      odds_ratio = unname(ft$estimate), p = ft$p.value
    )
  }
}
res <- bind_rows(rows)
res$padj <- p.adjust(res$p, method = "BH")
write_tsv(res, file.path(opt$outdir, "chip_promoter_summary.tsv"))
print(as.data.frame(res %>% select(stage, direction, n_deg, n_bound, pct_bound, bg_pct, odds_ratio, p, padj)), digits = 3)

lab <- c(de = "All DEGs", down = "Down in mutant", up = "Up in mutant")
res <- res %>%
  mutate(panel = factor(lab[direction], levels = unname(lab)),
         stage = factor(stage, levels = c("L1", "Day 0", "Day 3", "Day 6")),
         star  = ifelse(is.na(padj), "", ifelse(padj < 0.001, "***",
                 ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", "")))))

bg_line <- res %>% group_by(panel) %>% summarise(bg = mean(bg_pct), .groups = "drop")

p <- ggplot(res, aes(stage, pct_bound)) +
  geom_hline(data = bg_line, aes(yintercept = bg),
             linetype = "dashed", colour = "grey45", linewidth = 0.4) +
  geom_col(fill = "#3B6E8F", width = 0.62) +
  geom_text(aes(label = sprintf("%.0f%%%s", pct_bound, star)),
            vjust = -0.45, size = 3.1 * opt$font_scale, colour = "grey20") +
  facet_wrap(~ panel) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.16))) +
  labs(
    x = NULL, y = "DEGs with HLH-30 ChIP-seq binding (%)",
    title = sprintf("HLH-30 promoter binding among DEGs (ChIP-Atlas, +/-%s kb, score > %g)",
                    opt$dist, opt$score_min),
    subtitle = "Dashed line = expressed-gene background. Stars = BH-adjusted Fisher vs background."
  ) +
  theme_bw(base_size = 11 * opt$font_scale) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        strip.background = element_rect(fill = "grey96", colour = "grey80"),
        plot.title = element_text(size = 11 * opt$font_scale, face = "bold"),
        plot.subtitle = element_text(size = 8.5 * opt$font_scale, colour = "grey35"))

ggsave(file.path(opt$outdir, "chip_promoter.pdf"), p, width = 8.6, height = 4.0)
ggsave(file.path(opt$outdir, "chip_promoter.png"), p, width = 8.6, height = 4.0, dpi = 200)
message("wrote: ", file.path(opt$outdir, "chip_promoter.{pdf,png,tsv}"))
