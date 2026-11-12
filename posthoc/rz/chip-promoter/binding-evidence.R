#!/usr/bin/env Rscript
# Compare three lines of HLH-30 binding evidence over the same DEG lists:
#   ChIP-Atlas  - experimental peak near TSS (measurement)
#   E-box       - CACGTG in the -2000/+200 promoter (sequence prediction)
#   co-bound    - E-box with a DBE (DAF-16) or DAE (PQM-1) within 500 bp
# Percentages are of DEGs; the dashed line is the expressed-gene background.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(ggplot2); library(RSQLite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--chipatlas",  type = "character"),
  make_option("--dist",       type = "character", default = "1"),
  make_option("--motif-dir",  type = "character", dest = "motif_dir",
              help = "dir with genes_ebox.txt, genes_hlh30_dafl16_cobound.txt"),
  make_option("--deg-dir",    type = "character", dest = "deg_dir"),
  make_option("--de-dir",     type = "character", dest = "de_dir"),
  make_option("--contrasts",  type = "character"),
  make_option("--gene-db",    type = "character", dest = "gene_db"),
  make_option("--score-min",  type = "double", default = 50, dest = "score_min"),
  make_option("--outdir",     type = "character")
)))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
contrasts <- trimws(strsplit(opt$contrasts, ",")[[1]])

rd <- function(f) { p <- file.path(opt$motif_dir, f)
                    if (file.exists(p)) { x <- readLines(p); x[nzchar(x)] } else character(0) }
ebox    <- rd("genes_ebox.txt")
cobound <- rd("genes_hlh30_dafl16_cobound.txt")

ca <- read_tsv(opt$chipatlas, show_col_types = FALSE)
avg <- grep("^hlh-30\\|Average$", names(ca), value = TRUE)[1]
ca <- ca %>% rename(symbol = Target_genes, score = !!sym(avg)) %>% filter(!is.na(score))
bound_syms <- ca$symbol[ca$score > opt$score_min]

con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name, sequence_name FROM annotations")
ann$lookup <- ifelse(!is.na(ann$gene_name) & nzchar(ann$gene_name), ann$gene_name, ann$sequence_name)
sym2wb <- setNames(ann$wormbase_gene, ann$lookup)
chip_wb <- unique(na.omit(sym2wb[bound_syms]))   # ChIP evidence as WBGene IDs

evidence <- list(`ChIP-Atlas peak` = chip_wb, `E-box in promoter` = ebox,
                 `E-box + DAF-16/PQM-1` = cobound)
# Captured here: inside mutate() below, `evidence` resolves to the column of
# that name, so names(evidence) would be NULL and every level would go NA.
ev_levels <- names(evidence)

stage_of <- function(c) {
  s <- sub("^tm1978_", "", sub("_vs_WT_.*$", "", c))
  c(L1 = "L1", L4 = "Day 0", Day3 = "Day 3", Day6 = "Day 6")[s]
}

rows <- list()
for (ct in contrasts) {
  de <- read_csv(file.path(opt$de_dir, paste0(ct, ".csv")), show_col_types = FALSE)
  bg <- de$gene_id[!is.na(de$padj)]
  for (dir_ in c("de", "down", "up")) {
    f <- file.path(opt$deg_dir, sprintf("%s_genes_%s.txt", ct, dir_))
    if (!file.exists(f)) next
    dg <- intersect(readLines(f), bg)
    for (ev in names(evidence)) {
      set <- intersect(evidence[[ev]], bg)
      hit <- intersect(dg, set)
      m <- matrix(c(length(hit), length(dg) - length(hit),
                    length(set) - length(hit),
                    length(bg) - length(dg) - length(set) + length(hit)), nrow = 2)
      ft <- tryCatch(fisher.test(m), error = function(e) list(p.value = NA, estimate = NA))
      rows[[length(rows) + 1]] <- tibble(
        stage = stage_of(ct), direction = dir_, evidence = ev,
        n_deg = length(dg), n_hit = length(hit),
        pct = 100 * length(hit) / max(length(dg), 1),
        bg_pct = 100 * length(set) / length(bg),
        odds_ratio = unname(ft$estimate), p = ft$p.value)
    }
  }
}
res <- bind_rows(rows); res$padj <- p.adjust(res$p, "BH")
write_tsv(res, file.path(opt$outdir, "binding_evidence.tsv"))
print(as.data.frame(res %>% filter(direction != "de") %>%
      select(stage, direction, evidence, n_deg, pct, bg_pct, odds_ratio, padj)), digits = 3)

lab <- c(de = "All DEGs", down = "Down in mutant", up = "Up in mutant")
res <- res %>% mutate(
  panel = factor(lab[direction], levels = unname(lab)),
  stage = factor(stage, levels = c("L1", "Day 0", "Day 3", "Day 6")),
  evidence = factor(evidence, levels = ev_levels),
  star = ifelse(is.na(padj) | padj >= 0.05, "", ifelse(padj < 0.001, "***",
          ifelse(padj < 0.01, "**", "*"))))
bg <- res %>% distinct(panel, stage, evidence, bg_pct)

pal <- c("ChIP-Atlas peak" = "#1F77B4", "E-box in promoter" = "#E66101",
         "E-box + DAF-16/PQM-1" = "#7B3294")

p <- ggplot(res, aes(stage, pct, fill = evidence)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.72) +
  geom_errorbar(data = bg, inherit.aes = FALSE,
                aes(x = stage, ymin = bg_pct, ymax = bg_pct, group = evidence),
                position = position_dodge(width = 0.8), width = 0.72,
                colour = "grey25", linewidth = 0.45, linetype = "dashed") +
  geom_text(aes(label = star, group = evidence), position = position_dodge(width = 0.8),
            vjust = -0.2, size = 3.4, colour = "grey20") +
  facet_wrap(~ panel) +
  scale_fill_manual(values = pal, name = NULL) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.12))) +
  labs(x = NULL, y = "DEGs with evidence (%)",
       title = sprintf("HLH-30 binding evidence among DEGs: measured vs predicted (ChIP +/-%s kb)", opt$dist),
       subtitle = "Dashes = expressed-gene background for that evidence class. Stars = BH-adj. Fisher vs background.") +
  theme_bw(base_size = 11) +
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank(),
        legend.position = "top", legend.justification = "left",
        strip.background = element_rect(fill = "grey96", colour = "grey80"),
        plot.title = element_text(size = 11, face = "bold"),
        plot.subtitle = element_text(size = 8.5, colour = "grey35"))

ggsave(file.path(opt$outdir, "binding_evidence.pdf"), p, width = 9.6, height = 4.3)
ggsave(file.path(opt$outdir, "binding_evidence.png"), p, width = 9.6, height = 4.3, dpi = 200)
message("wrote: ", file.path(opt$outdir, "binding_evidence.{pdf,png,tsv}"))
