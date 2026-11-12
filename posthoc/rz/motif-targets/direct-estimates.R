#!/usr/bin/env Rscript
# Estimated share of direct HLH-30 targets per DEG list, from the motif scan
# (motif-targets.R) and from the L4 ChIP peaks (hlh30-integration bound-vs-dge),
# side by side. If every target carries the mark and the other genes carry it
# at the background rate: est = (list% - background%) / (100 - background%).

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(readr); library(tidyr) })

opt <- parse_args(OptionParser(option_list = list(
  make_option("--summary", type = "character", default = "out/hlh30-all-ruv/motif_targets_summary.tsv"),
  make_option("--chip",    type = "character", default = "../hlh30-integration/out/bound/fisher_bound_vs_deg.tsv"),
  make_option("--out",     type = "character", default = "out/hlh30-all-ruv/direct_target_estimates.csv")
)))

est  <- function(p, bg) pmax(0, 100 * (p - bg) / (100 - bg))
star <- function(q) ifelse(q < 0.001, "***", ifelse(q < 0.01, "**", ifelse(q < 0.05, "*", "ns")))
cell <- function(p, bg, q) sprintf("%.0f%% vs %.0f%%, %s (%s)", p, bg,
                                   ifelse(p < bg, "below chance", sprintf("est. %.0f%%", est(p, bg))), star(q))

stage_disp <- c(L1 = "L1", `Day 0` = "L4", `Day 3` = "D3", `Day 6` = "D6")
mot <- read_tsv(opt$summary, show_col_types = FALSE) %>%
  filter(family == "deg", direction != "all") %>%
  transmute(stage = stage_disp[group], direction, n, motif,
            v = cell(pct, bg_pct, padj)) %>%
  pivot_wider(names_from = motif, values_from = v)

chip <- read_tsv(opt$chip, show_col_types = FALSE) %>%
  transmute(stage, direction, `ChIP peak (L4)` = cell(pct_bound_deg, pct_bound_all, q))

tab <- mot %>% left_join(chip, by = c("stage", "direction")) %>%
  transmute(Stage = factor(stage, levels = c("L1", "L4", "D3", "D6")), Direction = direction,
            `DEGs` = n,
            `Any HLH-30 motif` = any, `Canonical E-box (JASPAR, fig S5C)` = hlh30, `ChIP-seq motif` = chip,
            `ChIP peak (L4)`) %>%
  arrange(Stage, desc(Direction))
write_csv(tab, opt$out)
print(as.data.frame(tab), right = FALSE)
message("wrote ", opt$out)
