#!/usr/bin/env Rscript
# Two extras for the "% direct HLH-30 targets" analysis, both reusing the cached
# genome-wide E-box scan (no re-scan):
#   (1) compare-DE across all developmental transitions -> a compact dot-grid
#       (category x transition), so the L1->L4 view extends to L4->D3/D6, D3->D6.
#   (2) a per-GO-term lookup table: for EVERY enriched GO term (k-means clusters
#       + compare-DE), the % of its genes that are direct HLH-30 targets.
# Enrichment via the hypergeometric tail (vectorised) so thousands of terms are
# fast. Background = the same expressed-gene E-box rate as the DEG dumbbell.

suppressPackageStartupMessages({
  library(optparse); library(dplyr); library(readr); library(tidyr); library(stringr)
  library(ggplot2); library(RSQLite)
})

opt <- parse_args(OptionParser(option_list = list(
  make_option("--scan", type = "character",
              default = "cache/has_cayrtg_MA1449.1_u2000_d200_p0.0001.tsv"),
  make_option("--gene-db", type = "character", dest = "gene_db", default = "../../../gene_db.sqlite"),
  make_option("--cde-l1l4", type = "character", dest = "cde_l1l4",
              default = "../compare-de/out/hlh30-all-ruv"),
  make_option("--cde-other", type = "character", dest = "cde_other",
              default = "/tmp/claude-1000/-home-rayandrew-code-uwisc-celegans-pipeline/c3af1e84-db37-4245-9648-54f2c6e121b8/scratchpad/cde_other"),
  make_option("--go-kmeans", type = "character", dest = "go_kmeans",
              default = "~/downloads/hlh30-all-ruv_go_kmeans"),
  make_option("--go-cde", type = "character", dest = "go_cde",
              default = "/tmp/claude-1000/-home-rayandrew-code-uwisc-celegans-pipeline/c3af1e84-db37-4245-9648-54f2c6e121b8/scratchpad/cde_go"),
  make_option("--outdir", type = "character", default = "out/hlh30-all-ruv"),
  make_option("--accent", type = "character", default = "#7B3294")
)))
ACC <- opt$accent; BLU <- "#3B6E8F"
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

# --- cached scan: background + per-gene E-box ---------------------------------
hm <- read_tsv(opt$scan, show_col_types = FALSE)
uni <- hm$gene_id; N <- length(uni); N_hit <- sum(hm$has_motif); bg <- 100 * N_hit / N
has_motif <- setNames(hm$has_motif, hm$gene_id)
message(sprintf("background: %d/%d = %.1f%% expressed genes carry ~CAYRTG", N_hit, N, bg))

con <- dbConnect(SQLite(), opt$gene_db); on.exit(dbDisconnect(con))
ann <- dbGetQuery(con, "SELECT wormbase_gene, gene_name, sequence_name FROM annotations")
sym2wb <- c(setNames(ann$wormbase_gene, ann$gene_name), setNames(ann$wormbase_gene, ann$sequence_name))
sym2wb <- sym2wb[!is.na(names(sym2wb)) & nzchar(names(sym2wb))]

# hypergeometric over-representation of E-box among a gene set vs the universe
enrich <- function(genes) {
  g <- intersect(unique(genes), uni); n <- length(g); k <- sum(has_motif[g])
  p <- if (n == 0) NA_real_ else phyper(k - 1, N_hit, N - N_hit, n, lower.tail = FALSE)
  data.frame(n = n, n_target = k, pct = ifelse(n > 0, 100 * k / n, NA),
             fold = ifelse(n > 0, (k / n) / (N_hit / N), NA), p = p)
}
star_of <- function(padj) ifelse(is.na(padj), "", ifelse(padj < 0.001, "***",
                          ifelse(padj < 0.01, "**", ifelse(padj < 0.05, "*", ""))))

# ============================================================ (1) compare-DE ==
STG <- c(L1 = "L1", L4 = "L4", Day3 = "D3", Day6 = "D6")
parse_cde <- function(dir) {
  fs <- list.files(dir, pattern = "_genes_.*\\.txt$", full.names = TRUE)
  bind_rows(lapply(fs, function(f) {
    b <- basename(f)
    comp <- sub("__X__.*", "", b)                       # mutant-side, e.g. hlh-30_Day3_vs_hlh-30_L4
    m <- str_match(comp, "_(L1|L4|Day3|Day6)_vs_.*_(L1|L4|Day3|Day6)$")
    trans <- paste0(STG[m[3]], " -> ", STG[m[2]])       # T1 -> T2
    cat_raw <- sub(".*_genes_", "", sub("\\.txt$", "", b))
    if (grepl("exclusive_(up|down)$", cat_raw)) {        # clean exclusive_up / _down only
      side <- if (grepl("hlh-30|tm1978", cat_raw)) "hlh-30(tm1978)-only" else "WT-only"
      dir_ <- if (grepl("_up$", cat_raw)) "up" else "down"
      cat <- paste(side, dir_)
    } else if (cat_raw %in% c("shared", "differential")) {
      cat <- cat_raw
    } else return(NULL)                                 # skip higher_in / bare exclusive
    wb <- readLines(f); wb <- wb[nzchar(wb)]
    cbind(transition = trans, category = cat, enrich(wb))
  }))
}
cde <- bind_rows(parse_cde(opt$cde_l1l4), parse_cde(opt$cde_other))
cde$padj <- p.adjust(cde$p, "BH"); cde$star <- star_of(cde$padj)
write_tsv(cde, file.path(opt$outdir, "targets_compareDE_alltrans.tsv"))

cat_lv <- c("hlh-30(tm1978)-only up","hlh-30(tm1978)-only down",
            "WT-only up","WT-only down","shared","differential")
trans_lv <- c("L1 -> L4","L4 -> D3","L4 -> D6","D3 -> D6")
cde <- cde %>% mutate(category = factor(category, levels = rev(cat_lv)),
                      transition = factor(transition, levels = trans_lv),
                      log2fold = log2(fold))
lim <- max(abs(cde$log2fold), na.rm = TRUE)
pg <- ggplot(cde, aes(transition, category)) +
  geom_point(aes(fill = log2fold, size = n), shape = 21, colour = "white", stroke = 0.5) +
  geom_text(aes(label = sprintf("%.0f%%%s", pct, star)), size = 2.5, colour = "grey15", nudge_y = 0.34) +
  scale_fill_gradient2(low = BLU, mid = "white", high = ACC, midpoint = 0, limits = c(-lim, lim),
                       name = expression(log[2]~"fold vs bg")) +
  scale_size_area(max_size = 11, name = "genes in set") +
  labs(x = "developmental transition", y = NULL,
       title = "Direct HLH-30 targets among compare-DE lists, all transitions",
       subtitle = sprintf("%% of each category's genes with a ~CAYRTG E-box; fill = fold vs the %.1f%% background.", bg),
       caption = "L1->L4 from the L1-included run (WT/tm1978); L4->D3, L4->D6, D3->D6 from the VZRJXQ run (N2/hlh-30). Stars = BH-adjusted hypergeometric.") +
  theme_minimal(base_size = 11) +
  theme(panel.grid.major = element_line(colour = "grey93"), panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 9), plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 8.5, colour = "grey40"),
        plot.caption = element_text(size = 7, colour = "grey55", hjust = 0))
ggsave(file.path(opt$outdir, "targets_compareDE_alltrans.pdf"), pg, width = 8.2, height = 4.6, device = cairo_pdf)
ggsave(file.path(opt$outdir, "targets_compareDE_alltrans.png"), pg, width = 8.2, height = 4.6, dpi = 200)
message("wrote targets_compareDE_alltrans.{pdf,png,tsv}")

# ============================================================ (2) GO term table
go_table <- function(dir, source_prefix) {
  fs <- list.files(path.expand(dir), pattern = "go_(BP|CC|MF)\\.csv$", recursive = TRUE, full.names = TRUE)
  bind_rows(lapply(fs, function(f) {
    go <- suppressWarnings(read_csv(f, show_col_types = FALSE))
    if (!all(c("Description","geneID","p.adjust") %in% names(go)) || !nrow(go)) return(NULL)
    listname <- basename(dirname(f)); ont <- str_match(basename(f), "go_(..)")[, 2]
    stat <- bind_rows(lapply(strsplit(go$geneID, "/"), function(sy)
      enrich(unique(na.omit(sym2wb[sy])))))
    tibble(source = paste0(source_prefix, "/", listname), ontology = ont,
           GO_ID = if ("ID" %in% names(go)) go$ID else NA, term = go$Description,
           go_padj = go$p.adjust, n_genes = stat$n, n_ebox = stat$n_target,
           pct_ebox = round(stat$pct, 1), fold = round(stat$fold, 2), ebox_p = stat$p)
  }))
}
gt <- bind_rows(go_table(opt$go_kmeans, "kmeans"), go_table(opt$go_cde, "compareDE"))
gt <- gt %>% filter(!is.na(pct_ebox)) %>%
  mutate(ebox_padj = p.adjust(ebox_p, "BH"),
         ebox_p = signif(ebox_p, 3), ebox_padj = signif(ebox_padj, 3),
         go_padj = signif(go_padj, 3)) %>%
  arrange(desc(fold))
write_tsv(gt, file.path(opt$outdir, "go_term_target_pct.tsv"))
message(sprintf("wrote go_term_target_pct.tsv : %d GO-term rows (%d unique terms)",
                nrow(gt), length(unique(gt$term))))
