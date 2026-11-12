#!/usr/bin/env Rscript
# Gene-list panels that the author first drew by hand (Fig 1F, Fig S1A), now
# made from the DESeq2 tables so the numbers and the type match the other
# panels. Each gene is a row of three cells, its shrunken log2FC in
# hlh-30(tm1978) vs WT at L4, D3 and D6, with the value printed in the cell.
# A cell that is not a DEG at that stage is drawn pale. Up and down genes are
# separate blocks, each sorted by the mean log2FC.
#
#   fig-1f   DEGs at all three stages
#   fig-s1a  DEGs at D3 and D6 but not at L4
#
# DEG = padj < 0.05 and |log2FC| >= 1, as in the volcanoes (Fig 1E).

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(ggplot2); library(patchwork) })
here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1]))
source(file.path(here, "..", "_common", "palette.R"))

opt <- parse_args(OptionParser(option_list = list(
  make_option("--de-dir", type = "character", default = "../volcano/cache/VZRJXQ-ruv-rbe", dest = "de_dir"),
  make_option("--gene-db", type = "character", default = "../../../gene_db.sqlite", dest = "gene_db",
              help = "repo gene database, for names of genes the DESeq2 table leaves unnamed"),
  make_option("--padj",   type = "double", default = 0.05),
  make_option("--lfc",    type = "double", default = 1),
  make_option("--width-1f",  type = "double", default = 4.2, dest = "width_1f",
              help = "fig 1F width (in); the panel is laid out 1:1 at 10pt"),
  make_option("--width-s1a", type = "double", default = 4.2, dest = "width_s1a",
              help = "fig S1A width (in)"),
  make_option("--outdir", type = "character", default = "out")
)))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
STAGES <- c(L4 = "L4", D3 = "Day3", D6 = "Day6")

d <- bind_rows(lapply(names(STAGES), function(st) {
  read.csv(file.path(opt$de_dir, sprintf("hlh-30_%s_vs_N2_%s.csv", STAGES[[st]], STAGES[[st]]))) %>%
    transmute(gene_id, gene_name, stage = st, lfc = log2FoldChange,
              deg = !is.na(padj) & padj < opt$padj & abs(log2FoldChange) >= opt$lfc)
})) %>% mutate(stage = factor(stage, levels = names(STAGES)))

# genes with no CGC name carry their WormBase id; use the WormBase name from
# the repo's gene database (it covers pseudogenes the GTF leaves unnamed)
if (file.exists(opt$gene_db)) {
  con <- DBI::dbConnect(RSQLite::SQLite(), opt$gene_db)
  ann <- DBI::dbGetQuery(con, "SELECT wormbase_gene, gene_name FROM annotations WHERE gene_name <> ''")
  DBI::dbDisconnect(con)
  nm <- setNames(ann$gene_name, ann$wormbase_gene)
  fix <- d$gene_name == d$gene_id & d$gene_id %in% names(nm)
  d$gene_name[fix] <- nm[d$gene_id[fix]]
}

degs_at <- function(st) unique(d$gene_id[d$stage == st & d$deg])
SETS <- list(
  "fig-1f"  = list(genes = Reduce(intersect, lapply(names(STAGES), degs_at)),
                   title = "DEGs at L4, D3 and D6"),
  "fig-s1a" = list(genes = setdiff(intersect(degs_at("D3"), degs_at("D6")), degs_at("L4")),
                   title = "DEGs at D3 and D6, not at L4")
)

LIM <- 4  # colour scale cap; the printed value is exact
block <- function(g, dir) {
  x <- d %>% filter(gene_id %in% g)
  # direction from the D6 change, the stage every list shares
  keep <- x %>% filter(stage == "D6") %>% filter(if (dir == "up") lfc > 0 else lfc < 0) %>% pull(gene_id)
  x <- x %>% filter(gene_id %in% keep)
  if (!nrow(x)) return(NULL)
  ord <- x %>% group_by(gene_id, gene_name) %>% summarise(m = mean(lfc), .groups = "drop") %>%
    arrange(if (dir == "up") -m else m)
  x <- x %>% mutate(row = factor(gene_id, levels = rev(ord$gene_id)),
                    fillv = pmax(pmin(lfc, LIM), -LIM))
  ggplot(x, aes(stage, row)) +
    geom_tile(aes(fill = fillv, alpha = deg), colour = "white", linewidth = 0.8) +
    # "-0.0" reads as a sign error, so round first
    geom_text(aes(label = sprintf("%.1f", round(lfc, 1) + 0), colour = abs(fillv) > 2.5 & deg), size = 10 / .pt) +
    scale_fill_gradient2(low = PAL$div[["low"]], mid = PAL$div[["mid"]], high = PAL$div[["high"]],
                         midpoint = 0, limits = c(-LIM, LIM), guide = "none") +
    scale_alpha_manual(values = c(`TRUE` = 1, `FALSE` = 0.35), guide = "none") +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "grey15"), guide = "none") +
    scale_y_discrete(labels = setNames(ord$gene_name, ord$gene_id)) +
    scale_x_discrete(position = "top") +
    labs(x = NULL, y = NULL,
         title = bquote(.(if (dir == "up") "up" else "down") ~ "in" ~ italic("hlh-30(tm1978)") ~ "(" * .(nrow(ord)) * ")")) +
    theme_minimal(base_size = 10) +
    theme(panel.grid = element_blank(),
          axis.text.y = element_text(face = "italic", colour = "grey10", size = 10),
          axis.text.x.top = element_text(colour = "grey10", size = 10),
          plot.title.position = "plot",
          plot.title = element_text(size = 10, face = "bold", hjust = 0,
                                    colour = if (dir == "up") PAL$up else PAL$down),
          plot.margin = margin(2, 6, 2, 2))
}

for (id in names(SETS)) {
  s <- SETS[[id]]
  up <- block(s$genes, "up"); dn <- block(s$genes, "down")
  n_up <- sum(d$stage == "D6" & d$gene_id %in% s$genes & d$lfc > 0); n_dn <- length(s$genes) - n_up
  # side by side; each block spans as many layout rows as it has genes, so
  # every gene row has the same height in both blocks
  rows <- max(n_up, n_dn) + 2  # + title and stage header
  blocks <- Filter(Negate(is.null), list(up, dn)); nr <- c(n_up, n_dn)[!sapply(list(up, dn), is.null)]
  design <- do.call(c, lapply(seq_along(blocks), function(i) area(1, i, nr[i] + 2, i)))
  # fig-1f holds only genes that are DEGs at every stage, so it has no pale cell
  any_pale <- any(!d$deg[d$gene_id %in% s$genes])
  caption <- paste0("log2FC, *hlh-30(tm1978)* vs WT", if (any_pale) "<br>pale = not a DEG at that stage")
  fig <- wrap_plots(blocks) + plot_layout(design = design, heights = rep(1, rows)) +
    plot_annotation(caption = caption,
                    theme = theme(plot.caption = ggtext::element_markdown(size = 10, colour = "grey35",
                                                                          hjust = 0, lineheight = 1.2)))
  h <- 0.6 + 0.26 * rows   # 0.26 in per gene row
  w <- if (id == "fig-s1a") opt$width_s1a else opt$width_1f
  ggsave(file.path(opt$outdir, paste0(id, ".pdf")), fig, width = w, height = h, device = cairo_pdf)
  ggsave(file.path(opt$outdir, paste0(id, ".png")), fig, width = w, height = h, dpi = 200)
  write.csv(d %>% filter(gene_id %in% s$genes) %>%
              select(gene_id, gene_name, stage, lfc, deg) %>%
              tidyr::pivot_wider(names_from = stage, values_from = c(lfc, deg)),
            file.path(opt$outdir, paste0(id, "_genes.csv")), row.names = FALSE)
  message("[wrote] ", id, ": ", length(s$genes), " genes (", n_up, " up, ", n_dn, " down)")
}
