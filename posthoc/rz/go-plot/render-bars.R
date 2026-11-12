#!/usr/bin/env Rscript
# Render filtered horizontal-bar GO plots from pipeline go_<ONT>.csv outputs.
#
# Reads targets.tsv (run, path, ontology, terms). For each row:
#   - downloads s3://<BUCKET>/results/<run>/enrichment/go/<path>/go_<ONT>.csv
#     to cache/<run>/<path>/go_<ONT>.csv (idempotent),
#   - filters rows by exact `Description` match against the terms list,
#   - applies optional "original => display" relabeling,
#   - emits a horizontal bar plot with x = -log10(p.adjust)
#     and gene count printed at the bar tip.

suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(stringr)
  library(ggforce)
})
source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                 "go-plot-common.R"))


opt <- parse_args(OptionParser(option_list = list(
  make_option("--config", type = "character"),
  make_option("--outdir", type = "character"),
  make_option("--bucket", type = "character",
              default = "s3://celegans-pipeline-rayandrew"),
  make_option("--cache",  type = "character", default = "cache"),
  make_option("--top-n",  type = "integer",   default = 20,
              dest = "top_n",
              help = "When terms column is empty, plot the top N most significant terms")
)))

stopifnot(!is.null(opt$config), !is.null(opt$outdir))
dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(opt$cache,  showWarnings = FALSE, recursive = TRUE)

targets <- read_tsv(opt$config, show_col_types = FALSE)


render_one <- function(row) {
  csv <- ensure_csv(row$run, row$path, row$ontology, opt$bucket, opt$cache)
  if (is.null(csv)) return(invisible(NULL))

  go  <- read_csv(csv, show_col_types = FALSE)
  hit <- select_terms(go, row$terms, opt$top_n,
                      sprintf("%s | %s | %s", row$run, row$path, row$ontology))

  hit <- hit %>%
    mutate(neglog10q = -log10(p.adjust),
           display_wrapped = str_wrap(display, width = 38)) %>%
    arrange(neglog10q) %>%
    mutate(display_wrapped = factor(display_wrapped, levels = unique(display_wrapped)))

  if (nrow(hit) == 0) {
    warning(sprintf("[%s | %s | %s] no terms matched, skipping figure",
                    row$run, row$path, row$ontology))
    return(invisible(NULL))
  }

  # Bars are polygons rather than geom_col so the data-end can carry the app's
  # corner radius; y therefore has to be numeric with the labels applied back.
  hit$ypos <- as.integer(hit$display_wrapped)
  bar_h <- 0.62
  bars <- do.call(rbind, lapply(seq_len(nrow(hit)), function(i) {
    data.frame(x = c(0, hit$neglog10q[i], hit$neglog10q[i], 0),
               y = hit$ypos[i] + c(-1, -1, 1, 1) * bar_h / 2,
               grp = i)
  }))

  p <- ggplot() +
    geom_shape(data = bars, aes(x = x, y = y, group = grp),
               fill = WB$primary, colour = NA, radius = WB$radius_sm) +
    geom_text(data = hit, aes(x = neglog10q, y = ypos, label = paste0("n=", Count)),
              hjust = -0.2, size = 3.2, colour = WB$muted, family = WB$font) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.18))) +
    scale_y_continuous(breaks = hit$ypos, labels = levels(hit$display_wrapped),
                       expand = expansion(add = 0.55)) +
    labs(x = expression(-log[10] ~ "(adj. " * italic(p) * ")"), y = NULL) +
    theme_minimal(base_size = 12, base_family = WB$font) +
    theme(
      plot.background    = element_rect(fill = WB$background, colour = NA),
      panel.background   = element_rect(fill = WB$background, colour = NA),
      panel.grid.major.y = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_line(colour = WB$border, linewidth = 0.4),
      axis.text.y        = element_text(colour = WB$foreground, lineheight = 0.9),
      axis.text.x        = element_text(colour = WB$muted),
      axis.title.x       = element_text(colour = WB$muted, size = 10),
      axis.ticks.x       = element_line(colour = WB$border),
      plot.margin        = margin(8, 12, 8, 8)
    )

  out_dir  <- file.path(opt$outdir, row$run, row$path)
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  base <- file.path(out_dir, sprintf("go_%s_bar", row$ontology))
  # height scales with row count + label wrapping
  n_rows <- nrow(hit)
  height <- max(2.2, 0.55 * n_rows + 1.4)
  ggsave(paste0(base, ".pdf"), p, width = 6.5, height = height)
  ggsave(paste0(base, ".png"), p, width = 6.5, height = height, dpi = 200)
  message("[wrote] ", base, ".{pdf,png}  (", n_rows, " terms)")
}

for (i in seq_len(nrow(targets))) {
  row <- targets[i, ]
  tryCatch(render_one(row),
           error = function(e) warning(sprintf("[row %d] %s", i, conditionMessage(e))))
}
