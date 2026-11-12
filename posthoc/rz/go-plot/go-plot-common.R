# Shared pieces for the GO figure renderers (render-bars.R, render-dots.R):
# design tokens ported from the wb.rs.ht app, the targets.tsv term parser, and
# the S3 fetch. Sourced, not run directly.

suppressPackageStartupMessages({
  library(dplyr); library(readr); library(stringr)
})

# --- design tokens, from the wb.rs.ht app (src/app.css, shadcn on a zinc base)
#   --primary 290 19% 46%   --foreground 240 10% 3.9%
#   --muted-foreground 240 3.8% 46.1%   --border 240 5.9% 90%
#   --radius 0.5rem  (sm = radius - 4px)
WB <- list(
  primary    = "#845F8C",   # brand purple
  foreground = "#09090B",
  muted      = "#71717A",
  border     = "#E4E4E7",
  background = "#FFFFFF",
  radius_sm  = grid::unit(2.5, "pt")
)
# Inter is the app's UI face; fall back to the default sans when unavailable.
WB$font <- if ("Inter" %in% systemfonts::system_fonts()$family) "Inter" else ""

# The app's dot-plot colour ramp: a 5-stop viridis approximation with evenly
# spaced stops (EnrichmentDotPlot.svelte, colorFor()).
WB$viridis5 <- c(
  rgb( 68,   1,  84, maxColorValue = 255),
  rgb( 59,  82, 139, maxColorValue = 255),
  rgb( 33, 145, 140, maxColorValue = 255),
  rgb( 94, 201,  98, maxColorValue = 255),
  rgb(253, 231,  37, maxColorValue = 255)
)

# Brand sequential ramp, light -> dark, all at the primary's hue (~290):
#   --accent (290 30% 95%) -> --primary (290 19% 46%) -> --accent-foreground
# Single-hue, so it works as a sequential scale; the app itself uses viridis.
WB$purple_seq <- c("#F5EEF6", "#C9AECF", "#845F8C", "#4C2D53")

# "original => display; other" -> data.frame(source, display)
parse_terms <- function(s) {
  parts <- str_trim(str_split(s, ";")[[1]])
  parts <- parts[nzchar(parts)]
  do.call(rbind, lapply(parts, function(p) {
    if (str_detect(p, "=>")) {
      kv <- str_trim(str_split(p, "=>", n = 2)[[1]])
      data.frame(source = kv[1], display = kv[2], stringsAsFactors = FALSE)
    } else {
      data.frame(source = p, display = p, stringsAsFactors = FALSE)
    }
  }))
}

# Fetch go_<ONT>.csv for one target into the cache (idempotent).
ensure_csv <- function(run, path, ont, bucket, cache) {
  local_dir <- file.path(cache, run, path)
  dir.create(local_dir, showWarnings = FALSE, recursive = TRUE)
  local_file <- file.path(local_dir, sprintf("go_%s.csv", ont))
  if (!file.exists(local_file)) {
    s3 <- sprintf("%s/results/%s/enrichment/go/%s/go_%s.csv", bucket, run, path, ont)
    message("[download] ", s3)
    rc <- system2("aws", c("s3", "cp", s3, local_file))
    if (rc != 0 || !file.exists(local_file)) {
      warning("download failed: ", s3); return(NULL)
    }
  }
  local_file
}

# Resolve targets.tsv rows to the term set actually plotted.
# A term prefixed with "!" in the terms cell is EXCLUDED: it is dropped from the
# plotted set (from the explicit include list, or from the top-N default when no
# includes are given). This lets a figure show "the usual top terms, minus one".
select_terms <- function(go, terms_cell, top_n, label) {
  terms_cell <- if (is.na(terms_cell)) "" else str_trim(terms_cell)
  excl <- character(0); incl <- NULL
  if (nzchar(terms_cell)) {
    spec <- parse_terms(terms_cell)
    is_excl <- str_detect(spec$source, "^!")
    excl <- str_trim(sub("^!", "", spec$source[is_excl]))
    if (any(!is_excl)) incl <- spec[!is_excl, , drop = FALSE]
  }
  if (is.null(incl)) {
    hit <- go %>% filter(!Description %in% excl) %>%
      arrange(p.adjust) %>% slice_head(n = top_n) %>% mutate(display = Description)
  } else {
    missing <- setdiff(incl$source, go$Description)
    if (length(missing))
      warning(sprintf("[%s] not found in go CSV: %s", label, paste(missing, collapse = "; ")))
    hit <- go %>% filter(Description %in% incl$source, !Description %in% excl) %>%
      left_join(incl, by = c("Description" = "source"))
  }
  hit
}

# --- paper labels (render-dots.R --label-in, full-tables.R) -------------------
# Greedy word wrap by measured string width, so every label fits a column of
# max_in inches at fontsize pt in the output font.
wrap_to <- function(x, max_in, fontsize) {
  grDevices::pdf(NULL); on.exit(grDevices::dev.off())
  gp <- grid::gpar(fontsize = fontsize, fontfamily = WB$font)
  w_in <- function(s) grid::convertWidth(grid::grobWidth(grid::textGrob(s, gp = gp)), "in", TRUE)
  vapply(x, function(s) {
    words <- strsplit(s, " ")[[1]]; lines <- character(0); cur <- ""
    for (w in words) {
      try <- if (nzchar(cur)) paste(cur, w) else w
      if (nzchar(cur) && w_in(try) > max_in) { lines <- c(lines, cur); cur <- w } else cur <- try
    }
    paste(c(lines, cur), collapse = "\n")
  }, character(1), USE.NAMES = FALSE)
}
# The full GO name, wrapped to max_in, when it needs at most max_lines lines;
# otherwise the short display name from targets.tsv, wrapped the same way.
paper_label <- function(full, short, max_in, fontsize, max_lines = 2) {
  wf <- wrap_to(full, max_in, fontsize); ws <- wrap_to(short, max_in, fontsize)
  lines <- function(x) stringr::str_count(x, "\n") + 1
  out <- ifelse(lines(wf) <= max_lines, wf, ws)
  over <- lines(out) > max_lines
  if (any(over)) warning("label over ", max_lines, " lines: ", paste(gsub("\n", " ", out[over]), collapse = "; "))
  out
}
# Label width and size of a paper panel placed placed_in wide, as render-dots.R
# sets them with --target-pt 10 --label-frac 0.457 on its 7 in canvas.
paper_label_for <- function(full, short, placed_in, width = 7, label_frac = 0.457, target_pt = 10) {
  fs <- ceiling(9.6 * target_pt / (9.6 * placed_in / width)) / 9.6
  paper_label(full, short, 0.92 * label_frac * width, 9.6 * fs)
}
