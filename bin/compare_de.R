#!/usr/bin/env Rscript

# compare_de.R - Compare DE results across conditions/studies
#
# Input:  Directory with DESeq2 result CSVs
# Output: Comparison report, Venn-style overlaps, classified gene lists

suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(ggplot2)
})

opts <- list(
  make_option("--input_dir",      type = "character", default = ".", help = "Directory with DESeq2 result CSVs"),
  make_option("--padj_threshold", type = "numeric",   default = 0.05, help = "Adj. p-value cutoff"),
  make_option("--lfc_threshold",  type = "numeric",   default = 1.0,  help = "Min |log2FC| for DE"),
  make_option("--lfc_diff",       type = "numeric",   default = 1.0,  help = "Min |lfc_diff| for differential"),
  make_option("--age_order",      type = "character", default = NULL,
              help = "Comma-separated age tokens, oldest LAST (e.g. 'L4,Day3,Day6'). If omitted, auto-derived alphabetically from result filenames."),
  make_option("--geno_order",     type = "character", default = NULL,
              help = "Comma-separated genotype tokens, baseline FIRST (e.g. 'N2,hlh-30'). If omitted, auto-derived alphabetically from result filenames."),
  make_option("--pairs",          type = "character", default = NULL,
              help = "Pair allow-list: comma-separated 'A_vs_B__X__C_vs_D' specs. If set, ONLY these contrast-pairs are compared. Useful when DESeq2 emits N contrasts but you only care about specific cross-comparisons. Pair canonicalization is applied before matching, so order within each side doesn't matter (you can write the pair in either orientation)."),
  make_option("--include_terms",  type = "character", default = NULL,
              help = "Optional comma-separated allow-list of <genotype>_<age> terms. DESeq2 result files whose numerator AND denominator both appear in this list are kept; others are skipped. Useful for trimming the input set before pair iteration (e.g. drop L2D contrasts entirely)."),
  make_option("--diff_method",    type = "character", default = "overlap",
              help = "How to compare two contrasts: 'overlap' (default, legacy: classify by independent DE calls + |lfc_diff| threshold) or 'ztest' (Wald test on the LFC difference using unshrunken MLE LFC+SE; emits concordant/divergent/reversal)."),
  make_option("--canonicalize",   type = "logical",   default = TRUE,
              help = "If TRUE (default; pairwise / legacy mode), parse contrast filenames as '<genotype>_<age>_vs_<genotype>_<age>' and canonicalize orientation (older numerator, non-baseline genotype numerator). If FALSE (multi-factor mode), skip parsing/canonicalization and use the raw filename basename as the comparison label; useful when filenames are DESeq2-coefficient-style such as 'treatment__high_vs_none.csv' that don't have a genotype/age structure."),
  make_option("--tested_in_both", type = "logical",   default = FALSE,
              help = "overlap method only. If TRUE, compare only genes with a non-NA padj in BOTH contrasts, so a gene is not 'exclusive' to one contrast just because DESeq2 filtered it out of the other (independent filtering, outliers). Default FALSE keeps the legacy full outer join."),
  make_option("--outdir",         type = "character",  default = ".",   help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)

# Inter-comparison separator (distinguishes from intra "_vs_")
SEP <- "__X__"

CANONICALIZE_ON <- isTRUE(args$canonicalize)
DIFF_METHOD     <- match.arg(args$diff_method, c("overlap", "ztest"))

# Grouped DESeq2 runs prefix contrast filenames with '<group>__'. Strip it so
# pairwise parsing/canonicalization sees the bare '<geno>_<age>_vs_<geno>_<age>'.
strip_group_prefix <- function(label) {
  if (CANONICALIZE_ON && grepl("__", label, fixed = TRUE))
    sub("^.+?__", "", label, perl = TRUE)
  else
    label
}

if (CANONICALIZE_ON) {
  # ---- LEGACY / PAIRWISE MODE ----
  #
  # parse_term: split a `<genotype>_<age>` token on the LAST underscore.
  #   "N2_L4"          -> genotype=N2,      age=L4
  #   "hlh-30_Day3"    -> genotype=hlh-30,  age=Day3
  # Convention: genotype may contain hyphens but no underscore. Age token is
  # the substring after the final underscore.
  parse_term <- function(t) {
    idx <- regexpr("_[^_]+$", t)
    if (idx < 1) return(NULL)
    list(genotype = substr(t, 1, idx - 1),
         age      = substr(t, idx + 1, nchar(t)))
  }

  parse_contrast <- function(label) {
    parts <- strsplit(label, "_vs_", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NULL)
    num <- parse_term(parts[1])
    den <- parse_term(parts[2])
    if (is.null(num) || is.null(den)) return(NULL)
    list(num = num, den = den)
  }

  term_label <- function(t) paste0(t$genotype, "_", t$age)
  contrast_label <- function(c) paste0(term_label(c$num), "_vs_", term_label(c$den))

  canonicalize <- function(c) {
    flipped <- FALSE
    if (c$num$genotype == c$den$genotype) {
      # age contrast: numerator should be older
      if (AGE_ORDER[[c$num$age]] < AGE_ORDER[[c$den$age]]) {
        tmp <- c$num; c$num <- c$den; c$den <- tmp
        flipped <- TRUE
      }
    } else if (c$num$age == c$den$age) {
      # genotype contrast: numerator should be hlh-30 (N2 = baseline)
      if (GENO_ORDER[[c$num$genotype]] < GENO_ORDER[[c$den$genotype]]) {
        tmp <- c$num; c$num <- c$den; c$den <- tmp
        flipped <- TRUE
      }
    }
    c$flipped <- flipped
    c
  }
} else {
  # ---- MULTI-FACTOR / RAW MODE ----
  #
  # Filenames are DESeq2-coefficient style: `<factor>__<level>_vs_<reference>`
  # (e.g. `treatment__high_vs_none`). No genotype/age decomposition; no
  # orientation flipping. The double-underscore separates the factor from
  # the level pair so we can group / filter by factor name.
  message("compare_de: canonicalize=FALSE (multi-factor / raw mode)")

  parse_contrast <- function(label) {
    # Detect double-underscore (multi-factor coefficient): "<factor>__<lvl>_vs_<ref>"
    dunder <- strsplit(label, "__", fixed = TRUE)[[1]]
    factor_name <- NA_character_
    body_label  <- label
    if (length(dunder) == 2) {
      factor_name <- dunder[1]
      body_label  <- dunder[2]
    }
    # Split level vs reference on _vs_
    parts <- strsplit(body_label, "_vs_", fixed = TRUE)[[1]]
    if (length(parts) != 2) return(NULL)
    list(
      raw    = label,
      factor = factor_name,
      num    = list(level = parts[1]),
      den    = list(level = parts[2])
    )
  }

  term_label     <- function(t) t$level
  contrast_label <- function(c) c$raw

  canonicalize <- function(c) {
    # No reorientation in raw mode -- the DESeq2 coefficient direction is the
    # canonical one (numerator is the non-reference level).
    c$flipped <- FALSE
    c
  }
}

# Find result files
result_files <- list.files(args$input_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(result_files) < 2) {
  message("Only ", length(result_files), " result file(s) - need at least 2 to compare. Skipping.")
  writeLines("Need at least 2 comparisons", file.path(args$outdir, "compare_de_report.csv"))
  quit(save = "no", status = 0)
}

# Pre-pass: derive AGE_ORDER and GENO_ORDER from the parseable files (or use
# user-supplied --age_order / --geno_order). Conventions:
#   - AGE_ORDER: oldest LAST (so a contrast num=older,den=younger satisfies
#                AGE_ORDER[num] > AGE_ORDER[den])
#   - GENO_ORDER: baseline FIRST (so num != baseline ranks higher)
parse_csv_list <- function(s) {
  if (is.null(s) || !nzchar(s)) return(NULL)
  trimws(strsplit(s, ",", fixed = TRUE)[[1]])
}
build_order <- function(values, override = NULL, sort_alpha = TRUE) {
  if (!is.null(override)) {
    missing <- setdiff(values, override)
    if (length(missing) > 0)
      warning("These tokens were in input but missing from order spec (will sort to end alphabetically): ",
              paste(missing, collapse = ", "))
    ord <- c(override, sort(missing))
  } else {
    ord <- if (sort_alpha) sort(unique(values)) else unique(values)
  }
  setNames(seq_along(ord), ord)
}

if (CANONICALIZE_ON) {
  discovered_genos <- character()
  discovered_ages  <- character()
  for (f in result_files) {
    raw_label <- strip_group_prefix(tools::file_path_sans_ext(basename(f)))
    parsed <- parse_contrast(raw_label)
    if (is.null(parsed)) next
    discovered_genos <- c(discovered_genos, parsed$num$genotype, parsed$den$genotype)
    discovered_ages  <- c(discovered_ages,  parsed$num$age,      parsed$den$age)
  }
  AGE_ORDER  <- build_order(discovered_ages,  parse_csv_list(args$age_order))
  GENO_ORDER <- build_order(discovered_genos, parse_csv_list(args$geno_order))
  message("AGE_ORDER:  ", paste(names(AGE_ORDER),  collapse = " < "))
  message("GENO_ORDER: ", paste(names(GENO_ORDER), collapse = " < "), "  (baseline first)")
} else {
  # Raw mode: no ordering needed (parse_contrast doesn't rely on AGE/GENO_ORDER)
  AGE_ORDER  <- setNames(integer(0), character(0))
  GENO_ORDER <- setNames(integer(0), character(0))
}

# Optional term allow-list (drop result files whose contrast touches a term not in the list)
include_terms_set <- if (!is.null(args$include_terms))
  trimws(strsplit(args$include_terms, ",", fixed = TRUE)[[1]]) else NULL

# Read all results, canonicalizing orientation (flip LFC sign in-memory if needed)
results_list  <- list()
contrasts     <- list()
for (f in result_files) {
  raw_label <- strip_group_prefix(tools::file_path_sans_ext(basename(f)))
  parsed <- parse_contrast(raw_label)
  if (is.null(parsed)) {
    message("  Skipping unparseable file: ", basename(f))
    next
  }
  canon <- canonicalize(parsed)
  canon_label <- contrast_label(canon)

  if (!is.null(include_terms_set)) {
    if (CANONICALIZE_ON) {
      terms_here <- c(term_label(canon$num), term_label(canon$den))
      if (!all(terms_here %in% include_terms_set)) {
        message("  Skipping (not in --include_terms): ", raw_label)
        next
      }
    } else {
      # In raw mode --include_terms is interpreted as factor-name allow-list
      # (e.g. include only the treatment and genotype coefficients).
      if (!is.na(canon$factor) && !(canon$factor %in% include_terms_set)) {
        message("  Skipping (factor not in --include_terms): ", raw_label)
        next
      }
    }
  }

  df <- read.csv(f, stringsAsFactors = FALSE)
  if (canon$flipped) {
    message("  Flipping orientation: ", raw_label, " -> ", canon_label)
    df$log2FoldChange <- -df$log2FoldChange
    if ("stat" %in% names(df))               df$stat               <- -df$stat
    if ("log2FoldChange_MLE" %in% names(df)) df$log2FoldChange_MLE <- -df$log2FoldChange_MLE
  }
  df$comparison <- canon_label
  results_list[[canon_label]] <- df
  contrasts[[canon_label]]    <- canon
}

comparisons <- names(results_list)
message("Comparing ", length(comparisons), " result sets (after canonicalization)")

# Optional pair allow-list. Specs look like "A_vs_B__X__C_vs_D"; we canonicalize
# both sides AND the pair (sorted by canonical contrast labels) so the user can
# write the pair in either orientation/order.
canon_pair_key <- function(label_a, label_b) {
  ks <- sort(c(label_a, label_b))
  paste(ks[1], ks[2], sep = SEP)
}
pair_allow_set <- NULL
if (!is.null(args$pairs)) {
  raw_pairs <- trimws(strsplit(args$pairs, ",", fixed = TRUE)[[1]])
  for (rp in raw_pairs) {
    parts <- strsplit(rp, SEP, fixed = TRUE)[[1]]
    if (length(parts) != 2) {
      warning("--pairs entry not in 'A_vs_B__X__C_vs_D' form, skipping: ", rp)
      next
    }
    ca <- parse_contrast(parts[1]); cb <- parse_contrast(parts[2])
    if (is.null(ca) || is.null(cb)) {
      warning("--pairs entry has unparseable side, skipping: ", rp); next
    }
    la <- contrast_label(canonicalize(ca))
    lb <- contrast_label(canonicalize(cb))
    pair_allow_set <- c(pair_allow_set, canon_pair_key(la, lb))
  }
  message("Pair allow-list (", length(pair_allow_set), " entries): ",
          paste(pair_allow_set, collapse = "; "))
}

if (CANONICALIZE_ON) {
  # PAIRWISE MODE: two canonical contrasts are comparable iff:
  #   (a) they share exactly one literal term (e.g. both anchored at N2_Day3), OR
  #   (b) parallel age axis: both are age contrasts on the same age pair, different genotypes
  #       (e.g. N2_Day3_vs_N2_L4 vs hlh-30_Day3_vs_hlh-30_L4), OR
  #   (c) parallel genotype axis: both are genotype contrasts on the same genotype pair, different ages
  #       (e.g. hlh-30_L4_vs_N2_L4 vs hlh-30_Day3_vs_N2_Day3)
  shares_axis <- function(ca, cb) {
    ta <- c(term_label(ca$num), term_label(ca$den))
    tb <- c(term_label(cb$num), term_label(cb$den))
    if (length(intersect(ta, tb)) == 1) return(TRUE)

    is_age_a  <- ca$num$genotype == ca$den$genotype
    is_age_b  <- cb$num$genotype == cb$den$genotype
    if (is_age_a && is_age_b &&
        ca$num$age      == cb$num$age &&
        ca$den$age      == cb$den$age &&
        ca$num$genotype != cb$num$genotype) return(TRUE)

    is_geno_a <- ca$num$age == ca$den$age
    is_geno_b <- cb$num$age == cb$den$age
    if (is_geno_a && is_geno_b &&
        ca$num$genotype == cb$num$genotype &&
        ca$den$genotype == cb$den$genotype &&
        ca$num$age      != cb$num$age) return(TRUE)

    FALSE
  }
} else {
  # MULTI-FACTOR / RAW MODE: any two coefficient comparisons are eligible
  # for cross-comparison. Two cases the user typically cares about:
  #   - WITHIN-FACTOR: e.g., treatment__high_vs_none vs treatment__low_vs_none
  #     ("is the high-dose effect similar to the low-dose effect?")
  #   - CROSS-FACTOR:  e.g., treatment__high_vs_none vs genotype__mut_vs_WT
  #     ("is the treatment effect similar to the genotype effect?")
  # Both are biologically meaningful. shares_axis returns TRUE for everything
  # here, so the user controls which pairs to actually compare via --pairs.
  # For small coefficient sets (typical 3-10 in multi-factor mode), N^2
  # pairwise is manageable; for larger sets the user should set --pairs.
  shares_axis <- function(ca, cb) TRUE
}

# Pairwise comparisons
all_comparisons <- list()

for (i in 1:(length(comparisons) - 1)) {
  for (j in (i + 1):length(comparisons)) {
    label_a <- comparisons[i]
    label_b <- comparisons[j]

    if (!is.null(pair_allow_set)) {
      if (!(canon_pair_key(label_a, label_b) %in% pair_allow_set)) {
        next
      }
    } else if (!shares_axis(contrasts[[label_a]], contrasts[[label_b]])) {
      message("  Skipping (no shared axis): ", label_a, " ", SEP, " ", label_b)
      next
    }

    res_a <- results_list[[label_a]]
    res_b <- results_list[[label_b]]

    if (DIFF_METHOD == "ztest") {
      pair_label <- paste(label_a, label_b, sep = SEP)
      need <- c("gene_id", "log2FoldChange", "padj", "log2FoldChange_MLE", "lfcSE_MLE")
      for (nm in c(label_a, label_b)) {
        miss <- setdiff(need, names(results_list[[nm]]))
        if (length(miss) > 0)
          stop("diff_method=ztest needs column(s) ", paste(miss, collapse = ", "),
               " in '", nm, "' - regenerate DESeq2 results with the updated deseq2.R")
      }
      m <- merge(res_a[, need], res_b[, need], by = "gene_id", suffixes = c(".a", ".b"))

      # Wald test on the LFC difference, using unshrunken MLE estimates
      m$lfc_diff    <- m$log2FoldChange_MLE.a - m$log2FoldChange_MLE.b
      se_diff       <- sqrt(m$lfcSE_MLE.a^2 + m$lfcSE_MLE.b^2)
      m$z           <- m$lfc_diff / se_diff
      m$pvalue_diff <- 2 * pnorm(-abs(m$z))
      m$padj_diff   <- p.adjust(m$pvalue_diff, method = "BH")

      de_a     <- !is.na(m$padj.a) & m$padj.a < args$padj_threshold
      de_b     <- !is.na(m$padj.b) & m$padj.b < args$padj_threshold
      diff_sig <- !is.na(m$padj_diff) & m$padj_diff < args$padj_threshold
      opp_sign <- sign(m$log2FoldChange_MLE.a) != sign(m$log2FoldChange_MLE.b)

      m$category <- case_when(
        diff_sig & de_a & de_b & opp_sign ~ "reversal",
        diff_sig                          ~ "divergent",
        de_a & de_b                       ~ "concordant",
        TRUE                              ~ "unchanged"
      )
      m$comparison_pair <- pair_label

      out <- m[, c("gene_id", "log2FoldChange.a", "padj.a", "log2FoldChange.b", "padj.b",
                   "lfc_diff", "z", "pvalue_diff", "padj_diff", "category", "comparison_pair")]
      names(out)[2:5] <- c(paste0("lfc.", label_a),  paste0("padj.", label_a),
                           paste0("lfc.", label_b),  paste0("padj.", label_b))
      out <- out[out$category != "unchanged", ]
      all_comparisons[[pair_label]] <- out

      write.csv(out, file.path(args$outdir, paste0(pair_label, "_compared.csv")), row.names = FALSE)
      # Per-category gene lists - consumed by the GO/KEGG/WormCat _COMPARE modules
      for (cat in c("concordant", "divergent", "reversal")) {
        g <- out$gene_id[out$category == cat]
        if (length(g) > 0)
          writeLines(g, file.path(args$outdir, paste0(pair_label, "_genes_", cat, ".txt")))
      }
      summ <- out %>% count(category) %>% mutate(line = paste0(category, ": ", n))
      writeLines(summ$line, file.path(args$outdir, paste0(pair_label, "_summary.txt")))
      message("  ", label_a, " ", SEP, " ", label_b, ": ", paste(summ$line, collapse = ", "))
      next
    }

    # Full outer join on gene_id
    merged <- merge(
      res_a[, c("gene_id", "log2FoldChange", "padj")],
      res_b[, c("gene_id", "log2FoldChange", "padj")],
      by = "gene_id", all = TRUE,
      suffixes = c(paste0(".", label_a), paste0(".", label_b))
    )

    if (isTRUE(args$tested_in_both)) {
      tested <- !is.na(merged[[paste0("padj.", label_a)]]) & !is.na(merged[[paste0("padj.", label_b)]])
      message("  tested_in_both: dropping ", sum(!tested), " genes without a padj in both contrasts")
      merged <- merged[tested, ]
    }

    # Fill NAs
    lfc_a <- paste0("log2FoldChange.", label_a)
    lfc_b <- paste0("log2FoldChange.", label_b)
    padj_a <- paste0("padj.", label_a)
    padj_b <- paste0("padj.", label_b)

    merged[[lfc_a]][is.na(merged[[lfc_a]])]   <- 0
    merged[[lfc_b]][is.na(merged[[lfc_b]])]   <- 0
    merged[[padj_a]][is.na(merged[[padj_a]])] <- 1
    merged[[padj_b]][is.na(merged[[padj_b]])] <- 1

    # Classify genes
    de_a <- merged[[padj_a]] < args$padj_threshold & abs(merged[[lfc_a]]) >= args$lfc_threshold
    de_b <- merged[[padj_b]] < args$padj_threshold & abs(merged[[lfc_b]]) >= args$lfc_threshold
    lfc_diff <- merged[[lfc_a]] - merged[[lfc_b]]

    merged$category <- case_when(
      de_a & !de_b  ~ paste0(label_a, "_exclusive"),
      !de_a & de_b  ~ paste0(label_b, "_exclusive"),
      de_a & de_b & abs(lfc_diff) >= args$lfc_diff ~ "differential",
      de_a & de_b   ~ "shared",
      TRUE           ~ "not_de"
    )

    merged$lfc_diff <- lfc_diff
    merged$comparison_pair <- paste(label_a, SEP, label_b)

    all_comparisons[[paste(label_a, label_b, sep = SEP)]] <- merged

    # Write per-pair outputs (flat names to avoid staging collisions)
    pair_label <- paste(label_a, label_b, sep = SEP)

    write.csv(merged, file.path(args$outdir, paste0(pair_label, "_compared_genes.csv")), row.names = FALSE)

    # Write comprehensive gene lists per category + direction
    # Exclusive genes
    for (lbl in c(label_a, label_b)) {
      excl_cat <- paste0(lbl, "_exclusive")
      excl <- merged[merged$category == excl_cat, ]
      lfc_col <- if (lbl == label_a) lfc_a else lfc_b

      if (nrow(excl) > 0) {
        writeLines(excl$gene_id, file.path(args$outdir, paste0(pair_label, "_genes_", lbl, "_exclusive.txt")))
        excl_up <- excl[excl[[lfc_col]] > 0, ]
        excl_down <- excl[excl[[lfc_col]] < 0, ]
        if (nrow(excl_up) > 0)
          writeLines(excl_up$gene_id, file.path(args$outdir, paste0(pair_label, "_genes_", lbl, "_exclusive_up.txt")))
        if (nrow(excl_down) > 0)
          writeLines(excl_down$gene_id, file.path(args$outdir, paste0(pair_label, "_genes_", lbl, "_exclusive_down.txt")))
      }
    }

    # Differential genes (DE in both, different magnitude)
    diff_genes <- merged[merged$category == "differential", ]
    if (nrow(diff_genes) > 0) {
      writeLines(diff_genes$gene_id, file.path(args$outdir, paste0(pair_label, "_genes_differential.txt")))
      higher_a <- diff_genes[diff_genes$lfc_diff > 0, ]
      higher_b <- diff_genes[diff_genes$lfc_diff < 0, ]
      if (nrow(higher_a) > 0)
        writeLines(higher_a$gene_id, file.path(args$outdir, paste0(pair_label, "_genes_differential_higher_in_", label_a, ".txt")))
      if (nrow(higher_b) > 0)
        writeLines(higher_b$gene_id, file.path(args$outdir, paste0(pair_label, "_genes_differential_higher_in_", label_b, ".txt")))
    }

    # Shared genes
    shared <- merged[merged$category == "shared", ]
    if (nrow(shared) > 0)
      writeLines(shared$gene_id, file.path(args$outdir, paste0(pair_label, "_genes_shared.txt")))

    # Composite: exclusive + higher in same direction
    for (lbl in c(label_a, label_b)) {
      excl_cat <- paste0(lbl, "_exclusive")
      excl_up <- merged[merged$category == excl_cat & merged[[if (lbl == label_a) lfc_a else lfc_b]] > 0, ]
      diff_higher <- merged[merged$category == "differential" &
                            ((lbl == label_a & merged$lfc_diff > 0) | (lbl == label_b & merged$lfc_diff < 0)), ]
      combo <- rbind(excl_up, diff_higher)
      if (nrow(combo) > 0)
        writeLines(combo$gene_id, file.path(args$outdir, paste0(pair_label, "_genes_exclusive_and_higher_in_", lbl, ".txt")))
    }

    # Summary
    summary_lines <- merged %>%
      filter(category != "not_de") %>%
      count(category) %>%
      mutate(line = paste0(category, ": ", n))
    writeLines(summary_lines$line, file.path(args$outdir, paste0(pair_label, "_summary.txt")))

    message("  ", label_a, " vs ", label_b, ":")
    message("    ", paste(summary_lines$line, collapse = ", "))
  }
}

# Combined report
combined <- bind_rows(all_comparisons)
write.csv(
  combined[combined$category != "not_de", ],
  file.path(args$outdir, "compare_de_report.csv"),
  row.names = FALSE
)

message("Done. Output: ", args$outdir)
