#!/usr/bin/env Rscript

# deseq2.R - DESeq2 differential expression analysis
#
# Input:  gene_counts.tsv (genes x samples) + sample_meta.tsv
# Output: results CSV, VST counts, normalized counts, size factors

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(dplyr)
  library(RSQLite)
})

.find_lib <- function() {
  for (p in c(Sys.which("gene_db_lib.R"), file.path("bin", "gene_db_lib.R")))
    if (nchar(p) > 0 && file.exists(p)) return(p)
  NULL
}
lib <- .find_lib()
if (!is.null(lib)) source(lib)

opts <- list(
  make_option("--counts",         type = "character", help = "Gene count matrix TSV"),
  make_option("--metadata",       type = "character", help = "Sample metadata TSV"),
  make_option("--condition_from", type = "character", default = "genotype_timepoint",
              help = "How to build condition: column name or 'col1_col2' to combine. Ignored when --design_formula is set."),
  make_option("--reference",      type = "character", help = "Reference level for condition (pairwise mode)"),
  make_option("--padj",           type = "numeric",   default = 0.05, help = "Adjusted p-value threshold"),
  make_option("--lfc",            type = "numeric",   default = 1.0,  help = "Log2 fold change threshold"),
  make_option("--shrinkage",      type = "character", default = "apeglm", help = "LFC shrinkage method"),
  make_option("--batch",          type = "character", default = NULL, help = "Column name for batch covariate"),
  make_option("--ruv_factors",    type = "character", default = NULL, help = "TSV of W factors from RUVSeq (sample_id + W_1..W_k)"),
  make_option("--apply_batch_correction", action = "store_true", default = FALSE,
              help = "Apply limma::removeBatchEffect to VST output for visualization/clustering"),
  make_option("--min_count",      type = "integer",   default = 10, help = "Minimum total count per gene"),
  make_option("--gene_db",        type = "character", default = NULL, help = "Path to gene_db.sqlite"),
  make_option("--contrasts_json", type = "character", default = NULL,
              help = "JSON file mapping contrast name -> { coef_name: weight, ... }. Computed via DESeq2::results(contrast = numeric_vector)."),
  make_option("--design_formula", type = "character", default = NULL,
              help = "MULTI-FACTOR MODE: full design formula override (e.g. '~ treatment + genotype + timepoint + replicate'). If set, replaces the auto-built ~[RUV]+[batch]+condition design. Output becomes one CSV per DESeq2 coefficient, named <factor>__<level>_vs_<reference>.csv. Disables the pairwise-contrast loop. Ignores --condition_from and --reference."),
  make_option("--reference_levels", type = "character", default = NULL,
              help = "Multi-factor mode only: comma-separated 'col1:level1,col2:level2' pairs to set the reference level for each factor (e.g. 'treatment:none,genotype:WT'). If unset, DESeq2 picks the alphabetically first level."),
  make_option("--subset_col",     type = "character", default = NULL,
              help = "Grouped mode: metadata column to subset samples by before fitting."),
  make_option("--subset_val",     type = "character", default = NULL,
              help = "Grouped mode: value of --subset_col to keep."),
  make_option("--label_prefix",   type = "character", default = "",
              help = "Prefix prepended to result CSV filenames (e.g. '<group>__')."),
  make_option("--outdir",         type = "character", default = ".", help = "Output directory")
)

args <- parse_args(OptionParser(option_list = opts))

# Open gene_db for annotation (optional)
gdb <- NULL
if (!is.null(args$gene_db) && file.exists(args$gene_db)) {
  gdb <- open_gene_db(args$gene_db)
  message("Using gene_db for annotation")
}

dir.create(args$outdir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(args$outdir, "results"), recursive = TRUE, showWarnings = FALSE)

# Read data
message("Reading count matrix...")
counts <- read.delim(args$counts, row.names = 1, check.names = FALSE)
meta   <- read.delim(args$metadata, stringsAsFactors = FALSE)
rownames(meta) <- meta$sample_id

# Grouped mode: keep only the samples for this group before anything else
if (!is.null(args$subset_col) && nzchar(args$subset_col)) {
  if (!(args$subset_col %in% colnames(meta)))
    stop("subset_col '", args$subset_col, "' not in metadata columns")
  before <- nrow(meta)
  meta <- meta[!is.na(meta[[args$subset_col]]) &
                 meta[[args$subset_col]] == args$subset_val, , drop = FALSE]
  message("Subset ", args$subset_col, " == ", args$subset_val, ": ",
          nrow(meta), "/", before, " samples")
  if (nrow(meta) == 0) stop("No samples left after subset")
}

# Ensure sample order matches
shared_samples <- intersect(colnames(counts), meta$sample_id)
stopifnot("No matching samples between counts and metadata" = length(shared_samples) > 0)
counts <- counts[, shared_samples, drop = FALSE]
meta   <- meta[shared_samples, , drop = FALSE]

# Optional: merge RUVSeq W factors (needed in both modes for design assembly)
ruv_cols <- character(0)
if (!is.null(args$ruv_factors) && file.exists(args$ruv_factors)) {
  W <- read.delim(args$ruv_factors, stringsAsFactors = FALSE)
  stopifnot("ruv_factors TSV must have sample_id column" = "sample_id" %in% colnames(W))
  ruv_cols <- grep("^W_", colnames(W), value = TRUE)
  if (length(ruv_cols) == 0) stop("No W_* columns found in ", args$ruv_factors)
  W <- W[match(rownames(meta), W$sample_id), , drop = FALSE]
  for (w in ruv_cols) meta[[w]] <- as.numeric(W[[w]])
  message("RUV factors merged: ", paste(ruv_cols, collapse = ", "))
}

# -------- DESIGN MODE BRANCH --------
# multifactor_mode = TRUE  ->  use --design_formula directly, output per-coefficient
#                              CSVs (multi-factor mode; bypasses condition/pairwise logic)
# multifactor_mode = FALSE ->  build categorical `condition` from --condition_from
#                              and run pairwise contrasts (legacy behavior)
multifactor_mode <- !is.null(args$design_formula) && nzchar(args$design_formula)

if (multifactor_mode) {
  message("MULTI-FACTOR MODE: --design_formula override active")
  user_formula <- as.formula(args$design_formula)
  user_vars    <- all.vars(user_formula)
  missing_cols <- setdiff(user_vars, c(colnames(meta), ruv_cols))
  if (length(missing_cols) > 0) {
    stop("design_formula references metadata columns not in the samplesheet: ",
         paste(missing_cols, collapse = ", "))
  }

  # Auto-prepend nuisance terms (mirrors pairwise mode behavior). The user
  # writes BIOLOGICAL terms only; the pipeline prepends [W_*] (RUV factors)
  # and [batch] so DE statistics are corrected for them. Skips prepending
  # if the user explicitly included those terms in their formula already.
  nuisance_terms <- character(0)
  for (rc in ruv_cols) {
    if (!(rc %in% user_vars)) nuisance_terms <- c(nuisance_terms, rc)
  }
  if (!is.null(args$batch) && args$batch %in% colnames(meta) &&
      !(args$batch %in% user_vars)) {
    meta[[args$batch]] <- factor(meta[[args$batch]])
    nuisance_terms <- c(nuisance_terms, args$batch)
  }
  # Bio terms in user-specified order (right-hand side of user_formula)
  bio_terms <- attr(terms(user_formula), "term.labels")
  if (length(bio_terms) == 0) {
    stop("design_formula has no RHS terms; cannot run multi-factor analysis")
  }
  all_terms <- c(nuisance_terms, bio_terms)
  design_formula <- as.formula(paste("~", paste(all_terms, collapse = " + ")))
  design_vars    <- all.vars(design_formula)

  # Factorize categorical columns (including any newly added batch)
  for (col in intersect(design_vars, colnames(meta))) {
    if (!is.numeric(meta[[col]])) meta[[col]] <- factor(meta[[col]])
  }

  # Apply per-factor reference levels (BEFORE DESeq2 sees the data)
  if (!is.null(args$reference_levels) && nzchar(args$reference_levels)) {
    for (kv_raw in trimws(strsplit(args$reference_levels, ",", fixed = TRUE)[[1]])) {
      kv <- strsplit(kv_raw, ":", fixed = TRUE)[[1]]
      if (length(kv) != 2) {
        warning("Ignoring malformed reference_levels entry: ", kv_raw)
        next
      }
      col <- kv[1]; ref_lvl <- kv[2]
      if (!(col %in% colnames(meta))) {
        warning("reference_levels: column '", col, "' not in metadata; skipping")
        next
      }
      if (!is.factor(meta[[col]])) meta[[col]] <- factor(meta[[col]])
      if (!(ref_lvl %in% levels(meta[[col]]))) {
        warning("reference_levels: level '", ref_lvl, "' not in column '", col, "' levels (",
                paste(levels(meta[[col]]), collapse = ","), "); skipping")
        next
      }
      meta[[col]] <- relevel(meta[[col]], ref = ref_lvl)
      message("  reference for ", col, ": ", ref_lvl)
    }
  }

  # Record bio_vars separately so coefficient output skips nuisance coefs
  bio_vars <- intersect(all.vars(as.formula(paste("~", paste(bio_terms, collapse = " + ")))),
                        c(colnames(meta), ruv_cols))

  message("User formula:        ", deparse(user_formula))
  if (length(nuisance_terms) > 0) {
    message("Nuisance auto-added: ", paste(nuisance_terms, collapse = " + "))
  }
  message("Final design:        ", deparse(design_formula))
} else {
  # Build condition column (legacy / pairwise mode)
  condition_parts <- strsplit(args$condition_from, "_")[[1]]
  if (length(condition_parts) == 1 && condition_parts[1] %in% colnames(meta)) {
    meta$condition <- as.character(meta[[condition_parts[1]]])
  } else {
    meta$condition <- apply(meta[, condition_parts, drop = FALSE], 1, paste, collapse = "_")
  }
  meta$condition <- factor(meta$condition)
  if (!is.null(args$reference)) {
    stopifnot("Reference not found in conditions" = args$reference %in% levels(meta$condition))
    meta$condition <- relevel(meta$condition, ref = args$reference)
  }
  message("Conditions: ", paste(levels(meta$condition), collapse = ", "))
  message("Reference:  ", levels(meta$condition)[1])

  # Build design formula: ~ [W_1+..+W_k] + [batch] + condition
  design_terms <- character(0)
  if (length(ruv_cols) > 0) design_terms <- c(design_terms, ruv_cols)
  if (!is.null(args$batch) && args$batch %in% colnames(meta)) {
    meta[[args$batch]] <- factor(meta[[args$batch]])
    design_terms <- c(design_terms, args$batch)
  }
  design_terms <- c(design_terms, "condition")
  design_formula <- as.formula(paste("~", paste(design_terms, collapse = " + ")))
  message("Design (pairwise): ", deparse(design_formula))
}

# Filter low-count genes
keep <- rowSums(counts) >= args$min_count
counts <- counts[keep, ]
message("Genes after filtering (min count ", args$min_count, "): ", nrow(counts))

# Run DESeq2
message("Running DESeq2...")
dds <- DESeqDataSetFromMatrix(
  countData = as.matrix(counts),
  colData   = meta,
  design    = design_formula
)
dds <- DESeq(dds)

# Persist the fitted object so downstream tasks (e.g. compound contrasts run
# separately, ad-hoc result extraction) can reuse it without refitting.
saveRDS(dds, file = file.path(args$outdir, "dds.rds"))

# -------- RESULTS OUTPUT BRANCH --------
if (multifactor_mode) {
  # MULTI-FACTOR MODE: one CSV per DESeq2 coefficient.
  #
  # Coefficient names from resultsNames(dds) look like "genotype_mut_vs_WT",
  # "timepoint_Day3_vs_L4", "W_1" (RUV), etc.
  # We rewrite filenames to the convention "<factor>__<level>_vs_<reference>"
  # so they're easier to read and group.
  coefs <- resultsNames(dds)
  message("Multi-factor coefficients: ", paste(coefs, collapse = ", "))

  format_output_name <- function(coef, vars) {
    if (coef == "(Intercept)") return(NULL)        # skip intercept
    # Find which formula variable this coef belongs to and split out the suffix.
    for (var in vars) {
      pat <- paste0("^", var, "_")
      if (grepl(pat, coef)) {
        rest <- sub(pat, "", coef)
        return(paste0(var, "__", rest))
      }
    }
    coef                                            # fallback (RUV W_*, etc.)
  }

  # Decide which variable a coefficient belongs to (used for nuisance filtering)
  coef_to_var <- function(coef, vars) {
    if (coef == "(Intercept)") return(NULL)
    for (var in vars) {
      pat <- paste0("^", var, "_")
      if (grepl(pat, coef)) return(var)
    }
    # RUV W_* names equal their var name exactly
    if (coef %in% vars) return(coef)
    NULL
  }

  for (coef in coefs) {
    output_name <- format_output_name(coef, design_vars)
    if (is.null(output_name)) next                  # skipped intercept

    # Skip nuisance coefficients (RUV W_*, batch_*) -- they're in the model
    # to correct DE statistics but aren't biological effects we want output for
    coef_var <- coef_to_var(coef, design_vars)
    if (!is.null(coef_var) && !(coef_var %in% bio_vars)) {
      message("  Skipping nuisance coef: ", coef, " (var=", coef_var, ")")
      next
    }

    message("Extracting coefficient: ", coef, " -> ", output_name)
    res <- results(dds, name = coef, alpha = args$padj)
    # Unshrunken MLE stats keyed by gene; joined back below (compare_de z-test
    # needs MLE LFC+SE, and gene_id join fixes a post-arrange() row misalignment).
    mle_df <- data.frame(
      gene_id            = rownames(res),
      stat               = res$stat,
      log2FoldChange_MLE = res$log2FoldChange,
      lfcSE_MLE          = res$lfcSE,
      stringsAsFactors   = FALSE
    )

    if (args$shrinkage != "none") {
      # apeglm works with `coef = ...`; for multi-factor we can use it directly
      # (unlike pairwise mode where contrast = c(name, A, B) forced ashr).
      shrink_type <- args$shrinkage
      res <- tryCatch(
        lfcShrink(dds, coef = coef, type = shrink_type, quiet = TRUE),
        error = function(e) {
          message("  ", shrink_type, " shrinkage failed (", e$message, "); using ashr")
          lfcShrink(dds, coef = coef, type = "ashr", quiet = TRUE)
        }
      )
    }

    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      left_join(mle_df, by = "gene_id") %>%
      arrange(padj, desc(abs(log2FoldChange)))
    if (!is.null(gdb)) res_df <- annotate_df(gdb, res_df)

    write.csv(
      res_df,
      file      = file.path(args$outdir, "results", paste0(args$label_prefix, output_name, ".csv")),
      row.names = FALSE
    )
    message("  Genes in results: ", nrow(res_df))
  }
} else {
  # PAIRWISE MODE (legacy): pairwise contrasts within / across genotype_timepoint.
  #   - Within genotype across timepoints (e.g., hlh-30_Day3 vs hlh-30_L4)
  #   - Same timepoint across genotypes (e.g., hlh-30_L4 vs N2_L4)
  all_conditions <- levels(meta$condition)

  # Parse genotype and timepoint from condition names
  cond_parts <- do.call(rbind, strsplit(all_conditions, "_"))
  cond_info <- data.frame(
    condition = all_conditions,
    genotype  = apply(cond_parts[, 1:(ncol(cond_parts)-1), drop = FALSE], 1, paste, collapse = "_"),
    timepoint = cond_parts[, ncol(cond_parts)],
    stringsAsFactors = FALSE
  )

  valid_pairs <- list()
  ref_genotype <- sub("_[^_]+$", "", levels(meta$condition)[1])

  for (i in 1:(nrow(cond_info) - 1)) {
    for (j in (i + 1):nrow(cond_info)) {
      same_genotype  <- cond_info$genotype[i] == cond_info$genotype[j]
      same_timepoint <- cond_info$timepoint[i] == cond_info$timepoint[j]
      if (same_genotype || same_timepoint) {
        ci <- cond_info$condition[i]
        cj <- cond_info$condition[j]
        gi <- cond_info$genotype[i]
        gj <- cond_info$genotype[j]

        if (same_timepoint) {
          if (gi == ref_genotype) {
            pair <- c(cj, ci)
          } else {
            pair <- c(ci, cj)
          }
        } else {
          pair <- c(cj, ci)
        }
        valid_pairs[[length(valid_pairs) + 1]] <- pair
      }
    }
  }

  message("Valid comparisons: ", length(valid_pairs))

  for (pair in valid_pairs) {
    contrast_label <- paste0(pair[1], "_vs_", pair[2])
    message("Extracting results: ", contrast_label)

    res <- results(dds, contrast = c("condition", pair[1], pair[2]), alpha = args$padj)
    # Unshrunken MLE stats keyed by gene; joined back below (compare_de z-test
    # needs MLE LFC+SE, and gene_id join fixes a post-arrange() row misalignment).
    mle_df <- data.frame(
      gene_id            = rownames(res),
      stat               = res$stat,
      log2FoldChange_MLE = res$log2FoldChange,
      lfcSE_MLE          = res$lfcSE,
      stringsAsFactors   = FALSE
    )

    if (args$shrinkage != "none") {
      shrink_type <- if (args$shrinkage == "apeglm") "ashr" else args$shrinkage
      res <- lfcShrink(dds, contrast = c("condition", pair[1], pair[2]),
                       type = shrink_type, quiet = TRUE)
    }

    res_df <- as.data.frame(res) %>%
      tibble::rownames_to_column("gene_id") %>%
      left_join(mle_df, by = "gene_id") %>%
      arrange(padj, desc(abs(log2FoldChange)))
    if (!is.null(gdb)) res_df <- annotate_df(gdb, res_df)

    write.csv(
      res_df,
      file      = file.path(args$outdir, "results", paste0(args$label_prefix, contrast_label, ".csv")),
      row.names = FALSE
    )
    message("  Genes in results: ", nrow(res_df))
  }
}

# -------- COMPOUND CONTRASTS BRANCH --------
# User-defined linear combinations of fitted coefficients. JSON format:
#   { "ContrastName": {"coef_name_1": 1, "coef_name_2": -1, ...}, ... }
# Each contrast is computed via DESeq2::results(dds, contrast = numeric_vector)
# and lands in contrast_results/<safe_name>.csv. Lives in its own subdir so it
# doesn't accidentally cascade into downstream enrichment that consumes results/.
if (!is.null(args$contrasts_json) && nzchar(args$contrasts_json) &&
    file.exists(args$contrasts_json)) {
  suppressPackageStartupMessages(library(jsonlite))
  contrasts_map <- jsonlite::fromJSON(args$contrasts_json, simplifyDataFrame = FALSE)
  if (length(contrasts_map)) {
    contrast_dir <- file.path(args$outdir, "contrast_results")
    dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)
    coef_names <- resultsNames(dds)
    message("Available coefficients for contrasts:")
    for (cn in coef_names) message("  ", cn)

    # Most lfcShrink types require coef= rather than contrast=, so for compound
    # contrasts we use ashr (the only shrinker that accepts arbitrary contrast
    # vectors). If user picked --shrinkage none, we honor that.
    shrink_compound <- if (args$shrinkage == "none") "none" else "ashr"

    for (cname in names(contrasts_map)) {
      weights_map <- contrasts_map[[cname]]
      cvec <- numeric(length(coef_names))
      names(cvec) <- coef_names
      missing_coefs <- setdiff(names(weights_map), coef_names)
      if (length(missing_coefs)) {
        message("[", cname, "] coefficients not in fit, skipping: ",
                paste(missing_coefs, collapse = ", "))
        next
      }
      for (k in names(weights_map)) cvec[k] <- as.numeric(weights_map[[k]])
      message("Compound contrast: ", cname, " = ",
              paste(sprintf("%g*%s", cvec[cvec != 0], names(cvec)[cvec != 0]),
                    collapse = " + "))

      res <- results(dds, contrast = cvec, alpha = args$padj)
      mle_df <- data.frame(
        gene_id            = rownames(res),
        stat               = res$stat,
        log2FoldChange_MLE = res$log2FoldChange,
        lfcSE_MLE          = res$lfcSE,
        stringsAsFactors   = FALSE
      )
      if (shrink_compound != "none") {
        res <- tryCatch(
          lfcShrink(dds, contrast = cvec, type = shrink_compound, quiet = TRUE),
          error = function(e) {
            message("[", cname, "] lfcShrink failed (", e$message, "); using unshrunken LFC")
            res
          })
      }

      res_df <- as.data.frame(res) %>%
        tibble::rownames_to_column("gene_id") %>%
        left_join(mle_df, by = "gene_id") %>%
        arrange(padj, desc(abs(log2FoldChange)))
      if (!is.null(gdb)) res_df <- annotate_df(gdb, res_df)

      safe <- gsub("[^A-Za-z0-9_.-]", "_", cname)
      write.csv(
        res_df,
        file      = file.path(contrast_dir, paste0(args$label_prefix, safe, ".csv")),
        row.names = FALSE
      )
      message("  Genes in results: ", nrow(res_df))
    }
  }
}

# VST-normalized counts
message("Computing VST counts...")
vst_obj <- tryCatch(
  vst(dds, blind = FALSE),
  error = function(e) {
    message("vst() failed (", e$message, "), using varianceStabilizingTransformation()")
    varianceStabilizingTransformation(dds, blind = FALSE)
  }
)
vst_counts <- assay(vst_obj)

# Always export the raw (uncorrected) VST so downstream consumers can compare
# corrected vs uncorrected (e.g. PCA before/after). The "active" matrix written
# below as vst_counts.csv may equal raw or be limma-corrected.
write.csv(vst_counts, file = file.path(args$outdir, "vst_counts_raw.csv"))

# Optionally apply limma::removeBatchEffect for downstream visualization/clustering.
# DE statistics already account for batch/W via the design formula; this only
# affects the exported VST matrix used by exploration/clustering.
if (isTRUE(args$apply_batch_correction)) {
  has_batch <- !is.null(args$batch) && args$batch %in% colnames(meta)
  has_ruv   <- length(ruv_cols) > 0
  if (!has_batch && !has_ruv) {
    message("apply_batch_correction set but no --batch or --ruv_factors - skipping VST correction")
  } else {
    suppressPackageStartupMessages(library(limma))
    sample_order <- colnames(vst_counts)
    meta_ord <- meta[sample_order, , drop = FALSE]
    # Build a "biology of interest" model matrix for removeBatchEffect to
    # preserve. In pairwise mode this is ~condition. In multi-factor mode
    # we preserve everything in the design formula EXCEPT the batch/RUV
    # nuisance variables.
    if (multifactor_mode) {
      # Preserve biological terms; remove nuisance (RUV/batch). `bio_vars`
      # was set during multi-factor design construction.
      if (length(bio_vars) == 0) {
        design_vst <- model.matrix(~ 1, data = meta_ord)
      } else {
        keep_formula <- as.formula(paste("~", paste(bio_vars, collapse = " + ")))
        design_vst <- model.matrix(keep_formula, data = meta_ord)
      }
    } else {
      design_vst <- model.matrix(~ condition, data = meta_ord)
    }
    batch_vec  <- if (has_batch) factor(meta_ord[[args$batch]]) else NULL
    cov_mat    <- if (has_ruv) as.matrix(meta_ord[, ruv_cols, drop = FALSE]) else NULL
    msg_parts <- c()
    if (has_batch) msg_parts <- c(msg_parts, paste0("batch=", args$batch))
    if (has_ruv)   msg_parts <- c(msg_parts, paste0("covariates=", paste(ruv_cols, collapse = ",")))
    message("VST correction (limma::removeBatchEffect): ", paste(msg_parts, collapse = "; "))
    vst_counts <- removeBatchEffect(vst_counts, batch = batch_vec,
                                    covariates = cov_mat, design = design_vst)
  }
}

write.csv(vst_counts, file = file.path(args$outdir, "vst_counts.csv"))

# Normalized counts
norm_counts <- counts(dds, normalized = TRUE)
write.csv(norm_counts, file = file.path(args$outdir, "normalized_counts.csv"))

# Size factors
sf <- sizeFactors(dds)
write.csv(
  data.frame(sample = names(sf), size_factor = sf),
  file      = file.path(args$outdir, "size_factors.csv"),
  row.names = FALSE
)

if (!is.null(gdb)) close_gene_db(gdb)
message("Done. Output: ", args$outdir)
