#!/usr/bin/env Rscript
# The DEGs behind each significant SEA motif (sea_significant.tsv): every DEG
# whose promoter has a match to the motif, one sheet per transcription factor.
#
# SEA picks each motif's score threshold on a random 10% of the genes (the
# hold-out set) and tests the other 90%, so its TP count leaves the hold-out
# genes out. The gene lists here keep them, since the threshold applies to
# every gene; in_holdout marks them, and TP = the rows with in_holdout FALSE.
#
#   Rscript sea-genes.R            # -> out/sea/sea_motif_genes.xlsx

suppressPackageStartupMessages({ library(optparse); library(dplyr); library(openxlsx) })
opt <- parse_args(OptionParser(option_list = list(
  make_option("--deg-dir", type = "character", default = "../volcano/cache/hlh30-all-ruv", dest = "deg_dir"),
  make_option("--padj",    type = "double", default = 0.05),
  make_option("--lfc",     type = "double", default = 1),
  make_option("--seadir",  type = "character", default = "out/sea")
)))
STAGES <- c(L1 = "L1", D0 = "L4", D3 = "Day3", D6 = "Day6")

# DEG log2FC per stage (the same DEG call as sea-dge.R)
deg <- bind_rows(lapply(names(STAGES), function(st) {
  d <- read.csv(file.path(opt$deg_dir, sprintf("tm1978_%s_vs_WT_%s.csv", STAGES[[st]], STAGES[[st]])),
                stringsAsFactors = FALSE)
  d %>% filter(!is.na(padj), padj < opt$padj, abs(log2FoldChange) >= opt$lfc) %>%
    transmute(gene_id, gene_name, code = st, stage = STAGES[[st]], log2FC = round(log2FoldChange, 2),
              padj = signif(padj, 3))
}))
names_of <- deg %>% distinct(gene_id, gene_name)

sig <- read.delim(file.path(opt$seadir, "sea_significant.tsv"), stringsAsFactors = FALSE, check.names = FALSE)
sig$run <- sprintf("sea_%s_u%d", gsub(" ", "_", sig$set), sig$window)

rows <- list()
for (run in unique(sig$run)) {
  s <- read.delim(file.path(opt$seadir, run, "sequences.tsv"), comment.char = "#",
                  stringsAsFactors = FALSE, check.names = FALSE)
  s <- s[s$seq_Class == "tp", ]
  for (i in which(sig$run == run)) {
    r <- sig[i, ]
    g <- s[s$motif_ID == r$ID, ]
    sc <- sub(" .*", "", r$set); dir <- sub(".* ", "", r$set)
    out <- tibble(tf = r$tf, motif_ID = r$ID, database = ifelse(grepl("JASPAR", r$DB), "JASPAR 2022", "CIS-BP 2.0"),
                  consensus = r$CONSENSUS, DEG_set = if (sc == "any") "any stage" else STAGES[[sc]],
                  direction = dir, window = paste0(r$window, " bp upstream"),
                  gene_id = g$seq_ID, motif_score = g$seq_Score, in_holdout = g[["is_holdout?"]] == 1,
                  SEA_q = r$QVALUE)
    # log2FC: that stage for a stage set; every stage where it is a DEG for "any stage"
    fc <- deg %>% filter(gene_id %in% out$gene_id, if (sc == "any") TRUE else code == sc,
                         if (dir == "up") log2FC > 0 else log2FC < 0) %>%
      group_by(gene_id) %>% summarise(stages = paste(stage, collapse = ", "),
                                      log2FC = paste(log2FC, collapse = ", "), .groups = "drop")
    rows[[length(rows) + 1]] <- out %>% left_join(names_of, by = "gene_id") %>% left_join(fc, by = "gene_id")
  }
}
# Most significant first: TFs by their smallest q, tests within a TF by q,
# genes within a test by motif score
all <- bind_rows(rows) %>%
  mutate(DEG_set = factor(DEG_set, c(unname(STAGES), "any stage"))) %>%
  group_by(tf) %>% mutate(tf_best_q = min(SEA_q)) %>% ungroup() %>%
  arrange(tf_best_q, tf, SEA_q, DEG_set, direction, motif_ID, desc(motif_score)) %>%
  select(tf, DEG_set, direction, window, motif_ID, database, consensus, SEA_q,
         gene_name, gene_id, stages, log2FC, motif_score, in_holdout, tf_best_q)

summ <- all %>% group_by(tf, DEG_set, direction, window, motif_ID, database, SEA_q, tf_best_q) %>%
  summarise(genes_with_motif = n(), SEA_TP = sum(!in_holdout), .groups = "drop") %>%
  arrange(tf_best_q, tf, SEA_q, DEG_set, direction, motif_ID) %>% select(-tf_best_q)
all <- select(all, -tf_best_q)
stopifnot(nrow(summ) == nrow(sig))

wb <- createWorkbook()
hs <- createStyle(textDecoration = "bold")
note <- c("DEGs (hlh-30(tm1978) vs wild type, adjusted P < 0.05 and |log2FC| >= 1) whose promoter has a match to a motif that SEA found enriched (q < 0.05).",
          "One sheet per transcription factor. A TF can have more than one motif (JASPAR and CIS-BP) and be enriched in more than one DEG set; the summary has one row per test (TF motif x DEG set x window).",
          "Order: most significant first. TF sheets and summary rows are sorted by the TF's smallest SEA q, then by the q of each test; genes within a test by motif score, highest first.",
          "DEG set 'any stage' = the DEGs of any stage; stages and log2FC then list every stage where the gene is a DEG in that direction.",
          "window = the promoter region scanned (500 or 2,000 bp upstream of the TSS). motif_score = SEA's best match score in that promoter.",
          "SEA picks the score cutoff on a random 10% of the genes (in_holdout TRUE) and tests the other 90%. Its TP count (SEA_TP) leaves the 10% out; the lists keep them, since the cutoff applies to every gene.")
addWorksheet(wb, "README"); writeData(wb, "README", data.frame(note = note), colNames = FALSE)
addWorksheet(wb, "summary"); writeData(wb, "summary", summ, headerStyle = hs); freezePane(wb, "summary", firstRow = TRUE)
used <- c("README", "summary")
for (t in unique(all$tf)) {
  sh <- substr(gsub("[][*?:/\\\\]", "_", t), 1, 31); k <- 1
  while (tolower(sh) %in% tolower(used)) { k <- k + 1; sh <- paste0(substr(sh, 1, 28), "_", k) }
  used <- c(used, sh)
  addWorksheet(wb, sh); writeData(wb, sh, all[all$tf == t, ], headerStyle = hs)
  freezePane(wb, sh, firstRow = TRUE); setColWidths(wb, sh, 1:14, "auto")
}
f <- file.path(opt$seadir, "sea_motif_genes.xlsx")
saveWorkbook(wb, f, overwrite = TRUE)
chk <- summ %>% mutate(window = as.integer(sub(" .*", "", window))) %>%
  inner_join(sig %>% transmute(motif_ID = ID, window, SEA_q = QVALUE, TP,
                               DEG_set = ifelse(grepl("^any", set), "any stage", STAGES[sub(" .*", "", set)]),
                               direction = sub(".* ", "", set)),
             by = c("motif_ID", "window", "SEA_q", "DEG_set", "direction"))
message("wrote ", f, ": ", length(unique(all$tf)), " TF sheets, ", nrow(all), " rows; ",
        sum(chk$SEA_TP == chk$TP), " of ", nrow(sig), " tests match SEA's TP")
