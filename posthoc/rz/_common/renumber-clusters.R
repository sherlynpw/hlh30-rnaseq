#!/usr/bin/env Rscript
# Renumber the pipeline's k-means clusters to the numbers the paper uses.
#
# k-means labels are arbitrary: a re-run (2026-10-12, forward-stranded counts)
# found the same clusters (Jaccard 0.84-0.94 with the old ones) under other
# numbers. kmeans-renumber.tsv maps each run's new number to the paper number,
# so cluster names, GO panels and colours keep their meaning.
#
#   Rscript renumber-clusters.R --run hlh30-all-ruv --csv cache/cluster_assignments.csv

suppressPackageStartupMessages(library(optparse))
opt <- parse_args(OptionParser(option_list = list(
  make_option("--run",  type = "character"),
  make_option("--csv",  type = "character", help = "cluster_assignments.csv, rewritten in place"),
  make_option("--map",  type = "character",
              default = file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE)[1])),
                                  "kmeans-renumber.tsv"))
)))
map <- read.delim(opt$map, comment.char = "#")
map <- map[map$run == opt$run, ]
if (!nrow(map)) stop("no renumbering for run ", opt$run, " in ", opt$map)
d <- read.csv(opt$csv, stringsAsFactors = FALSE, check.names = FALSE)
# idempotent: a file that was already renumbered carries a marker column
if ("pipeline_cluster" %in% names(d)) { message("already renumbered: ", opt$csv); quit(save = "no") }
if (!all(d$cluster %in% map$pipeline)) stop("cluster numbers in ", opt$csv, " not all in the map")
d$pipeline_cluster <- d$cluster
d$cluster <- map$paper[match(d$cluster, map$pipeline)]
write.csv(d, opt$csv, row.names = FALSE)
message("renumbered ", opt$csv, ": ", paste(sprintf("%d->%d", map$pipeline, map$paper), collapse = " "))
