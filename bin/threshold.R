#!/usr/bin/env Rscript
# Hamming distance threshold for clonal grouping from the distance to the
# nearest sequence with the same V, J and junction length (SHazaM, density
# method). When no threshold can be fitted, the fallback is used and said
# on stderr and in THRESHOLD_FILE.source.
# Usage: threshold.R DB_TSV THRESHOLD_FILE FALLBACK
suppressPackageStartupMessages({library(shazam); library(airr)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: threshold.R DB_TSV THRESHOLD_FILE FALLBACK")
db <- read_rearrangement(args[1])
db <- distToNearest(db, sequenceColumn = "junction", vCallColumn = "v_call", jCallColumn = "j_call",
                    model = "ham", normalize = "len", nproc = 1)
fit <- findThreshold(db$dist_nearest, method = "density")
thr <- if (is.null(fit)) NA else fit@threshold
source <- "density"
if (is.na(thr)) {
  thr <- as.numeric(args[3])
  source <- "fallback"
}
writeLines(sprintf("%.4f", thr), args[2])
writeLines(source, paste0(args[2], ".source"))
message(sprintf("threshold %.4f (%s)", thr, source))
