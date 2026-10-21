#!/usr/bin/env Rscript
# Clonal grouping with SCOPer: identical junctions, hierarchical single
# linkage at the SHazaM threshold, or spectral clustering of junctions (spectral) or of junctions with the
# mutations shared in the V and J segments (spectral_vj).
# Usage: clones_scoper.R DB_TSV METHOD THRESHOLD OUT_TSV
suppressPackageStartupMessages({library(scoper); library(airr)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: clones_scoper.R DB_TSV METHOD THRESHOLD OUT_TSV")
db <- read_rearrangement(args[1])
thr <- as.numeric(args[3])
common <- list(db = db, junction = "junction", v_call = "v_call", j_call = "j_call", summarize_clones = FALSE)
res <- switch(args[2],
  identical = do.call(identicalClones, c(common, method = "nt")),
  hierarchical = do.call(hierarchicalClones, c(common, threshold = thr, method = "nt", linkage = "single",
                                               normalize = "len")),
  # "novj" clusters on the junction alone, as the other methods here do
  spectral = do.call(spectralClones, c(common, method = "novj")),
  # "vj" also weighs mutations the reads share in V and J against the germline (sequence_alignment and
  # germline_alignment), so families with one V, one J and one junction length can still be told apart (run 19, F2)
  spectral_vj = do.call(spectralClones, c(common, method = "vj", sequence = "sequence_alignment",
                                          germline = "germline_alignment")),
  stop("unknown method: ", args[2]))
write_rearrangement(as.data.frame(res), args[4])
