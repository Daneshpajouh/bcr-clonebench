#!/usr/bin/env Rscript
# Observed mutation frequency over the IMGT V region (SHazaM), against
# germlines built per clone with the D region masked; written as mu_freq.
# Usage: mutation.R GERMLINES_TSV OUT_TSV
suppressPackageStartupMessages({library(shazam); library(airr)})
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) stop("usage: mutation.R GERMLINES_TSV OUT_TSV")
db <- read_rearrangement(args[1])
db <- observedMutations(db, sequenceColumn = "sequence_alignment", germlineColumn = "germline_alignment_d_mask",
                        regionDefinition = IMGT_V, frequency = TRUE, combine = TRUE, nproc = 1)
write_rearrangement(db, args[2])
