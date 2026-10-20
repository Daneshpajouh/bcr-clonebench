#!/usr/bin/env Rscript
# Naive human IGH rearrangements from immuneSIM (no hypermutation), and the
# germline V alleles they were drawn from.
# Usage: simulate_founders.R N SEED REPERTOIRE_TSV V_GERMLINE_FASTA
suppressPackageStartupMessages(library(immuneSIM))
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 4) stop("usage: simulate_founders.R N SEED REPERTOIRE_TSV V_GERMLINE_FASTA")
n <- as.integer(args[1])
set.seed(as.integer(args[2]))

rep <- immuneSIM(number_of_seqs = n, species = "hs", receptor = "ig", chain = "h",
                 shm.mode = "none", airr_compliant = TRUE, verbose = FALSE)
need <- c("sequence", "v_call", "d_call", "j_call")
missing <- setdiff(need, colnames(rep))
if (length(missing) > 0) stop("immuneSIM output lacks: ", paste(missing, collapse = ", "))
write.table(rep[, need], args[3], sep = "\t", quote = FALSE, row.names = FALSE)

# The germline V alleles: the table immuneSIM samples from for human IGH.
v <- list_germline_genes_allele_01[["hs"]][["ig"]][["h"]][["V"]]
if (!is.data.frame(v) || ncol(v) < 2) stop("unexpected germline table layout: ", paste(class(v), collapse = " "))
is_seq <- vapply(v, function(col) all(grepl("^[ACGTacgt.]+$", as.character(col))), logical(1))
is_name <- vapply(v, function(col) any(grepl("^IGHV", as.character(col))), logical(1))
if (!any(is_seq) || !any(is_name)) stop("cannot find the allele name and sequence columns: ", paste(colnames(v), collapse = ", "))
names_col <- as.character(v[[which(is_name)[1]]])
seq_col <- toupper(gsub(".", "", as.character(v[[which(is_seq)[1]]]), fixed = TRUE))
writeLines(paste0(">", names_col, "\n", seq_col), args[4])
