# Installs immuneSIM 0.8.7 and its imports. CRAN archived immuneSIM on
# 2025-06-12, so it is installed from its checked archive tarball, and its
# imports must be in place first: install.packages with repos = NULL does not
# resolve dependencies.
#
#   Rscript install.R <CRAN snapshot URL> <Bioconductor release> <tarball>
#
# install.packages only warns when a package fails, so every step checks that
# its packages load and stops the build if one does not.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 3) stop("usage: install.R <cran-url> <bioc-version> <tarball>")
cran <- args[1]
bioc <- args[2]
tarball <- args[3]

options(repos = c(CRAN = cran), Ncpus = parallel::detectCores())

need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("not installed: ", paste(missing, collapse = ", "))
}

# The Imports of immuneSIM 0.8.7, less grid (part of base R) and Biostrings
# (Bioconductor).
cran_imports <- c(
  "poweRlaw", "stringdist", "igraph", "stringr", "data.table", "plyr",
  "reshape2", "ggplot2", "ggthemes", "RColorBrewer", "Metrics", "repmis"
)
install.packages(c(cran_imports, "BiocManager"))
need(c(cran_imports, "BiocManager"))

BiocManager::install("Biostrings", version = bioc, ask = FALSE, update = FALSE)
need("Biostrings")
if (as.character(BiocManager::version()) != bioc) {
  stop("Bioconductor ", BiocManager::version(), " in use, expected ", bioc)
}

install.packages(tarball, repos = NULL, type = "source")
need("immuneSIM")
if (as.character(packageVersion("immuneSIM")) != "0.8.7") {
  stop("immuneSIM ", packageVersion("immuneSIM"), " installed, expected 0.8.7")
}
