# immuneSIM with R and Python, for the founder step. The base image is fixed
# by digest, CRAN by snapshot date and Bioconductor by release, so a rebuild
# gets the same packages; 'make pin' records the built image's digest in
# images.lock.
#
# CRAN removed immuneSIM on 2025-06-12 (e-mail to its maintainer was
# undeliverable), so install.packages('immuneSIM') no longer finds it. The
# last release, 0.8.7, comes from the CRAN archive and must match a fixed
# SHA-256 or the build stops. Its imports are installed first: the CRAN ones
# from a snapshot inside the R 4.4 period, Biostrings from Bioconductor 3.20,
# the last Bioconductor release for R 4.4.
FROM rocker/r-ver:4.4.1@sha256:78cb94ce2db23aaaf7b546450fcf70b5a3f2ace5a9b5fa1f87217da329211312
ARG CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2025-04-10
ARG BIOC_VERSION=3.20
ARG IMMUNESIM_URL=https://cran.r-project.org/src/contrib/Archive/immuneSIM/immuneSIM_0.8.7.tar.gz
ARG IMMUNESIM_SHA256=54da791105f77a8f066d68d0d479fd952383077669ac1d921db6fe69bb08322f
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 curl ca-certificates \
      libcurl4-openssl-dev libssl-dev libxml2-dev libglpk-dev libgmp3-dev \
      zlib1g-dev libicu-dev \
 && apt-get clean \
 && rm -rf /var/lib/apt/lists/*
# Fetch and check the tarball before the long package builds, so a changed
# archive file fails in seconds.
RUN curl -fsSL -o /opt/immuneSIM_0.8.7.tar.gz "${IMMUNESIM_URL}" \
 && echo "${IMMUNESIM_SHA256}  /opt/immuneSIM_0.8.7.tar.gz" | sha256sum -c -
COPY install.R /opt/install.R
RUN Rscript /opt/install.R "${CRAN_SNAPSHOT}" "${BIOC_VERSION}" /opt/immuneSIM_0.8.7.tar.gz \
 && Rscript -e "library(immuneSIM)"
