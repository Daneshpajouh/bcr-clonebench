# immuneSIM with R and Python, for the founder step. The base image and the
# CRAN snapshot date are fixed so a rebuild gets the same packages; 'make pin'
# records the built image's digest in images.lock.
FROM rocker/r-ver:4.4.1
ARG CRAN_SNAPSHOT=https://packagemanager.posit.co/cran/2026-10-01
RUN apt-get update \
 && apt-get install -y --no-install-recommends python3 \
 && apt-get clean
RUN Rscript -e "install.packages('immuneSIM', repos = '${CRAN_SNAPSHOT}')" \
 && Rscript -e "library(immuneSIM)"
