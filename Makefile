# bcr-clonebench. 'make test' and 'make stub' need only Python 3 and
# Nextflow; the other targets need Docker.
SUITE_TAG ?= immcantation/suite:4.5.0
IMMUNESIM_TAG ?= bcr-clonebench/immunesim:local
PARTIS_TAG ?= quay.io/matsengrp/partis:latest
PARTIS_R_TAG ?= bcr-clonebench/partis-r:local
SUITE = $(shell awk -F'\t' '$$1=="suite_image"{print $$3}' images.lock)
IMMUNESIM = $(shell awk -F'\t' '$$1=="immunesim_image"{print $$3}' images.lock)
PARTIS = $(shell awk -F'\t' '$$1=="partis_image"{print $$3}' images.lock)
PARTIS_R = $(shell awk -F'\t' '$$1=="partis_r_image"{print $$3}' images.lock)
REFERENCE ?= imgt

.PHONY: test stub pin check-partis-r reference known-answer

test:
	python3 -B -m unittest discover -s tests

stub:
	nextflow run main.nf -stub -profile test

# Build and pull the images, then write their digests into images.lock.
# A pulled image already locked by digest must still resolve to it, or this stops.
# The R build goes first, while the disk is emptiest, then one pull at a time;
# free space is printed after each image. partis-r is built on the digest of
# the partis image just pulled, so it shares that image's layers and matches
# the digest pin_images.py records (or refuses, if the tag moved), and is checked as the pipeline runs it (host user, R and TreeSim) before pinning.
pin:
	docker build -t $(IMMUNESIM_TAG) containers/immunesim
	df -h /
	docker pull $(SUITE_TAG)
	df -h /
	docker pull $(PARTIS_TAG)
	df -h /
	base=$$(docker image inspect --format '{{index .RepoDigests 0}}' $(PARTIS_TAG)) \
	  && { echo "$$base" | grep -Eq '@sha256:[0-9a-f]{64}$$' || { echo "pulled partis image has no registry digest: $$base" >&2; exit 1; }; } \
	  && docker build --build-arg PARTIS_BASE=$$base -t $(PARTIS_R_TAG) containers/partis-r
	df -h /
	$(MAKE) check-partis-r PARTIS_R_CHECK=$(PARTIS_R_TAG)
	python3 bin/pin_images.py images.lock suite_image $(SUITE_TAG) immunesim_image $(IMMUNESIM_TAG) partis_image $(PARTIS_TAG) partis_r_image $(PARTIS_R_TAG)

# R starts and TreeSim loads in the partis-r image, as the host user, the way
# Nextflow runs it (docker.runOptions). PARTIS_R_CHECK defaults to the pinned one.
PARTIS_R_CHECK ?= $(PARTIS_R)
check-partis-r:
	docker run --rm -u $$(id -u):$$(id -g) $(PARTIS_R_CHECK) sh -c 'R --version && Rscript -e "library(TreeSim)"'

# The germline reference, built inside the pinned suite image.
reference:
	docker run --rm -u $$(id -u):$$(id -g) -v $(CURDIR):/work -w /work $(SUITE) bin/build_reference.sh $(REFERENCE) reference

# The whole known-answer test in the pinned containers.
known-answer: reference
	nextflow run main.nf -profile docker,test --suite_image $(SUITE) --immunesim_image $(IMMUNESIM) --partis_image $(PARTIS_R)
