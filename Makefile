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
REFERENCE ?= ogrdb

.PHONY: test stub pin pin-reference check-partis-r reference known-answer mutation-tests full

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
	$(MAKE) pin-reference

# Download the OGRDB set named in reference.lock, check it against the MD5
# Zenodo publishes, and write its SHA-256 into reference.lock if not yet
# pinned. A pinned SHA-256 is never changed: a file that differs stops here.
pin-reference:
	python3 bin/fetch_ogrdb.py --pin reference.lock reference-src

# R starts and TreeSim loads in the partis-r image, as the host user, the way
# Nextflow runs it (docker.runOptions). PARTIS_R_CHECK defaults to the pinned one.
PARTIS_R_CHECK ?= $(PARTIS_R)
check-partis-r:
	docker run --rm -u $$(id -u):$$(id -g) $(PARTIS_R_CHECK) sh -c 'R --version && Rscript -e "library(TreeSim)"'

# The germline reference, built inside the pinned suite image. The default is
# the OGRDB set pinned in reference.lock (downloaded and checked on the host);
# REFERENCE=imgt uses the IMGT set the suite image ships instead.
reference:
ifeq ($(REFERENCE),ogrdb)
	python3 bin/fetch_ogrdb.py reference.lock reference-src
	docker run --rm -u $$(id -u):$$(id -g) -v $(CURDIR):/work -w /work $(SUITE) bin/build_reference.sh ogrdb reference reference-src/IGHV.fasta reference-src/IGHD.fasta reference-src/IGHJ.fasta
else
	docker run --rm -u $$(id -u):$$(id -g) -v $(CURDIR):/work -w /work $(SUITE) bin/build_reference.sh $(REFERENCE) reference
endif

# The whole known-answer test in the pinned containers, then the named checks
# (conf/known_answer_floors.json) on its summary table.
RUN = nextflow run main.nf -profile docker,test --suite_image $(SUITE) --immunesim_image $(IMMUNESIM) --partis_image $(PARTIS_R)
CHECK = python3 bin/check_known_answer.py
known-answer: reference
	$(RUN)
	$(CHECK) results/results.tsv conf/known_answer_floors.json

# Mutation tests, after known-answer in the same folder (-resume reuses the
# simulation): a deliberately wrong clone threshold and a deliberately wrong
# germline set must each make their named check fail.
mutation-tests:
	$(RUN) -resume --threshold_override 0.5 --outdir results-mutant-threshold
	$(CHECK) results-mutant-threshold/results.tsv conf/known_answer_floors.json --expect-fail KA-CHANGEO-PRECISION
	docker run --rm -u $$(id -u):$$(id -g) -v $(CURDIR):/work -w /work $(SUITE) bin/build_reference.sh shuffled-v reference-mutant-germline reference 20261021
	$(RUN) -resume --reference_dir $(CURDIR)/reference-mutant-germline --outdir results-mutant-germline
	$(CHECK) results-mutant-germline/results.tsv conf/known_answer_floors.json --expect-fail KA-V-INCLUDES

# The full-size run with the default parameters (2,000 founder families per level, the three levels, the three
# sequencing error rates), in the pinned containers, into results-full; then the same named checks on its table.
# Run by the manual 'full' workflow in CI (its time on a hosted runner is not measured yet).
FULL = nextflow run main.nf -profile docker --suite_image $(SUITE) --immunesim_image $(IMMUNESIM) --partis_image $(PARTIS_R) --outdir results-full
full: reference
	$(FULL)
	$(CHECK) results-full/results.tsv conf/known_answer_floors.json
