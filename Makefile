# bcr-clonebench. 'make test' and 'make stub' need only Python 3 and
# Nextflow; the other targets need Docker.
SUITE_TAG ?= immcantation/suite:4.5.0
IMMUNESIM_TAG ?= bcr-clonebench/immunesim:local
PARTIS_TAG ?= psathyrella/partis:latest
SUITE = $(shell awk -F'\t' '$$1=="suite_image"{print $$3}' images.lock)
IMMUNESIM = $(shell awk -F'\t' '$$1=="immunesim_image"{print $$3}' images.lock)
PARTIS = $(shell awk -F'\t' '$$1=="partis_image"{print $$3}' images.lock)
REFERENCE ?= imgt

.PHONY: test stub pin reference known-answer

test:
	python3 -B -m unittest discover -s tests

stub:
	nextflow run main.nf -stub -profile test

# Pull and build the images, then write their digests into images.lock.
pin:
	docker pull $(SUITE_TAG)
	docker pull $(PARTIS_TAG)
	docker build -t $(IMMUNESIM_TAG) containers/immunesim
	python3 bin/pin_images.py images.lock suite_image $(SUITE_TAG) immunesim_image $(IMMUNESIM_TAG) partis_image $(PARTIS_TAG)

# The germline reference, built inside the pinned suite image.
reference:
	docker run --rm -u $$(id -u):$$(id -g) -v $(CURDIR):/work -w /work $(SUITE) bin/build_reference.sh $(REFERENCE) reference

# The whole known-answer test in the pinned containers.
known-answer: reference
	nextflow run main.nf -profile docker,test --suite_image $(SUITE) --immunesim_image $(IMMUNESIM) --partis_image $(PARTIS)
