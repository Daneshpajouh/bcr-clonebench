# bcr-clonebench

bcr-clonebench checks how well B-cell receptor (BCR) repertoire tools recover the truth. It simulates human heavy-chain (IGH) repertoires in which every read's clonal family, germline genes and mutations are known, runs them through a standard Immcantation pipeline, and scores each clonal grouping method against those known answers. It is for anyone who groups BCR sequences into clones and wants to know which method to trust at which mutation level, and how much sequencing error changes the answer.

It is a benchmark harness, not a new grouping method.

## Status

- Built and unit-tested: the scorer, the truth converters, the second lineage generator, sequencing-error injection and the summary table (`make test`).
- The whole Nextflow pipeline runs end to end in stub mode (`make stub`), which checks the wiring without running any tool.
- Not yet done: the first run in the pinned containers, so there are no benchmark numbers in this repository yet. The real-data run (below) has not started. Results will be added here when they exist, with the exact command and image digests that produced them.

## How it works

1. **Simulate clonal families.** partis simulates human IGH clonal families from scratch at three hypermutation levels (low, medium, high), and writes the true family, genes and mutation state of every read.
2. **Add sequencing errors.** Each level is copied at each sequencing-error rate (0, 0.1 and 0.5 percent per base by default). The errors are kept out of the truth: they are noise the pipeline has to see through.
3. **Annotate.** IgBLAST assigns V, D and J genes against the OGRDB human IGH germline set, and Change-O writes an AIRR table.
4. **Pick a clustering threshold.** SHazaM estimates the distance threshold from the data. If it finds none, a fixed fallback is used and the run says so.
5. **Group into clones** with four methods: Change-O DefineClones at the SHazaM threshold, and SCOPer identical, hierarchical and spectral.
6. **Estimate hypermutation** per read with SHazaM.
7. **Score** every method against the known answers (`bin/score_clones.py`), then collect one table over all levels, error rates and methods.

Every tool runs in a container named by digest in `images.lock`. The pipeline refuses to start with an image that is not pinned. partis draws each family's tree with the R package TreeSim, which the published partis image lacks, so the lineage step runs in a small image built on the pinned partis digest that adds R 4.5.3 and TreeSim 2.4 from conda-forge, fixed by build string (`containers/partis-r`). CI rebuilds that small image on every run, so its own digest differs from run to run; the pinned part is the base partis digest in `images.lock` it is built on, and the R and TreeSim builds it adds.

## What is scored

- Clonal grouping: pairwise precision, recall and F1, and the adjusted Rand index (ARI).
- V and J gene calls: share of reads with the true gene.
- Hypermutation: mean absolute error of the per-read mutation frequency.
- Diversity: Hill numbers (q = 0, 1 and 2) of the predicted clone sizes against the true ones.
- Number of clones, true and predicted.

A read the pipeline drops (failed QC or annotation) still counts: as its own one-read clone and as a wrong V and J call. Dropping hard reads therefore cannot raise a score. The same metrics on the retained reads only are reported beside them.

## Two independent generators

partis (optionally with its bcr-phylo mode) is the primary source of truth. A second generator, `bin/expand_lineages.py`, grows naive immuneSIM rearrangements into families on random lineage trees, with mutations weighted toward WRC and GYW hot spots and no stop codons. It is simpler than partis and is reported in its own rows (`--lineage_simulator inrepo`).

The point of having both is that a ranking of methods that holds under two generators built in different ways is stronger evidence than one that holds under one. Where the two disagree, the disagreement is reported as it is. Neither generator is preferred after the fact.

## Prior work and what this adds

The closest published benchmark is Voss, Kaur, Banerjee, Breden and Pennell (2024). As we read the journal article, they simulated 1,200 heavy-chain repertoires with partis (`--simulate-from-scratch`) over 24 parameter settings, including six hypermutation rates from 0.001 to 0.2 per position, with 16 clonal families in the main setup. They compared MiXCR, Change-O, SCOPer hierarchical and SCOPer spectral with a phylogenetic species-delimitation method (mPTP), scored with precision, recall, specificity, F1, the number of families and the error of the median family size, and applied mPTP to a cattle repertoire.

This repository does not repeat their method comparison and does not include MiXCR or mPTP. What it adds:

- **Downstream accuracy, not only grouping.** V and J call accuracy, hypermutation error and diversity error are scored next to the grouping, on the same reads.
- **Sequencing error as its own axis.** Every hypermutation level is run at several error rates, so the effect of read errors on each method is measured directly.
- **ARI beside the pairwise scores.** ARI corrects for chance agreement and is less sensitive to how many families there are.
- **A second, independent generator**, as described above.
- **A human vaccination time course** as the real-data run (planned, below), where Voss et al. used a cattle repertoire.
- **One command, pinned images.** The whole known-answer test reruns from `make known-answer` with every tool in a container fixed by digest.

## Real data (planned)

One public study from the Observed Antibody Space (OAS), which shares its data under CC BY 4.0: Galson et al. (2015), heavy-chain repertoires from 9 people after hepatitis B vaccination, sampled on days 0, 7, 14, 21 and 28. The first run uses one person's five time points, all isotypes. There is no known answer for real data, so this run reports what the methods agree and disagree on, and how clone expansion after vaccination looks under each.

The subset is fixed by rule before any file is read (`bin/oas_subset.py`): the first subject, in plain sort order, whose units cover all five days; isotype units pooled per day; productive rows with a V call, a J call and a junction; one row per distinct nucleotide sequence; and, if a day still has more than 20,000 rows, the 20,000 with the smallest SHA-256 of a fixed salt and the sequence. The cap does not depend on file or row order. Every count before and after each step, and every input file's SHA-256, goes to `manifest.json`. Sampling within a day can split a clone across the cut, so clone sizes in this run are lower bounds.

## Running it

| Command | Needs | What it does |
|---|---|---|
| `make test` | Python 3 | Unit tests, standard library only |
| `make stub` | Nextflow, Java | Runs the pipeline wiring with no tools |
| `make pin` | Docker | Pulls and builds the images and writes their digests to `images.lock`; stops, writing nothing, if a locked image's tag now points at another digest |
| `make reference` | Docker | Builds the OGRDB human IGH germline reference in the suite image |
| `make known-answer` | Docker, Nextflow | The whole known-answer test on the small test profile |

The full run uses 2,000 clonal families per level (`params.n_founders`); the test profile uses 60. Settings are in `nextflow.config`.

## Data and tool licences

- partis and immuneSIM: GPL-3.0. They run in their own containers and are not copied into this repository.
- immuneSIM is no longer on CRAN: it was archived on 2025-06-12 because e-mail to its maintainer was undeliverable (https://CRAN.R-project.org/package=immuneSIM). The image installs its last release, 0.8.7 (packaged 2019-09-20), from the CRAN archive and stops if the file does not match a fixed SHA-256. Its imports come from a CRAN snapshot of 2025-04-10 and Bioconductor 3.20, the last Bioconductor release for R 4.4. See `containers/immunesim/`.
- OGRDB human IGH germline set: CC0.
- OAS data: CC BY 4.0, with credit to the original study.

## References

- Voss K, Kaur K, Banerjee R, Breden F, Pennell M. Applying phylogenetic methods for species delimitation to distinguish B-cell clonal families. Frontiers in Immunology 15:1505032 (2024). https://doi.org/10.3389/fimmu.2024.1505032
- Voss K, Kaur K, Banerjee R, Breden F, Pennell M. Evaluating methods for B-cell clonal family assignment. bioRxiv preprint (2024). https://doi.org/10.1101/2024.05.29.596491
- Galson JD, et al. Analysis of B cell repertoire dynamics following hepatitis B vaccination in humans, and enrichment of vaccine-specific antibody sequences. EBioMedicine 2:2070-2079 (2015). https://doi.org/10.1016/j.ebiom.2015.11.034
