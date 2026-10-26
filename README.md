# bcr-clonebench

bcr-clonebench checks how well B-cell receptor (BCR) repertoire tools recover the truth. It simulates human heavy-chain (IGH) repertoires in which every read's clonal family, germline genes and mutations are known, runs them through a standard Immcantation pipeline, and scores each clonal grouping method against those known answers. It is for anyone who groups BCR sequences into clones and wants to know which method to trust at which mutation level, and how much sequencing error changes the answer.

It is a benchmark harness, not a new grouping method.

## Status

What runs today:

- Unit tests for the scorer, the truth converters, the second lineage generator, sequencing-error injection and the summary table (`make test`), and the whole Nextflow pipeline in stub mode (`make stub`), which checks the wiring without running any tool.
- `make known-answer` in CI on every change, in the pinned containers, on the small `test` profile only: 60 founder families, the low, medium and high hypermutation levels, sequencing error rates 0 and 0.005, with shared and with separate V and J genes. It ends with the named checks and the mutation tests below. These runs check the harness; they are not the benchmark, and their numbers are not reported here as results. The first run with all three levels passed on 2026-10-25: [run 40251896057](https://github.com/Daneshpajouh/bcr-clonebench/actions/runs/40251896057) on commit bb6adb1, jobs unit, pipeline-stub and known-answer.

Not done yet, so this repository has no benchmark results:

- The full-size run with the default parameters has not finished. It is `make full`, run by hand in CI (the workflow `full`).
- The real-data run on Galson et al. (2015), below, has not started.
- There is no written report yet; the output is the summary table `results/results.tsv`. A Quarto report is planned.
- Diversity is scored by this repository's own Hill numbers (q = 0, 1 and 2, `bin/score_clones.py`). Alakazam's diversity curves are not in the pipeline.
- The CI runs up to 2026-10-21 annotated against the IMGT set that the suite image ships. The OGRDB set is now the default (`make reference`); no run on it has finished yet.
- Truth allele names come from the simulator's own germline set. A true allele that the OGRDB set lacks, or names differently, counts as a wrong V or J call, so part of a V or J miss may come from naming alone. This repository does not yet count how many true alleles are absent from the reference.

Results will be added here when they exist, with the exact command and image digests that produced them.

## How it works

1. **Simulate clonal families.** partis simulates human IGH clonal families from scratch at three hypermutation levels (low, medium, high), and writes the true family, genes and mutation state of every read. With `--shared_vj true` each level also runs in a harder setting, `<level>_shared`, where every family uses the same V allele, the same J allele and the same junction length (partis `--only-genes` and `--allowed-cdr3-lengths`, defaults `IGHV3-23*01:IGHJ4*02` and 48 nucleotides), so unrelated families look alike and a method can wrongly merge them.
2. **Add sequencing errors.** Each level is copied at each sequencing-error rate (0, 0.1 and 0.5 percent per base by default). The errors are kept out of the truth: they are noise the pipeline has to see through.
3. **Annotate.** IgBLAST assigns V, D and J genes against the OGRDB human IGH germline set (version 10, pinned in `reference.lock`), and Change-O writes an AIRR table.
4. **Pick a clustering threshold.** SHazaM estimates the distance threshold from the data. If it finds none, a fixed fallback is used and the run says so.
5. **Group into clones** with five methods: Change-O DefineClones at the SHazaM threshold, and SCOPer identical, hierarchical, spectral on the junction alone (`spectral`, SCOPer method `novj`) and spectral with the mutations shared in V and J (`spectral_vj`, method `vj`). Both spectral rows are kept because, in the shared V/J setting, the junction-only call puts nearly all reads in a few clones; the second row shows whether that is the method or the option. DefineClones and SCOPer hierarchical are one method in two implementations: single-linkage clustering of junctions by length-normalized Hamming distance at the same threshold (DefineClones' `--link` defaults to single in the pinned suite image, as its own `--help` says; the SCOPer call sets `linkage = "single"`). Their agreement is a check on the pipeline, not a second result.
6. **Estimate hypermutation** per read with SHazaM.
7. **Score** every method against the known answers (`bin/score_clones.py`), then collect one table over all levels, error rates and methods.

Each run keeps the `--help` text of DefineClones and of partis in `results/tool-help/`, read from the images it used, so every option named here can be checked against the installed version.

Every tool runs in a container named by digest in `images.lock`. The pipeline refuses to start with an image that is not pinned. partis draws each family's tree with the R package TreeSim, which the published partis image lacks, so the lineage step runs in a small image built on the pinned partis digest that adds R 4.5.3 and TreeSim 2.4 from conda-forge, fixed by build string (`containers/partis-r`). CI rebuilds that small image on every run, so its own digest differs from run to run; the pinned part is the base partis digest in `images.lock` it is built on, and the R and TreeSim builds it adds.

## What is scored

- Clonal grouping: pairwise precision, recall and F1, and the adjusted Rand index (ARI).
- V and J gene calls: share of reads whose call is exactly the true gene, and share whose call includes the true gene among tied genes (`v_gene_includes_true`, `j_gene_includes_true`). The two differ when the aligner cannot tell two genes apart. In the shared setting, the CI run of 2026-10-20 called IGHV3-23 together with IGHV3-23D on almost every read, so exact V was near 0 while the true gene was in every kept read's call.
- Hypermutation: mean absolute error of the per-read mutation frequency.
- Diversity: Hill numbers (q = 0, 1 and 2) of the predicted clone sizes against the true ones.
- Number of clones, true and predicted.

A read the pipeline drops (failed QC or annotation) still counts: as its own one-read clone and as a wrong V and J call. Dropping hard reads therefore cannot raise a score. The same metrics on the retained reads only are reported beside them.

Beside the summary table, `results/reads/` holds one table per level, error rate and method with a row per simulated read: its true and called clone, V and J, whether each call names exactly the true gene and whether the call set includes it, and whether the read was kept or, if dropped, why. The rows add up to the summary scores, so any figure in the table can be traced to the reads behind it.

The headline scores cover every simulated read. Where the simulator marks each read productive or not (partis does), the same clone and gene scores are also reported over the productive reads alone (n_truth_productive, pair_f1_productive, ari_productive). A productive read the pipeline drops still counts as its own clone there too.

## Named checks and mutation tests

`make known-answer` ends with eight named checks on the summary table (`bin/check_known_answer.py`, floors in `conf/known_answer_floors.json`): the true V and J genes are in the calls of the kept reads, exact V calls outside the shared setting, the pair precision of Change-O, hierarchical and identical grouping, Change-O's F1 on productive reads, and the hypermutation error. Each floor was set from CI run 20 of 2026-10-21 (test profile) with a margin, and the file records the value it was set from. The spectral rows have no floor: their collapse in the shared setting is a result, not a fault.

`make mutation-tests` then checks that the checks can fail. It reruns the pipeline twice from the same simulation, once with a deliberately wrong clone threshold (0.5) and once against a deliberately wrong germline set (every V gene renamed to another gene, seeded), and passes only if KA-CHANGEO-PRECISION and KA-V-INCLUDES fail in turn.

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
| `make pin-reference` | Python 3, network | Downloads the OGRDB set named in `reference.lock`, checks it against the MD5 Zenodo publishes and writes its SHA-256 there if not yet pinned (`make pin` runs it too); a pinned SHA-256 is never changed |
| `make reference` | Docker, network | Downloads the pinned OGRDB human IGH germline set (version 10), stops unless its MD5 and SHA-256 match `reference.lock`, and builds the reference in the suite image; `REFERENCE=imgt` uses the IMGT set the suite image ships instead |
| `make known-answer` | Docker, Nextflow | The whole known-answer test on the small test profile, then the named checks |
| `make mutation-tests` | Docker, Nextflow | After `make known-answer`: a wrong threshold and a wrong germline set must each fail their named check |
| `make full` | Docker, Nextflow | The full-size run with the default parameters into `results-full`, then the same named checks |

The full run uses 2,000 clonal families per level (`params.n_founders`); the test profile uses 60. Settings are in `nextflow.config`.

## Data and tool licences

- partis and immuneSIM: GPL-3.0. They run in their own containers and are not copied into this repository.
- immuneSIM is no longer on CRAN: it was archived on 2025-06-12 because e-mail to its maintainer was undeliverable (https://CRAN.R-project.org/package=immuneSIM). The image installs its last release, 0.8.7 (packaged 2019-09-20), from the CRAN archive and stops if the file does not match a fixed SHA-256. Its imports come from a CRAN snapshot of 2025-04-10 and Bioconductor 3.20, the last Bioconductor release for R 4.4. See `containers/immunesim/`.
- OGRDB human IGH germline set, IGH_VDJ version 10 (2026-05-27; Zenodo, https://doi.org/10.5281/zenodo.20409587): CC0 1.0, as the record's licence.txt states. It is downloaded at run time, not copied into this repository.
- OAS data: CC BY 4.0, with credit to the original study.

## References

- Voss K, Kaur K, Banerjee R, Breden F, Pennell M. Applying phylogenetic methods for species delimitation to distinguish B-cell clonal families. Frontiers in Immunology 15:1505032 (2024). https://doi.org/10.3389/fimmu.2024.1505032
- Voss K, Kaur K, Banerjee R, Breden F, Pennell M. Evaluating methods for B-cell clonal family assignment. bioRxiv preprint (2024). https://doi.org/10.1101/2024.05.29.596491
- Galson JD, et al. Analysis of B cell repertoire dynamics following hepatitis B vaccination in humans, and enrichment of vaccine-specific antibody sequences. EBioMedicine 2:2070-2079 (2015). https://doi.org/10.1016/j.ebiom.2015.11.034
