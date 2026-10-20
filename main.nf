#!/usr/bin/env nextflow
/*
 * bcr-clonebench: score B-cell receptor clonal grouping methods against
 * simulated human IGH repertoires whose clonal families are known.
 *
 * clonal families at three mutation levels (partis; or immuneSIM founders
 * grown by bin/expand_lineages.py) -> sequencing errors at each rate ->
 * IgBLAST annotation and Change-O
 * table -> clonal grouping by four methods -> hypermutation by SHazaM ->
 * scores against the known answers (bin/score_clones.py) -> one table.
 */
nextflow.enable.dsl = 2

include { FOUNDERS } from './modules/founders'
include { LINEAGES } from './modules/lineages'
include { LINEAGES_PARTIS } from './modules/lineages_partis'
include { SEQ_ERRORS } from './modules/seq_errors'
include { ANNOTATE } from './modules/annotate'
include { THRESHOLD } from './modules/threshold'
include { CLONES } from './modules/clones'
include { MUTATION } from './modules/mutation'
include { SCORE; SUMMARY } from './modules/score'

def requireDigest(name, image) {
    // a pulled image by registry digest, or a locally built one by image ID
    if (!workflow.stubRun && !(image ==~ /(.+@)?sha256:[0-9a-f]{64}/)) {
        error "params.${name} must name an image by digest (name@sha256:... or sha256:...); see images.lock and 'make pin'"
    }
}

workflow {
    requireDigest('suite_image', params.suite_image)
    levels = Channel.fromList(params.levels)
    if (params.lineage_simulator == 'partis') {
        // the clonal truth comes from partis (Director ruling 2026-10-19 19:52)
        requireDigest('partis_image', params.partis_image)
        lineages = LINEAGES_PARTIS(levels, file("${projectDir}/bin/airr_truth.py")).lineages   // [level, reads, truth]
    } else if (params.lineage_simulator == 'inrepo') {
        // immuneSIM founders grown by bin/expand_lineages.py: an independent second generator
        requireDigest('immunesim_image', params.immunesim_image)
        founders = FOUNDERS(params.n_founders, params.seed,
                           file("${projectDir}/bin/simulate_founders.R"), file("${projectDir}/bin/prepare_founders.py"))
        lineages = LINEAGES(founders.combine(levels))
    } else {
        error "params.lineage_simulator must be partis or inrepo"
    }
    // every level at every sequencing-error rate; labels become <level>_e<rate>
    noisy = SEQ_ERRORS(lineages.combine(Channel.fromList(params.error_rates)))
    annotated = ANNOTATE(noisy.map { label, reads, truth -> [label, reads] }, file(params.reference_dir),
                         file("${projectDir}/bin/drop_reasons.py"))
    thresholds = THRESHOLD(annotated.db, file("${projectDir}/bin/threshold.R"))   // [level, db, threshold file]
    methods = Channel.fromList(params.methods)
    clones = CLONES(thresholds.combine(methods), file("${projectDir}/bin/clones_scoper.R"))   // [level, method, clones]
    mutated = MUTATION(clones, file(params.reference_dir), file("${projectDir}/bin/mutation.R"))  // [level, method, clones with mu_freq]
    truths = noisy.map { label, reads, truth -> [label, truth] }
    scores = SCORE(mutated.combine(truths, by: 0).combine(annotated.drops, by: 0))   // json per level and method
    SUMMARY(scores.collect())
}
