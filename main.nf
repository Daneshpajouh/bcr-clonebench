#!/usr/bin/env nextflow
/*
 * bcr-clonebench: score B-cell receptor clonal grouping methods against
 * simulated human IGH repertoires whose clonal families are known.
 *
 * founders (immuneSIM, no hypermutation) -> lineages at three mutation
 * levels (bin/expand_lineages.py) -> IgBLAST annotation and Change-O
 * table -> clonal grouping by four methods -> hypermutation by SHazaM ->
 * scores against the known answers (bin/score_clones.py) -> one table.
 */
nextflow.enable.dsl = 2

include { FOUNDERS } from './modules/founders'
include { LINEAGES } from './modules/lineages'
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
    requireDigest('immunesim_image', params.immunesim_image)

    founders = FOUNDERS(params.n_founders, params.seed)
    levels = Channel.fromList(params.levels)
    lineages = LINEAGES(founders.combine(levels))          // [level, reads, truth]
    annotated = ANNOTATE(lineages.map { level, reads, truth -> [level, reads] }, file(params.reference_dir))
    thresholds = THRESHOLD(annotated)                       // [level, db, threshold file]
    methods = Channel.fromList(params.methods)
    clones = CLONES(thresholds.combine(methods))            // [level, method, clones]
    mutated = MUTATION(clones, file(params.reference_dir))  // [level, method, clones with mu_freq]
    truths = lineages.map { level, reads, truth -> [level, truth] }
    scores = SCORE(mutated.combine(truths, by: 0))          // json per level and method
    SUMMARY(scores.collect())
}
