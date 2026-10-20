// Scores against the known answers, then one table for all levels and methods.
// CI run 17 review, item (4): one row per truth read beside the scores, published under reads/.
process SCORE {
    tag "${level} ${method}"
    publishDir "${params.outdir}/reads", mode: 'copy', pattern: 'reads_*.tsv'

    input:
    tuple val(level), val(method), path(mutated), path(truth), path(drops)

    output:
    path "score_${level}_${method}.json", emit: json
    path "reads_${level}_${method}.tsv", emit: reads

    script:
    """
    score_clones.py --truth ${truth} --pred ${mutated} --mu-column mu_freq --label ${level}/${method} \\
        --drops ${drops} --out-json score_${level}_${method}.json --out-reads reads_${level}_${method}.tsv
    """
}

process SUMMARY {
    publishDir params.outdir, mode: 'copy'

    input:
    path scores

    output:
    path 'results.tsv'
    path 'results.md'

    script:
    """
    summarize_scores.py --out-tsv results.tsv --out-md results.md ${scores}
    """
}
