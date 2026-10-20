// Scores against the known answers, then one table for all levels and methods.
process SCORE {
    tag "${level} ${method}"

    input:
    tuple val(level), val(method), path(mutated), path(truth)

    output:
    path "score_${level}_${method}.json"

    script:
    """
    score_clones.py --truth ${truth} --pred ${mutated} --mu-column mu_freq --label ${level}/${method} \\
        --out-json score_${level}_${method}.json
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
