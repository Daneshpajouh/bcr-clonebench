// Clonal families with known hypermutation, one set per mutation level.
process LINEAGES {
    tag "${level}"

    input:
    tuple path(founders), val(level)

    output:
    tuple val(level), path("reads_${level}.fasta"), path("truth_${level}.tsv")

    script:
    """
    expand_lineages.py --founders ${founders} --level ${level} --seed ${params.seed} \\
        --mean-family-size ${params.mean_family_size} --error-rate ${params.error_rate} \\
        --reads reads_${level}.fasta --truth truth_${level}.tsv
    """

    stub:
    """
    expand_lineages.py --founders ${founders} --level ${level} --seed 1 \\
        --reads reads_${level}.fasta --truth truth_${level}.tsv
    """
}
