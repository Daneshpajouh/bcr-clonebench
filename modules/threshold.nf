// The Hamming distance threshold for Change-O and hierarchical SCOPer,
// from the nearest-neighbour distance distribution (SHazaM).
// threshold.R is staged into the task folder and run with Rscript, since the
// project's bin folder is not reachable inside a container (CI run 8).
process THRESHOLD {
    container params.suite_image
    tag "${level}"

    input:
    tuple val(level), path(db)
    path threshold_r

    output:
    tuple val(level), path(db), path("threshold_${level}.txt")

    script:
    """
    Rscript ${threshold_r} ${db} threshold_${level}.txt ${params.fallback_threshold}
    """

    stub:
    """
    echo ${params.fallback_threshold} > threshold_${level}.txt
    """
}
