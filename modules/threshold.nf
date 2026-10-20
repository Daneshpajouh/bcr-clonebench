// The Hamming distance threshold for Change-O and hierarchical SCOPer,
// from the nearest-neighbour distance distribution (SHazaM).
process THRESHOLD {
    container params.suite_image
    tag "${level}"

    input:
    tuple val(level), path(db)

    output:
    tuple val(level), path(db), path("threshold_${level}.txt")

    script:
    """
    threshold.R ${db} threshold_${level}.txt ${params.fallback_threshold}
    """

    stub:
    """
    echo ${params.fallback_threshold} > threshold_${level}.txt
    """
}
