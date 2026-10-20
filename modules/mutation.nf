// Germlines per clone with the D region masked, then the observed
// mutation frequency over the IMGT V region (SHazaM).
// mutation.R is staged into the task folder and run with Rscript, since the
// project's bin folder is not reachable inside a container (CI run 8).
process MUTATION {
    container params.suite_image
    tag "${level} ${method}"

    input:
    tuple val(level), val(method), path(clones)
    path reference
    path mutation_r

    output:
    tuple val(level), val(method), path("mutated_${level}_${method}.tsv")

    script:
    """
    CreateGermlines.py -d ${clones} -g dmask --cloned --format airr \\
        -r ${reference}/germlines/IGHV.fasta ${reference}/germlines/IGHD.fasta ${reference}/germlines/IGHJ.fasta \\
        --outname germ
    Rscript ${mutation_r} germ_germ-pass.tsv mutated_${level}_${method}.tsv
    """

    stub:
    """
    awk 'BEGIN{FS=OFS="\\t"} NR==1{print \$0,"mu_freq"; next} {print \$0,"0.05"}' ${clones} > mutated_${level}_${method}.tsv
    """
}
