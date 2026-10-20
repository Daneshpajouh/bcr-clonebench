// Clonal grouping by one method: changeo (DefineClones), or SCOPer's
// identical, hierarchical and spectral methods.
// clones_scoper.R is staged into the task folder and run with Rscript, since
// the project's bin folder is not reachable inside a container (CI run 8).
process CLONES {
    container params.suite_image
    tag "${level} ${method}"

    input:
    tuple val(level), path(db), path(threshold), val(method)
    path clones_scoper_r

    output:
    tuple val(level), val(method), path("clones_${level}_${method}.tsv")

    script:
    if (method == 'changeo')
        """
        DefineClones.py -d ${db} --act set --model ham --norm len --dist \$(cat ${threshold}) \\
            --format airr --nproc ${task.cpus} -o clones_${level}_${method}.tsv
        """
    else
        """
        Rscript ${clones_scoper_r} ${db} ${method} \$(cat ${threshold}) clones_${level}_${method}.tsv
        """

    stub:
    """
    awk 'BEGIN{FS=OFS="\\t"} NR==1{print \$0,"clone_id"; next} {split(\$1,a,"_m"); print \$0,a[1]}' ${db} \\
        > clones_${level}_${method}.tsv
    """
}
