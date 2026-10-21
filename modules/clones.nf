// Clonal grouping by one method: changeo (DefineClones), or SCOPer's
// identical, hierarchical, spectral (junction only) and spectral_vj methods.
// clones_scoper.R is staged into the task folder and run with Rscript, since
// the project's bin folder is not reachable inside a container (CI run 8).
// changeo also keeps DefineClones.py --help in the results (tool-help/), so
// the README's linkage label is read from the installed version (Director
// ruling 2026-10-20).
process CLONES {
    container params.suite_image
    tag "${level} ${method}"
    publishDir "${params.outdir}/tool-help", mode: 'copy', pattern: 'DefineClones_help_*.txt'

    input:
    tuple val(level), path(db), path(threshold), val(method)
    path clones_scoper_r

    output:
    tuple val(level), val(method), path("clones_${level}_${method}.tsv"), emit: clones
    path "DefineClones_help_${level}.txt", optional: true, emit: help

    script:
    if (method == 'changeo')
        """
        rc=0; DefineClones.py --help > DefineClones_help_${level}.txt 2>&1 || rc=\$?
        test \$rc -eq 0 || echo "DefineClones.py --help exited \$rc; help text kept in DefineClones_help_${level}.txt" >&2
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
    if [ "${method}" = changeo ]; then echo "stub: DefineClones not run" > DefineClones_help_${level}.txt; fi
    """
}
