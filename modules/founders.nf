// Naive human IGH rearrangements from immuneSIM, made into lineage founders.
// Both bin scripts are staged into the task folder and run with their
// interpreter, since the project's bin folder is not reachable inside a
// container (CI run 8, the partis step).
process FOUNDERS {
    container params.immunesim_image
    tag "n=${n}"

    input:
    val n
    val seed
    path simulate_founders_r
    path prepare_founders_py

    output:
    path 'founders.tsv'

    script:
    """
    Rscript ${simulate_founders_r} ${n * 2} ${seed} repertoire.tsv v_germline.fasta
    python3 ${prepare_founders_py} --immunesim repertoire.tsv --v-germline v_germline.fasta --out founders.tsv --max ${n}
    """

    stub:
    """
    printf 'founder_id\\tsequence\\tv_call\\td_call\\tj_call\\tv_end\\n' > founders.tsv
    for tail in GCGAGA GCGAAA GCCAGA ACGAGA GCGTGG; do
        printf 'F_%s\\tCAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTC%s\\tIGHV1-2*02\\tIGHD3-10*01\\tIGHJ4*02\\t54\\n' \$tail \$tail >> founders.tsv
    done
    """
}
