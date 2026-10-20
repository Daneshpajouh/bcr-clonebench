// Naive human IGH rearrangements from immuneSIM, made into lineage founders.
process FOUNDERS {
    container params.immunesim_image
    tag "n=${n}"

    input:
    val n
    val seed

    output:
    path 'founders.tsv'

    script:
    """
    simulate_founders.R ${n * 2} ${seed} repertoire.tsv v_germline.fasta
    prepare_founders.py --immunesim repertoire.tsv --v-germline v_germline.fasta --out founders.tsv --max ${n}
    """

    stub:
    """
    printf 'founder_id\\tsequence\\tv_call\\td_call\\tj_call\\tv_end\\n' > founders.tsv
    for tail in GCGAGA GCGAAA GCCAGA ACGAGA GCGTGG; do
        printf 'F_%s\\tCAGGTGCAGCTGGTGCAGTCTGGGGCTGAGGTGAAGAAGCCTGGGGCCTCAGTGAAGGTC%s\\tIGHV1-2*02\\tIGHD3-10*01\\tIGHJ4*02\\t54\\n' \$tail \$tail >> founders.tsv
    done
    """
}
