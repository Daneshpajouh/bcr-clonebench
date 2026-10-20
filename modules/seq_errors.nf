// Independent sequencing errors at one rate on one level's reads; the
// known answers are unchanged. The label becomes <level>_e<rate>.
process SEQ_ERRORS {
    tag "${level} e=${rate}"

    input:
    tuple val(level), path(reads), path(truth), val(rate)

    output:
    tuple val("${level}_e${rate}"), path("reads_${level}_e${rate}.fasta"), path(truth)

    script:
    """
    add_seq_errors.py --reads ${reads} --rate ${rate} --seed ${params.seed} --out reads_${level}_e${rate}.fasta
    """
}
