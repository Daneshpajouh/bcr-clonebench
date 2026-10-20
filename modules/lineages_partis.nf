// Clonal families simulated by partis, one set per mutation level, and
// their reads and known answers (bin/airr_truth.py). Mode 'partis' is
// partis simulate from scratch (bio++ mutates along each family's tree);
// mode 'bcr-phylo' runs partis's bcr-phylo-run.py, whose options go in
// params.bcr_phylo_args. The run keeps partis's own --help output beside
// the results, so every option used can be checked against the installed
// version.
process LINEAGES_PARTIS {
    container params.partis_image
    tag "${level}"

    input:
    val level

    output:
    tuple val(level), path("reads_${level}.fasta"), path("truth_${level}.tsv"), emit: lineages
    path "partis_help_${level}.txt", emit: help

    script:
    def mult = params.mutation_multiplier[level]
    if (params.lineage_mode == 'bcr-phylo')
        """
        bcr-phylo-run.py --help > partis_help_${level}.txt 2>&1 || true
        bcr-phylo-run.py --base-outdir bp_${level} --seed ${params.seed} ${params.bcr_phylo_args}
        simu=\$(find bp_${level} -name '*.tsv' | head -n 1)
        test -n "\$simu" || { echo "bcr-phylo-run.py wrote no AIRR table" >&2; exit 1; }
        airr_truth.py --airr "\$simu" --level ${level} --reads reads_${level}.fasta --truth truth_${level}.tsv
        """
    else
        """
        partis simulate --help > partis_help_${level}.txt 2>&1 || true
        partis simulate --simulate-from-scratch --n-sim-events ${params.n_founders} \\
            --n-leaves ${params.mean_family_size} --mutation-multiplier ${mult} --seed ${params.seed} \\
            --outfname simu_${level}.yaml --airr-output
        airr_truth.py --airr simu_${level}.tsv --level ${level} --reads reads_${level}.fasta --truth truth_${level}.tsv
        """

    stub:
    """
    echo "stub: partis not run" > partis_help_${level}.txt
    printf 'sequence_id\\tsequence\\tclone_id\\tv_call\\td_call\\tj_call\\tsequence_alignment\\tgermline_alignment\\tv_alignment_end\\n' > simu.tsv
    c=0
    for tail in GCA GCC GCT; do
        c=\$((c + 1))
        for m in 1 2 3; do
            seq=CAGGTGCAGCTGGTGCAGTCTGGGGCT\$tail
            printf 'c%s_m%s\\t%s\\tc%s\\tIGHV1-2*02\\tIGHD3-10*01\\tIGHJ4*02\\t%s\\tCAGGTGCAGCTGGTGCAGTCTGGGGCTGCG\\t27\\n' \$c \$m \$seq \$c \$seq >> simu.tsv
        done
    done
    airr_truth.py --airr simu.tsv --level ${level} --reads reads_${level}.fasta --truth truth_${level}.tsv
    """
}
