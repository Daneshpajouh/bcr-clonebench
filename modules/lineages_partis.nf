// Clonal families simulated by partis, one set per mutation level, and
// their reads and known answers (bin/airr_truth.py). Mode 'partis' is
// partis simulate from scratch (bio++ mutates along each family's tree);
// mode 'bcr-phylo' runs partis's bcr-phylo-run.py, whose options go in
// params.bcr_phylo_args. The run keeps partis's own --help output beside
// the results, so every option used can be checked against the installed
// version. The partis image keeps its tools in /partis/bin, which is not
// on the image's PATH, so each script adds it and stops at once (exit 127,
// with a message) if the tool is missing. The task's PATH is built by
// Nextflow, not taken from the image's ENV (CI run 7: R was in the image
// but not on the task's PATH), so the scripts also name partis-r's R
// folder and print the PATH they run with. bin/airr_truth.py is staged into
// the task folder as an input and run with python3 (CI run 8: the project's
// bin folder was not reachable inside the partis-r container), and each
// script lists the task folder so the run log shows what was staged.
// partis rearranges from the human IGH germline set it ships, not from a set
// it generates (CI run 9: a generated 55-allele set had alleles that rarely
// give productive rearrangements, and one sub-process kept retrying). This
// partis spells the switch --force-dont-generate-germline-set: bin/partis line
// 1796 adds it to subargs['simulate'] (store_true), and python/processargs.py
// line 472 turns generate_germline_set on whenever rearrange_from_scratch is
// set unless that switch is given (CI run 11 evidence). The script stops at
// once if the installed partis lacks it, after printing the help lines and
// partis source lines that name the germline set. It runs in one process
// (--n-procs 1): with several, partis merged only the .yaml and left the AIRR
// .tsv in its sub-process folders (CI run 12), so the step lists the task
// folder and stops at once if simu_<level>.tsv is missing. This partis build
// writes no AIRR table from simulate at all (CI run 13), so when it is missing
// the step converts the .yaml with partis's own bin/parse-output.py, keeping
// its help text and printing it if the conversion fails. matplotlib
// writes its cache in the task folder (MPLCONFIGDIR).
process LINEAGES_PARTIS {
    container params.partis_image
    tag "${level}"
    // The help text goes into the results too, so a README label can cite the installed version (Director ruling
    // 2026-10-20: options are read from --help in the image before they are named).
    publishDir "${params.outdir}/tool-help", mode: 'copy', pattern: 'partis_help_*.txt'

    input:
    val level
    path airr_truth_py

    output:
    tuple val(level), path("reads_${level}.fasta"), path("truth_${level}.tsv"), emit: lineages
    path "partis_help_${level}.txt", emit: help

    script:
    // <level>_shared (params.shared_vj): the level's mutation, with one V and one J allele and one junction length
    def shared = level.endsWith('_shared')
    def mult = params.mutation_multiplier[level.replaceAll(/_shared$/, '')]
    def shared_args = shared ? "--only-genes '${params.shared_vj_genes}' --allowed-cdr3-lengths ${params.shared_junction_length}" : ''
    def shared_genes = shared ? params.shared_vj_genes.tokenize(':').collect { "'${it}'" }.join(' ') : ''  // quoted: allele names hold '*'
    if (params.lineage_mode == 'bcr-phylo')
        """
        export PATH=/partis/bin:/opt/conda/bin:/opt/conda/envs/r/bin:\$PATH
        echo "PATH=\$PATH" >&2
        ls -la >&2
        test -x /partis/bin/bcr-phylo-run.py || { echo "bcr-phylo-run.py not found at /partis/bin/bcr-phylo-run.py in ${params.partis_image}" >&2; exit 127; }
        export MPLCONFIGDIR=\$PWD/mplconfig; mkdir -p \$MPLCONFIGDIR
        rc=0; bcr-phylo-run.py --help > partis_help_${level}.txt 2>&1 || rc=\$?
        test \$rc -eq 0 || echo "bcr-phylo-run.py --help exited \$rc; help text kept in partis_help_${level}.txt" >&2
        bcr-phylo-run.py --base-outdir bp_${level} --seed ${params.seed} ${params.bcr_phylo_args}
        simu=\$(find bp_${level} -name '*.tsv' | head -n 1)
        test -n "\$simu" || { echo "bcr-phylo-run.py wrote no AIRR table" >&2; exit 1; }
        python3 ${airr_truth_py} --airr "\$simu" --level ${level} --reads reads_${level}.fasta --truth truth_${level}.tsv
        """
    else
        """
        export PATH=/partis/bin:/opt/conda/bin:/opt/conda/envs/r/bin:\$PATH
        echo "PATH=\$PATH" >&2
        ls -la >&2
        test -x /partis/bin/partis || { echo "partis not found at /partis/bin/partis in ${params.partis_image}" >&2; exit 127; }
        command -v R > /dev/null || { echo "R not found on PATH in ${params.partis_image}: partis needs R and TreeSim to draw trees; use the partis-r image (images.lock partis_r_image)" >&2; exit 127; }
        export MPLCONFIGDIR=\$PWD/mplconfig; mkdir -p \$MPLCONFIGDIR
        rc=0; partis simulate --help > partis_help_${level}.txt 2>&1 || rc=\$?
        test \$rc -eq 0 || echo "partis simulate --help exited \$rc; help text kept in partis_help_${level}.txt" >&2
        echo "partis help lines naming germline:" >&2; grep -n -i germline partis_help_${level}.txt >&2 || echo "(none)" >&2
        echo "partis help, first 40 lines:" >&2; head -n 40 partis_help_${level}.txt >&2
        src=\$(grep -rn -e 'generate-germline-set' -e 'generate_germline_set' /partis/bin/partis /partis/python 2>/dev/null | head -n 40)
        echo "partis source lines naming the germline-set option:" >&2; if [ -n "\$src" ]; then printf '%s\\n' "\$src" >&2; else echo "(none)" >&2; fi
        for g in ${shared_genes}; do
            hit=\$(grep -rlF ">\$g" /partis/data/germlines/human/igh/ 2>/dev/null | head -n 1)
            test -n "\$hit" || { echo "\$g is not in partis's shipped human IGH germline set; first names there:" >&2; grep -h '^>' /partis/data/germlines/human/igh/*.fasta 2>/dev/null | head -n 20 >&2; exit 2; }
        done
        { grep -q -- '--force-dont-generate-germline-set' partis_help_${level}.txt || grep -qF "subargs['simulate'].append({'name' : '--force-dont-generate-germline-set'" /partis/bin/partis; } || { echo "partis in ${params.partis_image} has no --force-dont-generate-germline-set in partis_help_${level}.txt or in bin/partis's simulate options (lines above); cannot pin the shipped germline set" >&2; exit 2; }
        partis simulate --simulate-from-scratch --force-dont-generate-germline-set --n-sim-events ${params.n_founders} \\
            --n-leaves ${params.mean_family_size} --mutation-multiplier ${mult} --seed ${params.seed} ${shared_args} \\
            --n-procs 1 --outfname simu_${level}.yaml --airr-output
        if [ ! -s simu_${level}.tsv ]; then
            help_rc=0; parse-output.py --help > parse_output_help_${level}.txt 2>&1 || help_rc=\$?
            rc=0; parse-output.py simu_${level}.yaml simu_${level}.tsv --airr-output >&2 || rc=\$?
            test -s simu_${level}.tsv || { ls -la >&2; echo "parse-output.py (exit \$rc) wrote no simu_${level}.tsv (listing above); its help (exit \$help_rc) follows" >&2; cat parse_output_help_${level}.txt >&2; exit 1; }
        fi
        ls -la >&2
        test -s simu_${level}.tsv || { echo "partis simulate wrote no simu_${level}.tsv (AIRR) in the task folder (listing above)" >&2; exit 1; }
        python3 ${airr_truth_py} --airr simu_${level}.tsv --level ${level} --reads reads_${level}.fasta --truth truth_${level}.tsv
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
    python3 ${airr_truth_py} --airr simu.tsv --level ${level} --reads reads_${level}.fasta --truth truth_${level}.tsv
    """
}
