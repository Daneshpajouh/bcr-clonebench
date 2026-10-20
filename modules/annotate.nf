// IgBLAST against the germline reference, then a Change-O AIRR table of
// productive rearrangements.
process ANNOTATE {
    container params.suite_image
    tag "${level}"

    input:
    tuple val(level), path(reads)
    path reference

    output:
    tuple val(level), path("db_${level}.tsv")

    script:
    """
    AssignGenes.py igblast -s ${reads} -b ${reference}/igblast --organism human --loci ig \\
        --format blast --nproc ${task.cpus} -o igblast_${level}.fmt7
    MakeDb.py igblast -i igblast_${level}.fmt7 -s ${reads} \\
        -r ${reference}/germlines/IGHV.fasta ${reference}/germlines/IGHD.fasta ${reference}/germlines/IGHJ.fasta \\
        --extended --format airr -o all_${level}.tsv
    ParseDb.py select -d all_${level}.tsv -f productive -u T --outname productive_${level}
    mv productive_${level}_parse-select.tsv db_${level}.tsv
    """

    stub:
    """
    awk 'BEGIN{OFS="\\t"; print "sequence_id","sequence","v_call","j_call","junction","productive"}
         /^>/{id=substr(\$0,2); next} {print id,\$0,"IGHV1-2*02","IGHJ4*02","TGTGCGAGAGGG","T"}' ${reads} > db_${level}.tsv
    """
}
