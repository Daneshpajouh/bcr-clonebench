#!/usr/bin/env bash
# Build the germline reference the pipeline reads (reference/), inside the
# Immcantation suite container.
#
#   build_reference.sh imgt  OUT_DIR
#       the IMGT human IGH set the suite image ships (IMGT licence terms apply)
#   build_reference.sh ogrdb OUT_DIR GAPPED_V.fasta D.fasta J.fasta
#       an OGRDB human IGH germline set (CC0): V with IMGT gaps, D and J plain
#   build_reference.sh shuffled-v OUT_DIR SOURCE_DIR SEED
#       mutation test only: SOURCE_DIR's set with every V gene renamed to
#       another gene (shuffle_v_names.py), so V calls must come out wrong
#
# OUT_DIR gets germlines/IGH{V,D,J}.fasta (V with IMGT gaps, as MakeDb and
# CreateGermlines need) and igblast/, a copy of the suite's IgBLAST data
# folder whose human IG databases are rebuilt from those files (headers cut
# to the allele name by igblast_fasta.py), so AssignGenes finds them under
# the names it expects.
set -euo pipefail
SUITE_IGDATA=/usr/local/share/igblast
SUITE_GERM=/usr/local/share/germlines/imgt/human/vdj
kind=${1:?imgt, ogrdb or shuffled-v}
out=${2:?output folder}
mkdir -p "$out/germlines"
case "$kind" in
  imgt)
    cp "$SUITE_GERM/imgt_human_IGHV.fasta" "$out/germlines/IGHV.fasta"
    cp "$SUITE_GERM/imgt_human_IGHD.fasta" "$out/germlines/IGHD.fasta"
    cp "$SUITE_GERM/imgt_human_IGHJ.fasta" "$out/germlines/IGHJ.fasta"
    ;;
  ogrdb)
    cp "${3:?gapped V fasta}" "$out/germlines/IGHV.fasta"
    cp "${4:?D fasta}" "$out/germlines/IGHD.fasta"
    cp "${5:?J fasta}" "$out/germlines/IGHJ.fasta"
    ;;
  shuffled-v)
    src=${3:?source reference folder}
    python3 "$(dirname "$0")/shuffle_v_names.py" "$src/germlines/IGHV.fasta" "$out/germlines/IGHV.fasta" "${4:?seed}"
    cp "$src/germlines/IGHD.fasta" "$out/germlines/IGHD.fasta"
    cp "$src/germlines/IGHJ.fasta" "$out/germlines/IGHJ.fasta"
    ;;
  *) echo "unknown reference kind: $kind" >&2; exit 2 ;;
esac
cp -r "$SUITE_IGDATA" "$out/igblast"
here=$(dirname "$0")
for seg in v d j; do
  up=$(echo "$seg" | tr vdj VDJ)
  # headers cut to the allele name (makeblastdb ids are 50 characters at most)
  python3 "$here/igblast_fasta.py" "$out/germlines/IGH$up.fasta" "$out/igblast/fasta_ungapped_$seg.fasta"
  python3 "$here/igblast_fasta.py" --check "$out/igblast/fasta_ungapped_$seg.fasta"
  makeblastdb -parse_seqids -dbtype nucl -in "$out/igblast/fasta_ungapped_$seg.fasta" \
    -out "$out/igblast/database/imgt_human_ig_$seg" > /dev/null
done
sha256sum "$out"/germlines/*.fasta > "$out/SHA256SUMS"
echo "reference ($kind) in $out"
