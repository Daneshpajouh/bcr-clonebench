#!/usr/bin/env bash
# Build the germline reference the pipeline reads (reference/), inside the
# Immcantation suite container.
#
#   build_reference.sh imgt  OUT_DIR
#       the IMGT human IGH set the suite image ships (IMGT licence terms apply)
#   build_reference.sh ogrdb OUT_DIR GAPPED_V.fasta D.fasta J.fasta
#       an OGRDB human IGH germline set (CC0): V with IMGT gaps, D and J plain
#
# OUT_DIR gets germlines/IGH{V,D,J}.fasta (V with IMGT gaps, as MakeDb and
# CreateGermlines need) and igblast/, a copy of the suite's IgBLAST data
# folder whose human IG databases are rebuilt from those files, so
# AssignGenes finds them under the names it expects.
set -euo pipefail
SUITE_IGDATA=/usr/local/share/igblast
SUITE_GERM=/usr/local/share/germlines/imgt/human/vdj
kind=${1:?imgt or ogrdb}
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
  *) echo "unknown reference kind: $kind" >&2; exit 2 ;;
esac
cp -r "$SUITE_IGDATA" "$out/igblast"
for seg in v d j; do
  up=$(echo "$seg" | tr vdj VDJ)
  tr -d '.' < "$out/germlines/IGH$up.fasta" > "$out/igblast/fasta_ungapped_$seg.fasta"
  makeblastdb -parse_seqids -dbtype nucl -in "$out/igblast/fasta_ungapped_$seg.fasta" \
    -out "$out/igblast/database/imgt_human_ig_$seg" > /dev/null
done
sha256sum "$out"/germlines/*.fasta > "$out/SHA256SUMS"
echo "reference ($kind) in $out"
