#!/usr/bin/env python3
"""Write the ungapped germline FASTA that makeblastdb reads for IgBLAST.

makeblastdb -parse_seqids takes the first word of each header as the
sequence id and refuses ids over 50 characters. IMGT headers carry the
accession, allele, species and more ("M99641|IGHV1-18*01|Homo sapiens|F|..."),
so, as Immcantation's imgt2igblast.sh and IgBLAST's edit_imgt_file.pl do,
each header is cut to the allele name and the IMGT gap dots are removed.
The gapped germlines that MakeDb and CreateGermlines read are not changed.

    igblast_fasta.py IN_GAPPED.fasta OUT_UNGAPPED.fasta
    igblast_fasta.py --check OUT_UNGAPPED.fasta [...]

A name that repeats with the same sequence is kept once; a name that
repeats with a different sequence, an empty name or a name over 50
characters stops the build.
"""
import sys

MAX_ID = 50


def allele_name(header):
    """Allele name from a FASTA header line (without '>')."""
    fields = header.split("|")
    name = fields[1] if len(fields) > 1 else header
    words = name.split()
    return words[0] if words else ""


def read_fasta(path):
    records, name, seq = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    records.append((name, "".join(seq)))
                name, seq = line[1:], []
            elif line.strip():
                seq.append(line.strip())
    if name is not None:
        records.append((name, "".join(seq)))
    return records


def convert(records):
    """[(header, gapped seq)] -> [(allele name, ungapped seq)], checked."""
    out, seen = [], {}
    for header, gapped in records:
        name, seq = allele_name(header), gapped.replace(".", "")
        if not name:
            raise ValueError(f"empty allele name in header: {header!r}")
        if len(name) > MAX_ID:
            raise ValueError(f"allele name over {MAX_ID} characters: {name}")
        if name in seen:
            if seen[name].upper() != seq.upper():
                raise ValueError(f"allele name repeats with a different sequence: {name}")
            continue
        seen[name] = seq
        out.append((name, seq))
    return out


def check(path):
    """Every id in a written FASTA is 50 characters or fewer and unique."""
    names = [h.split()[0] if h.split() else "" for h, _ in read_fasta(path)]
    bad = [n for n in names if not n or len(n) > MAX_ID]
    dup = sorted({n for n in names if names.count(n) > 1})
    if bad or dup or not names:
        raise ValueError(f"{path}: {len(names)} ids, too long or empty {bad[:5]}, repeated {dup[:5]}")
    return len(names)


def main(argv):
    if argv[:1] == ["--check"] and len(argv) > 1:
        for p in argv[1:]:
            print(f"{p}: {check(p)} ids, all {MAX_ID} characters or fewer and unique")
        return 0
    if len(argv) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    with open(argv[1], "w") as fh:
        for name, seq in convert(read_fasta(argv[0])):
            fh.write(f">{name}\n{seq}\n")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except ValueError as err:
        print(f"igblast_fasta.py: {err}", file=sys.stderr)
        sys.exit(1)
