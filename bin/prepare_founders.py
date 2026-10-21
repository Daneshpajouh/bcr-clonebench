#!/usr/bin/env python3
"""Turn immuneSIM's naive rearrangements into lineage founders.

    prepare_founders.py --immunesim repertoire.tsv --v-germline v.fasta \
        --out founders.tsv [--max N]

repertoire.tsv is immuneSIM's AIRR table (run with shm.mode "none", so
each sequence is its germline V(D)J). v.fasta holds the germline V
alleles immuneSIM drew from, named as in its v_call column.

For each rearrangement this keeps the sequence in frame from its first
base (a trailing partial codon is cut), takes the first V call, and
sets v_end: the end of the conserved cysteine codon (IMGT 104) of the
germline V, the last in-frame TGT or TGC in the final 30 bases of the
allele. The V region scored for hypermutation is everything before
v_end, the stretch SHazaM's IMGT V region definition also scores.

Dropped, and counted on stderr: sequences with a stop codon in frame,
with characters other than ACGT, whose start does not match their
germline V up to v_end, or that repeat an earlier sequence.

Standard library only.
"""

import argparse
import csv
import sys

STOPS = {"TAA", "TAG", "TGA"}


def read_fasta(path):
    seqs, name = {}, None
    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line.upper().replace(".", ""))  # IMGT gaps, if any
    return {k: "".join(v) for k, v in seqs.items()}


def cys_end(v):
    """End (exclusive) of the conserved cysteine codon of a germline V, in
    frame from its first base; None when there is none near the end."""
    for start in range(len(v) - len(v) % 3 - 3, max(-1, len(v) - 33), -3):
        if v[start:start + 3] in ("TGT", "TGC"):
            return start + 3
    return None


def prepare(rows, germ, max_n=None):
    out, seen = [], set()
    dropped = {"stop": 0, "not ACGT": 0, "no germline": 0, "start differs from germline": 0, "repeat": 0}
    for r in rows:
        seq = r["sequence"].upper()
        seq = seq[:len(seq) - len(seq) % 3]
        v = r["v_call"].split(",")[0].strip()
        if set(seq) - set("ACGT"):
            dropped["not ACGT"] += 1
            continue
        if any(seq[i:i + 3] in STOPS for i in range(0, len(seq), 3)):
            dropped["stop"] += 1
            continue
        g = germ.get(v)
        end = cys_end(g) if g else None
        if end is None:
            dropped["no germline"] += 1
            continue
        if seq[:end] != g[:end]:
            dropped["start differs from germline"] += 1
            continue
        if seq in seen:
            dropped["repeat"] += 1
            continue
        seen.add(seq)
        out.append({"founder_id": f"F{len(out) + 1:05d}", "sequence": seq, "v_call": v,
                    "d_call": r.get("d_call", "").split(",")[0].strip(), "j_call": r["j_call"].split(",")[0].strip(),
                    "v_end": end})
        if max_n and len(out) >= max_n:
            break
    return out, dropped


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--immunesim", required=True)
    ap.add_argument("--v-germline", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--max", type=int, help="keep at most this many founders")
    a = ap.parse_args(argv)
    try:
        germ = read_fasta(a.v_germline)
        with open(a.immunesim, newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        need = {"sequence", "v_call", "j_call"} - set(rows[0] if rows else ())
        if not rows or need:
            raise ValueError(f"{a.immunesim}: " + ("no rows" if not rows else f"missing {', '.join(sorted(need))}"))
    except (OSError, ValueError) as e:
        print(f"prepare_founders: {e}", file=sys.stderr)
        return 2
    founders, dropped = prepare(rows, germ, a.max)
    if not founders:
        print(f"prepare_founders: no founder kept; dropped {dropped}", file=sys.stderr)
        return 1
    with open(a.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(founders[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(founders)
    print(f"{len(founders)} founders; dropped " + ", ".join(f"{k} {n}" for k, n in dropped.items() if n),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
