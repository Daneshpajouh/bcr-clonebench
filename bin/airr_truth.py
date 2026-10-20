#!/usr/bin/env python3
"""Reads and known answers from a simulator's true AIRR table.

    airr_truth.py --airr simu.tsv --level medium --reads reads.fasta --truth truth.tsv

simu.tsv is the AIRR rearrangement table a simulator writes for its own
true annotations (partis simulate with AIRR output): one row per
simulated sequence with sequence_id, sequence, clone_id, v_call, d_call,
j_call, sequence_alignment and germline_alignment, plus the V stretch:
either v_alignment_end (the AIRR position in sequence_alignment) or, as
partis writes it, v_sequence_alignment with v_germline_alignment.

true_mu_freq is the share of V positions (gaps and germline N skipped)
where the sequence differs from its germline: the stretch the pipeline's
mutation step scores. partis's AIRR table (partis/utils.py, airr_headers
and get_airr_line) has no v_alignment_end column; its
v_sequence_alignment and v_germline_alignment are the V slice of
sequence_alignment and germline_alignment, cut at the same place, so
they are read directly and checked to sit at one offset in both full
alignments. A row with neither form is refused, never given a guessed
answer, since score_clones.py scores against true_mu_freq.

true_productive copies the simulator's own productive flag (T or F;
empty when absent), so the score step can also report the reads the
simulator calls productive.

Standard library only.
"""

import argparse
import csv
import sys

NEED = ("sequence_id", "sequence", "clone_id", "v_call", "j_call", "sequence_alignment", "germline_alignment")
TRUTH = ("sequence_id", "true_clone_id", "true_v_call", "true_d_call", "true_j_call", "true_mu_freq", "level",
         "true_productive")
TRUE, FALSE = {"T", "TRUE", "1"}, {"F", "FALSE", "0"}


def v_stretch(row):
    """The V part of the alignment as (sequence, germline, route)."""
    seq, germ = row["sequence_alignment"].upper(), row["germline_alignment"].upper()
    if len(seq) != len(germ):
        raise ValueError(f"{row['sequence_id']}: sequence and germline alignments differ in length")
    end = (row.get("v_alignment_end") or "").strip()
    if end:
        if not end.isdigit():
            raise ValueError(f"{row['sequence_id']}: v_alignment_end {end!r} is not a position")
        return seq[:int(end)], germ[:int(end)], "v_alignment_end"
    v_seq = (row.get("v_sequence_alignment") or "").strip().upper()
    v_germ = (row.get("v_germline_alignment") or "").strip().upper()
    if not (v_seq and v_germ):
        raise ValueError(f"{row['sequence_id']}: no v_alignment_end and no v_sequence_alignment"
                         " with v_germline_alignment")
    if len(v_seq) != len(v_germ):
        raise ValueError(f"{row['sequence_id']}: v_sequence_alignment and v_germline_alignment differ in length")
    at = seq.find(v_seq)
    while at >= 0 and germ[at:at + len(v_germ)] != v_germ:
        at = seq.find(v_seq, at + 1)
    if at < 0:
        raise ValueError(f"{row['sequence_id']}: v_sequence_alignment and v_germline_alignment"
                         " are not one slice of sequence_alignment and germline_alignment")
    return v_seq, v_germ, "v_sequence_alignment"


def v_mu(row, routes=None):
    seq, germ, route = v_stretch(row)
    if routes is not None:
        routes[route] = routes.get(route, 0) + 1
    pairs = [(s, g) for s, g in zip(seq, germ) if s not in ".-" and g not in ".-N"]
    if not pairs:
        raise ValueError(f"{row['sequence_id']}: empty V alignment")
    return sum(s != g for s, g in pairs) / len(pairs)


def productive(row):
    """T or F from the simulator's own productive column (partis: is_functional, so no mutated conserved codon,
    in frame and no stop codon); empty when the column is absent or its value is neither, never guessed."""
    value = (row.get("productive") or "").strip().upper()
    return "T" if value in TRUE else "F" if value in FALSE else ""


def first_call(call):
    return (call or "").split(",")[0].strip()


def convert(rows, level, routes=None):
    seen = set()
    for r in rows:
        missing = [c for c in NEED if not r.get(c)]
        if missing:
            raise ValueError(f"{r.get('sequence_id', '?')}: missing {', '.join(missing)}")
        if r["sequence_id"] in seen:
            raise ValueError(f"{r['sequence_id']}: repeated sequence_id")
        seen.add(r["sequence_id"])
        yield r["sequence_id"], r["sequence"].upper(), {
            "sequence_id": r["sequence_id"], "true_clone_id": r["clone_id"], "true_v_call": first_call(r["v_call"]),
            "true_d_call": first_call(r.get("d_call")), "true_j_call": first_call(r["j_call"]),
            "true_mu_freq": f"{v_mu(r, routes):.6f}", "level": level, "true_productive": productive(r)}


def show_head(path, n=2):
    """On a failure, print the table's header line and first n data rows to stderr."""
    try:
        with open(path, encoding="utf-8", errors="replace") as f:
            lines = [f.readline() for _ in range(n + 1)]
    except OSError:
        return
    print(f"airr_truth: {path}: header and first {n} data rows follow", file=sys.stderr)
    for line in lines:
        if line:
            sys.stderr.write(line if line.endswith("\n") else line + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--airr", required=True)
    ap.add_argument("--level", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--truth", required=True)
    a = ap.parse_args(argv)
    routes = {}
    try:
        with open(a.airr, newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        if not rows:
            raise ValueError(f"{a.airr}: no rows")
        out = list(convert(rows, a.level, routes))
    except (OSError, ValueError) as e:
        print(f"airr_truth: {e}", file=sys.stderr)
        show_head(a.airr)
        return 2
    with open(a.reads, "w", encoding="utf-8") as fa, open(a.truth, "w", newline="", encoding="utf-8") as ft:
        w = csv.DictWriter(ft, fieldnames=TRUTH, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for sid, seq, truth in out:
            fa.write(f">{sid}\n{seq}\n")
            w.writerow(truth)
    print(f"{len(out)} reads, {len({t['true_clone_id'] for _, _, t in out})} families, level {a.level}")
    print("airr_truth: V stretch from " + ", ".join(f"{k} in {n} rows" for k, n in sorted(routes.items())),
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
