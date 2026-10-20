#!/usr/bin/env python3
"""Reads and known answers from a simulator's true AIRR table.

    airr_truth.py --airr simu.tsv --level medium --reads reads.fasta --truth truth.tsv

simu.tsv is the AIRR rearrangement table a simulator writes for its own
true annotations (partis simulate with AIRR output): one row per
simulated sequence with sequence_id, sequence, clone_id, v_call, d_call,
j_call, sequence_alignment, germline_alignment and v_alignment_end (or
v_sequence_end, counted on the sequence).

true_mu_freq is the share of V positions in the alignment (up to
v_alignment_end, gaps skipped) where the sequence differs from its
germline: the stretch the pipeline's mutation step scores. A row without
the alignments is refused, never given a guessed answer.

Standard library only.
"""

import argparse
import csv
import sys

NEED = ("sequence_id", "sequence", "clone_id", "v_call", "j_call", "sequence_alignment", "germline_alignment")
TRUTH = ("sequence_id", "true_clone_id", "true_v_call", "true_d_call", "true_j_call", "true_mu_freq", "level")


def v_mu(row):
    seq, germ = row["sequence_alignment"].upper(), row["germline_alignment"].upper()
    if len(seq) != len(germ):
        raise ValueError(f"{row['sequence_id']}: sequence and germline alignments differ in length")
    end = row.get("v_alignment_end") or ""
    end = int(end) if end.strip() else None
    if end is None:
        raise ValueError(f"{row['sequence_id']}: no v_alignment_end")
    pairs = [(s, g) for s, g in zip(seq[:end], germ[:end]) if s not in ".-" and g not in ".-N"]
    if not pairs:
        raise ValueError(f"{row['sequence_id']}: empty V alignment")
    return sum(s != g for s, g in pairs) / len(pairs)


def first_call(call):
    return (call or "").split(",")[0].strip()


def convert(rows, level):
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
            "true_mu_freq": f"{v_mu(r):.6f}", "level": level}


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--airr", required=True)
    ap.add_argument("--level", required=True)
    ap.add_argument("--reads", required=True)
    ap.add_argument("--truth", required=True)
    a = ap.parse_args(argv)
    try:
        with open(a.airr, newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f, delimiter="\t"))
        if not rows:
            raise ValueError(f"{a.airr}: no rows")
        out = list(convert(rows, a.level))
    except (OSError, ValueError) as e:
        print(f"airr_truth: {e}", file=sys.stderr)
        return 2
    with open(a.reads, "w", encoding="utf-8") as fa, open(a.truth, "w", newline="", encoding="utf-8") as ft:
        w = csv.DictWriter(ft, fieldnames=TRUTH, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for sid, seq, truth in out:
            fa.write(f">{sid}\n{seq}\n")
            w.writerow(truth)
    print(f"{len(out)} reads, {len({t['true_clone_id'] for _, _, t in out})} families, level {a.level}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
