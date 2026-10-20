#!/usr/bin/env python3
"""Why each read was dropped before clonal grouping.

    drop_reasons.py --reads reads.fasta --all all.tsv --kept db.tsv --out drops.tsv

reads.fasta is what the annotation step read; all.tsv is MakeDb's AIRR
table of every read IgBLAST annotated; db.tsv is the table kept after
the productive filter. Every read missing from db.tsv gets one row:

    not_annotated   no row in all.tsv, or a row without a V call, a J call
                    or a junction (IgBLAST or MakeDb could not place it)
    not_productive  annotated, but productive is not true (a stop codon or
                    a junction out of frame: often a result of mutation)
    other           annotated and productive, yet not kept (not expected)

so the score step can say where the dropped reads went (Director, CI run
15: 44 of 195 dropped at low and 121 to 129 at high mutation).

Standard library only.
"""

import argparse
import csv
import sys

TRUE = {"T", "TRUE", "1"}


def fasta_ids(path):
    with open(path, encoding="utf-8") as f:
        return [line[1:].split()[0] for line in f if line.startswith(">")]


def rows(path):
    with open(path, newline="", encoding="utf-8") as f:
        return {r["sequence_id"]: r for r in csv.DictReader(f, delimiter="\t")}


def reasons(read_ids, annotated, kept):
    out = []
    for sid in read_ids:
        if sid in kept:
            continue
        r = annotated.get(sid)
        if r is None or not all((r.get(c) or "").strip() for c in ("v_call", "j_call", "junction")):
            out.append((sid, "not_annotated"))
        elif (r.get("productive") or "").strip().upper() not in TRUE:
            out.append((sid, "not_productive"))
        else:
            out.append((sid, "other"))  # annotated and productive, yet not kept: named apart, never folded in
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--reads", required=True)
    ap.add_argument("--all", required=True)
    ap.add_argument("--kept", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    try:
        out = reasons(fasta_ids(a.reads), rows(a.all), rows(a.kept))
    except (OSError, KeyError) as e:
        print(f"drop_reasons: {e}", file=sys.stderr)
        return 2
    with open(a.out, "w", newline="", encoding="utf-8") as f:
        w = csv.writer(f, delimiter="\t", lineterminator="\n")
        w.writerow(("sequence_id", "reason"))
        w.writerows(out)
    counts = {}
    for _, reason in out:
        counts[reason] = counts.get(reason, 0) + 1
    print("drop_reasons: " + (", ".join(f"{k} {v}" for k, v in sorted(counts.items())) or "none dropped"))
    return 0


if __name__ == "__main__":
    sys.exit(main())
