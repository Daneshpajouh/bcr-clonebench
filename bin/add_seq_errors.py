#!/usr/bin/env python3
"""Add independent sequencing errors to simulated reads.

    add_seq_errors.py --reads reads.fasta --rate 0.005 --seed 1 --out noisy.fasta

Each base is replaced, with probability --rate, by one of the other three
bases chosen uniformly. The known answers do not change: errors are
noise the pipeline must see through, not hypermutation. Rate 0 copies
the reads unchanged. The benchmark runs every mutation level at each
rate in params.error_rates, so its tables show how each method holds up
as read quality falls.

Standard library only.
"""

import argparse
import random
import sys

BASES = "ACGT"


def read_fasta(path):
    name, seq = None, []
    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(seq)
                name, seq = line[1:], []
            elif line:
                seq.append(line)
    if name is not None:
        yield name, "".join(seq)


def noisy(seq, rate, rng):
    return "".join(rng.choice([x for x in BASES if x != b]) if b in BASES and rng.random() < rate else b
                   for b in seq.upper())


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--reads", required=True)
    ap.add_argument("--rate", type=float, required=True)
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)
    if not 0 <= a.rate < 0.75:
        print("add_seq_errors: --rate must be at least 0 and below 0.75", file=sys.stderr)
        return 2
    rng = random.Random(a.seed)
    n = changed = 0
    with open(a.out, "w", encoding="utf-8") as out:
        for name, seq in read_fasta(a.reads):
            new = noisy(seq, a.rate, rng)
            changed += sum(x != y for x, y in zip(new, seq.upper()))
            n += len(seq)
            out.write(f">{name}\n{new}\n")
    print(f"{changed} of {n} bases changed (rate {a.rate})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
