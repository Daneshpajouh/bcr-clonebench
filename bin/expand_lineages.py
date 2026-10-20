#!/usr/bin/env python3
"""Grow each simulated founder into a clonal family with known answers.

    expand_lineages.py --founders founders.tsv --level medium --seed 1 \
        --reads reads.fasta --truth truth.tsv

immuneSIM writes naive rearrangements (run with shm.mode "none") and
expands clones only as identical copies. A clonal grouping benchmark
needs families whose members differ, so this step builds them: each
founder becomes a family of related sequences on a random lineage tree,
with somatic mutations added along the tree's edges, then independent
sequencing errors on every read.

founders.tsv: founder_id, sequence (in frame from the first base),
v_call, d_call, j_call, and v_end (the length of the V region scored
for hypermutation, ending at the conserved cysteine codon).

Truth per read: the family, the V, D and J alleles, and true_mu_freq,
the share of V-region positions (before v_end) that differ from the
founder because of hypermutation. Sequencing errors are not counted
there: they are noise the pipeline has to see through.

Mutation levels set the expected V-region mutation frequency of family
members: low 0.01, medium 0.05, high 0.10. Mutations that would make a
stop codon are drawn again (selection keeps receptors productive), and
WRC and GYW hot spots are three times as likely to mutate as other
positions, a simple stand-in for AID targeting.

Standard library only.
"""

import argparse
import csv
import math
import random
import sys

LEVELS = {"low": 0.01, "medium": 0.05, "high": 0.10}
STOPS = {"TAA", "TAG", "TGA"}
BASES = "ACGT"
HOTSPOT_WEIGHT = 3.0
FOUNDER_COLUMNS = ("founder_id", "sequence", "v_call", "d_call", "j_call", "v_end")


def read_founders(path):
    with open(path, newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    if not rows:
        raise ValueError(f"{path}: no founders")
    missing = [c for c in FOUNDER_COLUMNS if c not in rows[0]]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)}")
    for r in rows:
        seq = r["sequence"].upper()
        if set(seq) - set(BASES):
            raise ValueError(f"founder {r['founder_id']}: sequence has characters other than ACGT")
        if has_stop(seq):
            raise ValueError(f"founder {r['founder_id']}: stop codon in frame")
        v_end = int(r["v_end"])
        if not 0 < v_end <= len(seq):
            raise ValueError(f"founder {r['founder_id']}: v_end {v_end} outside the sequence")
        r["sequence"], r["v_end"] = seq, v_end
    return rows


def has_stop(seq):
    return any(seq[i:i + 3] in STOPS for i in range(0, len(seq) - 2, 3))


def weights(seq):
    """Mutation weight per position: WRC (W = A or T, R = A or G) with the C
    at the position, and its reverse complement GYW with the G there."""
    w = [1.0] * len(seq)
    for i, b in enumerate(seq):
        if b == "C" and i >= 2 and seq[i - 2] in "AT" and seq[i - 1] in "AG":
            w[i] = HOTSPOT_WEIGHT
        if b == "G" and i + 2 < len(seq) and seq[i + 1] in "CT" and seq[i + 2] in "AT":
            w[i] = HOTSPOT_WEIGHT
    return w


def mutate(seq, n, rng, tries=50):
    """seq with n point mutations at weighted positions, none making a stop
    codon. A mutation that would is drawn again, up to `tries` times."""
    s = list(seq)
    for _ in range(n):
        for _ in range(tries):
            w = weights("".join(s))
            i = rng.choices(range(len(s)), weights=w)[0]
            new = rng.choice([b for b in BASES if b != s[i]])
            c = i - i % 3
            codon = "".join(s[c:c + 3])
            codon = codon[:i - c] + new + codon[i - c + 1:]
            if codon not in STOPS:
                s[i] = new
                break
    return "".join(s)


def tree(size, rng):
    """Parent index of each member of a random recursive tree; member 0
    is the founder's first descendant and hangs from the founder (-1)."""
    parents = [-1]
    for k in range(1, size):
        parents.append(rng.randrange(-1, k))
    return parents


def depths(parents):
    d = []
    for p in parents:
        d.append(1 if p < 0 else d[p] + 1)
    return d


def family_size(rng, mean):
    """At least 1; geometric, so a few families are large."""
    return 1 + int(math.log(1 - rng.random()) / math.log(1 - 1 / mean)) if mean > 1 else 1


def expand(founders, level, seed, mean_size=4.0, error_rate=0.001):
    """Yield (read_id, sequence, truth dict) for every read."""
    rng = random.Random(seed)
    target = LEVELS[level]
    for f in founders:
        size = family_size(rng, mean_size)
        parents = tree(size, rng)
        dep = depths(parents)
        # per-edge mutations so the members' mean V-region frequency is the target
        lam = target * f["v_end"] / (sum(dep) / len(dep))
        seqs = []
        for k, p in enumerate(parents):
            base = f["sequence"] if p < 0 else seqs[p]
            n = poisson(lam * len(base) / f["v_end"] if f["v_end"] else 0, rng)
            seqs.append(mutate(base, n, rng))
        for k, s in enumerate(seqs):
            mu = sum(a != b for a, b in zip(s[:f["v_end"]], f["sequence"][:f["v_end"]])) / f["v_end"]
            read = "".join(rng.choice([x for x in BASES if x != b]) if rng.random() < error_rate else b for b in s)
            rid = f"{f['founder_id']}_m{k}"
            yield rid, read, {"sequence_id": rid, "true_clone_id": f["founder_id"], "true_v_call": f["v_call"],
                              "true_d_call": f["d_call"], "true_j_call": f["j_call"],
                              "true_mu_freq": f"{mu:.6f}", "level": level}


def poisson(lam, rng):
    """Knuth's method; lam here is at most a few dozen."""
    if lam <= 0:
        return 0
    limit, k, p = math.exp(-lam), 0, 1.0
    while True:
        p *= rng.random()
        if p <= limit:
            return k
        k += 1


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--founders", required=True)
    ap.add_argument("--level", required=True, choices=sorted(LEVELS))
    ap.add_argument("--seed", type=int, required=True)
    ap.add_argument("--mean-family-size", type=float, default=4.0)
    ap.add_argument("--error-rate", type=float, default=0.001, help="sequencing error per base")
    ap.add_argument("--reads", required=True, help="FASTA out")
    ap.add_argument("--truth", required=True, help="TSV out")
    a = ap.parse_args(argv)
    try:
        founders = read_founders(a.founders)
    except (OSError, ValueError) as e:
        print(f"expand_lineages: {e}", file=sys.stderr)
        return 2
    cols = ("sequence_id", "true_clone_id", "true_v_call", "true_d_call", "true_j_call", "true_mu_freq", "level")
    n = 0
    with open(a.reads, "w", encoding="utf-8") as fa, open(a.truth, "w", newline="", encoding="utf-8") as ft:
        w = csv.DictWriter(ft, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        for rid, seq, truth in expand(founders, a.level, a.seed, a.mean_family_size, a.error_rate):
            fa.write(f">{rid}\n{seq}\n")
            w.writerow(truth)
            n += 1
    print(f"{n} reads from {len(founders)} families, level {a.level}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
