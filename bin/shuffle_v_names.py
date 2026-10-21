#!/usr/bin/env python3
"""Write a deliberately wrong V germline set: the same sequences under other genes' names.

    shuffle_v_names.py IN_IGHV.fasta OUT_IGHV.fasta SEED

Used only by the mutation test (make mutation-tests): the pipeline annotated
against this set must fail the named known-answer check KA-V-INCLUDES.
The gene part of every allele name (IGHV1-18 in IGHV1-18*01) is moved to
another gene by a seeded permutation with no fixed point; the allele suffix
stays. The permutation is one to one, so no two sequences share a name, and
the rest of each header (and every sequence) is unchanged.
"""

import random
import sys

from igblast_fasta import allele_name, read_fasta


def derangement(names, seed):
    """A seeded permutation of the sorted names in which no name maps to itself."""
    names = sorted(set(names))
    if len(names) < 2:
        raise ValueError("need at least two gene names to move")
    rng = random.Random(seed)
    order = names[:]
    rng.shuffle(order)
    # rotate the shuffled order by one: a cycle through every name, so no fixed point
    return {order[i]: order[(i + 1) % len(order)] for i in range(len(order))}


def shuffle(records, seed):
    genes = []
    for header, _ in records:
        name = allele_name(header)
        if "*" not in name:
            raise ValueError(f"no allele suffix in {name!r}")
        genes.append(name.split("*", 1)[0])
    move = derangement(genes, seed)
    out = []
    for (header, seq), gene in zip(records, genes):
        name = allele_name(header)
        new = move[gene] + name[len(gene):]
        if header.count(name) != 1:
            raise ValueError(f"allele name {name!r} is not once in its header")
        out.append((header.replace(name, new), seq))
    return out


def main(argv):
    if len(argv) != 3:
        print(__doc__, file=sys.stderr)
        return 2
    records = shuffle(read_fasta(argv[0]), int(argv[2]))
    with open(argv[1], "w") as fh:
        fh.writelines(f">{h}\n{s}\n" for h, s in records)
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
