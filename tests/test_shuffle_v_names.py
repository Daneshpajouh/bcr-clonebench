"""Tests for bin/shuffle_v_names.py, the wrong germline set of the mutation test."""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import shuffle_v_names as sv
from igblast_fasta import allele_name, read_fasta

RECORDS = [("M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|", "CAGGTT"),
           ("X62106|IGHV1-18*02|Homo sapiens|F|V-REGION|", "CAGGTC"),
           ("M99660|IGHV3-23*01|Homo sapiens|F|V-REGION|", "GAGGTG"),
           ("M99649|IGHV3-23D*01|Homo sapiens|F|V-REGION|", "GAGGTG"),
           ("L06612|IGHV5-51*01|Homo sapiens|F|V-REGION|", "GAGGTA")]


class Shuffle(unittest.TestCase):
    def test_every_gene_moves_and_names_stay_unique(self):
        out = sv.shuffle(RECORDS, 7)
        before = [allele_name(h) for h, _ in RECORDS]
        after = [allele_name(h) for h, _ in out]
        for b, a in zip(before, after):
            self.assertNotEqual(b.split("*")[0], a.split("*")[0])
            self.assertEqual(b.split("*")[1], a.split("*")[1])   # allele suffix kept
        self.assertEqual(len(set(after)), len(after))
        self.assertEqual([s for _, s in out], [s for _, s in RECORDS])   # sequences untouched
        self.assertEqual([h.split("|")[0] for h, _ in out], [h.split("|")[0] for h, _ in RECORDS])

    def test_seeded(self):
        self.assertEqual(sv.shuffle(RECORDS, 7), sv.shuffle(RECORDS, 7))

    def test_refusals(self):
        with self.assertRaises(ValueError):
            sv.shuffle(RECORDS[:2], 1)   # one gene only: nothing to move it to
        with self.assertRaises(ValueError):
            sv.shuffle([("x|IGHV1-2|y", "A"), ("z|IGHV1-3*01|w", "C")], 1)   # no allele suffix

    def test_file_round_trip(self):
        with tempfile.TemporaryDirectory() as d:
            src, dst = os.path.join(d, "in.fasta"), os.path.join(d, "out.fasta")
            with open(src, "w") as f:
                f.writelines(f">{h}\n{s}\n" for h, s in RECORDS)
            self.assertEqual(sv.main([src, dst, "3"]), 0)
            self.assertEqual(read_fasta(dst), sv.shuffle(RECORDS, 3))


if __name__ == "__main__":
    unittest.main()
