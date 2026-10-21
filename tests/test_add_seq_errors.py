"""Tests for bin/add_seq_errors.py."""

import os
import random
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import add_seq_errors as ase


class Errors(unittest.TestCase):
    def test_rate_zero_copies(self):
        self.assertEqual(ase.noisy("ACGTN", 0.0, random.Random(1)), "ACGTN")

    def test_rate_near_target_and_never_the_same_base(self):
        seq = "ACGT" * 25000
        new = ase.noisy(seq, 0.01, random.Random(2))
        diffs = [(a, b) for a, b in zip(seq, new) if a != b]
        # 100,000 bases at 1%: expect 1,000; three standard deviations is about 94
        self.assertLess(abs(len(diffs) - 1000), 100)
        self.assertTrue(all(b in "ACGT" for _, b in diffs))

    def test_cli_and_bad_rate(self):
        with tempfile.TemporaryDirectory() as d:
            src, out = os.path.join(d, "r.fa"), os.path.join(d, "o.fa")
            with open(src, "w") as f:
                f.write(">a\nACGT\nACGT\n>b\nTTTT\n")
            self.assertEqual(ase.main(["--reads", src, "--rate", "0", "--seed", "1", "--out", out]), 0)
            self.assertEqual(list(ase.read_fasta(out)), [("a", "ACGTACGT"), ("b", "TTTT")])
            self.assertEqual(ase.main(["--reads", src, "--rate", "0.9", "--seed", "1", "--out", out]), 2)


if __name__ == "__main__":
    unittest.main()
