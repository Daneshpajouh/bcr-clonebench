"""Tests for bin/expand_lineages.py on synthetic founders (random open
reading frames), so they need no simulator."""

import csv
import os
import random
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import expand_lineages as el

CODONS = [a + b + c for a in "ACGT" for b in "ACGT" for c in "ACGT" if a + b + c not in el.STOPS]


def founders(n, seed=7, length=360, v_end=290):
    rng = random.Random(seed)
    return [{"founder_id": f"F{i}", "sequence": "".join(rng.choice(CODONS) for _ in range(length // 3)),
             "v_call": "IGHV1-2*02", "d_call": "IGHD3-10*01", "j_call": "IGHJ4*02", "v_end": v_end}
            for i in range(n)]


class Expand(unittest.TestCase):
    def run_level(self, level, n=300, error_rate=0.0):
        return list(el.expand(founders(n), level, seed=11, error_rate=error_rate))

    def test_same_seed_same_reads(self):
        a = [(r, s) for r, s, _ in el.expand(founders(20), "medium", seed=3)]
        b = [(r, s) for r, s, _ in el.expand(founders(20), "medium", seed=3)]
        c = [(r, s) for r, s, _ in el.expand(founders(20), "medium", seed=4)]
        self.assertEqual(a, b)
        self.assertNotEqual(a, c)

    def test_mean_mutation_near_each_level(self):
        for level, target in el.LEVELS.items():
            with self.subTest(level):
                mus = [float(t["true_mu_freq"]) for _, _, t in self.run_level(level)]
                mean = sum(mus) / len(mus)
                # back mutations and hot-spot weighting move it a little; the
                # levels must stay apart and near their targets
                self.assertLess(abs(mean - target), 0.25 * target, (level, mean))

    def test_truth_matches_the_sequences(self):
        fs = {f["founder_id"]: f for f in founders(50)}
        for rid, seq, t in el.expand(list(fs.values()), "high", seed=5, error_rate=0.0):
            f = fs[t["true_clone_id"]]
            diff = sum(a != b for a, b in zip(seq[:290], f["sequence"][:290]))
            self.assertAlmostEqual(float(t["true_mu_freq"]), diff / 290, places=6)
            self.assertTrue(rid.startswith(t["true_clone_id"] + "_m"))

    def test_no_stop_codons_from_hypermutation(self):
        for _, seq, _ in self.run_level("high", n=100):
            self.assertFalse(el.has_stop(seq))

    def test_sequencing_errors_are_not_hypermutation(self):
        clean = {r: s for r, s, _ in el.expand(founders(100), "low", seed=9, error_rate=0.0)}
        noisy = list(el.expand(founders(100), "low", seed=9, error_rate=0.01))
        diffs = sum(sum(a != b for a, b in zip(s, clean[r])) for r, s, _ in noisy if r in clean)
        self.assertGreater(diffs, 0)

    def test_family_sizes(self):
        sizes = {}
        for _, _, t in self.run_level("low", n=1000):
            sizes[t["true_clone_id"]] = sizes.get(t["true_clone_id"], 0) + 1
        self.assertEqual(len(sizes), 1000)
        self.assertGreaterEqual(min(sizes.values()), 1)
        self.assertLess(abs(sum(sizes.values()) / 1000 - 4.0), 0.5)

    def test_hotspots(self):
        w = el.weights("AACGGCTA")
        # AAC: WRC with C at 2 -> hot; GCT: G at 4 followed by C,T -> GYW hot
        self.assertEqual(w, [1.0, 1.0, 3.0, 1.0, 3.0, 1.0, 1.0, 1.0])


class Cli(unittest.TestCase):
    def test_files_and_bad_input(self):
        with tempfile.TemporaryDirectory() as d:
            fp = os.path.join(d, "f.tsv")
            with open(fp, "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=el.FOUNDER_COLUMNS, delimiter="\t")
                w.writeheader()
                for r in founders(5):
                    w.writerow(r)
            reads, truth = os.path.join(d, "r.fa"), os.path.join(d, "t.tsv")
            self.assertEqual(el.main(["--founders", fp, "--level", "low", "--seed", "1",
                                      "--reads", reads, "--truth", truth]), 0)
            with open(reads) as f:
                n_reads = sum(1 for line in f if line.startswith(">"))
            with open(truth) as f:
                self.assertEqual(sum(1 for _ in f) - 1, n_reads)
            bad = os.path.join(d, "bad.tsv")
            with open(bad, "w") as f:
                f.write("founder_id\tsequence\tv_call\td_call\tj_call\tv_end\nF\tATGTAAGGG\tv\td\tj\t6\n")
            self.assertEqual(el.main(["--founders", bad, "--level", "low", "--seed", "1",
                                      "--reads", reads, "--truth", truth]), 2)


if __name__ == "__main__":
    unittest.main()
