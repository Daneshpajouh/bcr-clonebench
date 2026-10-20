"""Known-answer tests for bin/score_clones.py. Expected values are worked
by hand in the comments, not computed with the code under test."""

import json
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import score_clones as sc  # noqa: E402


class PairScores(unittest.TestCase):
    def test_hubert_arabie_example(self):
        # true {0,1}{2}{3}, predicted {0,1}{2,3}: shared pairs 1, true pairs 1,
        # predicted pairs 2, all pairs 6. Expected index 1*2/6 = 1/3, maximum
        # (1+2)/2 = 3/2, so ARI = (1 - 1/3) / (3/2 - 1/3) = 4/7.
        s = sc.pair_scores([0, 0, 1, 2], [0, 0, 1, 1])
        self.assertAlmostEqual(s["ari"], 4 / 7)
        self.assertEqual((s["pair_precision"], s["pair_recall"]), (0.5, 1.0))
        self.assertAlmostEqual(s["pair_f1"], 2 / 3)

    def test_relabelled_partition_is_perfect(self):
        s = sc.pair_scores(["a", "a", "b", "b", "b"], [7, 7, 3, 3, 3])
        self.assertEqual((s["ari"], s["pair_f1"]), (1.0, 1.0))

    def test_all_singletons_both_sides(self):
        s = sc.pair_scores([1, 2, 3], ["x", "y", "z"])
        self.assertEqual(s["ari"], 1.0)
        self.assertIsNone(s["pair_precision"])  # no predicted pairs: undefined, not 0 or 1

    def test_one_clone_against_singletons(self):
        # true one family of 4, predicted 4 singletons: recall 0, ARI 0
        s = sc.pair_scores([1, 1, 1, 1], [1, 2, 3, 4])
        self.assertEqual((s["pair_recall"], s["ari"]), (0.0, 0.0))


class Calls(unittest.TestCase):
    def test_gene_level(self):
        r = sc.call_scores(["IGHV1-2*02"] * 4,
                           ["IGHV1-2*04", "IGHV1-2*02,IGHV1-69*01", "IGHV3-23*01", ""])
        # exact: only the first (gene matches, one gene named) -> 1/4
        # includes: first and second -> 2/4
        self.assertEqual(r, {"exact": 0.25, "includes_true": 0.5})


class Hill(unittest.TestCase):
    def test_hand_values(self):
        # sizes 2,1,1: p = 1/2, 1/4, 1/4
        # q1: exp(0.5 ln 2 + 0.5 ln 4) = exp(1.5 ln 2) = 2 sqrt 2
        # q2: 1 / (1/4 + 1/16 + 1/16) = 8/3
        self.assertEqual(sc.hill([2, 1, 1], 0), 3.0)
        self.assertAlmostEqual(sc.hill([2, 1, 1], 1), 2 * math.sqrt(2))
        self.assertAlmostEqual(sc.hill([2, 1, 1], 2), 8 / 3)


def truth_row(i, clone, v="IGHV1-2*02", j="IGHJ4*02", mu="0.04"):
    return {"sequence_id": f"s{i}", "true_clone_id": clone, "true_v_call": v, "true_j_call": j, "true_mu_freq": mu}


def pred_row(i, clone, v="IGHV1-2*02", j="IGHJ4*02", mu="0.04"):
    return {"sequence_id": f"s{i}", "clone_id": clone, "v_call": v, "j_call": j, "mu_freq": mu}


class Score(unittest.TestCase):
    def test_dropped_sequence_is_a_singleton(self):
        truth = [truth_row(1, "A"), truth_row(2, "A"), truth_row(3, "A")]
        pred = [pred_row(1, "c1"), pred_row(2, "c1")]  # s3 dropped by the pipeline
        r = sc.score(truth, pred)
        self.assertEqual((r["n_scored"], r["n_dropped"]), (2, 1))
        # all: true pairs 3, predicted pairs 1 (s3 alone), shared 1
        self.assertEqual((r["all"]["pair_precision"], r["all"]["pair_recall"]), (1.0, 1 / 3))
        self.assertEqual(r["retained"]["pair_recall"], 1.0)
        self.assertEqual(r["v_gene"]["exact"], 2 / 3)  # the dropped read counts as a wrong call

    def test_empty_clone_id_counts_as_dropped(self):
        r = sc.score([truth_row(1, "A"), truth_row(2, "A")], [pred_row(1, "c"), pred_row(2, "")])
        self.assertEqual(r["n_dropped"], 1)

    def test_mutation_error(self):
        truth = [truth_row(1, "A", mu="0.04"), truth_row(2, "B", mu="0.04")]
        pred = [pred_row(1, "a", mu="0.05"), pred_row(2, "b", mu="0.02")]
        m = sc.score(truth, pred, "mu_freq")["mu_freq"]
        self.assertEqual(m["n"], 2)
        self.assertAlmostEqual(m["mean_abs_error"], 0.015)  # (0.01 + 0.02) / 2
        self.assertAlmostEqual(m["mean_error"], -0.005)  # (0.01 - 0.02) / 2

    def test_unknown_prediction_refused(self):
        with self.assertRaisesRegex(ValueError, "not in the truth"):
            sc.score([truth_row(1, "A")], [pred_row(1, "a"), pred_row(9, "a")])


class Cli(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.TemporaryDirectory()
        self.addCleanup(self.d.cleanup)

    def tsv(self, name, rows):
        p = os.path.join(self.d.name, name)
        with open(p, "w") as f:
            f.write("\t".join(rows[0]) + "\n")
            for r in rows[1:]:
                f.write("\t".join(r) + "\n")
        return p

    def test_run_and_errors(self):
        truth = self.tsv("t.tsv", [sc.TRUTH_COLUMNS, ("s1", "A", "IGHV1-2*02", "IGHJ4*02", "0.1"),
                                   ("s2", "A", "IGHV1-2*02", "IGHJ4*02", "0.1")])
        pred = self.tsv("p.tsv", [("sequence_id", "clone_id", "v_call", "j_call", "mu_freq"),
                                  ("s1", "1", "IGHV1-2*02", "IGHJ4*02", "0.1"),
                                  ("s2", "1", "IGHV1-2*02", "IGHJ4*02", "0.1")])
        out = os.path.join(self.d.name, "o.json")
        self.assertEqual(sc.main(["--truth", truth, "--pred", pred, "--mu-column", "mu_freq",
                                  "--label", "low", "--out-json", out]), 0)
        with open(out) as f:
            r = json.load(f)
        self.assertEqual((r["label"], r["all"]["ari"], r["mu_freq"]["mean_abs_error"]), ("low", 1.0, 0.0))
        bad = self.tsv("b.tsv", [("sequence_id", "clone_id"), ("s1", "1")])
        self.assertEqual(sc.main(["--truth", truth, "--pred", bad, "--out-json", out]), 2)
        dup = self.tsv("d.tsv", [sc.TRUTH_COLUMNS, ("s1", "A", "v", "j", "0"), ("s1", "A", "v", "j", "0")])
        self.assertEqual(sc.main(["--truth", dup, "--pred", pred, "--out-json", out]), 2)

    def test_drops_header_without_rows_is_zero_drops(self):
        # CI run 16: a drops table with a header and no rows (no read dropped) was refused as 'no rows'.
        truth = self.tsv("t.tsv", [sc.TRUTH_COLUMNS, ("s1", "A", "IGHV1-2*02", "IGHJ4*02", "0.1"),
                                   ("s2", "A", "IGHV1-2*02", "IGHJ4*02", "0.1")])
        pred = self.tsv("p.tsv", [("sequence_id", "clone_id", "v_call", "j_call"),
                                  ("s1", "1", "IGHV1-2*02", "IGHJ4*02"), ("s2", "1", "IGHV1-2*02", "IGHJ4*02")])
        drops = self.tsv("drops.tsv", [("sequence_id", "reason")])
        out = os.path.join(self.d.name, "o.json")
        self.assertEqual(sc.main(["--truth", truth, "--pred", pred, "--drops", drops, "--out-json", out]), 0)
        with open(out) as f:
            r = json.load(f)
        self.assertEqual((r["n_dropped"], r["dropped_by_reason"]), (0, {}))

    def test_drops_header_without_rows_still_names_reads_dropped_later(self):
        truth = self.tsv("t.tsv", [sc.TRUTH_COLUMNS, ("s1", "A", "v", "j", "0"), ("s2", "A", "v", "j", "0")])
        pred = self.tsv("p.tsv", [("sequence_id", "clone_id", "v_call", "j_call"), ("s1", "1", "v", "j")])
        drops = self.tsv("drops.tsv", [("sequence_id", "reason")])
        out = os.path.join(self.d.name, "o.json")
        self.assertEqual(sc.main(["--truth", truth, "--pred", pred, "--drops", drops, "--out-json", out]), 0)
        with open(out) as f:
            self.assertEqual(json.load(f)["dropped_by_reason"], {"dropped_later": 1})

    def test_drops_file_without_header_is_refused(self):
        truth = self.tsv("t.tsv", [sc.TRUTH_COLUMNS, ("s1", "A", "v", "j", "0")])
        pred = self.tsv("p.tsv", [("sequence_id", "clone_id", "v_call", "j_call"), ("s1", "1", "v", "j")])
        empty = os.path.join(self.d.name, "empty.tsv")
        open(empty, "w").close()
        out = os.path.join(self.d.name, "o.json")
        self.assertEqual(sc.main(["--truth", truth, "--pred", pred, "--drops", empty, "--out-json", out]), 2)
        with self.assertRaisesRegex(ValueError, "no header"):
            sc.read_tsv(empty, ("sequence_id", "reason"), allow_empty=True)

    def test_truth_and_prediction_with_a_header_and_no_rows_are_still_refused(self):
        head_only = self.tsv("h.tsv", [sc.TRUTH_COLUMNS])
        with self.assertRaisesRegex(ValueError, "no rows"):
            sc.read_tsv(head_only, sc.TRUTH_COLUMNS)


if __name__ == "__main__":
    unittest.main()
