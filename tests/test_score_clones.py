"""Known-answer tests for bin/score_clones.py. Expected values are worked
by hand in the comments, not computed with the code under test."""

import json
import math
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import score_clones as sc
from score_clones import PRED_COLUMNS


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
            f.writelines("\t".join(r) + "\n" for r in rows[1:])
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


class ProductiveView(unittest.TestCase):
    # Director ruling 2026-10-20 (truth view A): the all-reads score stays the headline; the same scores are added
    # over the reads the simulator marks productive. Values worked by hand.
    @staticmethod
    def truth(*rows):
        cols = (*sc.TRUTH_COLUMNS, "true_productive")
        return [dict(zip(cols, r)) for r in rows]

    @staticmethod
    def pred(*rows):
        return [dict(zip(PRED_COLUMNS, r)) for r in rows]

    def test_scores_over_productive_reads_only(self):
        # s1, s2, s3 productive in family A; s4 not productive, wrongly put with them. Over every read: 3 true
        # pairs, 6 predicted, 3 right, so precision 0.5. Over the productive reads: 3 of 3 pairs right.
        truth = self.truth(("s1", "A", "v", "j", "0", "T"), ("s2", "A", "v", "j", "0", "T"),
                           ("s3", "A", "v", "j", "0", "T"), ("s4", "B", "v", "j", "0", "F"))
        pred = self.pred(("s1", "1", "v", "j"), ("s2", "1", "v", "j"), ("s3", "1", "v", "j"), ("s4", "1", "v", "j"))
        view = sc.productive_view(truth, pred)
        self.assertEqual((view["n_truth"], view["n_not_productive"], view["n_unknown"]), (3, 1, 0))
        self.assertEqual((view["all"]["pair_precision"], view["all"]["pair_f1"]), (1.0, 1.0))
        self.assertEqual(sc.score(truth, pred)["all"]["pair_precision"], 0.5)

    def test_a_dropped_productive_read_still_counts_as_its_own_clone(self):
        # s3 is productive but dropped: 3 true pairs, 1 predicted (s1-s2), so recall 1/3 in the view too.
        truth = self.truth(("s1", "A", "v", "j", "0", "T"), ("s2", "A", "v", "j", "0", "T"),
                           ("s3", "A", "v", "j", "0", "T"))
        view = sc.productive_view(truth, self.pred(("s1", "1", "v", "j"), ("s2", "1", "v", "j")))
        self.assertEqual((view["n_scored"], view["n_dropped"]), (2, 1))
        self.assertAlmostEqual(view["all"]["pair_recall"], 1 / 3)

    def test_no_flag_gives_no_view_and_no_productive_read_gives_an_empty_one(self):
        plain = [dict(zip(sc.TRUTH_COLUMNS, ("s1", "A", "v", "j", "0")))]
        self.assertIsNone(sc.productive_view(plain, self.pred(("s1", "1", "v", "j"))))
        unknown = self.truth(("s1", "A", "v", "j", "0", ""), ("s2", "A", "v", "j", "0", "F"))
        view = sc.productive_view(unknown, self.pred(("s1", "1", "v", "j")))
        self.assertEqual((view["n_truth"], view["n_unknown"], view["all"]), (0, 1, None))

    def test_the_view_is_written_by_the_cli(self):
        with tempfile.TemporaryDirectory() as d:
            t, p, o = (os.path.join(d, n) for n in ("t.tsv", "p.tsv", "o.json"))
            with open(t, "w") as f:
                f.write("\t".join((*sc.TRUTH_COLUMNS, "true_productive")) + "\n")
                f.write("s1\tA\tv\tj\t0\tT\ns2\tA\tv\tj\t0\tF\n")
            with open(p, "w") as f:
                f.write("sequence_id\tclone_id\tv_call\tj_call\ns1\t1\tv\tj\ns2\t1\tv\tj\n")
            self.assertEqual(sc.main(["--truth", t, "--pred", p, "--out-json", o]), 0)
            with open(o) as f:
                r = json.load(f)
            self.assertEqual((r["n_truth"], r["productive_truth"]["n_truth"]), (2, 1))


if __name__ == "__main__":
    unittest.main()


class PerRead(unittest.TestCase):
    # CI run 17 review, item (4): one row per truth read, with its true and called V and J, whether the
    # call set includes the true gene, and whether the read was kept or why it was dropped. Values by hand.
    def setUp(self):
        cols = (*sc.TRUTH_COLUMNS, "true_productive")
        self.truth = [dict(zip(cols, r)) for r in (
            ("s1", "A", "IGHV1-2*02", "IGHJ4*02", "0.04", "T"),
            ("s2", "A", "IGHV1-2*02", "IGHJ4*02", "0.04", "T"),
            ("s3", "B", "IGHV3-23*01", "IGHJ6*02", "0.04", "F"),
            ("s4", "B", "IGHV3-23*01", "IGHJ6*02", "0.04", "T"))]
        # s1 exact; s2 tie that includes the true V; s3 dropped at annotation; s4 dropped later (not in drops)
        self.pred = [pred_row(1, "c1"), pred_row(2, "c1", v="IGHV1-2*04,IGHV1-69*01")]
        self.drops = {"s3": "not_productive"}

    def test_one_row_per_truth_read_in_truth_order(self):
        rows = sc.per_read(self.truth, self.pred, self.drops, "low", "changeo")
        self.assertEqual([r["sequence_id"] for r in rows], ["s1", "s2", "s3", "s4"])
        self.assertEqual({(r["level"], r["method"]) for r in rows}, {("low", "changeo")})
        s1, s2, s3, s4 = rows
        self.assertEqual((s1["status"], s1["drop_reason"], s1["v_exact"], s1["v_includes_true"]), ("kept", "", 1, 1))
        self.assertEqual((s2["v_call"], s2["v_exact"], s2["v_includes_true"], s2["j_exact"]),
                         ("IGHV1-2*04,IGHV1-69*01", 0, 1, 1))
        self.assertEqual((s3["status"], s3["drop_reason"], s3["true_productive"]), ("dropped", "not_productive", "F"))
        self.assertEqual((s3["v_call"], s3["v_exact"], s3["v_includes_true"], s3["clone_id"]), ("", 0, 0, ""))
        self.assertEqual((s4["status"], s4["drop_reason"]), ("dropped", "dropped_later"))

    def test_rows_add_up_to_the_summary_scores(self):
        rows = sc.per_read(self.truth, self.pred, self.drops, "low", "changeo")
        r = sc.score(self.truth, self.pred, drops=self.drops)
        n = len(rows)
        for g in ("v", "j"):
            self.assertEqual(sum(x[f"{g}_exact"] for x in rows) / n, r[f"{g}_gene"]["exact"])
            self.assertEqual(sum(x[f"{g}_includes_true"] for x in rows) / n, r[f"{g}_gene"]["includes_true"])
        reasons = {}
        for x in rows:
            if x["status"] == "dropped":
                reasons[x["drop_reason"]] = reasons.get(x["drop_reason"], 0) + 1
        self.assertEqual(reasons, r["dropped_by_reason"])

    def test_without_a_drops_table_a_dropped_read_says_dropped(self):
        rows = sc.per_read(self.truth, self.pred, None, "low", "changeo")
        self.assertEqual([x["drop_reason"] for x in rows], ["", "", "dropped", "dropped"])

    def test_cli_writes_the_table_and_splits_the_label(self):
        d = tempfile.TemporaryDirectory()
        self.addCleanup(d.cleanup)

        def tsv(name, header, rows):
            p = os.path.join(d.name, name)
            with open(p, "w") as f:
                f.write("\t".join(header) + "\n")
                f.writelines("\t".join(x[c] for c in header) + "\n" for x in rows)
            return p
        truth = tsv("t.tsv", (*sc.TRUTH_COLUMNS, "true_productive"), self.truth)
        pred = tsv("p.tsv", PRED_COLUMNS, self.pred)
        drops = tsv("d.tsv", ("sequence_id", "reason"), [{"sequence_id": "s3", "reason": "not_productive"}])
        out, reads = os.path.join(d.name, "o.json"), os.path.join(d.name, "reads.tsv")
        self.assertEqual(sc.main(["--truth", truth, "--pred", pred, "--drops", drops, "--label", "high/spectral",
                                  "--out-json", out, "--out-reads", reads]), 0)
        with open(reads) as f:
            lines = f.read().splitlines()
        self.assertEqual(lines[0].split("\t"), list(sc.READ_COLUMNS))
        self.assertEqual(len(lines), 5)
        first = dict(zip(sc.READ_COLUMNS, lines[1].split("\t")))
        self.assertEqual((first["level"], first["method"], first["v_exact"]), ("high", "spectral", "1"))
