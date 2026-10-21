"""Tests for bin/summarize_scores.py."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import summarize_scores as ss


def score(label, precision):
    return {"label": label, "n_truth": 10, "n_dropped": 1,
            "all": {"pair_precision": precision, "pair_recall": 0.5, "pair_f1": None, "ari": 0.25},
            "v_gene": {"exact": 0.9}, "j_gene": {"exact": 1.0}, "mu_freq": {"mean_abs_error": 0.01},
            "hill": {"q0": {"true": 4.0, "predicted": 5.0}}}


class Summary(unittest.TestCase):
    def test_order_and_na(self):
        with tempfile.TemporaryDirectory() as d:
            paths = []
            for i, (label, p) in enumerate((("high_e0.0/changeo", 0.8), ("low_e0.005/changeo", None),
                                            ("low_e0.0/changeo", 1.0))):
                paths.append(os.path.join(d, f"s{i}.json"))
                with open(paths[-1], "w") as f:
                    json.dump(score(label, p), f)
            tsv, md = os.path.join(d, "r.tsv"), os.path.join(d, "r.md")
            self.assertEqual(ss.main(["--out-tsv", tsv, "--out-md", md] + paths), 0)
            with open(tsv) as f:
                rows = [line.rstrip("\n").split("\t") for line in f]
            self.assertEqual(rows[0], list(ss.COLUMNS))
            self.assertEqual([r[:3] for r in rows[1:]], [["low", "0.0000", "changeo"], ["low", "0.0050", "changeo"],
                                                         ["high", "0.0000", "changeo"]])
            self.assertEqual(rows[2][5], "NA")  # undefined precision is NA, never 0
            self.assertEqual(rows[1][7], "NA")  # undefined F1
            with open(md) as f:
                self.assertEqual(len(f.read().splitlines()), 5)

    def test_productive_view_columns(self):
        with tempfile.TemporaryDirectory() as d:
            with_view = {**score("low_e0.0/changeo", 1.0),
                         "productive_truth": {"n_truth": 8, "all": {"pair_f1": 0.75, "ari": 0.5}}}
            paths = []
            for i, s in enumerate((with_view, score("low_e0.0/hierarchical", 1.0))):
                paths.append(os.path.join(d, f"s{i}.json"))
                with open(paths[-1], "w") as f:
                    json.dump(s, f)
            tsv, md = os.path.join(d, "r.tsv"), os.path.join(d, "r.md")
            self.assertEqual(ss.main(["--out-tsv", tsv, "--out-md", md] + paths), 0)
            with open(tsv) as f:
                rows = [dict(zip(ss.COLUMNS, line.rstrip("\n").split("\t"))) for line in list(f)[1:]]
            got = [(r["method"], r["n_truth_productive"], r["pair_f1_productive"], r["ari_productive"]) for r in rows]
            self.assertEqual(got, [("changeo", "8", "0.7500", "0.5000"), ("hierarchical", "NA", "NA", "NA")])

    def test_includes_true_columns(self):
        # A call that ties the true gene with an identical one (IGHV3-23 with IGHV3-23D) is never exact but does
        # include the true gene; the table shows both, and NA when a score file predates the includes share.
        with tempfile.TemporaryDirectory() as d:
            tied = {**score("low_shared_e0.0/changeo", 1.0),
                    "v_gene": {"exact": 0.0, "includes_true": 0.75}, "j_gene": {"exact": 0.7, "includes_true": 0.8}}
            paths = []
            for i, s in enumerate((tied, score("low_e0.0/changeo", 1.0))):
                paths.append(os.path.join(d, f"s{i}.json"))
                with open(paths[-1], "w") as f:
                    json.dump(s, f)
            tsv, md = os.path.join(d, "r.tsv"), os.path.join(d, "r.md")
            self.assertEqual(ss.main(["--out-tsv", tsv, "--out-md", md] + paths), 0)
            with open(tsv) as f:
                rows = [dict(zip(ss.COLUMNS, line.rstrip("\n").split("\t"))) for line in list(f)[1:]]
            got = [(r["level"], r["v_gene_exact"], r["v_gene_includes_true"], r["j_gene_includes_true"]) for r in rows]
            self.assertEqual(got, [("low", "0.9000", "NA", "NA"), ("low_shared", "0.0000", "0.7500", "0.8000")])


    def test_hill_q1_q2_columns(self):
        # Diversity is scored at q = 0, 1 and 2; q0 is the clone count, q1 and q2 get their own columns,
        # and NA when a score file has no q1 or q2 (or no predicted clones).
        with tempfile.TemporaryDirectory() as d:
            full = {**score("low_e0.0/changeo", 1.0),
                    "hill": {"q0": {"true": 4.0, "predicted": 5.0}, "q1": {"true": 3.5, "predicted": 4.25},
                             "q2": {"true": 3.0, "predicted": None}}}
            paths = []
            for i, s in enumerate((full, score("low_e0.0/hierarchical", 1.0))):
                paths.append(os.path.join(d, f"s{i}.json"))
                with open(paths[-1], "w") as f:
                    json.dump(s, f)
            tsv, md = os.path.join(d, "r.tsv"), os.path.join(d, "r.md")
            self.assertEqual(ss.main(["--out-tsv", tsv, "--out-md", md] + paths), 0)
            with open(tsv) as f:
                rows = [dict(zip(ss.COLUMNS, line.rstrip("\n").split("\t"))) for line in list(f)[1:]]
            cols = ("clones_true", "clones_predicted", "hill_q1_true", "hill_q1_predicted", "hill_q2_true",
                    "hill_q2_predicted")
            got = [(r["method"], *(r.get(c) for c in cols)) for r in rows]
            self.assertEqual(got, [("changeo", "4.0000", "5.0000", "3.5000", "4.2500", "3.0000", "NA"),
                                   ("hierarchical", "4.0000", "5.0000", "NA", "NA", "NA", "NA")])

if __name__ == "__main__":
    unittest.main()
