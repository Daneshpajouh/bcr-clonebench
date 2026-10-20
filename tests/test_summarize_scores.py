"""Tests for bin/summarize_scores.py."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import summarize_scores as ss  # noqa: E402


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


if __name__ == "__main__":
    unittest.main()
