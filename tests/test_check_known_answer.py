"""Tests for bin/check_known_answer.py: named checks pass on a good table, fail on a bad one, and the
--expect-fail mode passes only when a deliberately wrong setting is caught."""

import contextlib
import io
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import check_known_answer as ka

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
COLS = ("level", "error_rate", "method", "n_truth", "n_dropped", "pair_precision", "v_gene_includes_true", "mu_mae")


def table(rows):
    return "\t".join(COLS) + "\n" + "".join("\t".join(str(r[c]) for c in COLS) + "\n" for r in rows)


def good():
    return [
        {"level": "low", "error_rate": "0.0000", "method": "changeo", "n_truth": 200, "n_dropped": 50,
         "pair_precision": "1.0000", "v_gene_includes_true": "0.7500", "mu_mae": "0.0010"},
        {"level": "low_shared", "error_rate": "0.0000", "method": "changeo", "n_truth": 200, "n_dropped": 40,
         "pair_precision": "0.9900", "v_gene_includes_true": "0.8000", "mu_mae": "0.0020"},
    ]


CHECKS = {"checks": [
    {"name": "V", "metric": "v_gene_includes_true/kept", "min": 0.98, "where": {}},
    {"name": "P", "metric": "pair_precision", "min": 0.95, "where": {"method": "changeo"}},
    {"name": "P-SHARED", "metric": "pair_precision", "min": 0.995, "where": {"shared": True}},
    {"name": "MU", "metric": "mu_mae", "max": 0.02, "where": {}},
]}


class Check(unittest.TestCase):
    def run_on(self, rows, checks=CHECKS, *extra):
        with tempfile.TemporaryDirectory() as d:
            res, fl = os.path.join(d, "r.tsv"), os.path.join(d, "f.json")
            with open(res, "w") as f:
                f.write(table(rows))
            with open(fl, "w") as f:
                json.dump(checks, f)
            out = io.StringIO()
            with contextlib.redirect_stdout(out):
                rc = ka.main([res, fl, *extra])
            return rc, {line.split("\t")[1]: line.split("\t")[0] for line in out.getvalue().splitlines()
                        if "\t" in line}

    def test_kept_share_is_the_denominator(self):
        # 0.75 of all reads with 50 of 200 dropped is every kept read.
        self.assertAlmostEqual(ka.value(good()[0] | {"n_truth": "200", "n_dropped": "50"}, "v_gene_includes_true/kept"),
                               1.0)

    def test_where_shared_selects_the_shared_levels(self):
        rc, got = self.run_on(good())
        self.assertEqual(got, {"V": "PASS", "P": "PASS", "P-SHARED": "FAIL", "MU": "PASS"})
        self.assertEqual(rc, 1)

    def test_a_wrong_germline_shape_fails_the_v_check(self):
        rows = good()
        rows[1]["v_gene_includes_true"] = "0.1000"   # true V named for few kept reads
        rc, got = self.run_on(rows, {"checks": CHECKS["checks"][:2]})
        self.assertEqual((rc, got["V"], got["P"]), (1, "FAIL", "PASS"))

    def test_na_and_no_selected_row_fail(self):
        rows = good()
        rows[0]["mu_mae"] = "NA"
        rc, got = self.run_on(rows)
        self.assertEqual((rc, got["MU"]), (1, "NO DATA"))
        none = {"checks": [{"name": "NONE", "metric": "pair_precision", "min": 0.5, "where": {"method": "spectral"}}]}
        rc, got = self.run_on(good(), none)
        self.assertEqual((rc, got["NONE"]), (1, "NO DATA"))

    def test_no_data_is_not_a_caught_mutant(self):
        # A wrong setting that breaks the pipeline (no value) must not pass the mutation test.
        rows = good()
        for r in rows:
            r["v_gene_includes_true"] = "NA"
        self.assertEqual(self.run_on(rows, CHECKS, "--expect-fail", "V")[0], 1)
        rows = good()
        rows[1]["v_gene_includes_true"] = "0.1000"
        self.assertEqual(self.run_on(rows, CHECKS, "--expect-fail", "V")[0], 0)

    def test_expect_fail_passes_only_when_caught(self):
        self.assertEqual(self.run_on(good(), CHECKS, "--expect-fail", "P-SHARED")[0], 0)
        self.assertEqual(self.run_on(good(), CHECKS, "--expect-fail", "P")[0], 1)
        self.assertEqual(self.run_on(good(), CHECKS, "--expect-fail", "NO-SUCH-CHECK")[0], 2)

    def test_repository_floors_are_well_formed(self):
        with open(os.path.join(ROOT, "conf", "known_answer_floors.json")) as f:
            checks = json.load(f)["checks"]
        self.assertEqual(len({c["name"] for c in checks}), len(checks))
        for c in checks:
            with self.subTest(c["name"]):
                self.assertEqual(("min" in c) + ("max" in c), 1)
                bound = c.get("min", c.get("max"))
                # every floor holds on the value it was set from, with room to spare
                self.assertTrue(c["seen"] >= bound if "min" in c else c["seen"] <= bound)


    def test_make_targets_run_the_checks(self):
        # known-answer ends with the named checks; each mutation test names checks the floors file defines.
        with open(os.path.join(ROOT, "Makefile")) as f:
            make = f.read()
        with open(os.path.join(ROOT, "conf", "known_answer_floors.json")) as f:
            names = {c["name"] for c in json.load(f)["checks"]}
        ka_recipe = make[make.index("known-answer: reference"):make.index("mutation-tests:")]
        self.assertIn("$(CHECK) results/results.tsv conf/known_answer_floors.json", ka_recipe)
        mut = make[make.index("\nmutation-tests:"):]
        expected = [line.split("--expect-fail ")[1].strip() for line in mut.splitlines() if "--expect-fail" in line]
        self.assertEqual(len(expected), 2)
        for group in expected:
            for name in group.split(","):
                self.assertIn(name, names)
        self.assertIn("--threshold_override", mut)
        self.assertIn("shuffled-v", mut)
        with open(os.path.join(ROOT, "modules", "threshold.nf")) as f:
            self.assertIn("params.threshold_override != null", f.read())

if __name__ == "__main__":
    unittest.main()
