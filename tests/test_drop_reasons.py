"""Tests for bin/drop_reasons.py and the drop counts in the scores (CI run 15 review)."""

import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import drop_reasons as dr  # noqa: E402
import score_clones as sc  # noqa: E402
import summarize_scores as ss  # noqa: E402

HEAD = ("sequence_id", "v_call", "j_call", "junction", "productive")


def tsv(path, header, rows):
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n" + "".join("\t".join(r) + "\n" for r in rows))
    return path


class Reasons(unittest.TestCase):
    def test_each_dropped_read_gets_one_reason(self):
        annotated = {"a": dict(zip(HEAD, ("a", "IGHV1-2*02", "IGHJ4*02", "TGT", "T"))),
                     "b": dict(zip(HEAD, ("b", "IGHV1-2*02", "IGHJ4*02", "TGT", "F"))),
                     "c": dict(zip(HEAD, ("c", "IGHV1-2*02", "", "TGT", "T"))),
                     "e": dict(zip(HEAD, ("e", "IGHV1-2*02", "IGHJ4*02", "TGT", "TRUE")))}
        kept = {"a": {}}
        self.assertEqual(dr.reasons(["a", "b", "c", "d", "e"], annotated, kept),
                         [("b", "not_productive"), ("c", "not_annotated"), ("d", "not_annotated"), ("e", "other")])

    def test_cli_writes_the_table(self):
        with tempfile.TemporaryDirectory() as d:
            reads = os.path.join(d, "r.fa")
            with open(reads, "w") as f:
                f.write(">a\nACGT\n>b desc\nACGT\n>c\nACGT\n")
            allt = tsv(os.path.join(d, "all.tsv"), HEAD, [("a", "V", "J", "TGT", "T"), ("b", "V", "J", "TGT", "F")])
            kept = tsv(os.path.join(d, "db.tsv"), HEAD, [("a", "V", "J", "TGT", "T")])
            out = os.path.join(d, "drops.tsv")
            self.assertEqual(dr.main(["--reads", reads, "--all", allt, "--kept", kept, "--out", out]), 0)
            with open(out) as f:
                self.assertEqual(f.read(), "sequence_id\treason\nb\tnot_productive\nc\tnot_annotated\n")


def truth(n):
    return [{"sequence_id": f"s{i}", "true_clone_id": f"c{i // 2}", "true_v_call": "IGHV1-2*02",
             "true_j_call": "IGHJ4*02", "true_mu_freq": "0.01"} for i in range(n)]


class Scores(unittest.TestCase):
    def test_dropped_reads_are_split_by_reason_and_kept_calls_scored_apart(self):
        t = truth(6)
        pred = [{"sequence_id": f"s{i}", "clone_id": f"k{i // 2}", "v_call": "IGHV1-2*02", "j_call": "IGHJ4*02"}
                for i in (0, 1, 2)]
        res = sc.score(t, pred, drops={"s3": "not_productive", "s4": "not_annotated"})
        self.assertEqual(res["dropped_by_reason"], {"not_productive": 1, "not_annotated": 1, "dropped_later": 1})
        self.assertEqual(res["v_gene"]["exact"], 0.5)          # dropped reads count as wrong calls
        self.assertEqual(res["v_gene_kept"]["exact"], 1.0)     # the kept reads alone
        self.assertEqual(res["all"]["pair_precision"], res["retained"]["pair_precision"])  # drops add no pairs

    def test_summary_columns(self):
        t = truth(4)
        pred = [{"sequence_id": "s0", "clone_id": "k0", "v_call": "IGHV1-2*02", "j_call": "IGHJ4*02"}]
        res = sc.score(t, pred, drops={"s1": "not_productive"})
        res["label"] = "low_e0.0/changeo"
        r = ss.row(res)
        self.assertEqual((r["n_not_annotated"], r["n_not_productive"], r["n_dropped_later"]), (0, 1, 2))
        self.assertEqual((r["v_gene_exact"], r["v_gene_exact_kept"]), (0.25, 1.0))
        self.assertEqual(r["predicted_pairs"], 0)

    def test_old_score_files_still_summarize_with_na(self):
        res = sc.score(truth(2), [])
        res["label"] = "low_e0.0/changeo"
        del res["v_gene_kept"]
        r = ss.row(json.loads(json.dumps(res)))
        self.assertIsNone(r["n_not_productive"])
        self.assertIsNone(r["v_gene_exact_kept"])


if __name__ == "__main__":
    unittest.main()
