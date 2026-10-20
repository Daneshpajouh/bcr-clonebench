"""Tests for bin/airr_truth.py with hand-made AIRR rows."""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import airr_truth as at  # noqa: E402


def row(sid, clone, seq_aln, germ_aln, v_end, v="IGHV1-2*02,IGHV1-2*04"):
    return {"sequence_id": sid, "sequence": seq_aln.replace(".", ""), "clone_id": clone, "v_call": v,
            "d_call": "IGHD3-10*01", "j_call": "IGHJ4*02", "sequence_alignment": seq_aln,
            "germline_alignment": germ_aln, "v_alignment_end": str(v_end)}


class Convert(unittest.TestCase):
    def test_v_region_mutation(self):
        # V is the first 10 aligned positions; 2 differ there (positions 2 and 9),
        # one gap is skipped, and the difference after v_end does not count:
        # 2 of 9 compared positions
        r = row("s1", "c1", "ACTTACG.TAGGG", "AGTTACG.TTGCG", 10)
        self.assertAlmostEqual(at.v_mu(r), 2 / 9)

    def test_first_call_and_truth(self):
        out = list(at.convert([row("s1", "c1", "ACGTACGTAC", "ACGTACGTAC", 10)], "low"))
        sid, seq, t = out[0]
        self.assertEqual((t["true_clone_id"], t["true_v_call"], t["true_mu_freq"], t["level"]),
                         ("c1", "IGHV1-2*02", "0.000000", "low"))

    def test_refusals(self):
        bad = [row("s1", "c1", "ACGT", "ACG", 3), {**row("s2", "c1", "ACGT", "ACGT", 4), "v_alignment_end": ""},
               {**row("s3", "", "ACGT", "ACGT", 4)}]
        for r in bad:
            with self.subTest(r["sequence_id"]):
                with self.assertRaises(ValueError):
                    list(at.convert([r], "low"))
        with self.assertRaisesRegex(ValueError, "repeated"):
            list(at.convert([row("s1", "c", "ACGT", "ACGT", 4)] * 2, "low"))

    def test_cli(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "simu.tsv")
            rows = [row("s1", "c1", "ACGTACGTAC", "ACGTACGTAC", 10), row("s2", "c1", "ACCTACGTAC", "ACGTACGTAC", 10)]
            with open(p, "w") as f:
                f.write("\t".join(rows[0]) + "\n" + "".join("\t".join(r.values()) + "\n" for r in rows))
            reads, truth = os.path.join(d, "r.fa"), os.path.join(d, "t.tsv")
            self.assertEqual(at.main(["--airr", p, "--level", "high", "--reads", reads, "--truth", truth]), 0)
            with open(truth) as f:
                lines = f.read().splitlines()
            self.assertEqual(lines[2].split("\t")[5], "0.100000")
            self.assertEqual(at.main(["--airr", reads, "--level", "high", "--reads", reads + "x", "--truth", truth + "x"]), 2)


if __name__ == "__main__":
    unittest.main()
