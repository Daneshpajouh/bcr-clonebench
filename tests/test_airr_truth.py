"""Tests for bin/airr_truth.py with hand-made AIRR rows."""

import contextlib
import io
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



# Columns partis writes with --airr-output (partis/utils.py airr_headers, main 92c6aa00eba5): no
# v_alignment_end; the V slice comes as v_sequence_alignment and v_germline_alignment.
PARTIS_HEAD = ("sequence_id", "sequence", "rev_comp", "productive", "v_call", "d_call", "j_call",
               "sequence_alignment", "germline_alignment", "junction", "junction_aa", "v_cigar", "d_cigar",
               "j_cigar", "clone_id", "v_sequence_alignment", "v_germline_alignment", "v_sequence_start",
               "v_sequence_end")


def partis_row(sid, clone, seq_aln, germ_aln, v_start, v_end):
    """A partis-shaped row: full-length alignments, V given as the slice [v_start, v_end) of both."""
    r = dict.fromkeys(PARTIS_HEAD, "")
    r.update(sequence_id=sid, sequence=seq_aln.replace(".", "").replace("-", ""), rev_comp="False",
             productive="True", v_call="IGHV1-2*02", d_call="IGHD3-10*01", j_call="IGHJ4*02",
             sequence_alignment=seq_aln, germline_alignment=germ_aln, clone_id=clone,
             v_sequence_alignment=seq_aln[v_start:v_end], v_germline_alignment=germ_aln[v_start:v_end],
             v_sequence_start=str(v_start + 1), v_sequence_end=str(v_end))
    return r


class PartisShapedRows(unittest.TestCase):
    """CI run 14: partis's table has no v_alignment_end, so the V slice columns are read instead."""

    def test_a_row_without_v_alignment_end_is_scored_on_the_v_slice(self):
        # Two leading bases before V, V is 8 columns with one gap and 2 differences, a difference after V
        # does not count: 2 of 7 compared positions.
        r = partis_row("70d2a8c4c3", "c1", "GGACTTA.CGTAGGG", "GGAGTTA.CTTTGCG", 2, 10)
        self.assertNotIn("v_alignment_end", r)
        self.assertAlmostEqual(at.v_mu(r), 2 / 7)
        routes = {}
        out = list(at.convert([r], "low", routes))
        self.assertEqual(out[0][2]["true_mu_freq"], f"{2 / 7:.6f}")
        self.assertEqual(routes, {"v_sequence_alignment": 1})

    def test_v_alignment_end_still_wins_when_present(self):
        r = {**partis_row("s1", "c1", "ACGTACGTAC", "ACGTACGTAA", 0, 10), "v_alignment_end": "9"}
        self.assertEqual(at.v_mu(r), 0.0)

    def test_rows_with_no_usable_v_stretch_are_refused(self):
        good = partis_row("s1", "c1", "GGACTTACGTAGGG", "GGAGTTACTTTGCG", 2, 10)
        bad = {"neither": {**good, "v_sequence_alignment": ""}, "end": {**good, "v_alignment_end": "NA"},
               "lengths": {**good, "v_germline_alignment": good["v_germline_alignment"][:-1]},
               "not a slice": {**good, "v_germline_alignment": "A" * len(good["v_germline_alignment"])}}
        for name, r in bad.items():
            with self.subTest(name):
                with self.assertRaises(ValueError):
                    list(at.convert([r], "low"))

    def test_a_failed_read_prints_the_header_and_first_two_rows(self):
        rows = [partis_row(f"s{i}", "c1", "GGACTTACGTAGGG", "GGAGTTACTTTGCG", 2, 10) for i in range(3)]
        rows[0]["v_germline_alignment"] = ""
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "simu_low.tsv")
            with open(p, "w") as f:
                f.write("\t".join(PARTIS_HEAD) + "\n" + "".join("\t".join(r.values()) + "\n" for r in rows))
            err = io.StringIO()
            with contextlib.redirect_stderr(err):
                rc = at.main(["--airr", p, "--level", "low", "--reads", os.path.join(d, "r.fa"),
                              "--truth", os.path.join(d, "t.tsv")])
        self.assertEqual(rc, 2)
        lines = err.getvalue().splitlines()
        self.assertIn("s0: no v_alignment_end and no v_sequence_alignment", lines[0])
        self.assertEqual(lines[2], "\t".join(PARTIS_HEAD))
        self.assertTrue(lines[3].startswith("s0\t") and lines[4].startswith("s1\t"))
        self.assertEqual(len(lines), 5)

if __name__ == "__main__":
    unittest.main()
