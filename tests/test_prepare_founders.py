"""Tests for bin/prepare_founders.py on a synthetic germline V whose end
copies the usual IGHV ending (... GTG TAT TAC TGT GCG AGA GA)."""

import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import prepare_founders as pf  # noqa: E402

GERM_V = "CAG" * 30 + "GTGTATTACTGTGCGAGAGA"  # cysteine codon TGT at 99..101, so v_end 102
REST = "GGTTATTGGTTCGACCCCTGGGGCCAGGGAACC"  # an in-frame tail with no stop


def row(seq, v="IGHV1-2*02", j="IGHJ5*02"):
    return {"sequence": seq, "v_call": v, "d_call": "IGHD3-10*01", "j_call": j}


class Prepare(unittest.TestCase):
    def test_cys_end(self):
        self.assertEqual(pf.cys_end(GERM_V), 102)
        self.assertIsNone(pf.cys_end("CAG" * 40))

    def test_keeps_and_drops(self):
        good = GERM_V[:102] + "GC" + REST + "A"  # trailing partial codon is cut
        rows = [row(good), row(good),                                   # repeat
                row(GERM_V[:102] + "TAAGGG"),                           # stop in frame
                row("CTG" + GERM_V[3:102] + REST[:30]),                 # start differs from germline
                row(GERM_V[:102] + "NNN"),                              # not ACGT
                row(good + "", v="IGHV9-99*01")]                        # no germline
        out, dropped = pf.prepare(rows, {"IGHV1-2*02": GERM_V})
        self.assertEqual(len(out), 1)
        self.assertEqual(out[0]["v_end"], 102)
        self.assertEqual(len(out[0]["sequence"]) % 3, 0)
        self.assertEqual(dropped, {"stop": 1, "not ACGT": 1, "no germline": 1,
                                   "start differs from germline": 1, "repeat": 1})

    def test_first_call_only_and_max(self):
        seqs = [GERM_V[:102] + REST[i:i + 6] + "GGG" for i in range(0, 18, 3)]
        out, _ = pf.prepare([row(s, v="IGHV1-2*02,IGHV1-2*04") for s in seqs], {"IGHV1-2*02": GERM_V}, max_n=3)
        self.assertEqual([f["founder_id"] for f in out], ["F00001", "F00002", "F00003"])
        self.assertEqual(out[0]["v_call"], "IGHV1-2*02")

    def test_cli(self):
        with tempfile.TemporaryDirectory() as d:
            rep = os.path.join(d, "r.tsv")
            with open(rep, "w") as f:
                f.write("sequence\tv_call\td_call\tj_call\n" + GERM_V[:102] + REST + "\tIGHV1-2*02\tIGHD1*01\tIGHJ5*02\n")
            fa = os.path.join(d, "v.fasta")
            with open(fa, "w") as f:
                f.write(">IGHV1-2*02\n" + GERM_V[:60] + "\n" + GERM_V[60:] + "\n")
            out = os.path.join(d, "f.tsv")
            self.assertEqual(pf.main(["--immunesim", rep, "--v-germline", fa, "--out", out]), 0)
            with open(out) as f:
                self.assertEqual(f.readline().split(), ["founder_id", "sequence", "v_call", "d_call", "j_call", "v_end"])
            self.assertEqual(pf.main(["--immunesim", fa, "--v-germline", fa, "--out", out]), 2)


if __name__ == "__main__":
    unittest.main()
