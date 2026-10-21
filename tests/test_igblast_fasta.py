"""Tests for bin/igblast_fasta.py and its use in bin/build_reference.sh."""

import os
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import igblast_fasta as ig

# Header shape of the IMGT germline files the suite image ships.
IMGT = (
    ">M99641|IGHV1-18*01|Homo sapiens|F|V-REGION|188..483|296 nt|1| | | | |296+24=320|"
    "partial in 3'|nonfunctional|\n"
    "caggttcagctggtgcag...tctggagct\n"
    "gaggtgaagaagcct\n"
    ">X62106|IGHV1-2*02|Homo sapiens|F|V-REGION|..|296 nt|1| | | | | |\n"
    "caggtgcagctggtgcag...tctggg\n"
)


class ConvertTest(unittest.TestCase):
    def run_on(self, text):
        with tempfile.TemporaryDirectory() as d:
            src, dst = os.path.join(d, "in.fasta"), os.path.join(d, "out.fasta")
            with open(src, "w") as fh:
                fh.write(text)
            rc = ig.main([src, dst])
            with open(dst) as fh:
                return rc, fh.read(), ig.check(dst)

    def test_imgt_headers_cut_to_allele_and_gaps_removed(self):
        rc, out, n_ids = self.run_on(IMGT)
        self.assertEqual(rc, 0)
        self.assertEqual(out, ">IGHV1-18*01\ncaggttcagctggtgcagtctggagctgaggtgaagaagcct\n"
                              ">IGHV1-2*02\ncaggtgcagctggtgcagtctggg\n")
        self.assertEqual(n_ids, 2)

    def test_every_written_id_short_and_unique(self):
        # the run 4 failure: a 53-character first word in an IMGT header
        long_first_word = ">" + "A" * 53 + " rest\nacgt\n"
        self.assertGreater(len(long_first_word.split()[0]) - 1, ig.MAX_ID)
        rc, out, _ = self.run_on(IMGT + long_first_word.replace(">", ">X1|IGHV9-9*01|"))
        self.assertEqual(rc, 0)
        ids = [line[1:] for line in out.splitlines() if line.startswith(">")]
        self.assertTrue(all(0 < len(i) <= ig.MAX_ID for i in ids))
        self.assertEqual(len(ids), len(set(ids)))

    def test_plain_header_keeps_first_word(self):
        rc, out, _ = self.run_on(">IGHJ4*02 OGRDB\nACTACTTTGACTACTGG\n")
        self.assertEqual(rc, 0)
        self.assertEqual(out, ">IGHJ4*02\nACTACTTTGACTACTGG\n")

    def test_same_name_same_sequence_kept_once(self):
        rc, out, _ = self.run_on(">a|IGHD1-1*01|x\nGGTACAACTGGAACGAC\n"
                                 ">b|IGHD1-1*01|y\nggtacaactggaacgac\n")
        self.assertEqual(rc, 0)
        self.assertEqual(out.count(">"), 1)

    def test_same_name_other_sequence_stops(self):
        with self.assertRaises(ValueError):
            ig.convert([("a|IGHD1-1*01|x", "GGTACA"), ("b|IGHD1-1*01|y", "GGTACC")])

    def test_long_or_empty_name_stops(self):
        with self.assertRaises(ValueError):
            ig.convert([("x|" + "V" * 51 + "|y", "ACGT")])
        with self.assertRaises(ValueError):
            ig.convert([("x| |y", "ACGT")])

    def test_check_refuses_long_and_repeated_ids(self):
        with tempfile.TemporaryDirectory() as d:
            p = os.path.join(d, "f.fasta")
            for text in (">" + "A" * 51 + "\nAC\n", ">a\nAC\n>a\nAG\n", ""):
                with open(p, "w") as fh:
                    fh.write(text)
                with self.assertRaises(ValueError):
                    ig.check(p)

    def test_command_line_exit_code_on_bad_input(self):
        with tempfile.TemporaryDirectory() as d:
            src = os.path.join(d, "in.fasta")
            with open(src, "w") as fh:
                fh.write(">x|" + "V" * 51 + "|y\nACGT\n")
            r = subprocess.run([sys.executable, "-B", os.path.join(ROOT, "bin", "igblast_fasta.py"),
                                src, os.path.join(d, "out.fasta")], capture_output=True, text=True, check=False)
            self.assertEqual(r.returncode, 1)
            self.assertIn("over 50 characters", r.stderr)


class BuildReferenceTest(unittest.TestCase):
    def test_headers_cut_and_checked_before_makeblastdb(self):
        with open(os.path.join(ROOT, "bin", "build_reference.sh")) as fh:
            text = fh.read()
        loop = text[text.index("for seg in v d j"):]
        cut = loop.index("igblast_fasta.py\" \"$out/germlines")
        chk = loop.index("igblast_fasta.py\" --check")
        blast = loop.index("makeblastdb -parse_seqids")
        self.assertLess(cut, chk)
        self.assertLess(chk, blast)
        self.assertNotIn("tr -d '.'", loop)


if __name__ == "__main__":
    unittest.main()
