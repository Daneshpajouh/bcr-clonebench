"""Tests for bin/fetch_ogrdb.py (the download reads a local file:// URL)."""

import hashlib
import io
import os
import pathlib
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import fetch_ogrdb as fo

SET = (">IGHV1-2*02\nCAGGTGCAG...CTGGTG\n"
       ">IGHV4-30-4*i02\nCAGGTG...CAGCTG\n"
       ">IGHD3-10*01\nGTATTACTATGG\n"
       ">IGHJ4*02\nACTAC..TTTGACTAC\n").encode()


def lock_text(url, md5, sha256=fo.NOT_PINNED):
    return ("# comment line\n"
            f"ogrdb_igh_vdj\t10\t{url}\t{md5}\t{sha256}\n")


class FetchOgrdb(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.dir = self.tmp.name
        self.src = os.path.join(self.dir, "set_gapped.fasta")
        with open(self.src, "wb") as fh:
            fh.write(SET)
        self.url = pathlib.Path(self.src).as_uri()
        self.md5 = "md5:" + hashlib.md5(SET).hexdigest()
        self.sha = "sha256:" + hashlib.sha256(SET).hexdigest()
        self.lock = os.path.join(self.dir, "reference.lock")
        self.out = os.path.join(self.dir, "out")

    def tearDown(self):
        self.tmp.cleanup()

    def write_lock(self, text):
        with open(self.lock, "w", encoding="utf-8") as fh:
            fh.write(text)

    def run_main(self, argv):
        out, err = io.StringIO(), io.StringIO()
        with redirect_stdout(out), redirect_stderr(err):
            code = fo.main(argv)
        return code, out.getvalue(), err.getvalue()

    def read_lock_row(self):
        return fo.read_lock(self.lock)

    def test_pin_writes_the_sha256_and_splits_the_set(self):
        self.write_lock(lock_text(self.url, self.md5))
        code, out, _ = self.run_main(["--pin", self.lock, self.out])
        self.assertEqual(code, 0)
        self.assertEqual(self.read_lock_row()["sha256"], self.sha)
        self.assertIn("2 IGHV, 1 IGHD, 1 IGHJ", out)
        with open(os.path.join(self.out, "IGHV.fasta")) as fh:
            v = fh.read()
        with open(os.path.join(self.out, "IGHJ.fasta")) as fh:
            j = fh.read()
        self.assertIn("CAGGTGCAG...CTGGTG", v)  # V keeps its IMGT gaps
        self.assertEqual(j, ">IGHJ4*02\nACTACTTTGACTAC\n")  # J without dots

    def test_a_pinned_lock_downloads_and_checks_without_pin(self):
        self.write_lock(lock_text(self.url, self.md5, self.sha))
        code, _, _ = self.run_main([self.lock, self.out])
        self.assertEqual(code, 0)
        self.assertTrue(os.path.exists(os.path.join(self.out, "IGHD.fasta")))

    def test_an_unpinned_lock_stops_without_pin(self):
        self.write_lock(lock_text(self.url, self.md5))
        code, _, err = self.run_main([self.lock, self.out])
        self.assertEqual(code, 1)
        self.assertIn("make pin-reference", err)
        self.assertFalse(os.path.exists(self.out))

    def test_a_wrong_published_md5_stops_and_pins_nothing(self):
        self.write_lock(lock_text(self.url, "md5:" + "0" * 32))
        with self.assertRaisesRegex(ValueError, "published value"):
            self.run_main(["--pin", self.lock, self.out])
        self.assertEqual(self.read_lock_row()["sha256"], fo.NOT_PINNED)
        self.assertFalse(os.path.exists(self.out))

    def test_a_file_that_differs_from_the_pinned_sha256_stops(self):
        self.write_lock(lock_text(self.url, self.md5, "sha256:" + "a" * 64))
        for argv in ([self.lock, self.out], ["--pin", self.lock, self.out]):
            with self.assertRaisesRegex(ValueError, "pins sha256:a"):
                self.run_main(argv)
        self.assertEqual(self.read_lock_row()["sha256"], "sha256:" + "a" * 64)

    def test_a_pinned_sha256_is_never_rewritten(self):
        self.write_lock(lock_text(self.url, self.md5, self.sha))
        with self.assertRaisesRegex(ValueError, "already pinned"):
            fo.write_sha256(self.lock, "b" * 64)

    def test_the_lock_row_needs_five_fields(self):
        self.write_lock("ogrdb_igh_vdj\t10\turl\tmd5:x\n")
        with self.assertRaisesRegex(ValueError, "5 tab-separated fields"):
            fo.read_lock(self.lock)

    def test_a_record_outside_v_d_and_j_stops_the_split(self):
        with self.assertRaisesRegex(ValueError, "not IGHV, IGHD or IGHJ"):
            fo.split([("IGHV1-2*02", "CAG"), ("IGHD1-1*01", "GG"), ("IGHJ1*01", "AC"), ("IGHA1*01", "AC")])

    def test_a_missing_segment_stops_the_split(self):
        with self.assertRaisesRegex(ValueError, "no records for IGHD"):
            fo.split([("IGHV1-2*02", "CAG"), ("IGHJ1*01", "AC")])

    def test_the_repository_lock_names_version_10_and_its_published_md5(self):
        lock = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "reference.lock")
        row = fo.read_lock(lock)
        self.assertEqual(row["version"], "10")
        self.assertIn("/records/20409587/", row["url"])
        self.assertEqual(row["md5"], "md5:72f5e8eced4b7fe4e5536ce40b07ea99")
        self.assertRegex(row["sha256"], r"^sha256:([0-9a-f]{64}|not-yet-pinned)$")


if __name__ == "__main__":
    unittest.main()
