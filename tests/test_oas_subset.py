"""Tests for bin/oas_subset.py."""

import csv
import gzip
import io
import json
import os
import sys
import tempfile
import unittest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import oas_subset as osub  # noqa: E402

HEADER = ["sequence_id", "sequence", "productive", "v_call", "j_call", "junction"]


def unit(path, subject, day, rows, isotype="IGHM"):
    meta = json.dumps({"Subject": subject, "Longitudinal": day, "Isotype": isotype})
    buf = io.StringIO()
    w = csv.writer(buf)
    w.writerow([meta])
    w.writerow(HEADER)
    w.writerows(rows)
    with gzip.open(path, "wt", newline="") as fh:
        fh.write(buf.getvalue())


def row(i, seq=None, productive="T", v="IGHV1-2*02", j="IGHJ4*02", junction="TGTGCGAGAT"):
    return [f"s{i}", seq or f"ACGT{i:06d}", productive, v, j, junction]


class Subset(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.TemporaryDirectory()
        self.addCleanup(self.d.cleanup)
        self.p = lambda name: os.path.join(self.d.name, name)

    def write_subject(self, subject, days, n=10, extra=()):
        paths = []
        for day in days:
            path = self.p(f"{subject}_{day}.csv.gz")
            unit(path, subject, day, [row(k) for k in range(n)] + list(extra))
            paths.append(path)
        return paths

    def test_filters_dedupe_and_counts(self):
        bad = [row(900, productive="F"), row(901, v=""), row(902, junction=""), row(3)]  # row(3) repeats a sequence
        paths = self.write_subject("S1", osub.DAYS, n=10, extra=bad)
        out = self.p("out")
        osub.main([*paths, "--out", out, "--cap", "100"])
        m = json.load(open(os.path.join(out, "manifest.json")))
        c = m["days"]["Day-0"]
        self.assertEqual((c["rows_in"], c["not_productive"], c["no_v_j_or_junction"], c["duplicate_sequence"],
                          c["rows_out"], c["capped"]), (14, 1, 2, 1, 10, False))
        self.assertEqual(m["subject"], "S1")

    def test_cap_is_order_free_and_keeps_smallest_hashes(self):
        rows_ = [row(k) for k in range(50)]
        a, b = self.p("a.csv.gz"), self.p("b.csv.gz")
        unit(a, "S1", "Day-0", rows_)
        unit(b, "S1", "Day-0", list(reversed(rows_)))
        ka, _ = osub.subset_day([osub.Path(a)], 7, "salt")
        kb, _ = osub.subset_day([osub.Path(b)], 7, "salt")
        self.assertEqual([r["sequence"] for r in ka], [r["sequence"] for r in kb])
        want = sorted((osub.hashlib.sha256(f"salt\0{r[1]}".encode()).hexdigest(), r[1]) for r in rows_)[:7]
        self.assertEqual(sorted(r["sequence"] for r in ka), sorted(s for _, s in want))
        kc, cc = osub.subset_day([osub.Path(a)], 7, "other salt")
        self.assertTrue(cc["capped"])
        self.assertNotEqual(sorted(r["sequence"] for r in kc), sorted(r["sequence"] for r in ka))

    def test_first_full_subject_and_named_subject(self):
        paths = (self.write_subject("S2", osub.DAYS) + self.write_subject("S1", osub.DAYS[:4])
                 + self.write_subject("S3", osub.DAYS))
        out = self.p("out")
        osub.main([*paths, "--out", out])
        self.assertEqual(json.load(open(os.path.join(out, "manifest.json")))["subject"], "S2")
        with self.assertRaises(SystemExit):
            osub.main([*paths, "--out", self.p("out2"), "--subject", "S1"])  # misses Day-28

    def test_isotype_units_pool(self):
        paths = self.write_subject("S1", osub.DAYS, n=3)
        g = self.p("S1_Day-0_IGHG.csv.gz")
        unit(g, "S1", "Day-0", [row(k, seq=f"GGGG{k}") for k in range(4)], isotype="IGHG")
        out = self.p("out")
        osub.main([*paths, g, "--out", out])
        self.assertEqual(json.load(open(os.path.join(out, "manifest.json")))["days"]["Day-0"]["rows_out"], 7)

    def test_missing_metadata_key_and_column_stop(self):
        path = self.p("bad.csv.gz")
        with gzip.open(path, "wt") as fh:
            fh.write('"{""Subject"": ""S1""}"\n' + ",".join(HEADER) + "\n")
        with self.assertRaises(SystemExit):
            osub.read_meta(osub.Path(path))
        path2 = self.p("bad2.csv.gz")
        unit(path2, "S1", "Day-0", [])
        with gzip.open(path2, "wt") as fh:
            fh.write('"{""Subject"": ""S1"", ""Longitudinal"": ""Day-0""}"\nsequence_id,sequence\n')
        with self.assertRaises(SystemExit):
            list(osub.rows(osub.Path(path2)))

    def test_out_folder_never_reused(self):
        paths = self.write_subject("S1", osub.DAYS, n=2)
        out = self.p("out")
        osub.main([*paths, "--out", out])
        with self.assertRaises(FileExistsError):
            osub.main([*paths, "--out", out])


if __name__ == "__main__":
    unittest.main()
