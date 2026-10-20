"""Tests for bin/pin_images.py (docker is replaced by a fixed table)."""

import io
import os
import sys
import tempfile
import unittest
from contextlib import redirect_stderr, redirect_stdout

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "bin"))
import pin_images as pin  # noqa: E402

A = "a" * 64
B = "b" * 64
C = "c" * 64


def lock_text(partis_digest):
    return ("# comment line\n"
            "suite_image\timmcantation/suite:4.5.0\tdigest-not-yet-pinned\n"
            "immunesim_image\tcontainers/immunesim\tdigest-not-yet-pinned\n"
            f"partis_image\tquay.io/matsengrp/partis:latest\t{partis_digest}\n")


class PinImages(unittest.TestCase):
    def run_pin(self, text, found):
        with tempfile.TemporaryDirectory() as d:
            path = os.path.join(d, "images.lock")
            with open(path, "w", encoding="utf-8") as f:
                f.write(text)
            old = pin.digest
            pin.digest = lambda tag: found[tag]
            try:
                with redirect_stdout(io.StringIO()), redirect_stderr(io.StringIO()) as err:
                    rc = pin.main([path, "suite_image", "immcantation/suite:4.5.0",
                                   "immunesim_image", "bcr-clonebench/immunesim:local",
                                   "partis_image", "quay.io/matsengrp/partis:latest"])
            finally:
                pin.digest = old
            with open(path, encoding="utf-8") as f:
                return rc, f.read(), err.getvalue()

    def found(self, partis):
        return {"immcantation/suite:4.5.0": f"immcantation/suite@sha256:{A}",
                "bcr-clonebench/immunesim:local": f"sha256:{B}",
                "quay.io/matsengrp/partis:latest": f"quay.io/matsengrp/partis@sha256:{partis}"}

    def test_unpinned_rows_get_the_pulled_digest(self):
        rc, text, _ = self.run_pin(lock_text("digest-not-yet-pinned"), self.found(C))
        self.assertEqual(rc, 0)
        self.assertIn(f"suite_image\timmcantation/suite:4.5.0\timmcantation/suite@sha256:{A}\n", text)
        self.assertIn(f"immunesim_image\tbcr-clonebench/immunesim:local\tsha256:{B}\n", text)
        self.assertIn(f"partis_image\tquay.io/matsengrp/partis:latest\tquay.io/matsengrp/partis@sha256:{C}\n", text)

    def test_locked_digest_that_still_matches_is_kept(self):
        before = lock_text(f"quay.io/matsengrp/partis@sha256:{C}")
        rc, text, _ = self.run_pin(before, self.found(C))
        self.assertEqual(rc, 0)
        self.assertIn(f"\tquay.io/matsengrp/partis@sha256:{C}\n", text)

    def test_moved_tag_is_refused_and_the_lock_is_not_written(self):
        before = lock_text(f"quay.io/matsengrp/partis@sha256:{C}")
        rc, text, err = self.run_pin(before, self.found(A))
        self.assertNotEqual(rc, 0)
        self.assertEqual(text, before)
        self.assertIn("partis_image", err)
        self.assertIn("moved", err)

    def test_locked_digest_for_another_repository_is_refused(self):
        before = lock_text(f"psathyrella/partis@sha256:{C}")
        rc, text, err = self.run_pin(before, self.found(C))
        self.assertNotEqual(rc, 0)
        self.assertEqual(text, before)

    def test_a_new_tag_on_purpose_is_pinned_not_refused(self):
        before = lock_text(f"quay.io/matsengrp/partis@sha256:{C}").replace(
            "quay.io/matsengrp/partis:latest\t", "quay.io/matsengrp/partis:dev\t")
        rc, text, _ = self.run_pin(before, self.found(A))
        self.assertEqual(rc, 0)
        self.assertIn(f"partis_image\tquay.io/matsengrp/partis:latest\tquay.io/matsengrp/partis@sha256:{A}\n", text)


if __name__ == "__main__":
    unittest.main()
