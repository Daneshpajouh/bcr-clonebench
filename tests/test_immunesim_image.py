"""Static checks on the immuneSIM image recipe (no Docker needed)."""

import os
import re
import unittest

HERE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "containers", "immunesim")

# The Imports field of immuneSIM 0.8.7 (DESCRIPTION, packaged 2019-09-20).
IMPORTS = {
    "poweRlaw", "stringdist", "Biostrings", "igraph", "stringr", "data.table",
    "plyr", "reshape2", "ggplot2", "grid", "ggthemes", "RColorBrewer",
    "Metrics", "repmis",
}


def read(name):
    with open(os.path.join(HERE, name), encoding="utf-8") as f:
        return f.read()


def arg(dockerfile, name):
    m = re.search(rf"^ARG {name}=(\S+)$", dockerfile, re.M)
    return m.group(1) if m else None


class DockerfileTest(unittest.TestCase):
    def setUp(self):
        self.text = read("Dockerfile")

    def test_base_image_pinned_by_digest(self):
        self.assertRegex(self.text, r"(?m)^FROM \S+:\S+@sha256:[0-9a-f]{64}$")

    def test_tarball_checked_before_install(self):
        sha = arg(self.text, "IMMUNESIM_SHA256")
        self.assertRegex(sha or "", r"^[0-9a-f]{64}$")
        self.assertIn("sha256sum -c -", self.text)
        self.assertLess(self.text.index("sha256sum -c -"), self.text.index("Rscript /opt/install.R"))

    def test_tarball_from_cran_archive(self):
        self.assertEqual(
            arg(self.text, "IMMUNESIM_URL"),
            "https://cran.r-project.org/src/contrib/Archive/immuneSIM/immuneSIM_0.8.7.tar.gz",
        )

    def test_snapshot_and_bioc_fixed(self):
        self.assertRegex(arg(self.text, "CRAN_SNAPSHOT") or "", r"/cran/\d{4}-\d{2}-\d{2}$")
        self.assertRegex(arg(self.text, "BIOC_VERSION") or "", r"^\d+\.\d+$")


class InstallScriptTest(unittest.TestCase):
    def test_every_cran_import_listed(self):
        text = read("install.R")
        block = re.search(r"cran_imports <- c\((.*?)\)", text, re.S).group(1)
        listed = set(re.findall(r'"([^"]+)"', block))
        self.assertEqual(listed, IMPORTS - {"grid", "Biostrings"})
        self.assertIn('BiocManager::install("Biostrings"', text)


if __name__ == "__main__":
    unittest.main()
