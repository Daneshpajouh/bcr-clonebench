"""The harder setting (Director ruling 2026-10-20): with shared_vj each level also runs as <level>_shared, where
every partis family uses one V allele, one J allele and one junction length. The two partis options were read from
partis simulate --help in the pinned image (CI run 18 results/tool-help)."""

import os
import sys
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import summarize_scores as ss


def read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return f.read()


class SharedSetting(unittest.TestCase):
    def test_off_by_default_and_on_in_the_test_profile(self):
        config = read("nextflow.config")
        self.assertRegex(config, r"\n    shared_vj = false\n")
        self.assertRegex(config.split("test {", 1)[1], r"params\.shared_vj = true")
        self.assertIn("shared_vj_genes = 'IGHV3-23*01:IGHJ4*02'", config)
        self.assertIn("shared_junction_length = 48", config)

    def test_main_adds_shared_levels_for_partis_only(self):
        main = read("main.nf")
        self.assertIn('params.levels.collect { "${it}_shared" }', main)
        self.assertRegex(main, r"if \(params\.shared_vj && params\.lineage_simulator != 'partis'\)\s+error")
        self.assertRegex(main, r"if \(params\.shared_vj && params\.lineage_mode == 'bcr-phylo'\)\s+error")

    def test_partis_gets_both_options_and_checks_the_alleles_first(self):
        mod = read("modules", "lineages_partis.nf")
        self.assertIn("params.mutation_multiplier[level.replaceAll(/_shared$/, '')]", mod)
        self.assertIn("--only-genes '${params.shared_vj_genes}' --allowed-cdr3-lengths ${params.shared_junction_length}",
                      mod)
        partis = mod.split("else\n", 1)[1]
        self.assertLess(partis.index("for g in ${shared_genes}"), partis.index("partis simulate --simulate-from-scratch"))
        self.assertIn(""".collect { "'${it}'" }""", mod)   # each allele name quoted, so '*' never globs
        self.assertIn("${shared_args}", partis.split("partis simulate --simulate-from-scratch", 1)[1].split("--outfname")[0])

    def test_shared_rows_sort_after_the_plain_levels(self):
        order = sorted(["high_shared", "low", "high", "low_shared"], key=lambda l: ss.LEVEL_ORDER.get(l, 9))
        self.assertEqual(order, ["low", "high", "low_shared", "high_shared"])


if __name__ == "__main__":
    unittest.main()
