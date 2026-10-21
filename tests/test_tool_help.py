"""Each tool whose options a README label names keeps its --help text in
the results (tool-help/), read from the image the run used (Director ruling
2026-10-20: confirm single linkage from DefineClones --help in the suite
image, and the partis option from partis --help, before a label goes in)."""

import os
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return f.read()


class ToolHelp(unittest.TestCase):
    def test_changeo_keeps_defineclones_help_and_publishes_it(self):
        text = read("modules", "clones.nf")
        changeo = text.split("if (method == 'changeo')", 1)[1].split("else", 1)[0]
        self.assertIn("DefineClones.py --help > DefineClones_help_${level}.txt", changeo)
        self.assertLess(changeo.index("--help"), changeo.index("DefineClones.py -d"))  # kept even if grouping fails
        self.assertRegex(text, r"publishDir \"\$\{params.outdir\}/tool-help\".*pattern: 'DefineClones_help_\*\.txt'")
        self.assertRegex(text, r'path "DefineClones_help_\$\{level\}\.txt", optional: true, emit: help')

    def test_partis_help_is_published(self):
        text = read("modules", "lineages_partis.nf")
        self.assertRegex(text, r"publishDir \"\$\{params.outdir\}/tool-help\".*pattern: 'partis_help_\*\.txt'")
        self.assertIn('path "partis_help_${level}.txt", emit: help', text)

    def test_main_takes_the_named_clones_output(self):
        self.assertRegex(read("main.nf"), r"clones = CLONES\(.*\)\.clones\b")


if __name__ == "__main__":
    unittest.main()
