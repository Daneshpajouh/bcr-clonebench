"""Checks on modules/lineages_partis.nf: the partis image keeps its tools in
/partis/bin, which is not on PATH (CI run 5: 'partis: command not found'),
so each script must put that folder on PATH and stop with a clear message
when the tool is missing, before any other line runs."""

import os
import re
import unittest

MODULE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..",
                      "modules", "lineages_partis.nf")


def script_blocks():
    with open(MODULE, encoding="utf-8") as f:
        text = f.read()
    script = text.split("script:", 1)[1].split("stub:", 1)[0]
    return re.findall(r'"""(.*?)"""', script, re.S)


def lines(block):
    return [ln.strip() for ln in block.strip().splitlines() if ln.strip()]


class LineagesPartisScript(unittest.TestCase):
    def test_two_modes(self):
        self.assertEqual(len(script_blocks()), 2)

    def test_partis_bin_on_path_first(self):
        for block in script_blocks():
            first = lines(block)[0]
            self.assertRegex(first, r"^export PATH=\S*/partis/bin:")

    def test_r_folder_on_path_and_path_printed_before_checks(self):
        # CI run 7: Nextflow builds the task PATH itself, so the image's ENV PATH tail
        # (/opt/conda/envs/r/bin) was not seen; the scripts name it and print PATH.
        for block in script_blocks():
            first, second = lines(block)[:2]
            self.assertIn(":/opt/conda/envs/r/bin:", first)
            self.assertEqual(second, 'echo "PATH=\\$PATH" >&2')

    def test_task_folder_listed_after_path_print(self):
        # CI run 8: the run log must show what was staged into the task folder.
        for block in script_blocks():
            self.assertEqual(lines(block)[2], "ls -la >&2")

    def test_missing_tool_fails_fast_with_message(self):
        tools = ["bcr-phylo-run.py", "partis"]
        for block, tool in zip(script_blocks(), tools):
            check = lines(block)[3]
            self.assertIn(f"/partis/bin/{tool}", check)
            self.assertIn("test -x", check)
            self.assertIn("exit 127", check)
            self.assertIn(">&2", check)

    def test_help_capture_does_not_mask_failure(self):
        # A non-zero --help must not stop the run (bash -ue), but it must
        # not pass silently either: its exit code is kept and reported.
        for block in script_blocks():
            ls = lines(block)
            helps = [i for i, ln in enumerate(ls) if "--help >" in ln]
            self.assertEqual(len(helps), 1)
            ln = ls[helps[0]]
            self.assertNotIn("|| true", ln)
            self.assertIn("|| rc=\\$?", ln)
            self.assertRegex(ls[helps[0] + 1], r"test \\\$rc -eq 0 \|\| echo .*>&2$")

    def test_partis_uses_its_shipped_germline_set(self):
        # CI run 9: rearranging from scratch made partis generate a germline set whose
        # alleles rarely gave productive rearrangements; one sub-process kept retrying.
        block = script_blocks()[1]
        ls = lines(block)
        sims = [i for i, ln in enumerate(ls) if ln.startswith("partis simulate --simulate-from-scratch")]
        self.assertEqual(len(sims), 1)
        self.assertIn("--dont-generate-germline-set", ls[sims[0]])
        checks = [i for i, ln in enumerate(ls) if "grep -q -- '--dont-generate-germline-set'" in ln]
        self.assertEqual(len(checks), 1)
        self.assertLess(checks[0], sims[0])
        self.assertIn("partis_help_${level}.txt", ls[checks[0]])
        self.assertRegex(ls[checks[0]], r"exit [1-9]\d*; }$")
        self.assertIn(">&2", ls[checks[0]])

    def test_matplotlib_cache_in_the_task_folder_before_the_tool_runs(self):
        for block in script_blocks():
            ls = lines(block)
            mpl = [i for i, ln in enumerate(ls) if ln.startswith("export MPLCONFIGDIR=")]
            self.assertEqual(len(mpl), 1)
            self.assertIn("\\$PWD/", ls[mpl[0]])
            first_tool = next(i for i, ln in enumerate(ls) if "--help >" in ln)
            self.assertLess(mpl[0], first_tool)


if __name__ == "__main__":
    unittest.main()
