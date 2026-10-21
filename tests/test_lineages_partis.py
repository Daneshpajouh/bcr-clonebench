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
    return re.findall(r'"""(.*?)"""', script, re.DOTALL)


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
            helps = [i for i, ln in enumerate(ls) if "--help > partis_help_" in ln]   # parse-output's help is separate
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
        self.assertIn("--force-dont-generate-germline-set", ls[sims[0]])
        checks = [i for i, ln in enumerate(ls) if "grep -q -- '--force-dont-generate-germline-set'" in ln]
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

    def test_option_evidence_is_printed_before_the_guard(self):
        # CI run 10: the guard stopped the run and nothing showed why; print the help lines and
        # the partis source lines that name the option first, and accept the option from either.
        ls = lines(script_blocks()[1])
        guard = next(i for i, ln in enumerate(ls) if "grep -q -- '--force-dont-generate-germline-set'" in ln)
        before = "\n".join(ls[:guard])
        self.assertIn("grep -n -i germline partis_help_${level}.txt >&2", before)
        self.assertIn("head -n 40 partis_help_${level}.txt >&2", before)
        self.assertIn("/partis/bin/partis /partis/python", before)
        self.assertIn('if [ -n "\\$src" ]', before)       # '(none)' printed when the source has no match
        # from the source, only the option as bin/partis defines it for simulate counts (CI run 11:
        # partis builds options from dicts, so an add_argument pattern found nothing)
        self.assertIn("subargs['simulate'].append({'name' : '--force-dont-generate-germline-set'", ls[guard])
        self.assertIn("grep -qF", ls[guard])

    def test_one_process_and_the_airr_table_checked_before_airr_truth(self):
        # CI run 12: with several processes partis merged only the .yaml; the AIRR .tsv stayed in sub-process folders.
        ls = lines(script_blocks()[1])
        sim = next(i for i, ln in enumerate(ls) if ln.startswith("partis simulate --simulate-from-scratch"))
        self.assertIn("--n-procs 1", " ".join(ls[sim:sim + 3]))
        truth = next(i for i, ln in enumerate(ls) if "${airr_truth_py}" in ln)
        self.assertEqual(ls[truth - 2], "ls -la >&2")
        self.assertTrue(ls[truth - 1].startswith("test -s simu_${level}.tsv ||"))
        self.assertIn("exit 1", ls[truth - 1])

    def test_a_missing_airr_table_is_converted_from_the_yaml_with_partis_parse_output(self):
        # CI run 13: this partis writes only simu_<level>.yaml from simulate, whatever the process count.
        ls = lines(script_blocks()[1])
        sim = next(i for i, ln in enumerate(ls) if ln.startswith("partis simulate --simulate-from-scratch"))
        conv = next(i for i, ln in enumerate(ls) if "parse-output.py simu_${level}.yaml simu_${level}.tsv" in ln)
        listing = next(i for i, ln in enumerate(ls) if i > sim and ln == "ls -la >&2")
        self.assertLess(sim, conv)
        self.assertLess(conv, listing)
        self.assertIn("--airr-output", ls[conv])
        self.assertIn("|| rc=\\$?", ls[conv])
        self.assertTrue(any("parse-output.py --help > parse_output_help_${level}.txt" in ln for ln in ls[sim:conv]))
        stop = ls[conv + 1]
        self.assertIn("test -s simu_${level}.tsv", stop)
        self.assertIn("cat parse_output_help_${level}.txt >&2", stop)
        self.assertIn("exit 1", stop)
        self.assertEqual(ls[sim:listing].count("fi"), 1)
        # CI run 13 tidy: the help call keeps its own exit code, so the conversion's rc is not reset after it.
        self.assertTrue(any(ln.startswith("help_rc=0; parse-output.py --help") and ln.endswith("|| help_rc=\\$?")
                            for ln in ls[sim:conv]))
        self.assertIn("(exit \\$rc)", stop)
        self.assertIn("(exit \\$help_rc)", stop)


class Workflow(unittest.TestCase):
    def test_failed_task_logs_are_uploaded_and_results_always(self):
        with open(os.path.join(os.path.dirname(MODULE), "..", ".github", "workflows", "ci.yml"), encoding="utf-8") as f:
            ci = f.read()
        step = ci.split("name: Keep the task folders' scripts and logs when a step failed", 1)[1]
        self.assertIn("if: failure()", step.split("- ", 1)[0] + step[:400])
        for pattern in ("work/**/.command.err", "work/**/.command.sh", "work/**/partis_help_*.txt"):
            self.assertIn(pattern, step)
        results = ci.split("name: results", 1)[0].rsplit("- uses: actions/upload-artifact@v6", 1)[1]
        self.assertIn("if: always()", results)


if __name__ == "__main__":
    unittest.main()
