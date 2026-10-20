"""Every bin script that a container step runs is staged into the task folder
as an input and run with its interpreter (CI run 8: 'airr_truth.py: command
not found' inside the partis-r container, because the project's bin folder
was not reachable there). Steps without a container run on the host, where
Nextflow puts bin on PATH, and are not checked here."""

import os
import re
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
BIN = sorted(f for f in os.listdir(os.path.join(ROOT, "bin")) if f.endswith((".py", ".R")))
INTERPRETER = {".py": "python3", ".R": "Rscript"}


def read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return f.read()


def container_modules():
    for name in sorted(os.listdir(os.path.join(ROOT, "modules"))):
        text = read("modules", name)
        if re.search(r"^\s*container\s", text, re.M):
            yield name, text


def body(text):
    """script and stub sections, without comments."""
    code = text.split("script:", 1)[1]
    return "\n".join(ln for ln in code.splitlines() if not ln.strip().startswith("//"))


def input_paths(text):
    block = text.split("input:", 1)[1].split("output:", 1)[0]
    return re.findall(r"^\s*path\s+(\w+)\s*$", block, re.M)


class StagedBinScripts(unittest.TestCase):
    def test_no_container_step_calls_a_bin_script_by_bare_name(self):
        for name, text in container_modules():
            for script in BIN:
                self.assertNotRegex(body(text), r"(?<![\w$/{])" + re.escape(script),
                                    f"{name} calls {script} by name; stage it and run it with its interpreter")

    def test_each_staged_script_is_an_input_run_by_its_interpreter(self):
        found = {}
        for name, text in container_modules():
            for var in input_paths(text):
                m = re.fullmatch(r"(\w+)_(py|r)", var)
                if not m:
                    continue
                ext = ".py" if m.group(2) == "py" else ".R"
                self.assertIn(f"{INTERPRETER[ext]} ${{{var}}}", body(text), f"{name}: {var}")
                found[var] = name
        self.assertEqual(sorted(found), sorted([
            "airr_truth_py", "clones_scoper_r", "mutation_r", "prepare_founders_py",
            "simulate_founders_r", "threshold_r"]))

    def test_workflow_passes_each_script_from_the_project_bin(self):
        main = read("main.nf")
        for script in ["airr_truth.py", "simulate_founders.R", "prepare_founders.py",
                       "threshold.R", "clones_scoper.R", "mutation.R"]:
            self.assertIn(f'file("${{projectDir}}/bin/{script}")', main)
            self.assertTrue(os.path.isfile(os.path.join(ROOT, "bin", script)), script)

    def test_staged_python_scripts_need_only_the_standard_library(self):
        # The partis-r and immuneSIM images run them with their own python3; no
        # third-party module may be imported (pandas and friends are not there).
        allowed = {"argparse", "csv", "sys", "os", "re", "math", "random", "collections",
                   "itertools", "json", "gzip", "hashlib"}
        for script in ["airr_truth.py", "prepare_founders.py"]:
            mods = re.findall(r"^(?:import|from)\s+(\w+)", read("bin", script), re.M)
            self.assertTrue(set(mods) <= allowed, (script, sorted(set(mods) - allowed)))


if __name__ == "__main__":
    unittest.main()
