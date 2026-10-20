"""Static checks that the container job leaves room on a hosted runner's disk."""

import os
import re
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return f.read()


def recipe(makefile, target):
    m = re.search(rf"^{target}:.*\n((?:\t.*\n)+)", makefile, re.M)
    return [line.strip() for line in m.group(1).splitlines()]


class WorkflowTest(unittest.TestCase):
    def test_known_answer_frees_disk_first(self):
        text = read(".github", "workflows", "ci.yml")
        job = text[text.index("  known-answer:"):]
        first = re.search(r"^      - (.*)$", job, re.M).group(1)
        self.assertEqual(first, "name: Free runner disk")
        step = job[:job.index("actions/checkout")]
        self.assertEqual(step.count("df -h /"), 2)
        self.assertIn("docker system prune -af", step)
        for path in ("/usr/share/dotnet", "/usr/local/lib/android", "/opt/ghc",
                     "/opt/hostedtoolcache/CodeQL"):
            self.assertIn(path, step)


class PinTargetTest(unittest.TestCase):
    def setUp(self):
        self.lines = recipe(read("Makefile"), "pin")

    def test_build_before_any_pull(self):
        build = next(i for i, l in enumerate(self.lines) if l.startswith("docker build"))
        pulls = [i for i, l in enumerate(self.lines) if l.startswith("docker pull")]
        self.assertEqual(len(pulls), 2)
        self.assertLess(build, min(pulls))

    def test_free_space_after_each_image(self):
        for i, line in enumerate(self.lines):
            if line.startswith(("docker build", "docker pull")):
                self.assertEqual(self.lines[i + 1], "df -h /")


if __name__ == "__main__":
    unittest.main()
