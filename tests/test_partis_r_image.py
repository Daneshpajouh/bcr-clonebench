"""Static checks on the partis-r image and its wiring (no Docker needed).

CI run 6: partis simulate --simulate-from-scratch calls R for TreeSim, and the
partis image has no R. containers/partis-r adds R and TreeSim on top of the
locked partis digest; 'make check-partis-r' runs R and loads TreeSim in the
built image, as the host user, before 'make pin' records it."""

import os
import re
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")


def read(*parts):
    with open(os.path.join(ROOT, *parts), encoding="utf-8") as f:
        return f.read()


def arg(text, name):
    m = re.search(rf"^ARG {name}=(\S+)$", text, re.M)
    return m.group(1) if m else None


def make_target(text, name):
    m = re.search(rf"^{name}:.*?\n((?:\t.*\n)+)", text, re.M)
    return m.group(1) if m else ""


class DockerfileTest(unittest.TestCase):
    def setUp(self):
        self.text = read("containers", "partis-r", "Dockerfile")

    def test_base_is_the_build_argument_without_default(self):
        self.assertRegex(self.text, r"(?m)^ARG PARTIS_BASE$")
        self.assertRegex(self.text, r"(?m)^FROM \$\{PARTIS_BASE\}$")
        self.assertLess(self.text.index("ARG PARTIS_BASE"), self.text.index("FROM "))

    def test_r_and_treesim_fixed_by_version_and_build(self):
        self.assertRegex(arg(self.text, "R_BASE_SPEC") or "", r"^r-base=\d+\.\d+\.\d+=h[0-9a-f]{7}_\d+$")
        self.assertRegex(arg(self.text, "R_TREESIM_SPEC") or "", r"^r-treesim=\d+\.\d+=r\d\dh[0-9a-f]{7}_\d+$")

    def test_treesim_build_matches_r_minor_version(self):
        r_ver = arg(self.text, "R_BASE_SPEC").split("=")[1]
        build = arg(self.text, "R_TREESIM_SPEC").split("=")[2]
        major, minor = r_ver.split(".")[:2]
        self.assertTrue(build.startswith(f"r{major}{minor}h"), (r_ver, build))

    def test_conda_forge_only_in_own_environment(self):
        create = re.search(r"micromamba create .*", self.text).group(0)
        self.assertIn("-n r ", create)
        self.assertIn("--override-channels -c conda-forge", create)
        self.assertNotIn("-n base", self.text)
        self.assertNotIn("micromamba install", self.text)

    def test_solved_list_kept_and_cache_cleaned(self):
        self.assertIn("micromamba list -n r > /opt/conda/envs/r/packages.txt", self.text)
        self.assertIn("micromamba clean --all --yes", self.text)

    def test_r_environment_last_on_path(self):
        self.assertRegex(self.text, r"(?m)^ENV PATH=\$PATH:/opt/conda/envs/r/bin$")

    def test_build_stops_unless_r_runs_and_treesim_loads(self):
        check = self.text.index("RUN R --version")
        self.assertGreater(check, self.text.index("ENV PATH="))
        self.assertIn("library(TreeSim)", self.text[check:])


class WiringTest(unittest.TestCase):
    def setUp(self):
        self.make = read("Makefile")

    def test_pin_builds_on_the_digest_just_pulled(self):
        # not the digest already in images.lock: a tag changed on purpose, or a
        # first pin, must build on the image pin_images.py is about to record
        pin = make_target(self.make, "pin")
        inspect = pin.index("base=$$(docker image inspect --format '{{index .RepoDigests 0}}' $(PARTIS_TAG))")
        guard = pin.index("echo \"$$base\" | grep -Eq '@sha256:[0-9a-f]{64}$$'")
        build = pin.index("docker build --build-arg PARTIS_BASE=$$base -t $(PARTIS_R_TAG) containers/partis-r")
        self.assertLess(pin.index("docker pull $(PARTIS_TAG)"), inspect)
        self.assertLess(inspect, guard)
        self.assertLess(guard, build)
        self.assertNotIn("PARTIS_BASE=$(PARTIS)", pin)

    def test_pin_checks_the_image_before_recording_it(self):
        pin = make_target(self.make, "pin")
        self.assertLess(pin.index("$(MAKE) check-partis-r PARTIS_R_CHECK=$(PARTIS_R_TAG)"), pin.index("pin_images.py"))
        self.assertIn("partis_r_image $(PARTIS_R_TAG)", pin)

    def test_check_runs_r_and_treesim_as_host_user(self):
        check = make_target(self.make, "check-partis-r")
        self.assertIn("docker run --rm -u $$(id -u):$$(id -g) $(PARTIS_R_CHECK)", check)
        self.assertIn("PARTIS_R_CHECK ?= $(PARTIS_R)", self.make)
        self.assertIn("R --version", check)
        self.assertIn("library(TreeSim)", check)

    def test_known_answer_uses_partis_r(self):
        self.assertIn("--partis_image $(PARTIS_R)\n", make_target(self.make, "known-answer"))

    def test_lock_has_partis_r_row(self):
        rows = [ln.split("\t") for ln in read("images.lock").splitlines() if not ln.startswith("#")]
        self.assertIn("partis_r_image", [r[0] for r in rows if len(r) == 3])


class ModuleTest(unittest.TestCase):
    def setUp(self):
        text = read("modules", "lineages_partis.nf")
        script = text.split("script:", 1)[1].split("stub:", 1)[0]
        self.partis_block = re.findall(r'"""(.*?)"""', script, re.S)[1]

    def test_missing_r_fails_fast_before_simulation(self):
        lines = [ln.strip() for ln in self.partis_block.strip().splitlines() if ln.strip()]
        self.assertTrue(lines[2].startswith("command -v R > /dev/null ||"), lines[2])
        self.assertIn("exit 127", lines[2])
        self.assertLess(self.partis_block.index("command -v R"), self.partis_block.index("partis simulate --simulate"))

    def test_no_weibull_root_option(self):
        # with a Weibull root partis uses TreeSimGM, which the image does not carry
        self.assertNotIn("weibull", self.partis_block.lower())


if __name__ == "__main__":
    unittest.main()
