"""The ANNOTATE stub's own output goes through score_clones.py as the SCORE
step would read it (CI run 16: the stub's drops table had a header and no
rows, and score_clones refused it, so 'make stub' stopped). Nextflow cannot
run here, so the stub section is rendered as Nextflow would hand it to bash:
Groovy escapes resolved and the two inputs filled in."""

import json
import os
import re
import subprocess
import sys
import tempfile
import unittest

ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..")
sys.path.insert(0, os.path.join(ROOT, "bin"))
import score_clones as sc


def stub_script(level, reads):
    with open(os.path.join(ROOT, "modules", "annotate.nf"), encoding="utf-8") as f:
        text = f.read()
    body = re.search(r'^\s*stub:\s*\n\s*"""\n(.*?)\n\s*"""', text, re.DOTALL | re.MULTILINE).group(1)
    escapes = {"t": "\t", "n": "\n"}  # Groovy string escapes; \\ to \ and \$ to $ keep the second character
    body = re.sub(r"\\(.)", lambda m: escapes.get(m.group(1), m.group(1)), body)
    return body.replace("${level}", level).replace("${reads}", reads)


class AnnotateStub(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.TemporaryDirectory()
        self.addCleanup(self.d.cleanup)

    def run_stub(self):
        reads = os.path.join(self.d.name, "reads.fasta")
        with open(reads, "w") as f:
            f.write(">r1\nACGT\n>r2\nGGCC\n")
        subprocess.run(["bash", "-euc", stub_script("low", reads)], cwd=self.d.name, check=True)
        return os.path.join(self.d.name, "db_low.tsv"), os.path.join(self.d.name, "drops_low.tsv")

    def test_stub_writes_the_drops_header_and_no_rows(self):
        _, drops = self.run_stub()
        with open(drops, encoding="utf-8") as f:
            self.assertEqual(f.read(), "sequence_id\treason\n")

    def test_score_clones_reads_the_stub_output(self):
        db, drops = self.run_stub()
        rows = sc.read_tsv(db, ("sequence_id", "v_call", "j_call"))
        self.assertEqual([r["sequence_id"] for r in rows], ["r1", "r2"])
        truth = os.path.join(self.d.name, "truth.tsv")
        pred = os.path.join(self.d.name, "pred.tsv")
        with open(truth, "w") as t, open(pred, "w") as p:
            t.write("\t".join(sc.TRUTH_COLUMNS) + "\n")
            p.write("sequence_id\tclone_id\tv_call\tj_call\tmu_freq\n")
            for r in rows:
                t.write(f"{r['sequence_id']}\tA\t{r['v_call']}\t{r['j_call']}\t0\n")
                p.write(f"{r['sequence_id']}\t1\t{r['v_call']}\t{r['j_call']}\t0\n")
        out = os.path.join(self.d.name, "score.json")
        args = ["--truth", truth, "--pred", pred, "--mu-column", "mu_freq", "--label", "low/changeo",
                "--drops", drops, "--out-json", out]
        self.assertEqual(sc.main(args), 0)
        with open(out) as f:
            r = json.load(f)
        self.assertEqual((r["n_scored"], r["dropped_by_reason"]), (2, {}))


if __name__ == "__main__":
    unittest.main()
