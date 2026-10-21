#!/usr/bin/env python3
"""Named known-answer checks on the summary table.

    check_known_answer.py results.tsv conf/known_answer_floors.json
    check_known_answer.py results.tsv floors.json --expect-fail KA-CHANGEO-PRECISION

Each check names a metric, a floor ("min") or a ceiling ("max") and the rows it
applies to. "<column>/kept" divides the column by the kept share of reads, so a
gene call score is judged on the reads the pipeline kept. A check fails on its
worst row. A missing or NA value in a selected row, or a check that selects no
row, is reported as NO DATA and also stops the run. Exit 0 when every check passes.

With --expect-fail, the run is a mutation test: it exits 0 only when every named
check fails on a real value past its bound (a deliberately wrong germline set or
threshold must be caught). NO DATA does not count as caught: a pipeline that
broke is not a floor that worked.
"""

import argparse
import csv
import json
import sys


def value(row, metric):
    col, _, per = metric.partition("/")
    v = row.get(col, "NA")
    if v in ("", "NA"):
        return None
    v = float(v)
    if per == "kept":
        n = int(row["n_truth"])
        kept = (n - int(row["n_dropped"])) / n if n else 0.0
        return v / kept if kept else None
    if per:
        raise ValueError(f"unknown metric suffix in {metric}")
    return v


def selected(row, where):
    for key, want in where.items():
        if key == "shared":
            if row["level"].endswith("_shared") != want:
                return False
        elif row.get(key) != want:
            return False
    return True


def run_check(rows, check):
    """Return (PASS, FAIL or NO DATA, worst value or None, the row label of the worst or missing row)."""
    rows = [r for r in rows if selected(r, check.get("where", {}))]
    if not rows:
        return "NO DATA", None, "no row selected"
    lower = "min" in check
    bound = check["min"] if lower else check["max"]
    worst, worst_label = None, ""
    for r in rows:
        label = f"{r['level']} e{r['error_rate']} {r['method']}"
        v = value(r, check["metric"])
        if v is None:
            return "NO DATA", None, f"{label}: NA"
        if worst is None or (v < worst if lower else v > worst):
            worst, worst_label = v, label
    passed = worst >= bound if lower else worst <= bound
    return ("PASS" if passed else "FAIL"), worst, worst_label


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("results")
    ap.add_argument("floors")
    ap.add_argument("--expect-fail", default="", help="comma-separated check names that must fail")
    a = ap.parse_args(argv)
    with open(a.results, encoding="utf-8") as f:
        rows = list(csv.DictReader(f, delimiter="\t"))
    with open(a.floors, encoding="utf-8") as f:
        checks = json.load(f)["checks"]
    names = {c["name"] for c in checks}
    expect = [n for n in a.expect_fail.split(",") if n]
    unknown = [n for n in expect if n not in names]
    if unknown:
        print(f"unknown check name: {', '.join(unknown)}", file=sys.stderr)
        return 2
    failed, not_passed = set(), set()
    for c in checks:
        status, worst, where = run_check(rows, c)
        bound = f"min {c['min']}" if "min" in c else f"max {c['max']}"
        shown = "NA" if worst is None else f"{worst:.4f}"
        print(f"{status}\t{c['name']}\t{c['metric']}\t{bound}\tworst {shown}\t{where}")
        if status == "FAIL":
            failed.add(c["name"])
        if status != "PASS":
            not_passed.add(c["name"])
    if expect:
        missed = [n for n in expect if n not in failed]
        print(f"mutation test: expected to fail {', '.join(expect)}; "
              + ("all caught" if not missed else f"NOT caught: {', '.join(missed)}"))
        return 1 if missed else 0
    return 1 if not_passed else 0


if __name__ == "__main__":
    sys.exit(main())
