#!/usr/bin/env python3
"""One table from the score files of every level and method.

    summarize_scores.py --out-tsv results.tsv --out-md results.md score_*.json

Rows are sorted by mutation level (low, medium, high) and then method.
Undefined values (no predicted pairs) are written as NA, never as 0.
"""

import argparse
import json
import sys

LEVEL_ORDER = {"low": 0, "medium": 1, "high": 2}
COLUMNS = ("level", "method", "n_truth", "n_dropped", "pair_precision", "pair_recall", "pair_f1", "ari",
           "v_gene_exact", "j_gene_exact", "mu_mae", "clones_true", "clones_predicted")


def row(score):
    level, _, method = score["label"].partition("/")
    a = score["all"]
    return {"level": level, "method": method, "n_truth": score["n_truth"], "n_dropped": score["n_dropped"],
            "pair_precision": a["pair_precision"], "pair_recall": a["pair_recall"], "pair_f1": a["pair_f1"],
            "ari": a["ari"], "v_gene_exact": score["v_gene"]["exact"], "j_gene_exact": score["j_gene"]["exact"],
            "mu_mae": (score.get("mu_freq") or {}).get("mean_abs_error"),
            "clones_true": score["hill"]["q0"]["true"], "clones_predicted": score["hill"]["q0"]["predicted"]}


def fmt(v):
    if v is None:
        return "NA"
    if isinstance(v, float):
        return f"{v:.4f}"
    return str(v)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--out-tsv", required=True)
    ap.add_argument("--out-md", required=True)
    ap.add_argument("scores", nargs="+")
    a = ap.parse_args(argv)
    rows = []
    for p in a.scores:
        with open(p, encoding="utf-8") as f:
            rows.append(row(json.load(f)))
    rows.sort(key=lambda r: (LEVEL_ORDER.get(r["level"], 9), r["method"]))
    with open(a.out_tsv, "w", encoding="utf-8") as f:
        f.write("\t".join(COLUMNS) + "\n")
        for r in rows:
            f.write("\t".join(fmt(r[c]) for c in COLUMNS) + "\n")
    with open(a.out_md, "w", encoding="utf-8") as f:
        f.write("| " + " | ".join(COLUMNS) + " |\n|" + "---|" * len(COLUMNS) + "\n")
        for r in rows:
            f.write("| " + " | ".join(fmt(r[c]) for c in COLUMNS) + " |\n")
    print(f"{len(rows)} rows")
    return 0


if __name__ == "__main__":
    sys.exit(main())
