#!/usr/bin/env python3
"""One table from the score files of every level and method.

    summarize_scores.py --out-tsv results.tsv --out-md results.md score_*.json

Rows are sorted by mutation level (low, medium, high), then sequencing
error rate, then method.
Undefined values (no predicted pairs) are written as NA, never as 0.
"""

import argparse
import json
import sys

LEVEL_ORDER = {"low": 0, "medium": 1, "high": 2, "low_shared": 3, "medium_shared": 4, "high_shared": 5}
COLUMNS = ("level", "error_rate", "method", "n_truth", "n_dropped", "n_not_annotated", "n_not_productive",
           "n_dropped_later", "pair_precision", "predicted_pairs", "pair_recall", "pair_f1", "ari", "n_truth_productive",
           "pair_f1_productive", "ari_productive", "v_gene_exact",
           "v_gene_exact_kept", "j_gene_exact", "j_gene_exact_kept", "mu_mae", "clones_true", "clones_predicted")


def row(score):
    run, _, method = score["label"].partition("/")
    level, _, rate = run.partition("_e")
    a = score["all"]
    why = score.get("dropped_by_reason")
    kept = lambda k: (score.get(k) or {}).get("exact")  # noqa: E731
    prod = score.get("productive_truth") or {}
    prod_all = prod.get("all") or {}
    return {"n_not_annotated": why.get("not_annotated", 0) if why is not None else None,
            "n_not_productive": why.get("not_productive", 0) if why is not None else None,
            "n_dropped_later": (why.get("dropped_later", 0) + why.get("other", 0)) if why is not None else None,
            "predicted_pairs": a.get("predicted_pairs"), "v_gene_exact_kept": kept("v_gene_kept"),
            "j_gene_exact_kept": kept("j_gene_kept"),"level": level, "error_rate": float(rate) if rate else None, "method": method,
            "n_truth": score["n_truth"], "n_dropped": score["n_dropped"],
            "pair_precision": a["pair_precision"], "pair_recall": a["pair_recall"], "pair_f1": a["pair_f1"],
            "ari": a["ari"], "n_truth_productive": prod.get("n_truth"), "pair_f1_productive": prod_all.get("pair_f1"),
            "ari_productive": prod_all.get("ari"), "v_gene_exact": score["v_gene"]["exact"], "j_gene_exact": score["j_gene"]["exact"],
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
    rows.sort(key=lambda r: (LEVEL_ORDER.get(r["level"], 9), r["error_rate"] or 0.0, r["method"]))
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
