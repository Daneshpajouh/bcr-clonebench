#!/usr/bin/env python3
"""Score one pipeline result against the simulation's known answers.

    score_clones.py --truth truth.tsv --pred clones.tsv [--mu-column mu_freq]
                    [--label high] --out-json score.json

truth.tsv (written by the simulation step) has one row per simulated
read: sequence_id, true_clone_id, true_v_call, true_j_call and
true_mu_freq (the fraction of V positions mutated away from the true
germline). clones.tsv is the pipeline's AIRR rearrangement table with
sequence_id, clone_id, v_call and j_call, and the column named by
--mu-column for the hypermutation estimate.

Every truth sequence is scored. A sequence the pipeline dropped (failed
QC or annotation) counts as its own one-sequence clone and as a wrong V
and J call, so dropping hard reads cannot raise a score. The same
metrics on the retained sequences only are reported beside them, with
the V and J calls of the kept reads (v_gene_kept, j_gene_kept).

When the truth table carries true_productive (the simulator's own flag),
productive_truth repeats the clone and gene scores over the reads marked
productive (T), beside the headline over every read.

With --drops (the annotation step's drops table: sequence_id, reason),
n_dropped is split by reason; a truth read dropped after annotation, or
not named in the table, counts as dropped_later.

With --out-reads, one row per truth read is written as well (READ_COLUMNS):
its true and called clone, V and J, whether each call is exact and whether
the call set includes the true gene, and whether the read was kept or why
it was dropped. The rows add up to the summary scores. A --label of the
form level/method fills the level and method columns.

Standard library only, so CI needs nothing installed.
"""

import argparse
import csv
import json
import math
import sys

TRUTH_COLUMNS = ("sequence_id", "true_clone_id", "true_v_call", "true_j_call", "true_mu_freq")
PRED_COLUMNS = ("sequence_id", "clone_id", "v_call", "j_call")
READ_COLUMNS = ("sequence_id", "level", "method", "true_productive", "true_clone_id", "clone_id",
                "true_v_call", "v_call", "v_exact", "v_includes_true",
                "true_j_call", "j_call", "j_exact", "j_includes_true", "status", "drop_reason")


def read_tsv(path, required, allow_empty=False):
    """Rows of a tab-separated table with a header. A table with no header
    is always refused; a header with no rows only when allow_empty is
    false (a drops table with no rows is real: no read was dropped)."""
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f, delimiter="\t")
        header = reader.fieldnames or []
        rows = list(reader)
    if not header:
        raise ValueError(f"{path}: no header")
    missing = [c for c in required if c not in header]
    if missing:
        raise ValueError(f"{path}: missing column(s) {', '.join(missing)}")
    if not rows and not allow_empty:
        raise ValueError(f"{path}: no rows")
    ids = [r["sequence_id"] for r in rows]
    if len(set(ids)) != len(ids):
        raise ValueError(f"{path}: repeated sequence_id")
    return rows


def pairs(n):
    return n * (n - 1) // 2


def pair_scores(true_labels, pred_labels):
    """Pairwise precision, recall and F1, and the adjusted Rand index
    (Hubert and Arabie 1985), from the contingency table, so the cost is
    linear in the number of sequences."""
    cells, a, b = {}, {}, {}
    for t, p in zip(true_labels, pred_labels):
        cells[(t, p)] = cells.get((t, p), 0) + 1
        a[t] = a.get(t, 0) + 1
        b[p] = b.get(p, 0) + 1
    tp = sum(pairs(n) for n in cells.values())
    true_pairs = sum(pairs(n) for n in a.values())
    pred_pairs = sum(pairs(n) for n in b.values())
    all_pairs = pairs(len(true_labels))
    precision = tp / pred_pairs if pred_pairs else None
    recall = tp / true_pairs if true_pairs else None
    f1 = (2 * precision * recall / (precision + recall)
          if precision is not None and recall is not None and precision + recall > 0 else None)
    expected = true_pairs * pred_pairs / all_pairs if all_pairs else 0.0
    top = (true_pairs + pred_pairs) / 2 - expected
    # both partitions all singletons or one cluster: ARI is defined as 1 when they agree
    ari = (tp - expected) / top if top else (1.0 if true_pairs == pred_pairs == tp else 0.0)
    return {"pair_precision": precision, "pair_recall": recall, "pair_f1": f1, "ari": ari,
            "true_pairs": true_pairs, "predicted_pairs": pred_pairs, "shared_pairs": tp}


def gene(allele):
    return allele.split("*", 1)[0].strip()


def call_match(truth_call, pred_call):
    """(exact, includes_true) for one read, each 0 or 1: the call is exactly
    the true gene (one gene named); the call set includes it (ties allowed)."""
    genes = {gene(x) for x in (pred_call or "").split(",") if x.strip()}
    return int(genes == {gene(truth_call)}), int(gene(truth_call) in genes)


def call_scores(truth_calls, pred_calls):
    """Share of sequences whose call is exactly the true gene (one gene
    named), and share whose call set includes it (ties allowed)."""
    exact = anyof = 0
    for t, p in zip(truth_calls, pred_calls):
        e, a = call_match(t, p)
        exact += e
        anyof += a
    n = len(truth_calls)
    return {"exact": exact / n, "includes_true": anyof / n}


def hill(sizes, q):
    total = sum(sizes)
    p = [s / total for s in sizes if s]
    if q == 0:
        return float(len(p))
    if q == 1:
        return math.exp(-sum(x * math.log(x) for x in p))
    return sum(x ** q for x in p) ** (1 / (1 - q))


def clone_sizes(labels):
    sizes = {}
    for x in labels:
        sizes[x] = sizes.get(x, 0) + 1
    return list(sizes.values())


def score(truth, pred, mu_column=None, drops=None):
    by_id = {r["sequence_id"]: r for r in pred}
    unknown = [i for i in by_id if i not in {r["sequence_id"] for r in truth}]
    if unknown:
        raise ValueError(f"{len(unknown)} predicted sequence_id(s) not in the truth, e.g. {unknown[0]}")
    t_lab, p_lab, t_v, p_v, t_j, p_j, mu_err = [], [], [], [], [], [], []
    kept = []
    for r in truth:
        sid = r["sequence_id"]
        p = by_id.get(sid)
        t_lab.append(r["true_clone_id"])
        t_v.append(r["true_v_call"])
        t_j.append(r["true_j_call"])
        if p is None or not p.get("clone_id"):
            p_lab.append(("dropped", sid))  # its own clone
            p_v.append("")
            p_j.append("")
            continue
        kept.append(len(t_lab) - 1)
        p_lab.append(("clone", p["clone_id"]))
        p_v.append(p["v_call"])
        p_j.append(p["j_call"])
        if mu_column and p.get(mu_column) not in (None, "", "NA"):
            mu_err.append(float(p[mu_column]) - float(r["true_mu_freq"]))
    out = {"n_truth": len(truth), "n_scored": len(kept), "n_dropped": len(truth) - len(kept)}
    out["all"] = pair_scores(t_lab, p_lab)
    out["retained"] = pair_scores([t_lab[i] for i in kept], [p_lab[i] for i in kept]) if kept else None
    out["v_gene"] = call_scores(t_v, p_v)
    out["j_gene"] = call_scores(t_j, p_j)
    out["v_gene_kept"] = call_scores([t_v[i] for i in kept], [p_v[i] for i in kept]) if kept else None
    out["j_gene_kept"] = call_scores([t_j[i] for i in kept], [p_j[i] for i in kept]) if kept else None
    if drops is not None:
        kept_ids = {truth[i]["sequence_id"] for i in kept}
        by_reason = {}
        for r in truth:
            if r["sequence_id"] not in kept_ids:
                reason = drops.get(r["sequence_id"], "dropped_later")
                by_reason[reason] = by_reason.get(reason, 0) + 1
        out["dropped_by_reason"] = by_reason
    if mu_column:
        out["mu_freq"] = ({"n": len(mu_err), "mean_abs_error": sum(abs(e) for e in mu_err) / len(mu_err),
                           "mean_error": sum(mu_err) / len(mu_err)} if mu_err else {"n": 0})
    ts = clone_sizes(t_lab)
    ps = clone_sizes([p_lab[i] for i in kept])
    out["hill"] = {f"q{q}": {"true": hill(ts, q), "predicted": hill(ps, q) if ps else None} for q in (0, 1, 2)}
    return out


def productive_view(truth, pred, mu_column=None, drops=None):
    """The same scores over the truth reads the simulator marks productive (true_productive T). A productive read
    the pipeline dropped still counts as its own clone, so dropping reads cannot raise this score either. None when
    the truth table carries no productive flag."""
    flags = [(r.get("true_productive") or "").strip() for r in truth]
    if not any(flags):
        return None
    keep = {r["sequence_id"] for r, f in zip(truth, flags) if f == "T"}
    view = {"n_truth": len(keep), "n_not_productive": flags.count("F"), "n_unknown": flags.count("")}
    if not keep:
        return {**view, "n_scored": 0, "n_dropped": 0, "all": None, "v_gene": None, "j_gene": None}
    part = score([r for r in truth if r["sequence_id"] in keep],
                 [p for p in pred if p["sequence_id"] in keep], mu_column, drops)
    view.update({k: part[k] for k in ("n_scored", "n_dropped", "all", "v_gene", "j_gene")})
    if "dropped_by_reason" in part:
        view["dropped_by_reason"] = part["dropped_by_reason"]
    return view


def per_read(truth, pred, drops, level="", method=""):
    """One row per truth read, in truth order (READ_COLUMNS). A read is kept
    when the pipeline gave it a clone, as in score(); a dropped read has
    empty calls, so its matches are 0, and its reason comes from the drops
    table (dropped_later when not named there, 'dropped' with no table)."""
    by_id = {r["sequence_id"]: r for r in pred}
    rows = []
    for r in truth:
        sid = r["sequence_id"]
        p = by_id.get(sid)
        kept = p is not None and bool(p.get("clone_id"))
        v = p["v_call"] if kept else ""
        j = p["j_call"] if kept else ""
        ve, vi = call_match(r["true_v_call"], v)
        je, ji = call_match(r["true_j_call"], j)
        reason = "" if kept else ("dropped" if drops is None else drops.get(sid, "dropped_later"))
        rows.append({"sequence_id": sid, "level": level, "method": method,
                     "true_productive": (r.get("true_productive") or "").strip(),
                     "true_clone_id": r["true_clone_id"], "clone_id": p["clone_id"] if kept else "",
                     "true_v_call": r["true_v_call"], "v_call": v, "v_exact": ve, "v_includes_true": vi,
                     "true_j_call": r["true_j_call"], "j_call": j, "j_exact": je, "j_includes_true": ji,
                     "status": "kept" if kept else "dropped", "drop_reason": reason})
    return rows


def write_reads(path, rows):
    with open(path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=READ_COLUMNS, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--truth", required=True)
    ap.add_argument("--pred", required=True)
    ap.add_argument("--mu-column", help="column of the hypermutation estimate in --pred")
    ap.add_argument("--label", default="", help="written into the output, e.g. the mutation level")
    ap.add_argument("--drops", help="the annotation step's drops table (sequence_id, reason)")
    ap.add_argument("--out-json", required=True)
    ap.add_argument("--out-reads", help="also write one row per truth read here (tab-separated)")
    a = ap.parse_args(argv)
    try:
        drops = ({r["sequence_id"]: r["reason"]
                  for r in read_tsv(a.drops, ("sequence_id", "reason"), allow_empty=True)}
                 if a.drops else None)
        truth = read_tsv(a.truth, TRUTH_COLUMNS)
        pred = read_tsv(a.pred, PRED_COLUMNS + ((a.mu_column,) if a.mu_column else ()))
        res = score(truth, pred, a.mu_column, drops)
        res["productive_truth"] = productive_view(truth, pred, a.mu_column, drops)
        if a.out_reads:
            level, _, method = a.label.partition("/")
            write_reads(a.out_reads, per_read(truth, pred, drops, level, method))
    except (OSError, ValueError) as e:
        print(f"score_clones: {e}", file=sys.stderr)
        return 2
    res["label"] = a.label
    with open(a.out_json, "w", encoding="utf-8") as f:
        json.dump(res, f, indent=1, sort_keys=True)
        f.write("\n")
    s = res["all"]
    print(f"{a.label or 'result'}: {res['n_scored']}/{res['n_truth']} scored, pair F1 {s['pair_f1']}, ARI {s['ari']:.4f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
