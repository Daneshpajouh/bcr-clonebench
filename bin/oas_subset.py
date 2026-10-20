#!/usr/bin/env python3
"""Subset OAS unpaired data units for the real-data run, by a fixed rule.

Rule, fixed before any real file is read:
1. Subject: of the subjects whose units cover every time point asked for, the
   first in plain sort order of the OAS "Subject" field (or the one named with
   --subject). Isotype units of one subject and time point are pooled.
2. Rows kept: productive, with a V call, a J call and a junction; then one row
   per distinct nucleotide sequence within a time point.
3. Cap: when a time point still has more than --cap rows, the --cap rows with
   the smallest SHA-256 of (salt, sequence) are kept. The choice does not
   depend on file or row order, and the same salt and cap give the same rows.
4. Every count before and after each step, the inputs' SHA-256 and the rule's
   settings go to a manifest.

OAS unit files are gzip CSV whose first line is the unit's metadata (a JSON
object, written as one quoted CSV field) and whose second line is the AIRR
header. [CONFIRM the metadata key names and the productive flag values on the
first real unit; the script stops with a clear message if they differ.]

Standard library only.
"""
import argparse
import csv
import gzip
import hashlib
import heapq
import json
import sys
from pathlib import Path

DAYS = ("Day-0", "Day-7", "Day-14", "Day-21", "Day-28")
TRUE = {"T", "True", "TRUE", "true"}
NEEDED = ("sequence", "productive", "v_call", "j_call", "junction")
csv.field_size_limit(1 << 24)


def read_meta(path):
    """The unit's metadata from its first line."""
    with gzip.open(path, "rt", newline="") as fh:
        first = fh.readline()
    raw = next(csv.reader([first]))[0]
    meta = json.loads(raw)
    for key in ("Subject", "Longitudinal"):
        if key not in meta:
            sys.exit(f"STOP: {path.name}: metadata has no {key!r} (keys: {sorted(meta)})")
    return meta


def rows(path):
    with gzip.open(path, "rt", newline="") as fh:
        fh.readline()
        reader = csv.DictReader(fh)
        missing = [c for c in NEEDED if c not in (reader.fieldnames or [])]
        if missing:
            sys.exit(f"STOP: {path.name}: columns missing: {missing}")
        yield from reader


def sha256_file(path):
    h = hashlib.sha256()
    with open(path, "rb") as fh:
        for block in iter(lambda: fh.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def pick_subject(units, days, subject=None):
    cover = {}
    for path, meta in units:
        cover.setdefault(str(meta["Subject"]), set()).add(meta["Longitudinal"])
    full = sorted(s for s, d in cover.items() if set(days) <= d)
    if subject is not None:
        if subject not in full:
            sys.exit(f"STOP: subject {subject!r} does not cover {list(days)}")
        return subject, full
    if not full:
        sys.exit(f"STOP: no subject covers {list(days)}")
    return full[0], full


def subset_day(paths, cap, salt):
    """Filter, deduplicate and cap one time point; returns (kept rows, counts)."""
    counts = {"rows_in": 0, "not_productive": 0, "no_v_j_or_junction": 0, "duplicate_sequence": 0}
    seen, heap = set(), []  # heap of (-key, sequence, row): the cap smallest keys survive
    for path in paths:
        for row in rows(path):
            counts["rows_in"] += 1
            if row["productive"] not in TRUE:
                counts["not_productive"] += 1
                continue
            if not (row["v_call"] and row["j_call"] and row["junction"]):
                counts["no_v_j_or_junction"] += 1
                continue
            seq = row["sequence"]
            if seq in seen:
                counts["duplicate_sequence"] += 1
                continue
            seen.add(seq)
            key = hashlib.sha256(f"{salt}\0{seq}".encode()).hexdigest()
            item = (_neg(key), seq, row)
            if len(heap) < cap:
                heapq.heappush(heap, item)
            elif item > heap[0]:
                heapq.heapreplace(heap, item)
    counts["distinct_kept_before_cap"] = len(seen)
    kept = [row for _, _, row in sorted(heap, key=lambda t: t[1])]
    counts["rows_out"] = len(kept)
    counts["capped"] = len(seen) > cap
    return kept, counts


def _neg(hexkey):
    """Order key so that the heap's smallest item has the largest hash: popping it keeps the smallest hashes."""
    return -int(hexkey, 16)


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("units", nargs="+", type=Path, help="OAS unit files (.csv.gz)")
    ap.add_argument("--out", required=True, type=Path, help="new folder for the per-day TSVs and manifest.json")
    ap.add_argument("--cap", type=int, default=20000, help="rows kept per time point at most [CONFIRM by the first run's time]")
    ap.add_argument("--salt", default="bcr-clonebench-galson2015a", help="hash salt for the cap")
    ap.add_argument("--subject", help="OAS Subject to use instead of the first full one")
    ap.add_argument("--days", nargs="+", default=list(DAYS))
    a = ap.parse_args(argv)
    if a.cap < 1:
        sys.exit("STOP: --cap must be at least 1")
    units = [(p, read_meta(p)) for p in sorted(a.units)]
    subject, full = pick_subject(units, a.days, a.subject)
    a.out.mkdir(parents=False, exist_ok=False)
    manifest = {"rule": {"cap": a.cap, "salt": a.salt, "days": a.days, "subject_rule":
                         "named" if a.subject else "first full subject in sort order"},
                "subjects_covering_all_days": full, "subject": subject, "days": {}}
    for day in a.days:
        paths = [p for p, m in units if str(m["Subject"]) == subject and m["Longitudinal"] == day]
        kept, counts = subset_day(paths, a.cap, a.salt)
        out = a.out / f"{subject}_{day}.tsv"
        with open(out, "x", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(kept[0].keys()) if kept else list(NEEDED),
                               delimiter="\t", extrasaction="ignore")
            w.writeheader()
            w.writerows(kept)
        counts["inputs"] = {p.name: sha256_file(p) for p in paths}
        counts["output"] = {out.name: sha256_file(out)}
        manifest["days"][day] = counts
    with open(a.out / "manifest.json", "x") as fh:
        json.dump(manifest, fh, indent=1, sort_keys=True)
    print(json.dumps({"subject": subject, "rows_out": {d: c["rows_out"] for d, c in manifest["days"].items()}}))


if __name__ == "__main__":
    main()
