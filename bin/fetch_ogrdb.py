#!/usr/bin/env python3
"""Download the pinned OGRDB human IGH germline set and split it by segment.

reference.lock names one file of one published version (its Zenodo record),
the MD5 that Zenodo publishes for it and the SHA-256 written at the first
download. Every download is checked against both; a file that does not
match is never used.

    fetch_ogrdb.py reference.lock OUT_DIR
        download, check, and write OUT_DIR/IGH{V,D,J}.fasta
    fetch_ogrdb.py --pin reference.lock OUT_DIR
        the same, and if the SHA-256 is not yet pinned, write it into
        reference.lock (only after the published MD5 matches). A pinned
        SHA-256 is never changed: a file that differs from it stops here.

The OGRDB set is one gapped FASTA holding V, D and J. V keeps its IMGT gap
dots (MakeDb and CreateGermlines need them); D and J are written without
dots. A record that is not IGHV, IGHD or IGHJ, or a segment with no record,
stops the split.

A download that fails with HTTP 5xx, a connection error or a timeout is tried
again, five attempts in all, after 15, 30, 60 and 120 s, with one line on
standard error per retry. HTTP 4xx is not retried. The digest checks are the
same on every attempt and still stop the run.
"""
import hashlib
import os
import socket
import sys
import time
import urllib.error
import urllib.request

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from igblast_fasta import allele_name, read_fasta

KEY = "ogrdb_igh_vdj"
NOT_PINNED = "sha256:not-yet-pinned"
SEGMENTS = ("IGHV", "IGHD", "IGHJ")
BACKOFF_S = (15, 30, 60, 120)  # waits between the five attempts


def read_lock(path):
    """The lock row for KEY as a dict: version, url, md5, sha256."""
    with open(path, encoding="utf-8") as fh:
        for line in fh:
            if line.startswith("#") or not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if fields[0] == KEY:
                if len(fields) != 5:
                    raise ValueError(f"{path}: {KEY} needs 5 tab-separated fields, has {len(fields)}")
                return dict(zip(("version", "url", "md5", "sha256"), fields[1:]))
    raise ValueError(f"{path}: no {KEY} row")


def write_sha256(path, sha256):
    """Fill the not-yet-pinned SHA-256 of the KEY row; refuse to change a pinned one."""
    with open(path, encoding="utf-8") as fh:
        lines = fh.readlines()
    out, done = [], False
    for line in lines:
        fields = line.rstrip("\n").split("\t")
        if fields[0] == KEY and not line.startswith("#"):
            if fields[4] != NOT_PINNED:
                raise ValueError(f"{KEY} is already pinned to {fields[4]}; not changed")
            fields[4] = "sha256:" + sha256
            out.append("\t".join(fields) + "\n")
            done = True
        else:
            out.append(line)
    if not done:
        raise ValueError(f"{path}: no {KEY} row")
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(out)


def check_digests(data, row):
    """Raise unless data matches the published MD5 and, when pinned, the SHA-256."""
    md5 = "md5:" + hashlib.md5(data).hexdigest()  # the digest Zenodo publishes
    if md5 != row["md5"]:
        raise ValueError(f"downloaded file has {md5}, the published value is {row['md5']}")
    sha256 = hashlib.sha256(data).hexdigest()
    if row["sha256"] != NOT_PINNED and "sha256:" + sha256 != row["sha256"]:
        raise ValueError(f"downloaded file has sha256:{sha256}, reference.lock pins {row['sha256']}")
    return sha256


def download(url, backoff=BACKOFF_S, sleep=time.sleep, timeout=120):
    """The body at url. HTTP 5xx, connection errors and timeouts are retried after each wait in backoff; HTTP 4xx
    and the last failure raise ValueError."""
    for attempt in range(len(backoff) + 1):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:  # fixed URL from the lock
                return resp.read()
        except urllib.error.HTTPError as err:
            if err.code < 500:
                raise ValueError(f"download got HTTP {err.code} {err.reason}; not retried") from err
            why = f"HTTP {err.code} {err.reason}"
        except (urllib.error.URLError, TimeoutError, socket.timeout, ConnectionError) as err:
            why = f"{type(err).__name__}: {getattr(err, 'reason', err)}"
        if attempt == len(backoff):
            raise ValueError(f"download failed after {attempt + 1} attempts; last: {why}")
        print(f"fetch_ogrdb.py: attempt {attempt + 1} of {len(backoff) + 1} failed ({why}); "
              f"retrying in {backoff[attempt]} s", file=sys.stderr, flush=True)
        sleep(backoff[attempt])


def split(records):
    """[(header, seq)] -> {segment: [(allele, seq)]}; V gapped, D and J without dots."""
    parts = {s: [] for s in SEGMENTS}
    for header, seq in records:
        name = allele_name(header)
        seg = name[:4]
        if seg not in parts:
            raise ValueError(f"record is not IGHV, IGHD or IGHJ: {header!r}")
        parts[seg].append((name, seq if seg == "IGHV" else seq.replace(".", "")))
    empty = [s for s in SEGMENTS if not parts[s]]
    if empty:
        raise ValueError(f"no records for {', '.join(empty)}")
    return parts


def write_fasta(path, records):
    with open(path, "w", encoding="utf-8") as fh:
        fh.writelines(f">{name}\n{seq}\n" for name, seq in records)


def main(argv):
    pin = argv[:1] == ["--pin"]
    args = argv[1:] if pin else argv
    if len(args) != 2:
        print(__doc__, file=sys.stderr)
        return 2
    lock, out_dir = args
    row = read_lock(lock)
    if row["sha256"] == NOT_PINNED and not pin:
        print(f"{KEY}: SHA-256 not yet pinned; run 'make pin-reference' first", file=sys.stderr)
        return 1
    data = download(row["url"])
    sha256 = check_digests(data, row)
    if pin and row["sha256"] == NOT_PINNED:
        write_sha256(lock, sha256)
        print(f"{KEY}: pinned sha256:{sha256}")
    os.makedirs(out_dir, exist_ok=True)
    source = os.path.join(out_dir, os.path.basename(row["url"].split("?")[0]))
    with open(source, "wb") as fh:
        fh.write(data)
    parts = split(read_fasta(source))
    for seg in SEGMENTS:
        write_fasta(os.path.join(out_dir, f"{seg}.fasta"), parts[seg])
    print(f"OGRDB IGH_VDJ version {row['version']}: "
          + ", ".join(f"{len(parts[s])} {s}" for s in SEGMENTS) + f" in {out_dir}")
    return 0


if __name__ == "__main__":
    try:
        sys.exit(main(sys.argv[1:]))
    except ValueError as err:
        print(f"fetch_ogrdb.py: {err}", file=sys.stderr)
        sys.exit(1)
