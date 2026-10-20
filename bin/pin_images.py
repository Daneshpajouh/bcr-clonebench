#!/usr/bin/env python3
"""Write image digests into images.lock after 'docker pull' or 'docker build'.

    pin_images.py images.lock NAME TAG [NAME TAG ...]

For a pulled image the registry digest (name@sha256:...) is written; a
locally built image has none, so its image ID (sha256:...) is written,
which pins the exact build on that machine.

A row that already holds a registry digest is a lock: if the same tag now
resolves to any other digest, the tag has moved and nothing is written. A
different tag, set on purpose, is pinned afresh.
"""

import json
import subprocess
import sys


def digest(tag):
    out = subprocess.run(["docker", "image", "inspect", "--format", "{{json .RepoDigests}} {{.Id}}", tag],
                         capture_output=True, text=True, check=True).stdout.strip()
    repo_digests, image_id = out.rsplit(" ", 1)
    found = json.loads(repo_digests) or []
    return found[0] if found else image_id


def main(argv):
    if len(argv) < 3 or len(argv) % 2 == 0:
        print(__doc__.split("\n\n")[1], file=sys.stderr)
        return 2
    lock, pairs = argv[0], dict(zip(argv[1::2], argv[2::2]))
    with open(lock, encoding="utf-8") as f:
        lines = f.read().splitlines()
    out, moved = [], []
    for line in lines:
        parts = line.split("\t")
        if not line.startswith("#") and len(parts) == 3 and parts[0] in pairs:
            found = digest(pairs[parts[0]])
            if "@sha256:" in parts[2] and parts[1] == pairs[parts[0]] and found != parts[2]:
                moved.append(f"{parts[0]}: {pairs[parts[0]]} is locked to {parts[2]} but now resolves to {found}")
            parts = [parts[0], pairs[parts[0]], found]
            print("\t".join(parts))
        out.append("\t".join(parts))
    if moved:
        for m in moved:
            print(f"refused, tag moved: {m}", file=sys.stderr)
        print(f"{lock} not written; check the new image, then update the lock by hand", file=sys.stderr)
        return 1
    with open(lock, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
