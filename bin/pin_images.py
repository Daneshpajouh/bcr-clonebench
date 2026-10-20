#!/usr/bin/env python3
"""Write image digests into images.lock after 'docker pull' or 'docker build'.

    pin_images.py images.lock NAME TAG [NAME TAG ...]

For a pulled image the registry digest (name@sha256:...) is written; a
locally built image has none, so its image ID (sha256:...) is written,
which pins the exact build on that machine.
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
    out = []
    for line in lines:
        parts = line.split("\t")
        if not line.startswith("#") and len(parts) == 3 and parts[0] in pairs:
            parts = [parts[0], pairs[parts[0]], digest(pairs[parts[0]])]
            print("\t".join(parts))
        out.append("\t".join(parts))
    with open(lock, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
