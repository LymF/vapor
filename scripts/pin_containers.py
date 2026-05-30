#!/usr/bin/env python3
"""
pin_containers.py — resolve quay.io/biocontainers tags and generate containers.lock.yaml

Usage:
    python3 scripts/pin_containers.py               # generate / update lock file
    python3 scripts/pin_containers.py --check       # verify lock against containers.yaml
    python3 scripts/pin_containers.py --key fastp   # re-pin a single tool
    python3 scripts/pin_containers.py --dry-run     # print resolved URIs without writing

Reads:  containers.yaml  (tool key → image name + version)
Writes: containers.lock.yaml  (tool key → full docker:// URI)

Every bioconda package is automatically containerised by BioContainers and
pushed to quay.io/biocontainers/<name>:<version>--<build_hash>.
This script queries the Quay.io API to find the matching tag for each version
defined in containers.yaml and writes the full URI to containers.lock.yaml.
"""

import argparse
import json
import os
import sys
import time
import urllib.request
import urllib.error
from pathlib import Path

import yaml  # PyYAML — available in the snakemake environment


QUAY_API = "https://quay.io/api/v1/repository/biocontainers/{image}/tag/?onlyActiveTags=true&limit=100&page={page}"
REGISTRY  = "docker://quay.io/biocontainers"

ROOT = Path(__file__).resolve().parent.parent
CONTAINERS_YAML = ROOT / "containers.yaml"
LOCK_YAML       = ROOT / "containers.lock.yaml"


# ── helpers ───────────────────────────────────────────────────────────────


def fetch_tags(image: str) -> list[str]:
    """Return all active tag names for a BioContainers image."""
    tags = []
    page = 1
    while True:
        url = QUAY_API.format(image=image, page=page)
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "VAPOR/1.0"})
            with urllib.request.urlopen(req, timeout=20) as resp:
                data = json.loads(resp.read())
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return []
            raise
        batch = [t["name"] for t in data.get("tags", [])]
        tags.extend(batch)
        if not data.get("has_additional") or not batch:
            break
        page += 1
        time.sleep(0.2)
    return tags


def best_tag(tags: list[str], version: str) -> str | None:
    """
    Pick the best tag matching `version`.
    Prefers tags that start with exactly `<version>--` (the standard BioContainers
    pattern), then falls back to any tag that starts with the version string.
    Within each group, prefers `pyhdfd78af` (noarch/pure-Python) builds, then
    sorts lexicographically and returns the last entry (typically highest hash).
    """
    exact  = [t for t in tags if t.startswith(f"{version}--")]
    prefix = [t for t in tags if t.startswith(version)] if not exact else []
    pool   = exact or prefix
    if not pool:
        return None
    # prefer noarch pure-Python builds
    noarch = [t for t in pool if "pyhdfd78af" in t]
    return sorted(noarch or pool)[-1]


def resolve(key: str, cfg: dict) -> str | None:
    """Return full docker:// URI for one entry, or None on failure."""
    # Custom images (not from quay.io/biocontainers) pass through as-is.
    if cfg.get("custom"):
        return f"docker://{cfg['image']}:{cfg['version']}"

    image   = cfg["image"]
    version = str(cfg["version"])
    print(f"  [{key}] quay.io/biocontainers/{image}  version={version} ... ", end="", flush=True)

    try:
        tags = fetch_tags(image)
    except Exception as exc:
        print(f"ERROR ({exc})")
        return None

    tag = best_tag(tags, version)
    if tag:
        uri = f"{REGISTRY}/{image}:{tag}"
        print(f"→ {tag}")
        return f"docker://{uri.removeprefix('docker://')}"
    else:
        print(f"NOT FOUND (available: {tags[:5]}{'…' if len(tags) > 5 else ''})")
        return None


# ── main ──────────────────────────────────────────────────────────────────


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--check",   action="store_true", help="Verify lock file without writing")
    ap.add_argument("--dry-run", action="store_true", help="Print URIs without writing lock file")
    ap.add_argument("--key",     metavar="KEY", help="Re-pin only this tool key")
    args = ap.parse_args()

    if not CONTAINERS_YAML.exists():
        sys.exit(f"ERROR: {CONTAINERS_YAML} not found")

    with open(CONTAINERS_YAML) as f:
        spec = yaml.safe_load(f)

    # Load existing lock (may be partial)
    lock: dict = {}
    if LOCK_YAML.exists():
        with open(LOCK_YAML) as f:
            lock = yaml.safe_load(f) or {}

    keys = [args.key] if args.key else list(spec.keys())
    print(f"Resolving {len(keys)} container(s)…\n")

    errors = []
    for key in keys:
        if key not in spec:
            print(f"  WARNING: '{key}' not found in containers.yaml — skipped")
            continue
        cfg = spec[key]

        if args.check:
            existing = lock.get(key)
            status = "OK" if existing else "MISSING"
            print(f"  [{key}]  {status}  {existing or ''}")
            if not existing:
                errors.append(key)
            continue

        uri = resolve(key, cfg)
        if uri:
            lock[key] = uri
        else:
            errors.append(key)
            # keep old value if present
            if key in lock:
                print(f"    (keeping previous: {lock[key]})")

    print()

    if args.check:
        if errors:
            print(f"MISSING in lock: {errors}")
            sys.exit(1)
        print("All containers present in lock file.")
        return

    if args.dry_run:
        print("Resolved URIs (not written):")
        for k, v in lock.items():
            print(f"  {k}: {v}")
        return

    # Write lock file
    header = (
        "# containers.lock.yaml — AUTO-GENERATED by scripts/pin_containers.py\n"
        "# Do not edit manually. Re-run pin_containers.py to update.\n"
        "# Commit this file to pin container versions for reproducibility.\n\n"
    )
    with open(LOCK_YAML, "w") as f:
        f.write(header)
        yaml.dump(lock, f, default_flow_style=False, sort_keys=True)

    print(f"Written: {LOCK_YAML}")
    if errors:
        print(f"\nWARNING: {len(errors)} tool(s) could not be resolved: {errors}")
        print("Rules for those tools will fall back to --use-conda.")
        sys.exit(2)


if __name__ == "__main__":
    main()
