#!/usr/bin/env python3
"""
VAPOR — Viral And Prokaryotic genOme Recovery
CLI wrapper for the Snakemake metagenomics pipeline.

Usage:
  vapor [options]

Run 'vapor -h' for full help.
"""

import argparse
import subprocess
import sys
from pathlib import Path

__version__ = "1.1.0"

DESCRIPTION = """\
VAPOR — Viral And Prokaryotic genOme Recovery
A comprehensive pipeline for viral and prokaryotic genome recovery
from metagenomic short-read (Illumina) and long-read (ONT/HiFi) data.
"""

EPILOG = """\
Execution modes:
  vapor --threads 32             containers (default) — uses Apptainer/Singularity
  vapor --threads 32 --use-conda conda environments   — dev/local mode
  vapor --create-envs            create all conda environments without running

Examples:
  vapor --dry-run                        validate workflow without executing
  vapor --threads 32                     run with 32 cores (container mode)
  vapor --threads 32 --use-conda         run with conda environments
  vapor --config /path/to/config.yaml    use a custom config file
  vapor --rerun-incomplete               resume an interrupted run
  vapor --force                          re-run everything from scratch
  vapor --forcerun quast viral_detection force specific rules to re-run
  vapor --dag                            generate workflow DAG (dag.svg)
  vapor --unlock                         unlock directory after crash
  vapor --create-envs                    create all conda environments

Container setup (one-time, run from the VAPOR directory):
  python3 scripts/pin_containers.py      resolve exact quay.io tags → containers.lock.yaml
  apptainer --version                    verify Apptainer ≥ 1.1 is installed

Config file (config.yaml) must be edited before running.
See INSTALL.md for database and container setup instructions.
"""


def find_snakefile(cli_path):
    """Locate the Snakefile: CLI arg > next to vapor script > CWD."""
    if cli_path:
        p = Path(cli_path).resolve()
        if not p.exists():
            sys.exit(f"ERROR: Snakefile not found at {p}")
        return p
    candidates = [
        Path(__file__).resolve().parent / "Snakefile",
        Path.cwd() / "Snakefile",
    ]
    for p in candidates:
        if p.exists():
            return p
    sys.exit(
        "ERROR: Snakefile not found.\n"
        "Run vapor from the VAPOR directory, or use --snakefile to specify its path."
    )


def find_config(cli_path):
    """Locate config.yaml: CLI arg > CWD > next to Snakefile."""
    if cli_path:
        p = Path(cli_path).resolve()
        if not p.exists():
            sys.exit(f"ERROR: Config file not found: {p}")
        return p
    candidates = [
        Path.cwd() / "config.yaml",
        Path(__file__).resolve().parent / "config.yaml",
    ]
    for p in candidates:
        if p.exists():
            return p
    sys.exit(
        "ERROR: config.yaml not found.\n"
        "Copy config.yaml to your working directory and edit it before running."
    )


def _has_lock(snakefile: Path) -> bool:
    """True if containers.lock.yaml exists next to the Snakefile."""
    return (snakefile.parent / "containers.lock.yaml").exists()


def build_command(args, snakefile, config_path):
    use_conda = args.use_conda or not _has_lock(snakefile)

    if use_conda:
        mode_flags = ["--use-conda"]
    else:
        mode_flags = ["--use-apptainer", "--use-conda"]

    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--configfile", str(config_path),
        "--cores", str(args.threads) if args.threads else "all",
    ] + mode_flags

    if args.apptainer_args:
        cmd += ["--apptainer-args", args.apptainer_args]
    if args.dry_run:
        cmd.append("--dry-run")
    if args.rerun_incomplete:
        cmd.append("--rerun-incomplete")
    if args.force:
        cmd.append("--forceall")
    if args.forcerun:
        cmd += ["--forcerun"] + args.forcerun
    if args.unlock:
        cmd.append("--unlock")
    if args.target:
        cmd += args.target

    return cmd, "conda" if use_conda else "containers"


def main():
    parser = argparse.ArgumentParser(
        prog="vapor",
        description=DESCRIPTION,
        epilog=EPILOG,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--version", "-v",
        action="version",
        version=f"VAPOR {__version__}",
    )

    # ── Config / paths ───────────────────────────────────────────────
    parser.add_argument(
        "--config", "-c",
        metavar="FILE",
        default=None,
        help="Path to config YAML file (default: ./config.yaml)",
    )
    parser.add_argument(
        "--snakefile", "-s",
        metavar="FILE",
        default=None,
        help="Path to Snakefile (auto-detected if not specified)",
    )

    # ── Execution mode ────────────────────────────────────────────────
    parser.add_argument(
        "--use-conda",
        action="store_true",
        help=(
            "Use conda environments instead of containers. "
            "Default when containers.lock.yaml is absent."
        ),
    )
    parser.add_argument(
        "--apptainer-args",
        metavar="ARGS",
        default=None,
        help=(
            "Extra arguments forwarded to apptainer exec, e.g. "
            "'--nv' to enable GPU pass-through inside containers. "
            "Quote the whole string: --apptainer-args '--nv --bind /data'."
        ),
    )

    # ── Execution control ─────────────────────────────────────────────
    parser.add_argument(
        "--threads", "-t",
        type=int,
        metavar="N",
        default=None,
        help="Number of CPU cores to use (default: all available)",
    )
    parser.add_argument(
        "--dry-run", "-n",
        action="store_true",
        help="Validate workflow without executing any jobs",
    )
    parser.add_argument(
        "--rerun-incomplete",
        action="store_true",
        help="Re-run jobs that were left incomplete in a previous run",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Force re-execution of all rules from scratch",
    )
    parser.add_argument(
        "--forcerun", "-R",
        metavar="RULE",
        nargs="+",
        help="Force re-execution of specific rule(s), e.g. --forcerun quast",
    )
    parser.add_argument(
        "--unlock",
        action="store_true",
        help="Unlock the working directory (use after a crash)",
    )
    parser.add_argument(
        "--target",
        metavar="TARGET",
        nargs="+",
        help="Run pipeline only up to specific output file(s)",
    )
    parser.add_argument(
        "--create-envs",
        action="store_true",
        help="Create all conda environments without running the pipeline",
    )

    # ── Utilities ─────────────────────────────────────────────────────
    parser.add_argument(
        "--dag",
        action="store_true",
        help="Generate workflow DAG and save as dag.svg (requires graphviz)",
    )

    args = parser.parse_args()

    snakefile   = find_snakefile(args.snakefile)
    config_path = find_config(args.config)

    print(f"[VAPOR {__version__}] Snakefile : {snakefile}")
    print(f"[VAPOR {__version__}] Config    : {config_path}")

    # ── --create-envs ─────────────────────────────────────────────────
    if args.create_envs:
        cmd = [
            "snakemake",
            "--snakefile", str(snakefile),
            "--configfile", str(config_path),
            "--use-conda",
            "--cores", "1",
            "--conda-create-envs-only",
        ]
        print(f"[VAPOR] Creating conda environments: {' '.join(cmd)}\n")
        result = subprocess.run(cmd)
        sys.exit(result.returncode)

    # ── --dag ─────────────────────────────────────────────────────────
    if args.dag:
        dag_cmd = [
            "snakemake",
            "--snakefile", str(snakefile),
            "--configfile", str(config_path),
            "--dag",
        ]
        out_svg = Path.cwd() / "dag.svg"
        print(f"[VAPOR] Generating DAG → {out_svg}")
        try:
            p1 = subprocess.Popen(dag_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
            with open(out_svg, "w") as f:
                p2 = subprocess.Popen(
                    ["dot", "-Tsvg"], stdin=p1.stdout, stdout=f, stderr=subprocess.PIPE
                )
                p1.stdout.close()
                _, dot_err = p2.communicate()
            if p2.returncode != 0:
                print(f"[VAPOR] WARNING: dot exited with code {p2.returncode}", file=sys.stderr)
                if dot_err:
                    print(dot_err.decode(), file=sys.stderr)
        except FileNotFoundError:
            sys.exit("ERROR: 'dot' not found. Install graphviz: conda install -c conda-forge graphviz")
        print(f"[VAPOR] DAG saved to {out_svg}")
        return

    # ── normal run ────────────────────────────────────────────────────
    cmd, mode = build_command(args, snakefile, config_path)
    print(f"[VAPOR] Mode      : {mode}")
    print(f"[VAPOR] Running   : {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
