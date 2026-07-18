#!/usr/bin/env python3
"""
VAPOR — Viral And Prokaryotic genOme Recovery
CLI wrapper for the Snakemake metagenomics pipeline.

Usage:
  vapor [options]

Run 'vapor -h' for full help.
"""

import argparse
import os
import subprocess
import sys
from pathlib import Path

__version__ = "1.4.0"

DESCRIPTION = """\
VAPOR — Viral And Prokaryotic genOme Recovery
A comprehensive pipeline for viral and prokaryotic genome recovery
from metagenomic short-read (Illumina) and long-read (ONT/HiFi) data.
"""

EPILOG = """\
Examples:
  vapor --dry-run                             validate workflow without executing
  vapor --threads 32                          run with 32 CPU cores
  vapor --config /path/to/config.yaml         use a custom config file
  vapor --rerun-incomplete                    resume an interrupted run
  vapor --force                               re-run everything from scratch
  vapor --forcerun quast viral_detection      force specific rules to re-run
  vapor --dag                                 generate workflow DAG (dag.svg)
  vapor --unlock                              unlock directory after crash

Module switches (override config.yaml without editing):
  vapor --set use_spades=false                skip SPAdes (low-RAM server)
  vapor --set use_cobra=true                  enable viral contig extension
  vapor --set use_pharokka=false              skip phage annotation
  vapor --set use_phold=false                 skip structure-based annotation
  vapor --set use_bakta=false                 skip prokaryotic MAG annotation
  vapor --set use_eggnog=false                skip COG/KEGG functional annotation
  vapor --set use_comebin=false               skip COMEBin (no GPU available)
  vapor --set use_gunc=false                  skip chimera detection
  vapor --set use_mag_derep=false             skip MAG dereplication
  vapor --set use_amr_consensus=false         skip AMR consensus step
  vapor --set use_defense_viral=false         skip viral defense annotation
  vapor --set use_gpu=true                    enable GPU (VAMB/SemiBin2/COMEBin)

Performance and thresholds:
  vapor --set threads=40                      override thread count
  vapor --set min_contig=5000                 raise minimum contig length
  vapor --set votu_ani=95                     vOTU clustering ANI threshold
  vapor --set pharokka_min_completeness=70    lower phage annotation threshold
  vapor --set viral_consensus_mode=count      switch viral detection mode

Combined:
  vapor --threads 40 --set use_spades=false --set use_comebin=false
  vapor --dry-run --set use_pharokka=false --set use_phold=false

Config file (config.yaml) must be edited before running.
See INSTALL.md for database setup instructions.
"""

# Config keys that hold filesystem paths needing apptainer bind mounts.
_PATH_KEYS = [
    "fastq_dir", "outdir",
    "checkv_db", "vs2_db", "vibrant_base", "genomad_db",
    "checkm2_db", "inphared_db", "vcontact3_db", "gtdbtk_db",
    "gunc_db", "pharokka_db", "phold_db", "bakta_db", "eggnog_db",
    "card_db", "deeparg_db", "defense_finder_models_db", "apis_db",
    "custom_prok_mmseqs_db", "custom_viral_mmseqs_db",
    "host_genome", "host_index",
]


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


def _detect_executor():
    """Auto-detect available container runtime: apptainer > singularity > conda."""
    import shutil
    if shutil.which("apptainer"):
        return "apptainer"
    if shutil.which("singularity"):
        return "singularity"
    return "conda"


def _collect_bind_paths(config_path):
    """Return sorted list of existing directories to bind into the container.

    Resolves each config path to an existing directory (walks up if needed,
    converts file paths to their parent dir), then removes any path that is
    already covered by a higher-level entry to avoid redundant / conflicting
    apptainer bind mounts.
    """
    try:
        import yaml
        with open(config_path) as fh:
            cfg = yaml.safe_load(fh) or {}
    except Exception:
        return []

    raw: set[Path] = set()
    for key in _PATH_KEYS:
        val = cfg.get(key, "")
        if not val or not isinstance(val, str):
            continue
        p = Path(os.path.expanduser(val)).resolve()
        # Walk up to the first existing ancestor (outdir may not exist yet).
        while not p.exists() and p != p.parent:
            p = p.parent
        if not p.exists():
            continue
        # Always bind directories, not individual files.
        if p.is_file():
            p = p.parent
        raw.add(p)

    # Drop any path that is a sub-directory of another path already in the set
    # to prevent redundant / conflicting bind mounts (e.g. fastq_dir + outdir
    # both under the same parent).
    ordered = sorted(raw)
    deduped = [
        p for p in ordered
        if not any(
            p != other and str(p).startswith(str(other) + os.sep)
            for other in ordered
        )
    ]
    return [str(p) for p in deduped]


def build_command(args, snakefile, config_path):
    cmd = [
        "snakemake",
        "--snakefile", str(snakefile),
        "--configfile", str(config_path),
        "--cores", str(args.threads) if args.threads else "all",
    ]

    executor = args.executor if args.executor != "auto" else _detect_executor()
    args._resolved_executor = executor

    if executor == "conda":
        cmd.append("--use-conda")
    elif executor in ("singularity", "apptainer"):
        cmd.append("--use-apptainer" if executor == "apptainer" else "--use-singularity")

        # Auto-derive bind mounts from config paths and merge with user-supplied args.
        # --writable-tmpfs lets tools like VirSorter2 write to their install dir on
        # first run without modifying the read-only container image.
        auto_binds = _collect_bind_paths(config_path)
        # Always bind the pipeline directory so scripts/ is accessible inside containers
        # when vapor is invoked from outside the pipeline directory.
        pipeline_dir = str(Path(snakefile).resolve().parent)
        if not any(pipeline_dir == p or pipeline_dir.startswith(p + os.sep)
                   for p in auto_binds):
            auto_binds = [pipeline_dir] + auto_binds
        bind_flag  = "--bind " + ",".join(auto_binds) if auto_binds else ""
        extra      = args.singularity_args.strip()
        combined   = " ".join(filter(None, ["--writable-tmpfs", extra, bind_flag]))

        if combined:
            flag = "--apptainer-args" if executor == "apptainer" else "--singularity-args"
            cmd += [flag, combined]

    # --set KEY=VALUE overrides passed to Snakemake's --config mechanism.
    if args.set_config:
        cmd += ["--config"] + args.set_config

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

    return cmd


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
    parser.add_argument(
        "--set",
        metavar="KEY=VALUE",
        action="append",
        dest="set_config",
        default=[],
        help=(
            "Override a config.yaml value without editing the file. "
            "Repeatable: --set use_spades=false --set threads=32. "
            "Keys match config.yaml (e.g. use_spades, use_cobra, min_contig, "
            "pharokka_min_completeness, votu_ani)."
        ),
    )

    # ── Execution mode ────────────────────────────────────────────────
    parser.add_argument(
        "--executor", "-e",
        choices=["auto", "conda", "singularity", "apptainer"],
        default="auto",
        metavar="MODE",
        help="Software deployment: auto (default, detects apptainer > singularity > conda) | conda | singularity | apptainer",
    )
    parser.add_argument(
        "--singularity-args",
        metavar="ARGS",
        default="",
        dest="singularity_args",
        help='Extra args for Singularity/Apptainer, e.g. "--nv" for GPU pass-through',
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

    if args.dag:
        import shutil
        # Resolve dot binary: PATH first, then active conda env, then common prefixes.
        dot_bin = shutil.which("dot")
        if not dot_bin:
            for prefix_env in ("CONDA_PREFIX", "CONDA_EXE"):
                prefix = os.environ.get(prefix_env, "")
                if prefix:
                    # CONDA_EXE points to the conda binary; walk up to the env root.
                    root = Path(prefix) if prefix_env == "CONDA_PREFIX" else Path(prefix).parent.parent
                    candidate = root / "bin" / "dot"
                    if candidate.is_file():
                        dot_bin = str(candidate)
                        break
        if not dot_bin:
            sys.exit("ERROR: 'dot' not found. Install graphviz: conda install -c conda-forge graphviz")

        dag_cmd = [
            "snakemake",
            "--snakefile", str(snakefile),
            "--configfile", str(config_path),
            "--dag",
        ]
        out_svg = Path.cwd() / "dag.svg"
        print(f"[VAPOR] Generating DAG → {out_svg}")
        p1 = subprocess.Popen(dag_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        with open(out_svg, "w") as f:
            p2 = subprocess.Popen(
                [dot_bin, "-Tsvg"], stdin=p1.stdout, stdout=f, stderr=subprocess.PIPE
            )
            p1.stdout.close()
            _, dot_err = p2.communicate()
        if p2.returncode != 0:
            print(f"[VAPOR] WARNING: dot exited with code {p2.returncode}", file=sys.stderr)
            if dot_err:
                print(dot_err.decode(), file=sys.stderr)
        print(f"[VAPOR] DAG saved to {out_svg}")
        return

    cmd = build_command(args, snakefile, config_path)
    print(f"[VAPOR] Mode      : {args._resolved_executor}")
    if args.set_config:
        print(f"[VAPOR] Overrides : {' '.join(args.set_config)}")
    print(f"[VAPOR] Running   : {' '.join(cmd)}\n")
    result = subprocess.run(cmd)
    sys.exit(result.returncode)


if __name__ == "__main__":
    main()
