#!/usr/bin/env python3
"""
generate_report.py — VAPOR HTML report entry point for Snakemake.

This file is kept minimal; all logic lives in scripts/report/.
Snakemake adds scripts/ to sys.path when running via script: directive,
so `from report.renderer import build_report` resolves to scripts/report/renderer.py.
"""
import sys
import os

# Ensure scripts/ is on sys.path (fallback for direct execution)
_script_dir = os.path.dirname(os.path.abspath(__file__))
if _script_dir not in sys.path:
    sys.path.insert(0, _script_dir)

from report.renderer import build_report

build_report(snakemake)  # noqa: F821 — snakemake injected by Snakemake runtime
