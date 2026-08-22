"""Ponto de entrada Snakemake do report v2. Toda a logica esta no pacote."""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from report.renderer_v2 import build_data, write_report  # noqa: E402

write_report(build_data(snakemake), snakemake.output.html)  # noqa: F821
