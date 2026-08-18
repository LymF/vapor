#!/usr/bin/env python3
"""
filter_putative_amgs.py — Flag candidate auxiliary metabolic genes (AMGs)
from an eggNOG-mapper annotation table over viral ORFs.

Output is CANDIDATES, not confirmed AMGs. AMG calling from annotation alone
is a recognized source of false positives in the phage literature — a hit
on a metabolism KEGG pathway can come from contamination, mis-binning, or a
housekeeping gene the pathway map simply also lists — and confirming a real
AMG requires inspecting its genomic context (is it flanked by phage genes,
is it inside a clearly viral region, etc.). For that reason this script and
everything downstream of it names its output "putative AMGs", never "AMGs".

Criterion (the same operational definition VIBRANT and DRAM-v use): a
protein is a putative AMG if
  - KEGG_ko is non-empty (and not "-"), AND
  - KEGG_Pathway contains at least one metabolism map, i.e. an identifier
    matching ko00\\d{3} or map00\\d{3}.

Usage:
    python3 filter_putative_amgs.py <emapper.annotations> <putative_amgs.tsv>
"""
import csv
import re
import sys

METABOLISM_MAP_RE = re.compile(r"(?:ko|map)00\d{3}")

OUT_COLS = ["votu_id", "protein_id", "KEGG_ko", "KEGG_Pathway", "COG_category", "Description"]


def _votu_id_from_protein_id(protein_id):
    """Strip the trailing "_<N>" that prodigal appends to the contig ID."""
    return re.sub(r"_\d+$", "", protein_id)


def _read_emapper_annotations(path):
    """eggNOG-mapper writes a handful of '##' comment lines around the TSV
    body, and the header line itself starts with '#query' rather than
    'query'. Skip comments, normalize the header, and yield dict rows."""
    with open(path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            if line.startswith("##"):
                continue
            if line.startswith("#"):
                header = line.lstrip("#").split("\t")
                header[0] = header[0] if header[0] else "query"
                continue
            if header is None:
                continue
            values = line.split("\t")
            yield dict(zip(header, values))


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: filter_putative_amgs.py <emapper.annotations> <putative_amgs.tsv>")

    in_path, out_path = sys.argv[1], sys.argv[2]

    n_total = 0
    n_amg = 0
    with open(out_path, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, fieldnames=OUT_COLS, delimiter="\t")
        writer.writeheader()

        for row in _read_emapper_annotations(in_path):
            n_total += 1
            protein_id = row.get("query", "").strip()
            if not protein_id:
                continue
            kegg_ko = (row.get("KEGG_ko") or "").strip()
            kegg_pathway = (row.get("KEGG_Pathway") or "").strip()

            if not kegg_ko or kegg_ko == "-":
                continue
            if not kegg_pathway or kegg_pathway == "-":
                continue
            if not METABOLISM_MAP_RE.search(kegg_pathway):
                continue

            writer.writerow({
                "votu_id":      _votu_id_from_protein_id(protein_id),
                "protein_id":   protein_id,
                "KEGG_ko":      kegg_ko,
                "KEGG_Pathway": kegg_pathway,
                "COG_category": (row.get("COG_category") or "").strip(),
                "Description":  (row.get("Description") or "").strip(),
            })
            n_amg += 1

    print(f"[filter_putative_amgs] {n_amg} putative AMG(s) out of {n_total} annotated protein(s)",
          file=sys.stderr)


if __name__ == "__main__":
    main()
