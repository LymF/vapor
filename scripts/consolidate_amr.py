#!/usr/bin/env python3
"""
consolidate_amr.py — merge AMRFinderPlus + RGI/CARD + DeepARG AMR hits
by CDS locus and compute a cross-tool consensus score.

RGI is ARO-native (CARD).  AMRFinderPlus and DeepARG are normalized to ARO
by argNorm upstream.  Consensus score = n_tools_that_detected / N_TOOLS.

Output columns:
  locus, aro_accession, gene_name, drug_class, resistance_mechanism,
  n_tools, consensus_score, tools_detected
"""
import csv
import sys
from collections import defaultdict
from pathlib import Path

N_TOOLS = 3  # AMRFinderPlus, RGI, DeepARG


# ── helpers ─────────────────────────────────────────────────────────────────

def _norm_aro(raw: str) -> str:
    """Normalize ARO to 'ARO:NNNNNNN' format; return '' if not parseable."""
    s = (raw or "").strip()
    if not s or s in ("NA", "N/A", "-", ""):
        return ""
    if s.startswith("ARO:"):
        return s
    if s.isdigit():
        return f"ARO:{s}"
    if ":" in s:
        return s
    return ""


def _best(*candidates: str) -> str:
    """Return first non-empty candidate."""
    for c in candidates:
        v = (c or "").strip()
        if v and v not in ("NA", "N/A", "-"):
            return v
    return ""


def _has_data(path: str) -> bool:
    p = Path(path)
    if not p.exists() or p.stat().st_size == 0:
        return False
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#") and line.strip():
                return True
    return False


def _get(row: dict, *keys: str) -> str:
    for k in keys:
        v = row.get(k, "")
        if v and v.strip() not in ("NA", "N/A", "-", ""):
            return v.strip()
    return ""


# ── per-tool parsers ─────────────────────────────────────────────────────────

def _parse_amrfinder_normed(path: str) -> dict:
    """
    Returns: locus → {aro, gene, drug_class, resistance_mechanism}
    argNorm appends 'ARO Accession', 'Drug Class', 'Resistance Mechanism',
    'AMR Gene Family' to the original AMRFinderPlus TSV.
    """
    hits = {}
    if not _has_data(path):
        return hits
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            locus = _get(row, "Protein identifier")
            if not locus:
                continue
            hits[locus] = {
                "aro": _norm_aro(_get(row, "ARO Accession", "ARO")),
                "gene": _get(row, "Gene symbol", "gene_name", "AMR Gene Family"),
                "drug_class": _get(row, "Drug Class", "Class"),
                "resistance_mechanism": _get(row, "Resistance Mechanism"),
            }
    return hits


def _parse_rgi(path: str) -> dict:
    """
    RGI native output: locus → {aro, gene, drug_class, resistance_mechanism}
    RGI column 'ARO' is the numeric or 'ARO:XXXXXX' accession.
    'Best_Hit_ARO' is the human-readable ARO term name.
    """
    hits = {}
    if not _has_data(path):
        return hits
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            locus = _get(row, "ORF_ID")
            if not locus:
                continue
            hits[locus] = {
                "aro": _norm_aro(_get(row, "ARO")),
                "gene": _get(row, "Best_Hit_ARO", "AMR Gene Family"),
                "drug_class": _get(row, "Drug Class"),
                "resistance_mechanism": _get(row, "Resistance Mechanism"),
            }
    return hits


def _parse_deeparg_normed(path: str) -> dict:
    """
    DeepARG normed output: locus → {aro, gene, drug_class, resistance_mechanism}
    Header column is '#ARG' (or 'ARG' after stripping '#').
    argNorm adds same extra columns as for AMRFinderPlus.
    """
    hits = {}
    if not _has_data(path):
        return hits
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        # handle '#ARG' column name produced by DeepARG
        fieldnames = reader.fieldnames or []
        locus_col = next(
            (c for c in fieldnames if c.lstrip("#").strip().upper() == "ARG"),
            None
        )
        if locus_col is None:
            return hits
        for row in reader:
            locus = (row.get(locus_col) or "").strip()
            if not locus:
                continue
            hits[locus] = {
                "aro": _norm_aro(_get(row, "ARO Accession", "ARO")),
                "gene": _get(row, "gene_name", "AMR Gene Family", locus_col),
                "drug_class": _get(row, "Drug Class", "predicted_aro-drug_class"),
                "resistance_mechanism": _get(row, "Resistance Mechanism"),
            }
    return hits


# ── consolidation ────────────────────────────────────────────────────────────

def consolidate(
    amrfinder_normed: str,
    rgi_results: str,
    deeparg_normed: str,
) -> list[dict]:
    """
    Merge hits from all three tools by locus.
    Returns list of result dicts, sorted by descending n_tools then locus.
    """
    amrf = _parse_amrfinder_normed(amrfinder_normed)
    rgi  = _parse_rgi(rgi_results)
    deep = _parse_deeparg_normed(deeparg_normed)

    all_loci: set[str] = set(amrf) | set(rgi) | set(deep)

    rows = []
    for locus in sorted(all_loci):
        tool_hits = {}
        if locus in amrf:
            tool_hits["AMRFinderPlus"] = amrf[locus]
        if locus in rgi:
            tool_hits["RGI"] = rgi[locus]
        if locus in deep:
            tool_hits["DeepARG"] = deep[locus]

        n_tools = len(tool_hits)
        tools_detected = ",".join(sorted(tool_hits))

        # Prefer curated calls (RGI first, then AMRFinder, then DeepARG)
        # for the canonical annotation fields.
        ordered = (
            tool_hits.get("RGI")
            or tool_hits.get("AMRFinderPlus")
            or tool_hits.get("DeepARG")
            or {}
        )
        aro  = ordered.get("aro", "")
        gene = _best(
            ordered.get("gene", ""),
            *(h.get("gene", "") for h in tool_hits.values()),
        )
        drug = _best(
            ordered.get("drug_class", ""),
            *(h.get("drug_class", "") for h in tool_hits.values()),
        )
        mech = _best(
            ordered.get("resistance_mechanism", ""),
            *(h.get("resistance_mechanism", "") for h in tool_hits.values()),
        )

        rows.append({
            "locus":                locus,
            "aro_accession":        aro,
            "gene_name":            gene,
            "drug_class":           drug,
            "resistance_mechanism": mech,
            "n_tools":              n_tools,
            "consensus_score":      f"{n_tools / N_TOOLS:.3f}",
            "tools_detected":       tools_detected,
        })

    rows.sort(key=lambda r: (-int(r["n_tools"]), r["locus"]))
    return rows


COLS = [
    "locus", "aro_accession", "gene_name", "drug_class",
    "resistance_mechanism", "n_tools", "consensus_score", "tools_detected",
]


def write_output(rows: list[dict], out_path: str) -> None:
    Path(out_path).parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=COLS, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


# ── entry points ─────────────────────────────────────────────────────────────

def _snakemake_main() -> None:
    rows = consolidate(
        amrfinder_normed=str(snakemake.input.amrfinder_normed),  # noqa: F821
        rgi_results=str(snakemake.input.rgi_results),
        deeparg_normed=str(snakemake.input.deeparg_normed),
    )
    write_output(rows, str(snakemake.output.consensus))


def _cli_main() -> None:
    import argparse

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--amrfinder-normed", required=True)
    parser.add_argument("--rgi-results", required=True)
    parser.add_argument("--deeparg-normed", required=True)
    parser.add_argument("-o", "--output", required=True)
    args = parser.parse_args()

    rows = consolidate(
        amrfinder_normed=args.amrfinder_normed,
        rgi_results=args.rgi_results,
        deeparg_normed=args.deeparg_normed,
    )
    write_output(rows, args.output)
    print(f"Written {len(rows)} consolidated AMR hits → {args.output}", file=sys.stderr)


if __name__ == "__main__":
    if "snakemake" in dir():
        _snakemake_main()
    else:
        _cli_main()
