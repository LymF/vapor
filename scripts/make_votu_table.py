#!/usr/bin/env python3
"""
make_votu_table.py — Build the vOTU summary table for one sample.

Inputs (all wired via Snakemake):
  cluster_tsv       mmseqs/{sample}_cluster.tsv
  viral_nr_fasta    viral/consensus/{sample}_viral_nonredundant.fasta
  checkv_tsv        viral/checkv/quality_summary.tsv
  vibrant_dir       viral/vibrant/
  taxonomy_tsv      viral/taxonomy/viral_taxonomy_merged.tsv
  phist_csv         viral/phist/phist_results.csv

Output:
  viral/votu/{sample}_vOTU_table.tsv
"""

import csv
import glob
import os
import sys
from collections import defaultdict


# ──────────────────────────────────────────────────────────────────────────────
#  Helpers
# ──────────────────────────────────────────────────────────────────────────────

def read_fasta_ids(path):
    """Return ordered list of sequence IDs from a FASTA file."""
    ids = []
    if not os.path.exists(path):
        return ids
    with open(path) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].split()[0])
    return ids


def read_tsv(path, required=True):
    """Yield dicts from a TSV with a header row. Returns [] if missing."""
    if not os.path.exists(path):
        if required:
            print(f"[make_votu_table] WARNING: missing {path}", file=sys.stderr)
        return []
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            rows.append(dict(row))
    return rows


def read_csv(path, required=False):
    """Yield dicts from a CSV with a header row. Returns [] if missing."""
    if not os.path.exists(path):
        if required:
            print(f"[make_votu_table] WARNING: missing {path}", file=sys.stderr)
        return []
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f)
        for row in reader:
            rows.append(dict(row))
    return rows


# ──────────────────────────────────────────────────────────────────────────────
#  1. vOTU representatives from viral_nonredundant.fasta
# ──────────────────────────────────────────────────────────────────────────────

def load_votu_ids(viral_nr_fasta):
    ids = read_fasta_ids(viral_nr_fasta)
    print(f"[make_votu_table] vOTUs (viral_nonredundant): {len(ids)}", file=sys.stderr)
    return ids


# ──────────────────────────────────────────────────────────────────────────────
#  2. Cluster membership from MMseqs2 cluster.tsv
# ──────────────────────────────────────────────────────────────────────────────

def load_cluster_members(cluster_tsv, votu_ids):
    """
    cluster.tsv format: representative <TAB> member  (one pair per line).
    Count members per representative, filtered to viral representatives only.
    """
    votu_set = set(votu_ids)
    members = defaultdict(set)
    if not os.path.exists(cluster_tsv):
        print(f"[make_votu_table] WARNING: missing {cluster_tsv}", file=sys.stderr)
        return {}
    with open(cluster_tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            rep, mem = parts[0], parts[1]
            if rep in votu_set:
                members[rep].add(mem)
    result = {rep: len(mems) for rep, mems in members.items()}
    print(f"[make_votu_table] Cluster members found for {len(result)} vOTUs", file=sys.stderr)
    return result


# ──────────────────────────────────────────────────────────────────────────────
#  3. CheckV quality
# ──────────────────────────────────────────────────────────────────────────────

def load_checkv(checkv_tsv):
    """
    Returns {contig_id: {quality, completeness, length, genome_type}}.
    CheckV columns: contig_id, contig_length, proviral_length, provirus,
                    checkv_quality, miuvig_quality, completeness,
                    completeness_method, contamination, kmer_freq, warnings
    """
    data = {}
    for row in read_tsv(checkv_tsv, required=False):
        cid = row.get("contig_id", "").strip()
        if not cid:
            continue
        # genome_type: circular (provirus=No + circular in warnings) vs linear vs provirus
        provirus = row.get("provirus", "No").strip()
        warnings = row.get("warnings", "").lower()
        if provirus.lower() == "yes":
            gtype = "provirus"
        elif "circular" in warnings:
            gtype = "circular"
        else:
            gtype = "linear"
        data[cid] = {
            "checkv_quality":       row.get("checkv_quality", ""),
            "checkv_completeness":  row.get("completeness", ""),
            "checkv_length":        row.get("contig_length", ""),
            "genome_type":          gtype,
        }
    print(f"[make_votu_table] CheckV entries: {len(data)}", file=sys.stderr)
    return data


# ──────────────────────────────────────────────────────────────────────────────
#  4. Lifestyle from VIBRANT
#     - Scaffolds in VIBRANT_integrated_prophage_*.tsv → "lysogenic"
#     - All other VIBRANT scaffolds → "lytic"
#     - Not detected by VIBRANT → "unknown"
# ──────────────────────────────────────────────────────────────────────────────

def load_vibrant_lifestyle(vibrant_dir):
    """
    Returns {contig_id: "lytic" | "lysogenic" | "unknown"}.
    Also returns {contig_id: n_amgs} from VIBRANT AMG individuals TSV.
    """
    lysogenic = set()
    all_vibrant = set()
    amg_counts = defaultdict(int)

    # Integrated prophages → lysogenic
    for tsv in glob.glob(os.path.join(vibrant_dir, "**", "VIBRANT_integrated_prophage_*.tsv"),
                         recursive=True):
        for row in read_tsv(tsv, required=False):
            scaffold = row.get("scaffold", row.get("Scaffold", "")).strip()
            if scaffold:
                lysogenic.add(scaffold)
                all_vibrant.add(scaffold)

    # Summary results → all detected phages
    for tsv in glob.glob(os.path.join(vibrant_dir, "**", "VIBRANT_summary_results_*.tsv"),
                         recursive=True):
        for row in read_tsv(tsv, required=False):
            scaffold = row.get("scaffold", row.get("Scaffold", "")).strip()
            if scaffold:
                all_vibrant.add(scaffold)

    # AMG individuals → count per scaffold
    for tsv in glob.glob(os.path.join(vibrant_dir, "**", "VIBRANT_AMG_individuals_*.tsv"),
                         recursive=True):
        for row in read_tsv(tsv, required=False):
            scaffold = row.get("scaffold", row.get("Scaffold", "")).strip()
            if scaffold:
                amg_counts[scaffold] += 1

    # Build lifestyle dict
    lifestyle = {}
    for scaffold in all_vibrant:
        lifestyle[scaffold] = "lysogenic" if scaffold in lysogenic else "lytic"

    print(f"[make_votu_table] VIBRANT scaffolds: {len(all_vibrant)} "
          f"({len(lysogenic)} lysogenic)", file=sys.stderr)
    print(f"[make_votu_table] AMG counts: {len(amg_counts)} scaffolds with AMGs",
          file=sys.stderr)
    return lifestyle, dict(amg_counts)


# ──────────────────────────────────────────────────────────────────────────────
#  5. Taxonomy from viral_taxonomy_merged.tsv
# ──────────────────────────────────────────────────────────────────────────────

def load_taxonomy(taxonomy_tsv):
    """
    Returns {seq_name: {taxonomy_family, taxonomy_genus, taxonomy_source}}.
    Columns: seq_name, final_family, final_genus, best_taxonomy, source.
    """
    data = {}
    for row in read_tsv(taxonomy_tsv, required=False):
        name = row.get("seq_name", "").strip()
        if not name:
            continue
        data[name] = {
            "taxonomy_family": row.get("final_family", ""),
            "taxonomy_genus":  row.get("final_genus", ""),
            "taxonomy_order":  row.get("final_order", ""),
            "taxonomy_best":   row.get("best_taxonomy", ""),
            "taxonomy_source": row.get("source", ""),
        }
    print(f"[make_votu_table] Taxonomy entries: {len(data)}", file=sys.stderr)
    return data


# ──────────────────────────────────────────────────────────────────────────────
#  6. Host prediction from PHIST (soft-fail)
#     PHIST output: phist_results.csv
#     Columns: phage, host, #common-kmers, pvalue, adj-pvalue
# ──────────────────────────────────────────────────────────────────────────────

def load_phist(phist_csv):
    """
    Returns {virus_id: {host_bin, host_score}}.
    Extracts virus name from the phage path basename (without .fasta).
    Returns {} silently if file not found or empty.
    """
    data = {}
    if not os.path.exists(phist_csv):
        return data
    for row in read_csv(phist_csv, required=False):
        phage_path = row.get("phage", "").strip()
        host_path  = row.get("host",  "").strip()
        if not phage_path:
            continue
        virus = os.path.basename(phage_path)
        for ext in (".fasta", ".fa", ".fna"):
            if virus.endswith(ext):
                virus = virus[:-len(ext)]
                break
        host_bin = os.path.basename(host_path)
        for ext in (".fa", ".fasta", ".fna"):
            if host_bin.endswith(ext):
                host_bin = host_bin[:-len(ext)]
                break
        # Keep best (lowest adj-pvalue) hit per virus
        score = row.get("adj-pvalue", "")
        try:
            score_f = float(score)
        except (ValueError, TypeError):
            score_f = 1.0
        if virus not in data or score_f < data[virus]["_score"]:
            data[virus] = {
                "host_bin":   host_bin,
                "host_score": score,
                "_score":     score_f,
            }
    # Drop internal sort key
    for v in data:
        data[v].pop("_score", None)
    print(f"[make_votu_table] PHIST host predictions: {len(data)}", file=sys.stderr)
    return data


# ──────────────────────────────────────────────────────────────────────────────
#  MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    # Snakemake injects these names into globals when running as a script
    viral_nr_fasta = snakemake.input.viral_nr
    cluster_tsv    = snakemake.input.cluster
    checkv_tsv     = snakemake.input.checkv
    vibrant_dir    = snakemake.params.vibrant_dir
    taxonomy_tsv   = snakemake.input.taxonomy
    phist_csv      = snakemake.params.phist_csv
    out_tsv        = snakemake.output.tsv
    sample         = snakemake.wildcards.sample

    print(f"[make_votu_table] Sample: {sample}", file=sys.stderr)

    # Load all sources
    votu_ids  = load_votu_ids(viral_nr_fasta)
    members   = load_cluster_members(cluster_tsv, votu_ids)
    checkv    = load_checkv(checkv_tsv)
    lifestyle, amg_counts = load_vibrant_lifestyle(vibrant_dir)
    taxonomy  = load_taxonomy(taxonomy_tsv)
    hosts     = load_phist(phist_csv)

    # Build output rows
    COLS = [
        "vOTU_id", "sample",
        "n_members",
        "checkv_quality", "checkv_completeness", "checkv_length", "genome_type",
        "lifestyle",
        "taxonomy_family", "taxonomy_genus", "taxonomy_order",
        "taxonomy_best", "taxonomy_source",
        "n_AMGs",
        "host_bin", "host_score",
    ]

    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()
        for vid in votu_ids:
            cv   = checkv.get(vid, {})
            tax  = taxonomy.get(vid, {})
            host = hosts.get(vid, {})
            row = {
                "vOTU_id":             vid,
                "sample":              sample,
                "n_members":           members.get(vid, 1),
                "checkv_quality":      cv.get("checkv_quality", ""),
                "checkv_completeness": cv.get("checkv_completeness", ""),
                "checkv_length":       cv.get("checkv_length", ""),
                "genome_type":         cv.get("genome_type", ""),
                "lifestyle":           lifestyle.get(vid, "unknown"),
                "taxonomy_family":     tax.get("taxonomy_family", ""),
                "taxonomy_genus":      tax.get("taxonomy_genus", ""),
                "taxonomy_order":      tax.get("taxonomy_order", ""),
                "taxonomy_best":       tax.get("taxonomy_best", ""),
                "taxonomy_source":     tax.get("taxonomy_source", ""),
                "n_AMGs":              amg_counts.get(vid, 0),
                "host_bin":            host.get("host_bin", ""),
                "host_score":          host.get("host_score", ""),
            }
            writer.writerow(row)

    n_lytic = sum(1 for v in votu_ids if lifestyle.get(v, "unknown") == "lytic")
    n_lyso  = sum(1 for v in votu_ids if lifestyle.get(v, "unknown") == "lysogenic")
    n_amg   = sum(1 for v in votu_ids if amg_counts.get(v, 0) > 0)
    print(f"[make_votu_table] Written {len(votu_ids)} vOTUs to {out_tsv}", file=sys.stderr)
    print(f"[make_votu_table]   Lytic: {n_lytic}  Lysogenic: {n_lyso}  "
          f"With AMGs: {n_amg}", file=sys.stderr)


if __name__ == "__main__":
    main()
