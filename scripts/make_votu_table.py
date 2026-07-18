#!/usr/bin/env python3
"""
make_votu_table.py — Build the vOTU membership table for one sample.

One row per cluster member (including the representative itself), so the
full cluster structure is visible. Representative-level annotations
(CheckV quality, taxonomy, lifestyle, host) are propagated to all members.

Inputs (wired via Snakemake):
  votu_clusters   viral/votu/vOTU_clusters.tsv      (representative TAB member)
  votu_reps       viral/votu/votu_all_reps.fasta    (FASTA of all representatives)
  checkv_tsv      viral/checkv/quality_summary.tsv
  vibrant_dir     viral/vibrant/
  taxonomy_tsv    viral/taxonomy/viral_taxonomy_merged.tsv
  phist_csv       viral/phist/phist_results.csv

Output:
  viral/votu/{sample}_vOTU_table.tsv
"""

import csv
import glob
import os
import sys
from collections import defaultdict


def read_tsv(path, required=True):
    if not os.path.exists(path):
        if required:
            print(f"[make_votu_table] WARNING: missing {path}", file=sys.stderr)
        return []
    with open(path) as f:
        return list(csv.DictReader(f, delimiter="\t"))


def read_csv(path, required=False):
    if not os.path.exists(path):
        if required:
            print(f"[make_votu_table] WARNING: missing {path}", file=sys.stderr)
        return []
    with open(path) as f:
        return list(csv.DictReader(f))


def load_votu_clusters(clusters_tsv):
    """
    Returns:
      reps_ordered  — list of representative IDs in first-seen order
      rep_members   — {rep: [member, ...]} (members include the rep itself)
    """
    reps_ordered = []
    rep_members  = defaultdict(list)
    seen_reps    = set()

    for row in read_tsv(clusters_tsv, required=False):
        rep = row.get("representative", "").strip()
        mem = row.get("member", "").strip()
        if not rep:
            continue
        if rep not in seen_reps:
            reps_ordered.append(rep)
            seen_reps.add(rep)
        rep_members[rep].append(mem)

    print(f"[make_votu_table] vOTU clusters: {len(reps_ordered)} reps, "
          f"{sum(len(v) for v in rep_members.values())} total members", file=sys.stderr)
    return reps_ordered, dict(rep_members)


def load_rep_lengths(votu_reps_fasta):
    """Return {rep_id: length_bp} from the all-reps FASTA."""
    lengths = {}
    curr_id, curr_len = None, 0
    if not os.path.exists(votu_reps_fasta):
        return lengths
    with open(votu_reps_fasta) as f:
        for line in f:
            if line.startswith(">"):
                if curr_id:
                    lengths[curr_id] = curr_len
                curr_id = line[1:].split()[0]; curr_len = 0
            else:
                curr_len += len(line.strip())
    if curr_id:
        lengths[curr_id] = curr_len
    return lengths


def load_checkv(checkv_tsv):
    """Returns {contig_id: {quality, completeness, length, genome_type}}."""
    data = {}
    for row in read_tsv(checkv_tsv, required=False):
        cid = row.get("contig_id", "").strip()
        if not cid:
            continue
        provirus = row.get("provirus", "No").strip()
        warnings = row.get("warnings", "").lower()
        if provirus.lower() == "yes":
            gtype = "provirus"
        elif "circular" in warnings:
            gtype = "circular"
        else:
            gtype = "linear"
        data[cid] = {
            "checkv_quality":      row.get("checkv_quality", ""),
            "checkv_completeness": row.get("completeness", ""),
            "checkv_length":       row.get("contig_length", ""),
            "genome_type":         gtype,
        }
    print(f"[make_votu_table] CheckV entries: {len(data)}", file=sys.stderr)
    return data


def load_vibrant_lifestyle(vibrant_dir):
    """Returns ({contig: lifestyle}, {contig: n_amgs})."""
    lysogenic   = set()
    all_vibrant = set()
    amg_counts  = defaultdict(int)

    for tsv in glob.glob(os.path.join(vibrant_dir, "**",
                                      "VIBRANT_integrated_prophage_*.tsv"), recursive=True):
        for row in read_tsv(tsv, required=False):
            s = row.get("scaffold", row.get("Scaffold", "")).strip()
            if s:
                lysogenic.add(s); all_vibrant.add(s)

    for tsv in glob.glob(os.path.join(vibrant_dir, "**",
                                      "VIBRANT_summary_results_*.tsv"), recursive=True):
        for row in read_tsv(tsv, required=False):
            s = row.get("scaffold", row.get("Scaffold", "")).strip()
            if s:
                all_vibrant.add(s)

    for tsv in glob.glob(os.path.join(vibrant_dir, "**",
                                      "VIBRANT_AMG_individuals_*.tsv"), recursive=True):
        for row in read_tsv(tsv, required=False):
            s = row.get("scaffold", row.get("Scaffold", "")).strip()
            if s:
                amg_counts[s] += 1

    lifestyle = {s: ("lysogenic" if s in lysogenic else "lytic") for s in all_vibrant}
    print(f"[make_votu_table] VIBRANT: {len(all_vibrant)} scaffolds "
          f"({len(lysogenic)} lysogenic)", file=sys.stderr)
    return lifestyle, dict(amg_counts)


def load_taxonomy(taxonomy_tsv):
    """Returns {seq_name: {taxonomy_family, taxonomy_genus, taxonomy_order,
                           taxonomy_best, taxonomy_source}}."""
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


def load_phist(phist_csv):
    """Returns {virus_id: {host_bin, host_score}} — best (lowest adj-pvalue) hit."""
    data = {}
    for row in read_csv(phist_csv, required=False):
        phage_path = (row.get("phage") or "").strip()
        host_path  = (row.get("host")  or "").strip()
        if not phage_path:
            continue
        virus = os.path.basename(phage_path)
        for ext in (".fasta", ".fa", ".fna"):
            if virus.endswith(ext):
                virus = virus[:-len(ext)]; break
        host_bin = os.path.basename(host_path)
        for ext in (".fa", ".fasta", ".fna"):
            if host_bin.endswith(ext):
                host_bin = host_bin[:-len(ext)]; break
        score = row.get("adj-pvalue", "")
        try:
            score_f = float(score)
        except (ValueError, TypeError):
            score_f = 1.0
        if virus not in data or score_f < data[virus]["_score"]:
            data[virus] = {"host_bin": host_bin, "host_score": score, "_score": score_f}
    for v in data:
        data[v].pop("_score", None)
    print(f"[make_votu_table] PHIST host predictions: {len(data)}", file=sys.stderr)
    return data


# ──────────────────────────────────────────────────────────────────────────────
#  MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    clusters_tsv  = snakemake.input.votu_clusters
    votu_reps_fa  = snakemake.input.votu_reps
    checkv_tsv    = snakemake.input.checkv
    vibrant_dir   = snakemake.params.vibrant_dir
    taxonomy_tsv  = snakemake.input.taxonomy
    phist_csv     = snakemake.params.phist_csv
    out_tsv       = snakemake.output.tsv
    sample        = snakemake.wildcards.sample

    print(f"[make_votu_table] Sample: {sample}", file=sys.stderr)

    reps_ordered, rep_members = load_votu_clusters(clusters_tsv)
    rep_lengths               = load_rep_lengths(votu_reps_fa)
    checkv                    = load_checkv(checkv_tsv)
    lifestyle, amg_counts     = load_vibrant_lifestyle(vibrant_dir)
    taxonomy                  = load_taxonomy(taxonomy_tsv)
    hosts                     = load_phist(phist_csv)

    COLS = [
        "representative", "member", "is_rep", "cluster_size", "sample",
        "rep_length_bp",
        "checkv_quality", "checkv_completeness", "checkv_length", "genome_type",
        "lifestyle", "n_AMGs",
        "taxonomy_family", "taxonomy_genus", "taxonomy_order",
        "taxonomy_best", "taxonomy_source",
        "host_bin", "host_score",
    ]

    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    total_rows = 0
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        for rep in reps_ordered:
            members      = rep_members.get(rep, [rep])
            cluster_size = len(members)
            cv           = checkv.get(rep, {})
            tax          = taxonomy.get(rep, {})
            host         = hosts.get(rep, {})
            rep_lifestyle = lifestyle.get(rep, "unknown")
            rep_amgs      = amg_counts.get(rep, 0)
            rep_len       = rep_lengths.get(rep, "")

            for mem in members:
                writer.writerow({
                    "representative":      rep,
                    "member":              mem,
                    "is_rep":              "True" if mem == rep else "False",
                    "cluster_size":        cluster_size,
                    "sample":              sample,
                    "rep_length_bp":       rep_len,
                    "checkv_quality":      cv.get("checkv_quality", ""),
                    "checkv_completeness": cv.get("checkv_completeness", ""),
                    "checkv_length":       cv.get("checkv_length", ""),
                    "genome_type":         cv.get("genome_type", ""),
                    "lifestyle":           rep_lifestyle,
                    "n_AMGs":              rep_amgs,
                    "taxonomy_family":     tax.get("taxonomy_family", ""),
                    "taxonomy_genus":      tax.get("taxonomy_genus", ""),
                    "taxonomy_order":      tax.get("taxonomy_order", ""),
                    "taxonomy_best":       tax.get("taxonomy_best", ""),
                    "taxonomy_source":     tax.get("taxonomy_source", ""),
                    "host_bin":            host.get("host_bin", ""),
                    "host_score":          host.get("host_score", ""),
                })
                total_rows += 1

    n_clustered = sum(1 for rep in reps_ordered if len(rep_members.get(rep, [])) > 1)
    n_singletons = len(reps_ordered) - n_clustered
    n_lytic  = sum(1 for r in reps_ordered if lifestyle.get(r, "unknown") == "lytic")
    n_lyso   = sum(1 for r in reps_ordered if lifestyle.get(r, "unknown") == "lysogenic")
    n_amg    = sum(1 for r in reps_ordered if amg_counts.get(r, 0) > 0)

    print(f"[make_votu_table] Written {total_rows} rows "
          f"({len(reps_ordered)} vOTUs, {n_clustered} multi-member, "
          f"{n_singletons} singletons)", file=sys.stderr)
    print(f"[make_votu_table]   Lytic: {n_lytic}  Lysogenic: {n_lyso}  "
          f"With AMGs: {n_amg}", file=sys.stderr)


if __name__ == "__main__":
    main()
