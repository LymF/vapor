#!/usr/bin/env python3
"""
make_votu_table.py — Build the vOTU membership table for one sample.

One row per catalog vOTU that has at least one member from this sample.
Per-member annotations (CheckV quality, taxonomy, host) come
from this sample's own tables, keyed by that member's bare contig ID.

vOTU identity (votu_id, representative) is GLOBAL: the catalog in
rules/votu_catalog.smk clusters pooled viral sets across every sample
and co-assembly group, so a vOTU's representative may legitimately
belong to a different sample than the one this table is being built
for. `vOTU_clusters.tsv` therefore namespaces every ID as
"{source_id}|{contig_id}". This script:
  - keeps votu_id and representative namespaced (global identity, not
    stripped — stripping the representative would silently claim it
    belongs to this sample when it may not)
  - filters cluster members down to this sample's own contigs only
  - strips the "{sample}|" prefix off of SURVIVING members only, so the
    resulting bare IDs match the per-sample CheckV/taxonomy/
    PHIST tables, which were never namespaced
  - drops a vOTU entirely if it has no member in this sample

Inputs (wired via Snakemake):
  votu_clusters   votu_catalog/vOTU_clusters.tsv     (votu_id, representative, member — namespaced)
  votu_reps       votu_catalog/catalog_all_reps.fasta (FASTA of all representatives, namespaced headers)
  checkv_tsv      viral/checkv/quality_summary.tsv    (this sample, bare IDs)
  taxonomy_tsv    viral/taxonomy/viral_taxonomy_merged.tsv (this sample, bare IDs)
  phist_csv       viral/phist/phist_results.csv       (this sample, bare IDs)
  lifestyle_tsv   votu_catalog/bacphlip/votu_lifestyle.tsv (global, keyed by representative)
  amg_tsv         votu_catalog/eggnog_viral/putative_amgs.tsv (global, keyed by representative)

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


def load_votu_clusters(clusters_tsv, sample):
    """
    Read the global catalog's vOTU_clusters.tsv and filter it to this sample.

    The catalog's `member` column is namespaced ("{source_id}|{contig_id}")
    because it pools every sample and co-assembly group. Only members whose
    source prefix equals `sample` survive; their IDs are stripped back to
    bare contig IDs so they match this sample's own per-contig tables. The
    `representative` column is NEVER stripped: the representative of a vOTU
    may genuinely belong to a different sample, and flattening its prefix
    would silently misrepresent it as local.

    Returns:
      votu_order      — list of votu_id, in first-seen order, restricted to
                         vOTUs with >=1 member in this sample
      votu_rep        — {votu_id: representative}  (namespaced, as-is)
      votu_members    — {votu_id: [bare_member_id, ...]}  (this sample only)
      votu_total_size — {votu_id: n}  total member count across the WHOLE
                         catalog (all samples/groups), for context
    """
    prefix = f"{sample}|"
    votu_order      = []
    votu_rep        = {}
    votu_members    = defaultdict(list)
    votu_total_size = defaultdict(int)
    seen_votus      = set()
    all_votus       = set()
    n_rows          = 0

    for row in read_tsv(clusters_tsv, required=False):
        votu_id = row.get("votu_id", "").strip()
        rep     = row.get("representative", "").strip()
        mem     = row.get("member", "").strip()
        if not votu_id or not rep or not mem:
            continue
        n_rows += 1
        all_votus.add(votu_id)
        votu_total_size[votu_id] += 1

        if not mem.startswith(prefix):
            continue
        if votu_id not in seen_votus:
            votu_order.append(votu_id)
            seen_votus.add(votu_id)
            votu_rep[votu_id] = rep
        votu_members[votu_id].append(mem[len(prefix):])

    n_members_kept = sum(len(v) for v in votu_members.values())
    print(f"[make_votu_table] Catalog: {len(all_votus)} vOTUs, {n_rows} total members "
          f"(all samples/groups)", file=sys.stderr)
    print(f"[make_votu_table] Sample '{sample}' after filtering: "
          f"{len(votu_order)} vOTUs, {n_members_kept} members survive", file=sys.stderr)
    if n_rows and not votu_order:
        print(f"[make_votu_table] WARNING: 0 vOTUs matched prefix '{prefix}' — "
              f"check that sample names in vOTU_clusters.tsv match wildcards.sample",
              file=sys.stderr)
    return votu_order, votu_rep, dict(votu_members), dict(votu_total_size)


def load_rep_lengths(votu_reps_fasta):
    """Return {rep_id: length_bp} from the all-reps FASTA (namespaced headers)."""
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


def load_lifestyle(lifestyle_tsv):
    """Returns {representative_id: lifestyle}. Keyed by the namespaced vOTU
    representative ID (BACPHLIP runs once per representative, not per
    member) — see rules/votu_catalog.smk:bacphlip_votu."""
    data = {}
    for row in read_tsv(lifestyle_tsv, required=False):
        rid = (row.get("votu_id") or "").strip()
        if not rid:
            continue
        data[rid] = row.get("lifestyle", "unknown")
    print(f"[make_votu_table] BACPHLIP lifestyle calls: {len(data)}", file=sys.stderr)
    return data


def load_amg_counts(amg_tsv):
    """Returns {representative_id: n_putative_AMGs}. Keyed by the namespaced
    vOTU representative ID (eggNOG runs once per MQ+ representative) —
    see rules/votu_catalog.smk:eggnog_viral."""
    counts = defaultdict(int)
    for row in read_tsv(amg_tsv, required=False):
        rid = (row.get("votu_id") or "").strip()
        if not rid:
            continue
        counts[rid] += 1
    print(f"[make_votu_table] Putative AMG-bearing representatives: {len(counts)}", file=sys.stderr)
    return counts


def load_phist(phist_csv):
    """Returns {virus_id: {host_bin, host_score}} — best (lowest adj-pvalue) hit.

    A chave e o ID do genoma viral como o PHIST o conhece: o nome do ARQUIVO
    que `split_viral_fastas.py` escreveu, sem o prefixo `contig_` que ele
    poe e sem a extensao. Esse prefixo NAO era removido aqui ate 2026-08-19
    (o loader do relatorio ja o removia), entao a chave era
    "contig_S1|k141_10" e nada casava: as colunas host_bin/host_score da
    tabela de vOTU saiam vazias em toda amostra, sem erro nenhum.

    Os IDs sao NAMESPACED ("{source_id}|{contig_id}"), porque desde
    2026-08-13 o PHIST roda sobre as representantes do catalogo global.
    """
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
        if virus.startswith("contig_"):
            virus = virus[len("contig_"):]
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
    taxonomy_tsv  = snakemake.input.taxonomy
    phist_csv     = snakemake.params.phist_csv
    lifestyle_tsv = snakemake.params.lifestyle_tsv
    amg_tsv       = snakemake.params.amg_tsv
    out_tsv       = snakemake.output.tsv
    sample        = snakemake.wildcards.sample

    print(f"[make_votu_table] Sample: {sample}", file=sys.stderr)

    votu_order, votu_rep, votu_members, votu_total_size = load_votu_clusters(
        clusters_tsv, sample)
    rep_lengths            = load_rep_lengths(votu_reps_fa)
    checkv                 = load_checkv(checkv_tsv)
    # Lifestyle (BACPHLIP) and putative AMGs (eggNOG) are computed once per
    # vOTU representative by the global catalog rules, not per member —
    # keyed here by the representative's namespaced ID, not the bare member
    # ID used for CheckV/taxonomy/PHIST below.
    lifestyle                = load_lifestyle(lifestyle_tsv)
    amg_counts               = load_amg_counts(amg_tsv)
    taxonomy                 = load_taxonomy(taxonomy_tsv)
    hosts                   = load_phist(phist_csv)

    COLS = [
        "votu_id", "representative", "member", "is_rep", "cluster_size", "sample",
        "rep_length_bp",
        "checkv_quality", "checkv_completeness", "checkv_length", "genome_type",
        "lifestyle", "n_putative_AMGs",
        "taxonomy_family", "taxonomy_genus", "taxonomy_order",
        "taxonomy_best", "taxonomy_source",
        "host_bin", "host_score",
    ]

    os.makedirs(os.path.dirname(out_tsv), exist_ok=True)
    total_rows = 0
    n_lytic = n_lyso = n_amg = 0
    with open(out_tsv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=COLS, delimiter="\t",
                                extrasaction="ignore")
        writer.writeheader()

        for votu_id in votu_order:
            rep          = votu_rep[votu_id]
            members      = votu_members.get(votu_id, [])
            cluster_size = votu_total_size.get(votu_id, len(members))
            rep_len      = rep_lengths.get(rep, "")

            for mem in members:
                # Annotations are looked up by the MEMBER's own bare ID, not
                # the representative's — the representative may belong to a
                # different sample and this script only has this sample's
                # per-contig tables. Using the member ID is also more
                # accurate: each row gets its own contig's data instead of
                # inheriting the representative's.
                cv   = checkv.get(mem, {})
                tax  = taxonomy.get(mem, {})
                # PHIST so pontuou as REPRESENTANTES (roda sobre
                # votu_catalog_reps), entao o hospedeiro se herda do
                # representante, como lifestyle e AMG -- procurar pelo ID do
                # membro deixaria sem hospedeiro todo membro que nao e
                # representante, que e a maioria. O MAG hospedeiro continua
                # sendo desta amostra: o PHIST roda contra os bins locais.
                host = (hosts.get(rep)
                        or hosts.get(f"{sample}|{mem}")
                        or hosts.get(mem)
                        or {})
                # Lifestyle/AMGs are per-representative, so every member of
                # a given vOTU reports the same value (the representative's).
                mem_lifestyle = lifestyle.get(rep, "unknown")
                mem_amgs      = amg_counts.get(rep, 0)
                is_rep = "True" if f"{sample}|{mem}" == rep else "False"

                writer.writerow({
                    "votu_id":             votu_id,
                    "representative":      rep,
                    "member":              mem,
                    "is_rep":              is_rep,
                    "cluster_size":        cluster_size,
                    "sample":              sample,
                    "rep_length_bp":       rep_len,
                    "checkv_quality":      cv.get("checkv_quality", ""),
                    "checkv_completeness": cv.get("checkv_completeness", ""),
                    "checkv_length":       cv.get("checkv_length", ""),
                    "genome_type":         cv.get("genome_type", ""),
                    "lifestyle":           mem_lifestyle,
                    "n_putative_AMGs":     mem_amgs,
                    "taxonomy_family":     tax.get("taxonomy_family", ""),
                    "taxonomy_genus":      tax.get("taxonomy_genus", ""),
                    "taxonomy_order":      tax.get("taxonomy_order", ""),
                    "taxonomy_best":       tax.get("taxonomy_best", ""),
                    "taxonomy_source":     tax.get("taxonomy_source", ""),
                    "host_bin":            host.get("host_bin", ""),
                    "host_score":          host.get("host_score", ""),
                })
                total_rows += 1
                if mem_lifestyle == "lytic":
                    n_lytic += 1
                elif mem_lifestyle == "lysogenic":
                    n_lyso += 1
                if mem_amgs > 0:
                    n_amg += 1

    n_multi_sample_votus = sum(1 for v in votu_order if votu_total_size.get(v, 1) > 1)

    print(f"[make_votu_table] Written {total_rows} rows "
          f"({len(votu_order)} vOTUs with a member in this sample, "
          f"{n_multi_sample_votus} belong to catalog clusters with >1 member overall)",
          file=sys.stderr)
    print(f"[make_votu_table]   Lytic: {n_lytic}  Lysogenic: {n_lyso}  "
          f"With AMGs: {n_amg}", file=sys.stderr)


if __name__ == "__main__":
    main()
