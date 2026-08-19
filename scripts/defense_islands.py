#!/usr/bin/env python3
"""Deteccao de ilha de defesa -- nucleo puro, compartilhado.

Ate 2026-08-19 isto vivia so em scripts/report/data_loaders.py, ou seja, na
camada de RELATORIO: a pipeline nao alcancava. O portao do pangenoma
(rules/pangenome.smk) precisa da mesma definicao, e duas copias divergiriam
em silencio -- a familia de defeito que docs/ROADMAP_SIMPLIFICACAO.md
persegue. Uma definicao, dois consumidores.
"""

import os
from collections import defaultdict


def _contig_from_protein_id(pid):
    """'k141_219139_5' -> 'k141_219139'. O Prodigal numera '{seqid}_{n}'."""
    if not pid or "_" not in pid:
        return ""
    return pid.rsplit("_", 1)[0]


def genes_by_contig(faa_path):
    """Genes por contig em ordem genomica, com coordenadas e fita.

    O Prodigal codifica o locus no proprio header do FASTA:
        >{seqid}_{n} # {start} # {end} # {strand} # ID=...;partial=...
    entao coordenadas reais saem do .faa que a pipeline ja escreve.
    """
    by_contig = defaultdict(list)
    if not faa_path or not os.path.exists(faa_path):
        return by_contig
    with open(faa_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\n")
            pid = header.split()[0].strip()
            contig = _contig_from_protein_id(pid)
            if not contig:
                continue
            start = end = strand = None
            parts = [p.strip() for p in header.split("#")]
            if len(parts) >= 4:
                try:
                    start, end, strand = int(parts[1]), int(parts[2]), int(parts[3])
                except (ValueError, IndexError):
                    start = end = strand = None
            by_contig[contig].append({"Protein": pid, "Start": start,
                                      "End": end, "Strand": strand})
    return by_contig


def find_islands(genes_by_contig_map, prot_to_sys,
                 min_genes=5, min_systems=3, window=10):
    """Uma entrada por ilha: corrida de genes de defesa no mesmo contig onde
    posicoes consecutivas nunca distam mais que `window` genes, com
    >= min_genes genes de >= min_systems sistemas distintos.

    `prot_to_sys`: {protein_id: (bin, system, system_id)}.
    """
    islands = []
    for contig, gene_recs in genes_by_contig_map.items():
        ordered = [g["Protein"] for g in gene_recs]
        hits = [(i, p, prot_to_sys[p]) for i, p in enumerate(ordered)
                if p in prot_to_sys]
        if len(hits) < min_genes:
            continue

        def flush(cluster):
            if len(cluster) < min_genes:
                return
            systems = {h[2][1] for h in cluster if h[2][1]}
            if len(systems) < min_systems:
                return
            start_idx, end_idx = cluster[0][0], cluster[-1][0]
            defense_idx = {h[0]: h[2][1] for h in cluster}
            sysid_idx = {h[0]: h[2][2] for h in cluster}
            window_genes = []
            for i in range(start_idx, end_idx + 1):
                rec = gene_recs[i]
                window_genes.append({
                    "Protein": rec["Protein"], "Index": i,
                    "System": defense_idx.get(i, ""),
                    "System_id": sysid_idx.get(i, ""),
                    "Start": rec.get("Start"), "End": rec.get("End"),
                    "Strand": rec.get("Strand"),
                })
            coords = [g for g in window_genes
                      if g["Start"] is not None and g["End"] is not None]
            islands.append({
                "Contig": contig,
                "n_genes": len(cluster), "n_systems": len(systems),
                "Systems": sorted(systems), "window_genes": window_genes,
                "start_idx": start_idx, "end_idx": end_idx,
                "Start_bp": min(g["Start"] for g in coords) if coords else None,
                "End_bp": max(g["End"] for g in coords) if coords else None,
            })

        cluster = []
        for h in hits:
            if cluster and h[0] - cluster[-1][0] > window:
                flush(cluster)
                cluster = []
            cluster.append(h)
        flush(cluster)
    return islands
