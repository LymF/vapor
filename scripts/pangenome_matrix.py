#!/usr/bin/env python3
"""Matriz gene x membro dos clusters elegiveis, em TRES estados.

Ausencia em MAG nao e ausencia no organismo. Um "." pode significar "o
organismo nao tem o gene" ou "o MAG esta 74% completo e a regiao nao
montou"; tratar os dois como o mesmo zero produz conclusao errada com
aparencia de dado -- a mesma familia do done.txt vazio lido como zero
biologico (docs/ROADMAP_SIMPLIFICACAO.md).

Por isso: 'x' presente, '.' ausente, '?' NAO AVALIAVEL. E o '?' fica fora
do denominador da frequencia.
"""

MIN_COMPLETENESS = 70.0   # mesmo piso do mag_bakta
CORE_FRACTION = 0.90      # 99% zera o core com MAG (metaFun)


def _evaluable(members, completeness, min_completeness):
    return [m for m in members
            if completeness.get(m, 0.0) >= min_completeness]


def build_matrix(clusters, members_by_rep, gene_hits, completeness,
                 min_completeness=MIN_COMPLETENESS):
    """Uma linha por (cluster, gene)."""
    rows = []
    for rep in clusters:
        members = members_by_rep.get(rep, [])
        evaluable = set(_evaluable(members, completeness, min_completeness))
        genes = sorted({g for m in members for g in gene_hits.get(m, set())})
        for gene in genes:
            states, n_present = {}, 0
            for m in members:
                if m not in evaluable:
                    states[m] = "?"
                elif gene in gene_hits.get(m, set()):
                    states[m] = "x"
                    n_present += 1
                else:
                    states[m] = "."
            rows.append({
                "representative_id": rep,
                "gene": gene,
                "states": states,
                "n_present": n_present,
                "n_evaluable": len(evaluable),
                "freq": f"{n_present}/{len(evaluable)}" if evaluable else "0/0",
            })
    return rows


def summarize_clusters(matrix_rows, members_by_rep, completeness,
                       min_completeness=MIN_COMPLETENESS):
    """Uma linha por cluster: e isto que decide se a fase 2 se justifica."""
    by_rep = {}
    for row in matrix_rows:
        by_rep.setdefault(row["representative_id"], []).append(row)

    summaries = []
    for rep, rows in sorted(by_rep.items()):
        members = members_by_rep.get(rep, [])
        n_eval = len(_evaluable(members, completeness, min_completeness))
        core = variable = singleton = 0
        for row in rows:
            if not n_eval or row["n_present"] == 0:
                continue
            if row["n_present"] >= CORE_FRACTION * n_eval:
                core += 1
            else:
                variable += 1
            if row["n_present"] == 1:
                singleton += 1
        summaries.append({
            "representative_id":    rep,
            "n_members":            len(members),
            "n_members_avaliaveis": n_eval,
            "n_genes_core":         core,
            "n_genes_variaveis":    variable,
            "n_genes_singleton":    singleton,
        })
    return summaries
