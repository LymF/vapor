#!/usr/bin/env python3
"""Portao do pangenoma: quais clusters de MAG merecem anotacao por membro.

O principio (h) (docs/ROADMAP_SIMPLIFICACAO.md) computa na representante e
herda no membro. Comparar conteudo genico DENTRO de um cluster exige o
contrario, e por isso este portao existe: reintroduz anotacao por membro
apenas onde a pergunta de core/acessorio tem sentido -- dezenas de
genomas, nao todos.
"""

import csv

MIN_MEMBERS = 3
MIN_SYSTEMS = 3

# Ordem = precedencia na hora de nomear o criterio disparado. A ilha e a
# evidencia mais forte, o AMR entra por direito proprio, e o plasmidio NAO
# aparece aqui: e sinal de mobilidade que reforca, nunca criterio isolado.
_CRITERIA = (
    ("ilha",      lambda e: e.get("n_islands", 0) >= 1),
    ("sistemas",  lambda e: e.get("n_systems", 0) >= MIN_SYSTEMS),
    ("amr",       lambda e: e.get("n_args", 0) >= 1),
)


def load_membership(path):
    """{representative_id: [{'source_id','original_bin_id','member_id'}, ...]}"""
    by_rep = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rep = (row.get("representative_id") or "").strip()
            member = (row.get("member_id") or "").strip()
            if not rep or not member:
                continue
            by_rep.setdefault(rep, []).append({
                "source_id": (row.get("source_id") or "").strip(),
                "original_bin_id": (row.get("original_bin_id") or "").strip(),
                "member_id": member,
            })
    return by_rep


def load_completeness(path):
    """{genome: completeness} do checkm2_quality_report.tsv do catalogo."""
    out = {}
    try:
        with open(path, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                name = (row.get("Name") or "").strip()
                if not name:
                    continue
                try:
                    out[name] = float(row.get("Completeness") or 0)
                except ValueError:
                    continue
    except OSError:
        pass
    return out


def select_clusters(membership, evidence, min_members=MIN_MEMBERS):
    """Uma linha por cluster, elegivel ou nao -- a selecao tem de ser
    auditavel. `criterio` nomeia o primeiro criterio satisfeito, ou o
    motivo da recusa."""
    rows = []
    for rep, members in sorted(membership.items()):
        ev = evidence.get(rep, {})
        n_members = len(members)
        criterio, eligible = "", False
        if n_members < min_members:
            criterio = f"poucos membros ({n_members} < {min_members})"
        else:
            for name, test in _CRITERIA:
                if test(ev):
                    criterio, eligible = name, True
                    break
            if not eligible:
                criterio = "sem evidencia de defesa/amr"
        rows.append({
            "representative_id": rep,
            "n_members":  n_members,
            "n_islands":  ev.get("n_islands", 0),
            "n_systems":  ev.get("n_systems", 0),
            "n_args":     ev.get("n_args", 0),
            "n_plasmid":  ev.get("n_plasmid", 0),
            "criterio":   criterio,
            "eligible":   eligible,
        })
    return rows
