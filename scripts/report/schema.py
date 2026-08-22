"""Contrato de dados do report: um grafico so consome campo declarado aqui.

Existe por causa de um bug real: um loader copiava a linha crua do
classify_wf do GTDB-Tk (`base = dict(gtdb_bins[key])`), embarcando trinta
campos para consumir sete -- entre eles other_related_references, que e uma
string enorme por genoma. A projecao explicita mata essa classe de bug.
"""
from dataclasses import dataclass


class UndeclaredField(Exception):
    """Linha traz campo que o bloco nao declara."""


@dataclass(frozen=True)
class Block:
    name: str
    fields: tuple
    key: str = None


def project(block, rows):
    return [{f: row.get(f, '') for f in block.fields} for row in rows]


def project_strict(block, rows):
    declared = set(block.fields)
    for row in rows:
        extra = set(row) - declared
        if extra:
            raise UndeclaredField(
                f"bloco '{block.name}' recebeu campo(s) nao declarado(s): "
                f"{sorted(extra)}"
            )
    return project(block, rows)
