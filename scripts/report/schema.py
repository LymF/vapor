"""Contrato de dados do report: um grafico so consome campo declarado aqui.

Existe por causa de um bug real: um loader copiava a linha crua do
classify_wf do GTDB-Tk (`base = dict(gtdb_bins[key])`), embarcando trinta
campos para consumir sete -- entre eles other_related_references, que e uma
string enorme por genoma. A projecao explicita mata essa classe de bug.
"""
from dataclasses import dataclass
import json


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


class PayloadOverBudget(Exception):
    """O JSON embarcado passou do orcamento."""


def payload_report(data):
    # separators compactos: e o mesmo dumps que _data_script usa para
    # escrever o HTML de verdade (renderer_v2.py) -- medir com os espacos
    # default do json.dumps superestimaria o orcamento em relacao ao que
    # realmente vai para o disco.
    sizes = [(name, len(json.dumps(obj, ensure_ascii=False,
                                    separators=(',', ':')).encode('utf-8')))
             for name, obj in data.items()]
    return sorted(sizes, key=lambda item: item[1], reverse=True)


def check_budget(data, limit_mb=25.0):
    sizes = payload_report(data)
    total = sum(size for _, size in sizes)
    if total > limit_mb * 1024 * 1024:
        piores = ", ".join(f"{n} ({s / 1024 / 1024:.1f} MB)" for n, s in sizes[:3])
        raise PayloadOverBudget(
            f"payload de {total / 1024 / 1024:.1f} MB excede o orcamento de "
            f"{limit_mb} MB. Maiores blocos: {piores}. "
            f"Projete os campos com schema.project antes de embarcar."
        )
    return sizes
