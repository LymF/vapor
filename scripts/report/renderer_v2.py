"""renderer_v2.py — costura JSON + bundle React + CSS num HTML standalone.

Quatro passos e nada mais: montar o JSON, ler o bundle, ler o CSS, escrever.
Toda a logica de grafico vive em src/report-ui/, compilada em
scripts/report/assets/report-ui.js.
"""
import csv
import json
import os

from .data_loaders import (
    load_tool_status, parse_fasta_lengths, parse_quast_all,
    parse_total_reads, safe_int,
)
from .schema import check_budget

_HERE = os.path.dirname(__file__)
_ASSETS = os.path.join(_HERE, "assets")
_COMP = os.path.join(_HERE, "components")

TODAS = "__all__"

# checkv_quality vazio significa "CheckV nunca avaliou este contig", que NAO e o
# mesmo que "avaliado e ruim". viral_length_gate.format_discard_row preserva
# essa diferenca de proposito; aqui ela vira um rotulo proprio.
SEM_AVALIACAO = "sem avaliação CheckV"


def _read(path):
    with open(path, encoding='utf-8') as f:
        return f.read()


def _data_script(data):
    # O escape de "</" impede que uma string do dado feche o <script> que a
    # carrega -- e o mesmo cuidado do _jsstr do renderer antigo.
    payload = json.dumps(data, ensure_ascii=False).replace("</", "<\\/")
    return f"<script>window.VAPOR_DATA = {payload};</script>"


def render_html(data, assets_dir=_ASSETS, comp_dir=_COMP):
    shell = _read(os.path.join(comp_dir, "shell_v2.html"))
    return (shell
            .replace("{{CSS}}", _read(os.path.join(assets_dir, "report-ui.css")))
            .replace("{{DATA_JSON}}", _data_script(data))
            .replace("{{APP_JS}}", _read(os.path.join(assets_dir, "report-ui.js"))))


def write_report(data, out_html, assets_dir=_ASSETS, comp_dir=_COMP, limit_mb=25.0):
    check_budget(data, limit_mb=limit_mb)
    html = render_html(data, assets_dir, comp_dir)
    os.makedirs(os.path.dirname(out_html) or '.', exist_ok=True)
    with open(out_html, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"[VAPOR] Report (v2) escrito em {out_html}")
    return out_html


def _quebra_por_tier(tsv_path):
    # Descartes do portao composto (item (e)) agrupados por tier do CheckV.
    #
    # NAO e uma quebra por "motivo": toda linha desse arquivo falhou as tres
    # armas do portao ao mesmo tempo (sem bin do vRhyme, tier abaixo de MQ, e
    # comprimento abaixo de VIRAL_MIN_CONTIG). Eleger uma das armas como causa
    # seria inventar informacao que o dado nao tem.
    if not os.path.exists(tsv_path):
        return []
    contagem = {}
    with open(tsv_path, encoding='utf-8', newline='') as fh:
        for linha in csv.DictReader(fh, delimiter='\t'):
            tier = (linha.get("checkv_quality") or "").strip() or SEM_AVALIACAO
            contagem[tier] = contagem.get(tier, 0) + 1
    return [{"reason": t, "count": n}
            for t, n in sorted(contagem.items(), key=lambda kv: -kv[1])]


def _etapas(contagens):
    # Etapas na ordem do funil, omitindo aquela cuja fonte nao existe. Zero aqui
    # significa "nao consegui ler a fonte", nao "zero biologico" -- desenhar uma
    # barra zerada afirmaria o segundo. A etapa some.
    ordem = ["reads", "contigs", "candidatos virais", "vOTUs retidos"]
    return [{"name": nome, "value": contagens[nome]}
            for nome in ordem
            if contagens.get(nome)]


def _conta_fasta(caminho):
    return len(parse_fasta_lengths(caminho)) if os.path.exists(caminho) else 0


def _funil_da_amostra(outdir, sample):
    quast = parse_quast_all(os.path.join(outdir, sample, "quast", "report.tsv"))
    n_contigs = safe_int((quast or {}).get("# contigs", 0))
    descartado = os.path.join(outdir, sample, "final", "viral", "viral_discarded.tsv")
    return {
        "stages": _etapas({
            "reads": parse_total_reads(outdir, sample),
            "contigs": n_contigs,
            "candidatos virais": _conta_fasta(os.path.join(
                outdir, sample, "viral", "consensus",
                f"{sample}_viral_consensus.fasta")),
            "vOTUs retidos": _conta_fasta(os.path.join(
                outdir, sample, "final", "viral", "viral_nonredundant.fasta")),
        }),
        "losses": {"vOTUs retidos": _quebra_por_tier(descartado)},
    }


def _funil_agregado(outdir, samples):
    por_amostra = [_funil_da_amostra(outdir, s) for s in samples]
    soma_etapas, perdas = {}, {}
    for f in por_amostra:
        for etapa in f["stages"]:
            soma_etapas[etapa["name"]] = soma_etapas.get(etapa["name"], 0) + etapa["value"]
        for motivo in f["losses"].get("vOTUs retidos", []):
            perdas[motivo["reason"]] = perdas.get(motivo["reason"], 0) + motivo["count"]
    return {
        "stages": _etapas(soma_etapas),
        "losses": {"vOTUs retidos": [
            {"reason": r, "count": n}
            for r, n in sorted(perdas.items(), key=lambda kv: -kv[1])]},
    }


def build_data(snakemake):
    """Monta o dicionario do report a partir do que ja existe em disco.

    So o bloco 'overview' nesta fase; as demais abas entram no plano 2.
    """
    outdir = snakemake.params.outdir
    samples = list(snakemake.params.samples)
    grupos = list(getattr(snakemake.params, 'coassembly_groups', []) or [])

    # load_tool_status devolve {sample: {tool: {state, reason, raw}}}, com as
    # pseudo-amostras "(global)" e "(coassembly) <grupo>". A StatusMatrix quer
    # linhas, entao achatamos aqui -- e 'raw' NAO atravessa: e o conteudo bruto
    # do done.txt, que nenhum componente le.
    status = [
        {"rule": tool, "sample": unidade,
         "status": entrada.get("state", "unknown"),
         "reason": entrada.get("reason", "")}
        for unidade, ferramentas in load_tool_status(outdir, samples, grupos).items()
        for tool, entrada in sorted(ferramentas.items())
    ]

    contigs = {}
    for s in samples:
        caminho = os.path.join(outdir, s, "final", "viral", "viral_nonredundant.fasta")
        contigs[s] = len(parse_fasta_lengths(caminho)) if os.path.exists(caminho) else 0

    kpis = [
        {"label": "Amostras", "value": len(samples)},
        {"label": "Grupos", "value": len(grupos)},
        {"label": "vOTUs retidos", "value": sum(contigs.values())},
    ]

    funil = {TODAS: _funil_agregado(outdir, samples)}
    for s in samples:
        funil[s] = _funil_da_amostra(outdir, s)

    return {
        "run": {"title": "VAPOR", "samples": samples, "groups": grupos},
        "overview": {"kpis": kpis, "status": status, "funnel": funil},
    }
