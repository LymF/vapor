"""renderer_v2.py — costura JSON + bundle React + CSS num HTML standalone.

Quatro passos e nada mais: montar o JSON, ler o bundle, ler o CSS, escrever.
Toda a logica de grafico vive em src/report-ui/, compilada em
scripts/report/assets/report-ui.js.
"""
import csv
import json
import os

from .data_loaders import (
    load_tool_status, parse_fasta_lengths, parse_quast_all, safe_int,
    parse_fastp_json, parse_mapping_rate, collect_depth_data, parse_tsv,
    parse_support_combos, load_viral_taxonomy, load_votu_catalog,
    load_votu_presence, load_votu_lifestyle,
)
from .schema import Block, project, check_budget

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


ETAPA_FINAL = "sequências virais retidas"

# Unidade de cada etapa. A barra de perda so faz sentido entre etapas na MESMA
# unidade -- o vRhyme colapsa N contigs num vMAG, entao "candidatos virais" ->
# "sequências virais retidas" muda de unidade (contig -> sequencia) e a
# diferenca embute consolidacao, nao descarte.
_UNIDADES = {
    "contigs": "contig",
    "candidatos virais": "contig",
    ETAPA_FINAL: "sequência",
}


def _etapas(contagens):
    # Etapas na ordem do funil. contagens[nome] == None significa "fonte nao
    # existe" e a etapa some; contagens[nome] == 0 significa "fonte existe e a
    # contagem e genuinamente zero" e a etapa aparece com valor 0. Um dict
    # tradicional (sem a chave) tambem e tratado como fonte ausente.
    ordem = ["contigs", "candidatos virais", ETAPA_FINAL]
    return [{"name": nome, "value": contagens[nome], "unit": _UNIDADES[nome]}
            for nome in ordem
            if contagens.get(nome) is not None]


def _conta_fasta(caminho):
    return len(parse_fasta_lengths(caminho)) if os.path.exists(caminho) else None


def _contagem_votus_catalogo(outdir):
    # vOTU e cluster do catalogo global (95% ANI / 85% AF) -- uma linha por
    # MEMBRO em vOTU_clusters.tsv, entao a contagem certa e votu_id DISTINTOS,
    # nao o numero de linhas nem a soma de viral_nonredundant.fasta entre
    # amostras (isso conta o mesmo virus varias vezes).
    caminho = os.path.join(outdir, "votu_catalog", "vOTU_clusters.tsv")
    if not os.path.exists(caminho):
        return None
    votus = set()
    with open(caminho, encoding='utf-8', newline='') as fh:
        for linha in csv.DictReader(fh, delimiter='\t'):
            vid = (linha.get("votu_id") or "").strip()
            if vid:
                votus.add(vid)
    return len(votus)


def _funil_da_amostra(outdir, sample):
    quast = parse_quast_all(os.path.join(outdir, sample, "quast", "report.tsv"))
    # parse_quast_all devolve {rotulo_do_assembly: {metrica: valor}} -- a
    # metrica fica um nivel abaixo do dict que report.tsv indexa por rotulo.
    # "assembly" e o rotulo desde o item (d) (um so assembly por amostra); o
    # fallback para o primeiro rotulo tolera um report.tsv de rodada antiga.
    qd = quast.get("assembly") or next(iter(quast.values()), {}) if quast else {}
    bruto = qd.get("# contigs")
    n_contigs = safe_int(bruto) if bruto is not None else None
    descartado = os.path.join(outdir, sample, "final", "viral", "viral_discarded.tsv")
    return {
        "stages": _etapas({
            "contigs": n_contigs,
            "candidatos virais": _conta_fasta(os.path.join(
                outdir, sample, "viral", "consensus",
                f"{sample}_viral_consensus.fasta")),
            ETAPA_FINAL: _conta_fasta(os.path.join(
                outdir, sample, "final", "viral", "viral_nonredundant.fasta")),
        }),
        "losses": {ETAPA_FINAL: _quebra_por_tier(descartado)},
    }


def _funil_agregado(outdir, samples):
    por_amostra = [_funil_da_amostra(outdir, s) for s in samples]
    soma_etapas, perdas = {}, {}
    for f in por_amostra:
        for etapa in f["stages"]:
            soma_etapas[etapa["name"]] = soma_etapas.get(etapa["name"], 0) + etapa["value"]
        for motivo in f["losses"].get(ETAPA_FINAL, []):
            perdas[motivo["reason"]] = perdas.get(motivo["reason"], 0) + motivo["count"]
    return {
        "stages": _etapas(soma_etapas),
        "losses": {ETAPA_FINAL: [
            {"reason": r, "count": n}
            for r, n in sorted(perdas.items(), key=lambda kv: -kv[1])]},
    }


# ── Aba Sequenciamento ───────────────────────────────────────────────────────

QC_BLOCK = Block(name="qc", fields=("sample", "reads_before", "reads_after", "q30"))
QUAST_BLOCK = Block(name="quast", fields=(
    "# contigs", "Largest contig", "Total length", "GC (%)",
    "N50", "N75", "L50", "L75",
))
MAPPING_BLOCK = Block(name="mapping", fields=("sample", "rate"))


def _q30_pct(reads):
    # parse_fastp_json nao devolve q30 isolado: ele ja o embutiu na formula
    # mean_quality = 10 + 25*q30 (q30 como fracao 0-1) para cada bloco
    # before/after -- ver o comentario em parse_fastp_json. Invertemos essa
    # mesma formula em vez de reabrir o JSON do fastp aqui.
    for r in reads:
        if r.get("stage") == "trimmed":
            return round((r.get("mean_quality", 10.0) - 10.0) / 25.0 * 100, 2)
    return None


def build_sequencing(outdir, samples):
    saida = {}

    qc_rows = []
    for s in samples:
        fastp_path = os.path.join(outdir, s, "qc_raw", f"{s}_fastp.json")
        if not os.path.exists(fastp_path):
            continue
        parsed = parse_fastp_json(outdir, s)
        trim = parsed["trim"]
        qc_rows.append({
            "sample": s,
            "reads_before": trim["reads_in"],
            "reads_after": trim["reads_written"],
            "q30": _q30_pct(parsed["reads"]),
        })
    if qc_rows:
        saida["qc"] = project(QC_BLOCK, qc_rows)

    quast = {}
    for s in samples:
        report_path = os.path.join(outdir, s, "quast", "report.tsv")
        if not os.path.exists(report_path):
            continue
        parsed = parse_quast_all(report_path)
        # ver comentario em _funil_da_amostra: a metrica fica um nivel abaixo
        # do rotulo do assembly.
        qd = parsed.get("assembly") or next(iter(parsed.values()), {}) if parsed else {}
        quast[s] = project(QUAST_BLOCK, [qd])[0]
    if quast:
        saida["quast"] = quast

    mapping_rows = []
    for s in samples:
        flagstat_path = os.path.join(outdir, s, "mapping", "flagstat.txt")
        if not os.path.exists(flagstat_path):
            continue
        mapping_rows.append({"sample": s, "rate": parse_mapping_rate(outdir, s)})
    if mapping_rows:
        saida["mapping"] = project(MAPPING_BLOCK, mapping_rows)

    lengths = {}
    depth = {}
    for s in samples:
        depth_path = os.path.join(outdir, s, "mapping", f"{s}_depth.txt")
        if not os.path.exists(depth_path):
            continue
        comprimentos, profundidades = collect_depth_data(depth_path)
        if not comprimentos:
            continue
        lengths[s] = comprimentos
        depth[s] = [list(par) for par in zip(comprimentos, profundidades)]
    if lengths:
        saida["lengths"] = lengths
    if depth:
        saida["depth"] = depth

    return saida


# ── Aba Catálogo viral ───────────────────────────────────────────────────────

TAXONOMY_BLOCK = Block(name="viral_taxonomy", fields=(
    "sample", "Phylum", "Class", "Order", "Family", "Genus", "count",
))
EXPLORER_FEATURE_BLOCK = Block(
    name="explorer_feature", fields=("start", "end", "strand", "label", "kind"))


def _checkv_tier_counts(quality_summary_path):
    # A chave vazia ('') e uma categoria de verdade -- "CheckV nao avaliou
    # este contig" -- e precisa sobreviver distinta de qualquer tier baixo.
    # Ao contrario de _quebra_por_tier (funil), aqui NAO trocamos por
    # SEM_AVALIACAO: o contrato desta aba pede a chave vazia literal.
    if not os.path.exists(quality_summary_path):
        return None
    contagem = {}
    for row in parse_tsv(quality_summary_path):
        tier = (row.get("checkv_quality") or "").strip()
        contagem[tier] = contagem.get(tier, 0) + 1
    return contagem


def _build_taxonomy(outdir, samples):
    caminhos = [os.path.join(outdir, s, "viral", "taxonomy",
                              "viral_taxonomy_merged.tsv") for s in samples]
    if not any(os.path.exists(p) for p in caminhos):
        return None
    registros = load_viral_taxonomy(caminhos, samples)
    if not registros:
        return []
    contagem = {}
    for r in registros:
        chave = (r.get("sample", ""), r.get("Order", ""),
                 r.get("Family", ""), r.get("Genus", ""))
        contagem[chave] = contagem.get(chave, 0) + 1
    linhas = [
        {"sample": s, "Phylum": "", "Class": "", "Order": o,
         "Family": f, "Genus": g, "count": n}
        for (s, o, f, g), n in contagem.items()
    ]
    return project(TAXONOMY_BLOCK, linhas)


def _build_detectors(outdir, samples):
    combos_totais = {}
    algum_arquivo = False
    for s in samples:
        caminho = os.path.join(outdir, s, "viral", "consensus",
                                f"{s}_tool_support.tsv")
        if not os.path.exists(caminho):
            continue
        algum_arquivo = True
        for chave, n in parse_support_combos(caminho).items():
            combos_totais[chave] = combos_totais.get(chave, 0) + n
    if not algum_arquivo:
        return None

    sets = {}
    combos = []
    for chave, n in sorted(combos_totais.items(), key=lambda kv: -kv[1]):
        ferramentas = [t for t in chave.split(",") if t]
        combos.append({"tools": ferramentas, "count": n})
        for t in ferramentas:
            sets[t] = sets.get(t, 0) + n
    return {"sets": sets, "combos": combos}


def _pharokka_rows(outdir):
    caminho = os.path.join(outdir, "votu_catalog", "annotation", "pharokka",
                            "pharokka_cds_final_merged_output.tsv")
    if not os.path.exists(caminho):
        return None
    return parse_tsv(caminho)


def _limita_explorer(linhas):
    return sorted(linhas, key=lambda r: r["length"], reverse=True)[:50]


def _build_explorer(outdir):
    linhas_pharokka = _pharokka_rows(outdir)
    if linhas_pharokka is None:
        return None

    por_contig = {}
    for row in linhas_pharokka:
        contig = row.get("contig") or row.get("contig_id") or ""
        if not contig:
            continue
        try:
            start = int(float(row.get("start", "") or 0))
            end = int(float(row.get("stop", row.get("end", "")) or 0))
        except (TypeError, ValueError):
            continue
        strand_bruta = (row.get("strand", "") or "").strip()
        strand = -1 if strand_bruta in ("-", "-1") else 1
        feat = {
            "start": start,
            "end": end,
            "strand": strand,
            "label": row.get("gene") or row.get("product") or row.get("top_hit") or "",
            "kind": row.get("phrog_category") or row.get("category") or "",
        }
        por_contig.setdefault(contig, []).append(feat)

    linhas = [
        {"votu_id": contig,
         "length": max((f["end"] for f in feats), default=0),
         "features": project(EXPLORER_FEATURE_BLOCK, feats)}
        for contig, feats in por_contig.items()
    ]
    return _limita_explorer(linhas)


def build_viral(outdir, samples):
    saida = {}

    taxonomia = _build_taxonomy(outdir, samples)
    if taxonomia is not None:
        saida["taxonomy"] = taxonomia

    checkv_tiers = {}
    for s in samples:
        caminho = os.path.join(outdir, s, "viral", "checkv", "quality_summary.tsv")
        contagem = _checkv_tier_counts(caminho)
        if contagem is not None:
            checkv_tiers[s] = contagem
    if checkv_tiers:
        saida["checkv_tiers"] = checkv_tiers

    detectores = _build_detectors(outdir, samples)
    if detectores is not None:
        saida["detectors"] = detectores

    if os.path.exists(os.path.join(outdir, "votu_catalog", "vOTU_clusters.tsv")):
        saida["catalog"] = load_votu_catalog(outdir)

    if os.path.exists(os.path.join(outdir, "votu_catalog", "presence_matrix.tsv")):
        saida["presence"] = load_votu_presence(outdir, samples)

    if os.path.exists(os.path.join(
            outdir, "votu_catalog", "bacphlip", "votu_lifestyle.tsv")):
        saida["lifestyle"] = load_votu_lifestyle(outdir)

    explorer = _build_explorer(outdir)
    if explorer is not None:
        saida["explorer"] = explorer

    return saida


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

    kpis = [
        {"label": "Amostras", "value": len(samples)},
        {"label": "Grupos", "value": len(grupos)},
    ]
    n_votus = _contagem_votus_catalogo(outdir)
    if n_votus is not None:
        kpis.append({"label": "vOTUs no catálogo", "value": n_votus})

    funil = {TODAS: _funil_agregado(outdir, samples)}
    for s in samples:
        funil[s] = _funil_da_amostra(outdir, s)

    dados = {
        "run": {"title": "VAPOR", "samples": samples, "groups": grupos},
        "overview": {"kpis": kpis, "status": status, "funnel": funil},
    }

    sequencing = build_sequencing(outdir, samples)
    if sequencing:
        dados["sequencing"] = sequencing

    viral = build_viral(outdir, samples)
    if viral:
        dados["viral"] = viral

    return dados
