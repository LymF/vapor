"""renderer_v2.py — costura JSON + bundle React + CSS num HTML standalone.

Quatro passos e nada mais: montar o JSON, ler o bundle, ler o CSS, escrever.
Toda a logica de grafico vive em src/report-ui/, compilada em
scripts/report/assets/report-ui.js.
"""
import csv
import json
import math
import os
import re

from .data_loaders import (
    load_tool_status, parse_fasta_lengths, parse_quast_all, safe_int,
    parse_fastp_json, parse_mapping_rate, collect_depth_data, parse_tsv,
    parse_support_combos, load_viral_taxonomy, load_votu_catalog,
    load_votu_presence, load_votu_lifestyle, load_alpha_diversity,
    load_pcoord, load_reads_classify,
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
    # separators compactos: mesma escolha de report.schema.payload_report,
    # para que o orcamento medido corresponda ao que de fato vai pro disco.
    payload = json.dumps(data, ensure_ascii=False,
                          separators=(',', ':')).replace("</", "<\\/")
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


# O navegador so consome AGREGADOS de depth/lengths (DistPlot e uma grade 2D
# de densidade) -- nunca os pontos crus por contig. Numa rodada de 788 mil
# contigs isso sozinho era 24 MB de 25 MB de orcamento. Binamos no Python
# porque binagem usa TODO o dado (amostrar descartaria parte da medida).
#
# Abaixo de LIMIAR_BRUTO contigs por amostra os pontos continuam crus: com n
# pequeno a forma certa e um strip plot dos valores reais (REPORT_VIZ_GUIDE.md
# secao 4), e binar destruiria a unica informacao honesta que existe.
LIMIAR_BRUTO = 20
N_BINS_LENGTH = 60
N_BINS_DEPTH = 50


def _mediana(xs):
    ys = sorted(xs)
    n = len(ys)
    if n == 0:
        return None
    meio = n // 2
    return ys[meio] if n % 2 else (ys[meio - 1] + ys[meio]) / 2


def _bin_log1d(valores, n_bins):
    # Bins iguais em escala log10. Todo ponto cai em algum bin (o clamp nos
    # extremos e so para o arredondamento de ponto flutuante do ultimo valor
    # de cada amostra) -- e por isso soma(count) == len(valores) sempre, sem
    # excecao: e o teste que prova que binar nao descarta dado.
    lo, hi = min(valores), max(valores)
    if lo == hi:
        return [{"x0": lo, "x1": hi, "count": len(valores)}]
    log_lo, log_hi = math.log10(lo), math.log10(hi)
    largura = (log_hi - log_lo) / n_bins
    contagens = [0] * n_bins
    for v in valores:
        idx = int((math.log10(v) - log_lo) / largura)
        idx = max(0, min(idx, n_bins - 1))
        contagens[idx] += 1
    # Arredondar as bordas para inteiro (bp): sub-bp nao tem sentido biologico
    # e cada casa decimal a mais nos 60 bins e puro peso de JSON.
    return [
        {"x0": round(10 ** (log_lo + i * largura)),
         "x1": round(10 ** (log_lo + (i + 1) * largura)),
         "count": c}
        for i, c in enumerate(contagens)
    ]


def _sem_comprimento_invalido(comprimentos, profundidades=None):
    # Defesa contra fonte corrompida: comprimento <= 0 nao deveria existir
    # (um contig tem pelo menos 1 base), mas um depth.txt truncado ou
    # corrompido acontece, e math.log10(0)/math.log10(negativo) levanta
    # ValueError -- sem este filtro, UM contig ruim derrubava o build do
    # report inteiro em vez de so ficar de fora da distribuicao. Filtra os
    # dois arrays em paralelo (mesmo indice) quando profundidades e dado, para
    # nao desalinhar comprimento x profundidade em _build_depth_block.
    if profundidades is None:
        return [c for c in comprimentos if c > 0]
    pares = [(c, d) for c, d in zip(comprimentos, profundidades) if c > 0]
    return ([c for c, _ in pares], [d for _, d in pares])


def _build_length_block(comprimentos):
    comprimentos = _sem_comprimento_invalido(comprimentos)
    n = len(comprimentos)
    if n < LIMIAR_BRUTO:
        return {"values": list(comprimentos), "n": n}
    return {
        "bins": _bin_log1d(comprimentos, N_BINS_LENGTH),
        "n": n,
        "min": min(comprimentos),
        "max": max(comprimentos),
        "median": _mediana(comprimentos),
    }


def _build_depth_block(comprimentos, profundidades):
    # x = comprimento em log10 (sempre > 0); y = profundidade em log1p (a
    # jgi_summarize_bam_contig_depths emite 0.0 para contigs sem leitura
    # mapeada, e log10(0) nao existe -- log1p aceita zero sem inventar piso).
    comprimentos, profundidades = _sem_comprimento_invalido(comprimentos, profundidades)
    n = len(comprimentos)
    if n < LIMIAR_BRUTO:
        return {"values": [[l, d] for l, d in zip(comprimentos, profundidades)], "n": n}

    xs = [math.log10(l) for l in comprimentos]
    ys = [math.log1p(d) for d in profundidades]
    lo_x, hi_x = min(xs), max(xs)
    lo_y, hi_y = min(ys), max(ys)
    passo_x = ((hi_x - lo_x) or 1e-9) / N_BINS_DEPTH
    passo_y = ((hi_y - lo_y) or 1e-9) / N_BINS_DEPTH

    grade = {}
    for x, y in zip(xs, ys):
        ix = max(0, min(int((x - lo_x) / passo_x), N_BINS_DEPTH - 1))
        iy = max(0, min(int((y - lo_y) / passo_y), N_BINS_DEPTH - 1))
        grade[(ix, iy)] = grade.get((ix, iy), 0) + 1

    # Bins de contagem zero nao entram no dict acima -- e o que mantem a
    # grade 2D pequena mesmo com 50x50 = 2500 celulas possiveis. Cada bin
    # carrega so o INDICE da celula + contagem: repetir x0/x1/y0/y1 (4 bordas
    # de ponto flutuante) em cada uma das ate 2500 celulas era, medido na
    # rodada real, o proprio peso que a binagem deveria eliminar -- o cliente
    # reconstroi as bordas a partir de "grid" (compartilhado por toda a
    # amostra) com a MESMA transformacao (log10 / log1p) e o mesmo n_bins.
    bins2d = [{"ix": ix, "iy": iy, "count": c} for (ix, iy), c in grade.items()]
    return {
        "bins2d": bins2d,
        "grid": {
            "x0": round(10 ** lo_x), "x1": round(10 ** hi_x),
            "y0": round(math.expm1(lo_y), 2), "y1": round(math.expm1(hi_y), 2),
            "n_bins": N_BINS_DEPTH, "x_scale": "log10", "y_scale": "log1p",
        },
        "n": n,
    }


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
        lengths[s] = _build_length_block(comprimentos)
        depth[s] = _build_depth_block(comprimentos, profundidades)
    if lengths:
        saida["lengths"] = lengths
    if depth:
        saida["depth"] = depth

    return saida


# ── Aba Catálogo viral ───────────────────────────────────────────────────────

TAXONOMY_BLOCK = Block(name="viral_taxonomy", fields=(
    "sample", "Phylum", "Class", "Order", "Family", "Genus", "count",
))

# A nomenclatura viral do ICTV e padronizada por SUFIXO, nao por posicao no
# campo `lineage` -- o numero de campos e o prefixo de cada token variam entre
# ferramentas (geNomad e mmseqs2 nao concordam), entao contar campos (como o
# _deepest_level de data_loaders.py faz, so como fallback de family/genus/
# order) erra silenciosamente quando um rank intermediario falta. Sufixos mais
# longos primeiro: nao ha overlap real entre eles (nenhum e sufixo de outro),
# mas testar do mais longo pro mais curto e a defesa correta mesmo assim.
_ICTV_RANK_SUFFIXES = (
    ("Class", "viricetes"),
    ("Phylum", "viricota"),
    ("Order", "virales"),
    ("Family", "viridae"),
    ("Kingdom", "virae"),
    ("Genus", "virus"),
    ("Realm", "viria"),
)
_LINEAGE_NULL = {"unclassified", ""}


def _ictv_rank_of(token):
    # Prefixos observados na rodada real: "-_" (mmseqs2/custom) e estilos
    # tipo GTDB "d__"/"p__". Ambos sao ruido antes do nome do taxon --
    # descartamos qualquer prefixo nao alfabetico.
    limpo = re.sub(r'^[^A-Za-z]+', '', (token or '').strip())
    if not limpo or limpo.lower() in _LINEAGE_NULL:
        return None, None
    baixo = limpo.lower()
    for rank, sufixo in _ICTV_RANK_SUFFIXES:
        if baixo.endswith(sufixo):
            return rank, limpo
    return None, None


def _ranks_from_lineage(lineage):
    """Classifica cada token de `lineage` (separado por ';') pelo sufixo
    ICTV do seu nome, nao pela posicao. Devolve um dict rank -> nome, so com
    os ranks efetivamente encontrados."""
    ranks = {}
    for token in (lineage or "").split(";"):
        rank, nome = _ictv_rank_of(token)
        if rank:
            ranks[rank] = nome
    return ranks
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
        # final_family/final_genus/final_order sao a decisao final da
        # pipeline (load_viral_taxonomy ja aplicou seu proprio fallback
        # posicional) -- tem precedencia. `lineage` so preenche o que eles
        # nao trazem, e e a UNICA fonte de Phylum/Class (nunca existiram como
        # campo explicito).
        da_linhagem = _ranks_from_lineage(r.get("Lineage", ""))
        chave = (
            r.get("sample", ""),
            da_linhagem.get("Phylum", ""),
            da_linhagem.get("Class", ""),
            r.get("Order") or da_linhagem.get("Order", ""),
            r.get("Family") or da_linhagem.get("Family", ""),
            r.get("Genus") or da_linhagem.get("Genus", ""),
        )
        contagem[chave] = contagem.get(chave, 0) + 1
    linhas = [
        {"sample": s, "Phylum": p, "Class": c, "Order": o,
         "Family": f, "Genus": g, "count": n}
        for (s, p, c, o, f, g), n in contagem.items()
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
        # O consumidor e o GenomeTrack.jsx, que compara f.strand === '-'
        # (string) -- emitir -1|1 (numero) aqui faz a comparacao falhar
        # sempre, e todo gene em fita reversa desenharia apontando pro lado
        # errado, silenciosamente. String e o contrato, nao numero.
        strand_bruta = (row.get("strand", "") or "").strip()
        strand = "-" if strand_bruta in ("-", "-1") else "+"
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
         # APROXIMACAO: "length" e o maior "end" entre as features do
         # pharokka para este contig, NAO o comprimento real do contig (que
         # este loader nao le -- so tem o TSV de anotacao). So afeta a
         # ordenacao dos "50 mais longos" em _limita_explorer; quem ligar
         # isto ao catalogo global de vOTUs pode trocar por um comprimento
         # de verdade (ex.: FASTA de representantes) se a ordenacao exata
         # importar.
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


# ── Aba Catálogo de MAGs ─────────────────────────────────────────────────────
#
# A fonte e SEMPRE `mag_catalog/`, nunca as vistas `{sample}/bins/...`. As
# vistas existem para as ferramentas a jusante e implementam heranca: o mesmo
# MAG representante aparece sob cada amostra cujo cluster ele representa. Ler
# delas aqui contaria o mesmo organismo uma vez por amostra.

MAG_QUALITY_BLOCK = Block(name="mag_quality", fields=(
    "genome", "source", "completeness", "contamination", "css", "gunc_pass",
    "is_representative", "representative",
))
MAG_TAXONOMY_BLOCK = Block(name="mag_taxonomy", fields=(
    "genome", "source", "Phylum", "Class", "Order", "Family", "Genus",
    "count", "inherited", "representative",
))

# Prefixo do GTDB -> rank do RankSelector. `d__`/`s__` ficam de fora: o
# seletor vai de filo a genero, e embarcar campo que nenhum grafico consome e
# exatamente o padrao que o schema existe para impedir.
_GTDB_PREFIXOS = (
    ("p__", "Phylum"), ("c__", "Class"), ("o__", "Order"),
    ("f__", "Family"), ("g__", "Genus"),
)

GUNC_TSV = "GUNC.progenomes_2.1.maxCSS_level.tsv"


def _gtdb_ranks(classification):
    """'d__Bacteria;p__Bacillota;...' -> {rank: nome}, so os ranks do seletor.

    Um token sem nome depois do prefixo (`o__`) e "nao classificado neste
    rank" -- entra como string vazia, distinta de rank ausente na linhagem.
    """
    ranks = {}
    for token in (classification or "").split(";"):
        token = token.strip()
        for prefixo, rank in _GTDB_PREFIXOS:
            if token.startswith(prefixo):
                ranks[rank] = token[len(prefixo):].strip()
    return ranks


def _le_membership(caminho):
    """[(source_id, original_bin_id, member_id, representative_id)]."""
    linhas = []
    for row in parse_tsv(caminho):
        member = (row.get("member_id") or "").strip()
        if not member:
            continue
        linhas.append((
            (row.get("source_id") or "").strip(),
            (row.get("original_bin_id") or "").strip(),
            member,
            (row.get("representative_id") or "").strip() or member,
        ))
    return linhas


def _float_ou_none(valor):
    try:
        return float(valor)
    except (TypeError, ValueError):
        return None


def _le_gunc(outdir, membership):
    """{member_id: (css, pass)} a partir dos TSVs por fonte.

    A coluna `genome` do GUNC traz o nome ORIGINAL do bin, e esse nome colide
    entre fontes: o Binette emite `binette_bin1` em toda amostra e o VAMB
    emite inteiros nus. Casar sem prefixar a fonte daria o CSS de um
    organismo a outro, em silencio -- a mesma familia de bug que o namespace
    do catalogo existe para matar.
    """
    fontes = {}
    for source_id, bin_name, member_id, _rep in membership:
        fontes.setdefault(source_id, {})[bin_name] = member_id

    saida = {}
    for source_id, bins in fontes.items():
        candidatos = [
            os.path.join(outdir, source_id, "bins", "gunc", GUNC_TSV),
            os.path.join(outdir, "coassembly", source_id, "bins", "gunc", GUNC_TSV),
        ]
        caminho = next((c for c in candidatos if os.path.exists(c)), None)
        if caminho is None:
            continue
        for row in parse_tsv(caminho):
            member_id = bins.get((row.get("genome") or "").strip())
            if not member_id:
                continue
            passou = (row.get("pass.GUNC") or "").strip().lower()
            saida[member_id] = (
                _float_ou_none(row.get("clade_separation_score")),
                True if passou == "true" else (False if passou == "false" else None),
            )
    return saida


def _build_mag_quality(outdir, membership, catalog_dir):
    checkm2 = {}
    for row in parse_tsv(os.path.join(catalog_dir, "checkm2_quality_report.tsv")):
        nome = (row.get("Name") or "").strip()
        if nome:
            checkm2[nome] = (_float_ou_none(row.get("Completeness")),
                             _float_ou_none(row.get("Contamination")))
    gunc = _le_gunc(outdir, membership)

    linhas = []
    for source_id, _bin_name, member_id, rep in membership:
        completude, contaminacao = checkm2.get(member_id, (None, None))
        css, passou = gunc.get(member_id, (None, None))
        linhas.append({
            "genome": member_id,
            "source": source_id,
            "completeness": completude,
            "contamination": contaminacao,
            "css": css,
            "gunc_pass": passou,
            "is_representative": member_id == rep,
            "representative": rep,
        })
    return project(MAG_QUALITY_BLOCK, linhas)


def _build_mag_clusters(membership):
    clusters = {}
    for source_id, _bin_name, _member, rep in membership:
        d = clusters.setdefault(rep, {"n_members": 0, "sources": set()})
        d["n_members"] += 1
        d["sources"].add(source_id)
    tamanhos = [
        {"representative": rep, "n_members": d["n_members"],
         "n_sources": len(d["sources"])}
        for rep, d in clusters.items()
    ]
    tamanhos.sort(key=lambda c: (-c["n_members"], c["representative"]))
    return {
        "n_mags": len(membership),
        "n_clusters": len(clusters),
        "sizes": tamanhos,
    }


def _build_mag_taxonomy(membership, catalog_dir):
    """Uma linha por MAG, com a linhagem do seu REPRESENTANTE.

    Herdar e o ponto: o GTDB-Tk rodou uma vez, sobre as representantes. Um
    membro que nao e representante nunca aparece na tabela do classify_wf, e
    filtrar essa tabela pelo prefixo da fonte devolveria quase nada -- foi
    exatamente esse o bug do `viral_taxonomy` em 2026-08-18.
    """
    classificacao = {}
    for arquivo in ("gtdbtk.bac120.summary.tsv", "gtdbtk.ar53.summary.tsv"):
        caminho = os.path.join(catalog_dir, "gtdbtk", "classify", arquivo)
        for row in parse_tsv(caminho):
            genoma = (row.get("user_genome") or "").strip()
            if genoma:
                classificacao[genoma] = _gtdb_ranks(row.get("classification"))
    if not classificacao:
        return None

    linhas = []
    for source_id, _bin_name, member_id, rep in membership:
        ranks = classificacao.get(rep)
        if ranks is None:
            continue
        linhas.append({
            "genome": member_id,
            "source": source_id,
            "Phylum": ranks.get("Phylum", ""),
            "Class": ranks.get("Class", ""),
            "Order": ranks.get("Order", ""),
            "Family": ranks.get("Family", ""),
            "Genus": ranks.get("Genus", ""),
            # O sunburst e a barra somam `count`; aqui a unidade e o MAG.
            "count": 1,
            "inherited": member_id != rep,
            "representative": rep,
        })
    return project(MAG_TAXONOMY_BLOCK, linhas) if linhas else None


def _build_mag_kegg(catalog_dir):
    caminho = os.path.join(catalog_dir, "kegg_modules", "module_completeness.tsv")
    if not os.path.exists(caminho):
        return None
    valores, modulos, genomas = {}, {}, []
    for row in parse_tsv(caminho):
        mag = (row.get("mag") or "").strip()
        modulo = (row.get("module_accession") or "").strip()
        completude = _float_ou_none(row.get("completeness"))
        if not mag or not modulo or completude is None:
            continue
        if mag not in valores:
            valores[mag] = {}
            genomas.append(mag)
        valores[mag][modulo] = completude
        d = modulos.setdefault(modulo, {
            "module": modulo,
            "name": (row.get("pathway_name") or "").strip(),
            # missing_ko por GENOMA, nao um so por modulo: o passo que falta
            # e diferente em cada MAG, e e ele que torna a via interpretavel.
            "missing": {},
        })
        faltando = (row.get("missing_ko") or "").strip()
        if faltando:
            d["missing"][mag] = faltando
    if not valores:
        return None
    ordenados = sorted(modulos.values(), key=lambda m: m["module"])
    return {"genomes": genomas, "modules": ordenados, "values": valores}


def _build_mag_cazy(catalog_dir):
    caminho = os.path.join(catalog_dir, "kegg", "cazy_per_mag.tsv")
    if not os.path.exists(caminho):
        return None
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    try:
        from annotation_tables import cazy_class
    except ImportError:
        return None

    por_genoma = {}
    for row in parse_tsv(caminho):
        mag = (row.get("mag") or "").strip()
        familia = (row.get("cazy_family") or "").strip()
        if not mag or not familia:
            continue
        classes = por_genoma.setdefault(mag, {})
        cls = cazy_class(familia)
        classes[cls] = classes.get(cls, 0) + 1
    if not por_genoma:
        return None
    return [{"genome": g, "parts": partes} for g, partes in por_genoma.items()]


def build_prokaryotic(outdir, samples, groups):
    catalog_dir = os.path.join(outdir, "mag_catalog")
    membership = _le_membership(os.path.join(catalog_dir, "mag_membership.tsv"))
    if not membership:
        return {}

    saida = {
        "quality": _build_mag_quality(outdir, membership, catalog_dir),
        "clusters": _build_mag_clusters(membership),
    }
    for chave, valor in (
        ("taxonomy", _build_mag_taxonomy(membership, catalog_dir)),
        ("kegg", _build_mag_kegg(catalog_dir)),
        ("cazy", _build_mag_cazy(catalog_dir)),
    ):
        if valor is not None:
            saida[chave] = valor
    return saida


# ── Aba Defesa, AMR e plasmídeos ─────────────────────────────────────────────

DEFENSE_BLOCK = Block(name="defense", fields=("genome", "system", "count"))
AMR_BLOCK = Block(name="amr", fields=(
    "genome", "contig", "gene", "drug_class", "n_tools", "tools"))
PLASMID_BLOCK = Block(name="plasmid", fields=(
    "genome", "contig", "gene", "start", "end"))

# Rotulos das evidencias do UpSet. Ficam num so lugar porque aparecem em
# `sets`, em `combos` e no painel -- tres copias divergiriam.
EV_REPLICON = "replicon"
EV_ARG = "ARG de consenso"
EV_DEFESA = "sistema de defesa"

# O consenso de AMR so vale com duas ferramentas concordando. Uma so e
# exploratorio, e misturar os dois niveis na mesma barra apagaria a diferenca.
MIN_TOOLS_AMR = 2


def _contig_do_orf(orf_id):
    """'k141_1_5' -> 'k141_1'. Corte no ULTIMO '_', que e a convencao do
    Prodigal ('{contig}_{n}') -- o nome do contig contem '_', entao cortar no
    primeiro devolveria 'k141' e nenhuma colocalizacao casaria."""
    return (orf_id or "").rsplit("_", 1)[0] if "_" in (orf_id or "") else ""


def _representantes(outdir):
    membership = _le_membership(
        os.path.join(outdir, "mag_catalog", "mag_membership.tsv"))
    return membership, {rep for _s, _b, _m, rep in membership}


def _genoma_da_proteina(protein_id, representantes):
    """('{rep}__{orf}', reps) -> (genoma, orf).

    NUNCA cortar no primeiro '__': um ID do catalogo ja e '{source}__{bin}',
    entao 'S1__binette_bin1__k141_1_5' cortado ali devolve 'S1' e credita o
    achado a AMOSTRA. Casa-se contra os representantes conhecidos, do prefixo
    mais longo para o mais curto, como em scripts/checkv_provirus.py.
    """
    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from mag_catalog import resolve_prefixed_id
    return resolve_prefixed_id(protein_id, representantes)


def _build_defense(catalog_dir):
    caminho = os.path.join(catalog_dir, "defensefinder", "defensefinder_systems.tsv")
    if not os.path.exists(caminho):
        return None
    contagem = {}
    for row in parse_tsv(caminho):
        # `genome` ja e o ID do representante: o DefenseFinder roda uma vez
        # por genoma (o MacSyFinder precisa de ordem genica no repliconte),
        # entao a regra escreve a coluna e nao ha prefixo a resolver aqui.
        genome = (row.get("genome") or "").strip()
        sistema = (row.get("type") or row.get("subtype") or "").strip()
        if not genome or not sistema:
            continue
        contagem[(genome, sistema)] = contagem.get((genome, sistema), 0) + 1
    if not contagem:
        return None
    linhas = [{"genome": g, "system": s, "count": n}
              for (g, s), n in sorted(contagem.items())]
    return project(DEFENSE_BLOCK, linhas)


def _build_amr(catalog_dir, representantes):
    caminho = os.path.join(catalog_dir, "amr_consensus", "amr_consensus.tsv")
    if not os.path.exists(caminho):
        return None
    linhas = []
    for row in parse_tsv(caminho):
        locus = (row.get("locus") or "").strip()
        if not locus or safe_int(row.get("n_tools", 0)) < MIN_TOOLS_AMR:
            continue
        genome, orf = _genoma_da_proteina(locus, representantes)
        if not genome:
            continue
        linhas.append({
            "genome": genome,
            "contig": _contig_do_orf(orf),
            "gene": (row.get("gene_name") or "").strip(),
            "drug_class": (row.get("drug_class") or "").strip(),
            "n_tools": safe_int(row.get("n_tools", 0)),
            "tools": (row.get("tools_detected") or "").strip(),
        })
    return project(AMR_BLOCK, linhas) if linhas else None


def _build_plasmids(catalog_dir):
    caminho = os.path.join(catalog_dir, "abricate", "plasmidfinder_results.tsv")
    if not os.path.exists(caminho):
        return None
    linhas = []
    for row in parse_tsv(caminho):
        # O ABRicate identifica o genoma pelo CAMINHO do arquivo que recebeu.
        arquivo = (row.get("#FILE") or row.get("FILE") or "").strip()
        contig = (row.get("SEQUENCE") or "").strip()
        if not arquivo or not contig:
            continue
        genome = os.path.basename(arquivo)
        for ext in (".fa", ".fna", ".fasta"):
            if genome.endswith(ext):
                genome = genome[: -len(ext)]
                break
        linhas.append({
            "genome": genome,
            "contig": contig,
            "gene": (row.get("GENE") or "").strip(),
            "start": safe_int(row.get("START", 0)),
            "end": safe_int(row.get("END", 0)),
        })
    return project(PLASMID_BLOCK, linhas) if linhas else None


ISLAND_GENE_BLOCK = Block(name="island_gene", fields=(
    "start", "end", "strand", "label", "kind"))
MAX_ILHAS = 20


def _build_islands(catalog_dir):
    """Ilhas de defesa sobre os representantes, com coordenadas em bp.

    O nucleo e `scripts/defense_islands.py`, o MESMO que o portao do
    pangenoma usa -- duas definicoes de "ilha" divergiriam em silencio.

    Os IDs de `protein_in_syst` sao NUS (`k141_1_2`): o DefenseFinder roda
    uma vez por genoma, sobre o .faa daquele genoma. O prefixo `{genome}__`
    so existe no .faa concatenado que as ferramentas de AMR consomem, e
    confundir os dois faria nenhuma ilha casar com gene nenhum.
    """
    manifest = os.path.join(catalog_dir, "proteins", "manifest.txt")
    systems = os.path.join(catalog_dir, "defensefinder", "defensefinder_systems.tsv")
    if not (os.path.exists(manifest) and os.path.exists(systems)):
        return None

    import sys
    sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
    from defense_islands import genes_by_contig, find_islands

    por_genoma = {}
    for row in parse_tsv(systems):
        genome = (row.get("genome") or "").strip()
        sistema = (row.get("type") or row.get("subtype") or "").strip()
        for prot in (row.get("protein_in_syst") or "").split(","):
            prot = prot.strip()
            if genome and sistema and prot:
                por_genoma.setdefault(genome, {})[prot] = (
                    genome, sistema, (row.get("sys_id") or sistema).strip())

    ilhas = []
    for linha in _read_manifest_lines(manifest):
        genome, faa = linha[0], linha[3]
        prot_to_sys = por_genoma.get(genome)
        if not prot_to_sys:
            continue
        for isl in find_islands(genes_by_contig(faa), prot_to_sys):
            genes = [
                {"start": g.get("Start"), "end": g.get("End"),
                 "strand": "-" if g.get("Strand") == -1 else "+",
                 "label": g.get("System") or "",
                 # `kind` separa gene de defesa de gene qualquer dentro da
                 # janela: sem isso a trilha pintaria a ilha inteira como
                 # defesa, que e uma afirmacao mais forte que o dado.
                 "kind": "defense" if g.get("System") else "outro"}
                for g in isl.get("window_genes", [])
                if g.get("Start") is not None and g.get("End") is not None
            ]
            ilhas.append({
                "genome": genome,
                "contig": isl["Contig"],
                "start": isl.get("Start_bp"),
                "end": isl.get("End_bp"),
                "n_genes": isl["n_genes"],
                "n_systems": isl["n_systems"],
                "systems": isl["Systems"],
                "genes": project(ISLAND_GENE_BLOCK, genes),
            })
    if not ilhas:
        return None
    # Densidade de sistemas e o criterio de ranking: uma ilha de 3 sistemas em
    # 6 genes diz mais que uma de 3 sistemas em 40.
    ilhas.sort(key=lambda i: (-i["n_systems"], i["n_genes"]))
    return ilhas[:MAX_ILHAS]


def _read_manifest_lines(path):
    linhas = []
    with open(path, encoding='utf-8') as fh:
        for linha in fh:
            partes = linha.rstrip("\n").split("\t")
            if len(partes) >= 5:
                linhas.append(partes)
    return linhas


def _build_colocalizacao(amr, plasmids):
    """ARGs num contig que carrega replicon plasmidial, no MESMO genoma.

    A chave e o PAR (genoma, contig), nunca o contig sozinho: nomes de contig
    colidem entre MAGs -- todo assembly do MEGAHIT emite `k141_1` -- e casar
    so pelo nome afirmaria colocalizacao entre organismos diferentes.

    Isto e evidencia de ORIGEM plasmidial do ARG, nao prova de plasmidio
    intacto: a montagem pode ter quebrado o elemento. A ressalva vai escrita
    na aba, nao so aqui.
    """
    contigs_com_replicon = {(p["genome"], p["contig"]) for p in plasmids}
    marcados = [f"{a['genome']}|{a['contig']}|{a['gene']}" for a in amr
                if (a["genome"], a["contig"]) in contigs_com_replicon]
    return {
        "n_args": len(amr),
        "n_args_on_replicon": len(marcados),
        "args_on_replicon": sorted(marcados),
    }


def _build_upset(defense, amr, plasmids, islands):
    """MAGs por combinacao de evidencias. A unidade e o MAG, nunca o hit."""
    evidencias = {}
    for rows, rotulo in ((plasmids, EV_REPLICON), (amr, EV_ARG),
                         (defense, EV_DEFESA)):
        for r in rows or []:
            evidencias.setdefault(r["genome"], set()).add(rotulo)
    for i in islands or []:
        evidencias.setdefault(i["genome"], set()).add("ilha de defesa")

    sets, combos = {}, {}
    for genome, marcas in evidencias.items():
        chave = tuple(sorted(marcas))
        combos[chave] = combos.get(chave, 0) + 1
        for m in marcas:
            sets[m] = sets.get(m, 0) + 1
    return {
        "sets": sets,
        "combos": [{"tools": list(k), "count": n}
                   for k, n in sorted(combos.items(), key=lambda kv: -kv[1])],
    }


def build_defense_amr(outdir):
    catalog_dir = os.path.join(outdir, "mag_catalog")
    _membership, representantes = _representantes(outdir)
    if not representantes:
        return {}

    defense = _build_defense(catalog_dir)
    amr = _build_amr(catalog_dir, representantes)
    plasmids = _build_plasmids(catalog_dir)
    islands = _build_islands(catalog_dir)

    saida = {}
    for chave, valor in (("defense", defense), ("amr", amr),
                          ("plasmids", plasmids), ("islands", islands)):
        if valor is not None:
            saida[chave] = valor
    if not saida:
        return {}

    # Sem replicon nao existe pergunta de colocalizacao. Emitir zero aqui
    # seria afirmar "nenhum ARG em plasmidio", que e diferente de "a trilha
    # do PlasmidFinder nao rodou".
    if amr and plasmids:
        saida["colocalization"] = _build_colocalizacao(amr, plasmids)

    saida["upset"] = _build_upset(defense, amr, plasmids, islands)
    return saida


# ── Aba Pangenoma ────────────────────────────────────────────────────────────

CANDIDATE_BLOCK = Block(name="pangenome_candidate", fields=(
    "representative", "n_members", "n_islands", "n_systems", "n_args",
    "n_plasmid", "criterio", "eligible"))
CLUSTER_BLOCK = Block(name="pangenome_cluster", fields=(
    "representative", "n_members", "n_evaluable", "n_core", "n_variable",
    "n_singleton", "completeness_median", "size_median", "taxonomy"))

# Estados da matriz. Os TRES primeiros sao biologia; o quarto e formato: as
# colunas do TSV sao TODOS os membros de TODOS os clusters, e '-' marca a
# celula que nao pertence aquele cluster. Ele nunca pode virar ausencia.
NAO_AVALIAVEL = "?"
FORA_DO_CLUSTER = "-"


def _bool_tsv(valor):
    return (valor or "").strip().lower() in ("true", "1", "yes")


def _build_pangenome_candidates(pg_dir):
    caminho = os.path.join(pg_dir, "candidates.tsv")
    if not os.path.exists(caminho):
        return None
    linhas = []
    for row in parse_tsv(caminho):
        rep = (row.get("representative_id") or "").strip()
        if not rep:
            continue
        linhas.append({
            "representative": rep,
            "n_members": safe_int(row.get("n_members", 0)),
            "n_islands": safe_int(row.get("n_islands", 0)),
            "n_systems": safe_int(row.get("n_systems", 0)),
            "n_args": safe_int(row.get("n_args", 0)),
            # Sinal de MOBILIDADE, nunca criterio: o portao da regra nao
            # elege cluster por plasmidio sozinho, e o report nao pode
            # sugerir que elege.
            "n_plasmid": safe_int(row.get("n_plasmid", 0)),
            "criterio": (row.get("criterio") or "").strip(),
            "eligible": _bool_tsv(row.get("eligible")),
        })
    return project(CANDIDATE_BLOCK, linhas) if linhas else None


def _build_pangenome_clusters(pg_dir):
    caminho = os.path.join(pg_dir, "cluster_summary.tsv")
    if not os.path.exists(caminho):
        return None
    linhas = []
    for row in parse_tsv(caminho):
        rep = (row.get("representative_id") or "").strip()
        if not rep:
            continue
        linhas.append({
            "representative": rep,
            "n_members": safe_int(row.get("n_members", 0)),
            "n_evaluable": safe_int(row.get("n_members_avaliaveis", 0)),
            "n_core": safe_int(row.get("n_genes_core", 0)),
            "n_variable": safe_int(row.get("n_genes_variaveis", 0)),
            "n_singleton": safe_int(row.get("n_genes_singleton", 0)),
            "completeness_median": _float_ou_none(row.get("completude_mediana")),
            "size_median": _float_ou_none(row.get("tamanho_mediana_bp")),
            "taxonomy": (row.get("gtdb_taxonomy") or "").strip(),
        })
    return project(CLUSTER_BLOCK, linhas) if linhas else None


def _parse_completude_cabecalho(linha):
    """'# completude: S1__bin1=98.0, S1__bin8=42.0' -> {membro: float}.

    A completude viaja no cabecalho da matriz porque "1/2" so e
    interpretavel com ela a vista -- e o report nao pode separar as duas.
    """
    saida = {}
    _rotulo, _sep, corpo = linha.partition(":")
    for item in corpo.split(","):
        nome, _sep2, valor = item.strip().partition("=")
        v = _float_ou_none(valor)
        if nome and v is not None:
            saida[nome.strip()] = v
    return saida


def _build_pangenome_matrix(pg_dir):
    caminho = os.path.join(pg_dir, "gene_by_member.tsv")
    if not os.path.exists(caminho):
        return None, {}

    completude = {}
    cabecalho = None
    linhas_dados = []
    with open(caminho, encoding='utf-8') as fh:
        for linha in fh:
            linha = linha.rstrip("\n")
            if linha.startswith("#"):
                if "completude" in linha:
                    completude = _parse_completude_cabecalho(linha)
                continue
            if cabecalho is None:
                cabecalho = linha.split("\t")
                continue
            if linha.strip():
                linhas_dados.append(linha.split("\t"))
    if not cabecalho or not linhas_dados:
        return None, completude

    N_FIXAS = 6      # cluster, tipo, gene, freq, n_present, n_evaluable
    membros_globais = cabecalho[N_FIXAS:]

    por_cluster = {}
    for campos in linhas_dados:
        rep = campos[0]
        estados_crus = dict(zip(membros_globais, campos[N_FIXAS:]))
        # Os membros DESTE cluster sao os que nao trazem '-'. Filtrar aqui,
        # e nao no navegador, e o que impede a celula de formato de ser
        # desenhada como ausencia biologica.
        membros = [m for m in membros_globais
                   if estados_crus.get(m) != FORA_DO_CLUSTER]
        d = por_cluster.setdefault(rep, {"members": membros, "rows": []})
        d["rows"].append({
            "tipo": campos[1],
            "gene": campos[2],
            "freq": campos[3],
            "n_present": safe_int(campos[4]),
            # Denominador ja calculado pela regra, excluindo os '?'. Recontar
            # aqui criaria uma segunda fonte de verdade que poderia discordar.
            "n_evaluable": safe_int(campos[5]),
            "states": {m: estados_crus[m] for m in membros},
        })
    return por_cluster, completude


def build_pangenome(outdir):
    pg_dir = os.path.join(outdir, "mag_catalog", "pangenome")
    candidates = _build_pangenome_candidates(pg_dir)
    clusters = _build_pangenome_clusters(pg_dir)
    matrix, completude = _build_pangenome_matrix(pg_dir)

    saida = {}
    for chave, valor in (("candidates", candidates), ("clusters", clusters),
                          ("matrix", matrix)):
        if valor:
            saida[chave] = valor
    if saida and completude:
        saida["completeness"] = completude
    return saida


# ── Abas Diversidade e Leituras ──────────────────────────────────────────────

ALPHA_BLOCK = Block(name="alpha", fields=("sample", "domain", "index", "value"))
PCOA_BLOCK = Block(name="pcoa", fields=(
    "sample", "pc1", "pc2", "var_pc1", "var_pc2"))
PROCRUSTES_BLOCK = Block(name="procrustes_pair", fields=(
    "sample", "viral_pc1", "viral_pc2", "prok_pc1", "prok_pc2"))

# Os quatro indices que compute_diversity.py pode escrever. Simpson e Chao1
# saem VAZIOS quando nao ha contagens de reads -- os dois sao estimadores de
# contagem (f1/f2 e a*(a-1)) e nao se calculam sobre RPKM. Vazio vira lacuna
# declarada, nunca 0.0: zero e um valor, e Simpson 0 diria "uma unica especie
# domina completamente".
INDICES_ALFA = ("observed", "shannon", "simpson", "chao1")

PCOA_TRILHAS = (("viral", "beta_pcoord_viral.tsv"),
                ("prok", "beta_pcoord_prok.tsv"),
                ("combined", "beta_pcoord_combined.tsv"))

# Quantos taxa por dominio viajam para o navegador. A tabela do sylph tem uma
# linha por clado x uma coluna por amostra; embarca-la inteira e o padrao que
# o orcamento de payload existe para impedir. O corte e por abundancia total,
# e o painel diz quantos ficaram de fora.
MAX_TAXA_READS = 60


def build_diversity(outdir, samples):
    div_dir = os.path.join(outdir, "diversity")
    saida = {}

    alpha_path = os.path.join(div_dir, "alpha_diversity.tsv")
    if os.path.exists(alpha_path):
        linhas = load_alpha_diversity(alpha_path)
        if linhas:
            saida["alpha"] = project(ALPHA_BLOCK, linhas)
            presentes = {r["index"] for r in linhas}
            # O que a rodada NAO calculou tem de ser nomeado; senao o indice
            # simplesmente some do painel e o leitor conclui que a pipeline
            # nao roda Chao1, em vez de "esta rodada nao tinha contagens".
            saida["alpha_missing"] = [i for i in INDICES_ALFA if i not in presentes]

    pcoa = {}
    for nome, arquivo in PCOA_TRILHAS:
        caminho = os.path.join(div_dir, arquivo)
        if not os.path.exists(caminho):
            continue
        linhas = load_pcoord(caminho)
        if linhas:
            pcoa[nome] = project(PCOA_BLOCK, linhas)
    if pcoa:
        saida["pcoa"] = pcoa

    proc_path = os.path.join(div_dir, "procrustes_coords.tsv")
    if os.path.exists(proc_path):
        pares, disparidade = [], None
        for row in parse_tsv(proc_path):
            amostra = (row.get("sample") or "").strip()
            if not amostra:
                continue
            pares.append({
                "sample": amostra,
                "viral_pc1": _float_ou_none(row.get("viral_PC1")),
                "viral_pc2": _float_ou_none(row.get("viral_PC2")),
                "prok_pc1": _float_ou_none(row.get("prok_PC1")),
                "prok_pc2": _float_ou_none(row.get("prok_PC2")),
            })
            disparidade = _float_ou_none(row.get("disparity"))
        if pares:
            saida["procrustes"] = {
                "pairs": project(PROCRUSTES_BLOCK, pares),
                "disparity": disparidade,
            }
    return saida


# Escrito por extenso porque e uma restricao de INTERPRETACAO, nao um detalhe
# de implementacao: o sylph identifica genomas de referencia do IMG/VR
# ('t__IMGVR_UViG_...') e a trilha de montagem identifica contigs do MEGAHIT
# ('k141_...'). Nenhum join entre as duas e valido, e o report diz isso na
# cara do usuario em vez de deixar a tentacao de cruzá-los.
AVISO_ESPACO_IDS = (
    "Os identificadores desta aba são genomas de referência do banco do sylph "
    "(t__IMGVR_UViG_…) e não têm relação com os contigs montados (k141_…) das "
    "demais abas. Nenhum cruzamento entre as duas trilhas é válido."
)


def _corta_por_abundancia(linhas, samples, limite=MAX_TAXA_READS):
    ordenadas = sorted(linhas,
                       key=lambda r: sum(r.get(s, 0.0) or 0.0 for s in samples),
                       reverse=True)
    return ordenadas[:limite], max(len(ordenadas) - limite, 0)


def build_reads(outdir, samples):
    rc_dir = os.path.join(outdir, "reads_classify")
    otu = os.path.join(rc_dir, "otu_table.tsv")
    host = os.path.join(rc_dir, "viral_abundance_by_host.tsv")
    if not os.path.exists(otu) and not os.path.exists(host):
        return {}

    dados = load_reads_classify(otu if os.path.exists(otu) else "",
                                host if os.path.exists(host) else "",
                                samples)
    if not (dados.get("has_data") or dados.get("host")):
        return {}

    saida = {"id_space_warning": AVISO_ESPACO_IDS}
    for chave in ("viral", "prok", "archaea"):
        linhas = dados.get(chave) or []
        if not linhas:
            continue
        cortadas, omitidos = _corta_por_abundancia(linhas, samples)
        saida[chave] = cortadas
        if omitidos:
            saida.setdefault("truncated", {})[chave] = omitidos
    if dados.get("host"):
        saida["host"] = dados["host"][:MAX_TAXA_READS]
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

    prokaryotic = build_prokaryotic(outdir, samples, grupos)
    if prokaryotic:
        dados["prokaryotic"] = prokaryotic

    defesa = build_defense_amr(outdir)
    if defesa:
        dados["defense_amr"] = defesa

    pangenoma = build_pangenome(outdir)
    if pangenoma:
        dados["pangenome"] = pangenoma

    diversidade = build_diversity(outdir, samples)
    if diversidade:
        dados["diversity"] = diversidade

    leituras = build_reads(outdir, samples)
    if leituras:
        dados["reads"] = leituras

    return dados
