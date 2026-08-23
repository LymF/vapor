"""Bloco `prokaryotic` do report v2 — catálogo global de MAGs.

O que estes testes protegem, e por quê:

1. A fonte é `mag_catalog/`, nunca as vistas por amostra. As vistas existem
   para as ferramentas a jusante; ler delas aqui contaria o mesmo MAG uma vez
   por amostra em que a espécie apareceu.
2. Um ID do catálogo é `{source}__{bin}` e um ID de proteína é
   `{source}__{bin}__{orf}`. Cortar no primeiro `__` atribui o achado à
   AMOSTRA. Aqui o corte correto é contra os IDs conhecidos.
3. Taxonomia GTDB, KEGG e CAZy foram computados no REPRESENTANTE. Todo MAG
   que não é representante recebe a linha do seu representante marcada como
   herdada — filtrar a tabela global pelo prefixo da fonte devolveria quase
   nada (o bug do `viral_taxonomy` de 2026-08-18).
"""
import os

import pytest

from report.renderer_v2 import build_prokaryotic, _gtdb_ranks


def _escreve(caminho, texto):
    os.makedirs(os.path.dirname(caminho), exist_ok=True)
    with open(caminho, 'w', encoding='utf-8') as fh:
        fh.write(texto)


@pytest.fixture
def outdir(tmp_path):
    """Catálogo com dois clusters: um com dois membros (S1 representa S2) e
    um só de S2. É a forma mínima que distingue herança de identidade."""
    o = str(tmp_path / "results")
    cat = os.path.join(o, "mag_catalog")

    _escreve(os.path.join(cat, "mag_membership.tsv"),
             "source_id\toriginal_bin_id\tmember_id\trepresentative_id\n"
             "S1\tbinette_bin1\tS1__binette_bin1\tS1__binette_bin1\n"
             "S2\tbinette_bin1\tS2__binette_bin1\tS1__binette_bin1\n"
             "S2\tbinette_bin2\tS2__binette_bin2\tS2__binette_bin2\n")

    _escreve(os.path.join(cat, "checkm2_quality_report.tsv"),
             "Name\tCompleteness\tContamination\n"
             "S1__binette_bin1\t95.2\t1.4\n"
             "S2__binette_bin1\t72.0\t3.1\n"
             "S2__binette_bin2\t55.5\t9.8\n")

    _escreve(os.path.join(o, "S1", "bins", "gunc",
                          "GUNC.progenomes_2.1.maxCSS_level.tsv"),
             "genome\tclade_separation_score\tpass.GUNC\n"
             "binette_bin1\t0.02\tTrue\n")
    _escreve(os.path.join(o, "S2", "bins", "gunc",
                          "GUNC.progenomes_2.1.maxCSS_level.tsv"),
             "genome\tclade_separation_score\tpass.GUNC\n"
             "binette_bin1\t0.11\tTrue\n"
             "binette_bin2\t0.87\tFalse\n")

    _escreve(os.path.join(cat, "gtdbtk", "classify", "gtdbtk.bac120.summary.tsv"),
             "user_genome\tclassification\n"
             "S1__binette_bin1\td__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria;"
             "o__Enterobacterales;f__Enterobacteriaceae;g__Escherichia;s__Escherichia coli\n"
             "S2__binette_bin2\td__Bacteria;p__Bacillota;c__Bacilli;o__;f__;g__;s__\n")
    _escreve(os.path.join(cat, "gtdbtk", "classify", "gtdbtk.ar53.summary.tsv"),
             "user_genome\tclassification\n")

    _escreve(os.path.join(cat, "kegg_modules", "module_completeness.tsv"),
             "mag\tmodule_accession\tcompleteness\tpathway_name\tmissing_ko\n"
             "S1__binette_bin1\tM00001\t100.0\tGlycolysis\t\n"
             "S1__binette_bin1\tM00002\t66.7\tPentose phosphate\tK00615\n"
             "S2__binette_bin2\tM00001\t80.0\tGlycolysis\tK01810\n")

    _escreve(os.path.join(cat, "kegg", "cazy_per_mag.tsv"),
             "mag\tcazy_family\n"
             "S1__binette_bin1\tGH13\n"
             "S1__binette_bin1\tGH5\n"
             "S1__binette_bin1\tGT51\n"
             "S2__binette_bin2\tPL1\n")
    return o


def test_sem_catalogo_o_bloco_nao_existe(tmp_path):
    # Rodada sem a trilha procariótica: aba some, nunca aparece vazia com
    # eixo quebrado.
    assert build_prokaryotic(str(tmp_path), ["S1"], []) == {}


def test_qualidade_traz_um_mag_por_membro_com_gunc_casado(outdir):
    bloco = build_prokaryotic(outdir, ["S1", "S2"], [])
    por_genoma = {r["genome"]: r for r in bloco["quality"]}
    assert set(por_genoma) == {"S1__binette_bin1", "S2__binette_bin1",
                               "S2__binette_bin2"}
    # O GUNC é medido POR AMOSTRA e sua coluna `genome` traz o nome ORIGINAL
    # do bin, que colide entre amostras: `binette_bin1` existe em S1 e em S2 e
    # são organismos diferentes. Casar sem namespace daria o CSS de um ao
    # outro em silêncio.
    assert por_genoma["S1__binette_bin1"]["css"] == 0.02
    assert por_genoma["S2__binette_bin1"]["css"] == 0.11
    assert por_genoma["S2__binette_bin2"]["gunc_pass"] is False
    assert por_genoma["S1__binette_bin1"]["completeness"] == 95.2
    assert por_genoma["S2__binette_bin2"]["contamination"] == 9.8


def test_representante_marcado_e_membro_aponta_para_ele(outdir):
    bloco = build_prokaryotic(outdir, ["S1", "S2"], [])
    por_genoma = {r["genome"]: r for r in bloco["quality"]}
    assert por_genoma["S1__binette_bin1"]["is_representative"] is True
    assert por_genoma["S2__binette_bin1"]["is_representative"] is False
    assert por_genoma["S2__binette_bin1"]["representative"] == "S1__binette_bin1"
    assert por_genoma["S2__binette_bin1"]["source"] == "S2"


def test_estrutura_do_catalogo_conta_clusters_e_fontes(outdir):
    clusters = build_prokaryotic(outdir, ["S1", "S2"], [])["clusters"]
    assert clusters["n_mags"] == 3
    assert clusters["n_clusters"] == 2
    tamanhos = {c["representative"]: c for c in clusters["sizes"]}
    assert tamanhos["S1__binette_bin1"]["n_members"] == 2
    # Duas amostras no mesmo cluster é a informação que justifica o catálogo
    # existir; um cluster de dois membros da MESMA amostra não é a mesma coisa.
    assert tamanhos["S1__binette_bin1"]["n_sources"] == 2
    assert tamanhos["S2__binette_bin2"]["n_members"] == 1


def test_gtdb_e_herdado_pelo_membro_e_marcado_como_herdado(outdir):
    taxonomy = build_prokaryotic(outdir, ["S1", "S2"], [])["taxonomy"]
    por_genoma = {r["genome"]: r for r in taxonomy}
    # S2__binette_bin1 nunca passou pelo GTDB-Tk (só as representantes
    # passaram). Sem herança ele sumiria da taxonomia; com herança errada
    # (filtro por prefixo) sumiria também.
    assert por_genoma["S2__binette_bin1"]["Genus"] == "Escherichia"
    assert por_genoma["S2__binette_bin1"]["inherited"] is True
    assert por_genoma["S2__binette_bin1"]["representative"] == "S1__binette_bin1"
    assert por_genoma["S1__binette_bin1"]["inherited"] is False


def test_rank_vazio_do_gtdb_nao_vira_nome_vazio(outdir):
    taxonomy = build_prokaryotic(outdir, ["S1", "S2"], [])["taxonomy"]
    bacillota = next(r for r in taxonomy if r["genome"] == "S2__binette_bin2")
    assert bacillota["Phylum"] == "Bacillota"
    # `o__` sem nome é "não classificado neste rank", e o report precisa
    # distinguir isso de um nome literal vazio.
    assert bacillota["Order"] == ""


def test_gtdb_ranks_ignora_prefixo_e_rank_vazio():
    r = _gtdb_ranks("d__Bacteria;p__Bacillota;c__Bacilli;o__;f__;g__;s__")
    assert r["Phylum"] == "Bacillota"
    assert r["Class"] == "Bacilli"
    assert r["Order"] == ""
    assert _gtdb_ranks("") == {}


def test_kegg_preserva_missing_ko(outdir):
    kegg = build_prokaryotic(outdir, ["S1", "S2"], [])["kegg"]
    assert kegg["genomes"] == ["S1__binette_bin1", "S2__binette_bin2"]
    m1 = next(m for m in kegg["modules"] if m["module"] == "M00002")
    # Uma via a 66,7% só é interpretável com o passo que falta — é a razão de
    # a regra do Snakemake ter guardado `missing_ko`.
    assert m1["missing"]["S1__binette_bin1"] == "K00615"
    assert kegg["values"]["S1__binette_bin1"]["M00001"] == 100.0


def test_cazy_agrega_por_classe_preservando_a_familia(outdir):
    cazy = build_prokaryotic(outdir, ["S1", "S2"], [])["cazy"]
    por_genoma = {r["genome"]: r for r in cazy}
    assert por_genoma["S1__binette_bin1"]["parts"] == {"GH": 2, "GT": 1}
    assert por_genoma["S2__binette_bin2"]["parts"] == {"PL": 1}


def test_trilha_parcial_nao_derruba_o_bloco(outdir):
    # KEGG desligado em config.yaml: o resto da aba continua de pé, e a
    # chave some em vez de aparecer como zero.
    os.remove(os.path.join(outdir, "mag_catalog", "kegg_modules",
                           "module_completeness.tsv"))
    bloco = build_prokaryotic(outdir, ["S1", "S2"], [])
    assert "kegg" not in bloco
    assert len(bloco["quality"]) == 3
