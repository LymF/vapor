"""Bloco `pangenome` do report v2 — fase 1 dos clusters do catálogo.

O que estes testes protegem:

1. **O `?` não é ausência.** A regra do Snakemake se deu ao trabalho de
   codificar três estados porque um membro abaixo de 70% de completude (ou
   com falha de anotação) não é evidência de que o organismo não tem o gene.
   Ele sai do denominador da frequência, e o report não pode reintroduzi-lo.
2. **O `-` é um quarto estado do arquivo, não um terceiro do biológico.** As
   colunas do TSV são TODOS os membros de TODOS os clusters; `-` significa
   "esta coluna não pertence a este cluster" e nunca pode virar ausência.
3. **O PlasmidFinder nunca elege um cluster sozinho.** Ele é sinal de
   mobilidade em `candidates.tsv`; tratá-lo como critério mudaria a seleção
   que a pipeline fez.
"""
import os

import pytest

from report.renderer_v2 import build_pangenome


def _escreve(caminho, texto):
    os.makedirs(os.path.dirname(caminho), exist_ok=True)
    with open(caminho, 'w', encoding='utf-8') as fh:
        fh.write(texto)


@pytest.fixture
def outdir(tmp_path):
    o = str(tmp_path / "results")
    pg = os.path.join(o, "mag_catalog", "pangenome")

    _escreve(os.path.join(pg, "candidates.tsv"),
             "representative_id\tn_members\tn_islands\tn_systems\tn_args\t"
             "n_plasmid\tcriterio\teligible\n"
             "S1__bin1\t3\t1\t4\t2\t1\tilha de defesa\tTrue\n"
             "S2__bin9\t4\t0\t0\t0\t3\tsem evidencia de defesa/amr\tFalse\n"
             "S3__bin2\t2\t0\t5\t1\t0\tpoucos membros (2 < 3)\tFalse\n")

    # Colunas = TODOS os membros de TODOS os clusters. S2__bin9 nao pertence
    # ao cluster de S1__bin1, entao sua celula e '-'.
    _escreve(os.path.join(pg, "gene_by_member.tsv"),
             "# completude: S1__bin1=98.0, S1__bin7=71.5, S1__bin8=42.0, S2__bin9=90.0\n"
             "# estados: x=presente .=ausente ?=nao avaliavel -=membro nao "
             "pertence a este cluster\n"
             "cluster\ttipo\tgene\tfreq\tn_present\tn_evaluable\t"
             "S1__bin1\tS1__bin7\tS1__bin8\tS2__bin9\n"
             "S1__bin1\tdefesa\tCBASS\t2/2\t2\t2\tx\tx\t?\t-\n"
             "S1__bin1\tdefesa\ttetA\t1/2\t1\t2\tx\t.\t?\t-\n"
             "S1__bin1\tamr\ttetA\t1/2\t1\t2\t.\tx\t?\t-\n")

    _escreve(os.path.join(pg, "cluster_summary.tsv"),
             "representative_id\tn_members\tn_members_avaliaveis\tn_genes_core\t"
             "n_genes_variaveis\tn_genes_singleton\tcompletude_mediana\t"
             "tamanho_mediana_bp\tgtdb_taxonomy\n"
             "S1__bin1\t3\t2\t1\t2\t2\t84.8\t4200000\t"
             "d__Bacteria;p__Pseudomonadota;c__Gammaproteobacteria\n")
    return o


def test_sem_pangenoma_o_bloco_nao_existe(tmp_path):
    assert build_pangenome(str(tmp_path)) == {}


def test_a_matriz_preserva_os_tres_estados(outdir):
    matriz = build_pangenome(outdir)["matrix"]["S1__bin1"]
    linha = next(r for r in matriz["rows"]
                 if r["gene"] == "CBASS" and r["tipo"] == "defesa")
    assert linha["states"]["S1__bin1"] == "x"
    assert linha["states"]["S1__bin7"] == "x"
    # S1__bin8 tem 42% de completude: nao e avaliavel. Virar '.' aqui
    # afirmaria que o organismo nao tem CBASS, que o dado nao sustenta.
    assert linha["states"]["S1__bin8"] == "?"


def test_membro_de_outro_cluster_nao_entra_na_matriz(outdir):
    matriz = build_pangenome(outdir)["matrix"]["S1__bin1"]
    # '-' e uma propriedade do FORMATO (colunas globais), nao do biologico.
    assert matriz["members"] == ["S1__bin1", "S1__bin7", "S1__bin8"]
    assert all("S2__bin9" not in r["states"] for r in matriz["rows"])


def test_o_denominador_da_frequencia_exclui_o_nao_avaliavel(outdir):
    matriz = build_pangenome(outdir)["matrix"]["S1__bin1"]
    linha = next(r for r in matriz["rows"]
                 if r["gene"] == "CBASS" and r["tipo"] == "defesa")
    assert linha["n_evaluable"] == 2
    assert linha["freq"] == "2/2"
    # Tres membros no cluster, dois no denominador: a diferenca e o '?'.
    assert len(matriz["members"]) == 3
    assert sum(1 for e in linha["states"].values() if e == "?") == 1


def test_defesa_e_amr_com_o_mesmo_nome_nao_se_fundem(outdir):
    linhas = build_pangenome(outdir)["matrix"]["S1__bin1"]["rows"]
    teta = [r for r in linhas if r["gene"] == "tetA"]
    assert {r["tipo"] for r in teta} == {"defesa", "amr"}
    # Fundir as duas daria um perfil de presenca que nenhuma das duas tem.
    assert next(r for r in teta if r["tipo"] == "defesa")["states"]["S1__bin1"] == "x"
    assert next(r for r in teta if r["tipo"] == "amr")["states"]["S1__bin1"] == "."


def test_completude_do_cabecalho_viaja_com_a_matriz(outdir):
    bloco = build_pangenome(outdir)
    # "1/2" so e interpretavel com a completude a vista -- e por isso que a
    # regra a escreve no cabecalho em vez de deixa-la noutro arquivo.
    assert bloco["completeness"]["S1__bin8"] == 42.0
    assert bloco["completeness"]["S1__bin1"] == 98.0


def test_candidatos_trazem_o_criterio_e_os_recusados(outdir):
    candidatos = build_pangenome(outdir)["candidates"]
    por_rep = {c["representative"]: c for c in candidatos}
    # A selecao tem de ser auditavel: o recusado aparece COM o motivo, nao
    # some da tabela.
    assert por_rep["S3__bin2"]["criterio"] == "poucos membros (2 < 3)"
    assert por_rep["S3__bin2"]["eligible"] is False
    assert por_rep["S1__bin1"]["eligible"] is True


def test_plasmidfinder_sozinho_nao_elege_cluster(outdir):
    candidatos = build_pangenome(outdir)["candidates"]
    so_plasmidio = next(c for c in candidatos if c["representative"] == "S2__bin9")
    assert so_plasmidio["n_plasmid"] == 3
    # Plasmidio sem defesa nem ARG nao motiva um pangenoma: e sinal de
    # mobilidade, nunca criterio.
    assert so_plasmidio["eligible"] is False
    assert so_plasmidio["criterio"] == "sem evidencia de defesa/amr"


def test_sumario_por_cluster_traz_core_e_variaveis(outdir):
    clusters = build_pangenome(outdir)["clusters"]
    c = clusters[0]
    assert c["representative"] == "S1__bin1"
    assert c["n_core"] == 1
    assert c["n_variable"] == 2
    assert c["n_evaluable"] == 2
    assert c["taxonomy"].endswith("c__Gammaproteobacteria")


def test_matriz_ausente_nao_derruba_os_candidatos(outdir):
    os.remove(os.path.join(outdir, "mag_catalog", "pangenome",
                           "gene_by_member.tsv"))
    bloco = build_pangenome(outdir)
    assert "matrix" not in bloco
    assert len(bloco["candidates"]) == 3
