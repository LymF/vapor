"""Bloco `defense_amr` do report v2 — defesa, AMR e plasmídeos no catálogo.

A armadilha central destes testes é uma só, e ela já causou bug real na
pipeline: um ID de proteína do catálogo é `{source}__{bin}__{orf}`, ou seja
`S1__binette_bin1__k141_1_5`. Cortar no PRIMEIRO `__` devolve `S1` e atribui
todo achado de AMR à AMOSTRA em vez do MAG. As vistas do Snakemake reescrevem
o prefixo antes das ferramentas a jusante lerem; aqui, que lê o catálogo
direto, o corte tem de ser feito contra os representantes conhecidos.

A segunda armadilha é a colocalização: nomes de contig colidem entre MAGs
(todo assembly do MEGAHIT emite `k141_1`). Cruzar ARG com replicon pelo nome
do contig sem exigir o MESMO genoma afirmaria colocalização entre organismos
diferentes.
"""
import os

import pytest

from report.renderer_v2 import build_defense_amr, _contig_do_orf


def _escreve(caminho, texto):
    os.makedirs(os.path.dirname(caminho), exist_ok=True)
    with open(caminho, 'w', encoding='utf-8') as fh:
        fh.write(texto)


@pytest.fixture
def outdir(tmp_path):
    o = str(tmp_path / "results")
    cat = os.path.join(o, "mag_catalog")

    _escreve(os.path.join(cat, "mag_membership.tsv"),
             "source_id\toriginal_bin_id\tmember_id\trepresentative_id\n"
             "S1\tbinette_bin1\tS1__binette_bin1\tS1__binette_bin1\n"
             "S2\tbinette_bin1\tS2__binette_bin1\tS1__binette_bin1\n"
             "S2\tbinette_bin2\tS2__binette_bin2\tS2__binette_bin2\n")

    # `protein_in_syst` traz IDs NUS (`k141_1_2`), nao prefixados: o
    # DefenseFinder roda uma vez por genoma, sobre o .faa daquele genoma, e o
    # prefixo `{genome}__` so existe no .faa CONCATENADO que as ferramentas
    # de AMR (nivel de gene) consomem. Confundir os dois faria a ilha de
    # defesa nunca casar com nenhum gene.
    _escreve(os.path.join(cat, "defensefinder", "defensefinder_systems.tsv"),
             "genome\tsys_id\ttype\tsubtype\tgenes_count\tprotein_in_syst\n"
             "S1__binette_bin1\tCBASS_1\tCBASS\tCBASS_I\t2\tk141_1_1,k141_1_2\n"
             "S1__binette_bin1\tRM_1\tRM\tRM_I\t1\tk141_1_3\n"
             "S1__binette_bin1\tWadjet_1\tWadjet\tWadjet_I\t2\tk141_1_4,k141_1_5\n"
             "S2__binette_bin2\tGabija_1\tGabija\tGabija\t1\tk141_1_3\n")

    faa = os.path.join(cat, "proteins", "S1__binette_bin1.faa")
    _escreve(faa, "".join(
        f">k141_1_{i} # {i * 1000} # {i * 1000 + 800} # 1 # ID={i}\nMA\n"
        for i in range(1, 7)))
    _escreve(os.path.join(cat, "proteins", "manifest.txt"),
             f"S1__binette_bin1\tprodigal\t/x/S1.fna\t{faa}\t/x/S1.gff\n")

    _escreve(os.path.join(cat, "amr_consensus", "amr_consensus.tsv"),
             "locus\taro_accession\tgene_name\tdrug_class\t"
             "resistance_mechanism\tn_tools\tconsensus_score\ttools_detected\n"
             "S1__binette_bin1__k141_1_5\tARO:1\ttetA\ttetracycline\tefflux\t"
             "3\t1.0\tamrfinder,rgi,deeparg\n"
             "S2__binette_bin2__k141_1_9\tARO:2\tblaTEM\tbeta-lactam\t"
             "inactivation\t2\t0.67\tamrfinder,rgi\n"
             "S2__binette_bin2__k141_7_1\tARO:3\tsul1\tsulfonamide\ttarget\t"
             "1\t0.33\trgi\n")

    _escreve(os.path.join(cat, "abricate", "plasmidfinder_results.tsv"),
             "#FILE\tSEQUENCE\tSTART\tEND\tSTRAND\tGENE\t%IDENTITY\n"
             "/x/genomes/S1__binette_bin1.fa\tk141_1\t100\t900\t+\tIncFII\t99.1\n"
             "/x/genomes/S2__binette_bin2.fa\tk141_3\t50\t400\t+\tIncQ\t98.0\n")
    return o


def test_sem_catalogo_o_bloco_nao_existe(tmp_path):
    assert build_defense_amr(str(tmp_path)) == {}


def test_amr_atribui_o_hit_ao_MAG_e_nunca_a_amostra(outdir):
    amr = build_defense_amr(outdir)["amr"]
    genomas = {r["genome"] for r in amr}
    # 'S1' aqui significaria "todo ARG da amostra", que e a leitura errada e
    # silenciosa: o corte no primeiro '__' devolveria exatamente isso.
    assert "S1" not in genomas
    assert genomas == {"S1__binette_bin1", "S2__binette_bin2"}


def test_amr_exige_consenso_de_duas_ferramentas(outdir):
    amr = build_defense_amr(outdir)["amr"]
    # sul1 saiu de uma ferramenta so: fica de fora do consenso, e o painel
    # nunca o mostra como ARG confirmado.
    assert {r["gene"] for r in amr} == {"tetA", "blaTEM"}
    assert all(r["n_tools"] >= 2 for r in amr)


def test_contig_do_orf_corta_no_ultimo_underscore(outdir):
    # Convencao do Prodigal: '{contig}_{n}'. O nome do contig contem '_'
    # ('k141_1'), entao cortar no PRIMEIRO devolveria 'k141'.
    assert _contig_do_orf("k141_1_5") == "k141_1"
    assert _contig_do_orf("k141_1") == "k141"
    assert _contig_do_orf("") == ""


def test_colocalizacao_exige_o_mesmo_genoma_e_o_mesmo_contig(outdir):
    coloc = build_defense_amr(outdir)["colocalization"]
    # tetA esta em S1__binette_bin1/k141_1, que carrega o replicon IncFII:
    # colocalizado. blaTEM esta em S2__binette_bin2/k141_1 -- MESMO nome de
    # contig, genoma diferente, e o replicon desse genoma esta em k141_3.
    # Contar blaTEM afirmaria colocalizacao entre organismos diferentes.
    assert coloc["n_args"] == 2
    assert coloc["n_args_on_replicon"] == 1
    assert coloc["args_on_replicon"] == ["S1__binette_bin1|k141_1|tetA"]


def test_defesa_conta_sistemas_por_genoma(outdir):
    defense = build_defense_amr(outdir)["defense"]
    por_genoma = {}
    for r in defense:
        por_genoma.setdefault(r["genome"], {})[r["system"]] = r["count"]
    assert por_genoma["S1__binette_bin1"] == {"CBASS": 1, "RM": 1, "Wadjet": 1}
    assert por_genoma["S2__binette_bin2"] == {"Gabija": 1}


def test_ilha_de_defesa_sai_com_coordenadas_reais(outdir):
    ilhas = build_defense_amr(outdir)["islands"]
    assert len(ilhas) == 1
    ilha = ilhas[0]
    assert ilha["genome"] == "S1__binette_bin1"
    assert ilha["contig"] == "k141_1"
    assert ilha["n_systems"] == 3
    # Coordenadas em bp vindas do cabeçalho do Prodigal: a trilha genômica
    # desenha em base, não em ordem de gene.
    assert ilha["start"] == 1000
    assert ilha["end"] == 5800
    assert {g["label"] for g in ilha["genes"]} >= {"CBASS", "RM", "Wadjet"}


def test_ilha_entra_no_upset_como_evidencia_propria(outdir):
    upset = build_defense_amr(outdir)["upset"]
    assert upset["sets"]["ilha de defesa"] == 1


def test_upset_agrupa_MAGs_por_evidencia(outdir):
    upset = build_defense_amr(outdir)["upset"]
    # S1__binette_bin1 tem os tres sinais possiveis menos ilha; S2__binette_bin2
    # tambem tem replicon e ARG. A pergunta do UpSet e "quantos MAGs carregam
    # esta combinacao", nunca "quantos hits".
    assert upset["sets"]["replicon"] == 2
    assert upset["sets"]["ARG de consenso"] == 2
    combos = {tuple(sorted(c["tools"])): c["count"] for c in upset["combos"]}
    # S2__binette_bin2 carrega os tres sinais e nenhuma ilha; S1__binette_bin1
    # carrega os tres MAIS a ilha, e por isso cai numa combinacao propria --
    # e a leitura que o UpSet existe para dar.
    assert combos[("ARG de consenso", "replicon", "sistema de defesa")] == 1
    assert combos[("ARG de consenso", "ilha de defesa", "replicon",
                   "sistema de defesa")] == 1


def test_plasmidfinder_ausente_nao_derruba_o_bloco(outdir):
    os.remove(os.path.join(outdir, "mag_catalog", "abricate",
                           "plasmidfinder_results.tsv"))
    bloco = build_defense_amr(outdir)
    assert "plasmids" not in bloco
    # Sem replicon nao existe pergunta de colocalizacao -- a chave some, e o
    # painel mostra lacuna em vez de "0% dos ARGs em plasmidio", que seria
    # uma afirmacao que o dado nao sustenta.
    assert "colocalization" not in bloco
    assert len(bloco["amr"]) == 2
