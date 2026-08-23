"""Blocos `diversity` e `reads` do report v2.

Três coisas que já falharam em silêncio e que estes testes fixam:

1. **Simpson e Chao1 saem VAZIOS quando não há contagens de reads** (ambos
   são estimadores de contagem — f1/f2 e a*(a-1) — e `compute_diversity.py`
   se recusa a calculá-los sobre RPKM). Vazio tem de virar lacuna declarada,
   nunca 0.0 plotado.
2. **O cabeçalho da tabela do sylph é o CAMINHO do arquivo de reads**, que em
   paired-end é o `-1` (`{sample}_R1_fastp.fq.gz`). Resolver por extensão
   zerava a trilha inteira com `has_data` verdadeiro.
3. **O espaço de IDs do sylph não conversa com o dos contigs montados.** O
   bloco carrega esse aviso, porque nenhum join entre as duas trilhas é
   válido.
"""
import os

import pytest

from report.renderer_v2 import build_diversity, build_reads


def _escreve(caminho, texto):
    os.makedirs(os.path.dirname(caminho), exist_ok=True)
    with open(caminho, 'w', encoding='utf-8') as fh:
        fh.write(texto)


@pytest.fixture
def outdir(tmp_path):
    o = str(tmp_path / "results")
    div = os.path.join(o, "diversity")

    # simpson/chao1 vazios: rodada sem contagens de reads.
    _escreve(os.path.join(div, "alpha_diversity.tsv"),
             "sample\tdomain\trichness\tshannon\tsimpson\tchao1\n"
             "S1\tviral\t120\t3.41\t\t\n"
             "S2\tviral\t95\t3.02\t\t\n"
             "S1\tprokaryotic\t340\t4.80\t\t\n")

    _escreve(os.path.join(div, "beta_pcoord_viral.tsv"),
             "sample\tPC1\tPC2\tPC1_var\tPC2_var\n"
             "S1\t0.31\t-0.12\t42.5\t18.3\n"
             "S2\t-0.29\t0.15\t42.5\t18.3\n")

    _escreve(os.path.join(div, "procrustes_coords.tsv"),
             "sample\tviral_PC1\tviral_PC2\tprok_PC1\tprok_PC2\tdisparity\n"
             "S1\t0.30\t-0.10\t0.28\t-0.14\t0.214\n"
             "S2\t-0.28\t0.13\t-0.26\t0.16\t0.214\n")
    return o


@pytest.fixture
def outdir_reads(tmp_path):
    o = str(tmp_path / "reads_run")
    rc = os.path.join(o, "reads_classify")
    # Cabecalho = caminho do arquivo de reads que o sylph recebeu. Em
    # paired-end e o `-1`, ja aparado pelo fastp.
    _escreve(os.path.join(rc, "otu_table.tsv"),
             "#OTU_ID\tS1_R1_fastp.fq.gz\tS2_R1_fastp.fq.gz\n"
             "r__Duplodnaviria|k__Heunggongvirae|p__Uroviricota\t12.5\t8.0\n"
             "d__Bacteria|p__Pseudomonadota\t40.0\t35.5\n")
    _escreve(os.path.join(rc, "viral_abundance_by_host.tsv"),
             "host_genus\thost_source\tn_viral_taxa\tS1_R1_fastp.fq.gz\t"
             "S2_R1_fastp.fq.gz\n"
             "Escherichia\tdb\t4\t9.0\t5.0\n"
             "desconhecido\tnone\t2\t3.5\t3.0\n")
    return o


def test_sem_diversidade_o_bloco_nao_existe(tmp_path):
    assert build_diversity(str(tmp_path), ["S1"]) == {}


def test_simpson_e_chao1_vazios_viram_lacuna_declarada(outdir):
    bloco = build_diversity(outdir, ["S1", "S2"])
    indices = {r["index"] for r in bloco["alpha"]}
    # Nunca 0.0: zero e um valor de diversidade, e Simpson 0 significaria
    # "uma unica especie domina completamente".
    assert "simpson" not in indices
    assert "chao1" not in indices
    assert set(bloco["alpha_missing"]) == {"simpson", "chao1"}
    assert {"observed", "shannon"} <= indices


def test_pcoa_carrega_a_variancia_dos_dois_eixos(outdir):
    pcoa = build_diversity(outdir, ["S1", "S2"])["pcoa"]["viral"]
    # Sem a % de variancia, a distancia entre dois pontos num PCoA nao tem
    # interpretacao -- e a unica coisa que diz quanto do sinal o plano
    # explica.
    assert pcoa[0]["var_pc1"] == pytest.approx(0.425)
    assert pcoa[0]["var_pc2"] == pytest.approx(0.183)


def test_procrustes_pareia_as_duas_ordenacoes(outdir):
    proc = build_diversity(outdir, ["S1", "S2"])["procrustes"]
    assert proc["disparity"] == pytest.approx(0.214)
    par = next(p for p in proc["pairs"] if p["sample"] == "S1")
    # A seta liga o ponto viral ao procariotico da MESMA amostra -- e o par
    # que e a unidade, nao o ponto.
    assert (par["viral_pc1"], par["prok_pc1"]) == (0.30, 0.28)


def test_pcoa_de_uma_trilha_desligada_nao_inventa_eixo(outdir):
    bloco = build_diversity(outdir, ["S1", "S2"])
    # Só a ordenacao viral existe nesta rodada.
    assert set(bloco["pcoa"]) == {"viral"}


def test_reads_resolve_a_coluna_de_paired_end_para_a_amostra(outdir_reads):
    bloco = build_reads(outdir_reads, ["S1", "S2"])
    viral = bloco["viral"]
    assert viral, "a trilha de reads nao pode sair vazia com o arquivo presente"
    # '{sample}_R1_fastp' resolvido para '{sample}'. Sem isto a coluna nao
    # casa com amostra nenhuma e TODA a trilha sai zerada -- com has_data
    # verdadeiro, porque as linhas existem.
    assert viral[0]["S1"] == pytest.approx(12.5)
    assert viral[0]["S2"] == pytest.approx(8.0)


def test_reads_separa_viral_de_procariotico_pelo_prefixo_de_reino(outdir_reads):
    bloco = build_reads(outdir_reads, ["S1", "S2"])
    assert all("r__" in r["clade"] for r in bloco["viral"])
    assert all("r__" not in r["clade"] for r in bloco["prok"])


def test_host_source_viaja_com_o_hospedeiro(outdir_reads):
    hosts = build_reads(outdir_reads, ["S1", "S2"])["host"]
    por_genero = {h["host_genus"]: h for h in hosts}
    # Sem host_source, "hospedeiro desconhecido" e "hospedeiro anotado pelo
    # banco" ficam indistinguiveis no grafico.
    assert por_genero["Escherichia"]["host_source"] == "db"
    assert por_genero["desconhecido"]["host_source"] == "none"


def test_reads_declara_que_seus_ids_nao_casam_com_os_contigs(outdir_reads):
    bloco = build_reads(outdir_reads, ["S1", "S2"])
    # O aviso e dado, nao decoracao: 't__IMGVR_UViG_...' e 'k141_...' sao
    # espacos de ID distintos, e nenhum join entre as trilhas e valido.
    assert bloco["id_space_warning"]
