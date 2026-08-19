"""Catalogo global de MAGs procarioticos (2026-08-19)."""
import os
import pytest
from mag_catalog import (
    namespaced_id, split_namespaced_id, strip_ext, collect_bins, build_pool,
    merge_checkm2, parse_galah_clusters, representative_view, read_provenance,
    NamespaceCollision,
)


# ── namespace ────────────────────────────────────────────────────────────

def test_nomes_de_bin_colidem_entre_amostras_por_isso_o_prefixo():
    """Binette emite "binette_bin1" em TODA amostra e o VAMB emite numeros
    nus. Sem prefixo, juntar num diretorio so sobrescreveria organismos
    diferentes."""
    a = namespaced_id("P01_RNG_08_947", "binette_bin1")
    b = namespaced_id("P02_SOL_01_949", "binette_bin1")
    assert a != b
    assert split_namespaced_id(a) == ("P01_RNG_08_947", "binette_bin1")


def test_source_id_com_separador_falha_alto():
    """Um source_id com '__' faria a divisao devolver o pedaco errado, e a
    anotacao seria atribuida a amostra errada em silencio."""
    with pytest.raises(NamespaceCollision):
        namespaced_id("amostra__estranha", "bin1")


def test_split_usa_o_primeiro_separador():
    """Nome de bin pode conter '__'; source_id nao pode."""
    assert split_namespaced_id("S1__bin__x") == ("S1", "bin__x")


def test_split_de_id_sem_namespace_nao_inventa_fonte():
    assert split_namespaced_id("bin1") == ("", "bin1")


def test_strip_ext_cobre_fa_fna_fasta():
    assert strip_ext("/x/binette_bin1.fa") == "binette_bin1"
    assert strip_ext("/x/1.fna") == "1"
    assert strip_ext("/x/g.fasta") == "g"
    assert strip_ext("/x/sem_ext") == "sem_ext"


# ── pool ─────────────────────────────────────────────────────────────────

def _mk_bins(tmp_path, name, files, ext=".fa"):
    d = tmp_path / name
    d.mkdir(parents=True)
    for f in files:
        (d / (f + ext)).write_text(">c1\nACGT\n")
    return str(d)


def test_pool_namespaceia_e_registra_procedencia(tmp_path):
    s1 = _mk_bins(tmp_path, "S1", ["binette_bin1", "binette_bin2"])
    g1 = _mk_bins(tmp_path, "G1", ["1", "136"], ext=".fna")
    sources = [("sample", "S1", s1, ".fa"), ("group", "G1", g1, ".fna")]

    out = tmp_path / "genomes"
    prov = tmp_path / "provenance.tsv"
    stats = build_pool(sources, str(out), str(prov))

    assert stats == {"n_genomes": 4, "n_sources": 2}
    assert sorted(os.listdir(out)) == [
        "G1__1.fa", "G1__136.fa", "S1__binette_bin1.fa", "S1__binette_bin2.fa"]
    rows = read_provenance(str(prov))
    assert ("G1__136", "group", "G1", "136") in rows


def test_pool_normaliza_a_extensao(tmp_path):
    """VAMB emite .fna e Binette .fa; o pool uniformiza para .fa, senao todo
    consumidor a jusante precisaria de dois globs."""
    g1 = _mk_bins(tmp_path, "G1", ["1"], ext=".fna")
    build_pool([("group", "G1", g1, ".fna")], str(tmp_path / "g"),
               str(tmp_path / "p.tsv"))
    assert os.listdir(tmp_path / "g") == ["G1__1.fa"]


def test_pool_e_reexecutavel(tmp_path):
    """Symlink pre-existente nao pode abortar um resume do Snakemake."""
    s1 = _mk_bins(tmp_path, "S1", ["b1"])
    src = [("sample", "S1", s1, ".fa")]
    build_pool(src, str(tmp_path / "g"), str(tmp_path / "p.tsv"))
    stats = build_pool(src, str(tmp_path / "g"), str(tmp_path / "p.tsv"))
    assert stats["n_genomes"] == 1


def test_diretorio_ausente_nao_e_erro(tmp_path):
    """Uma amostra sem MAG e resultado real, nao falha."""
    stats = build_pool([("sample", "X", str(tmp_path / "naoexiste"), ".fa")],
                       str(tmp_path / "g"), str(tmp_path / "p.tsv"))
    assert stats["n_genomes"] == 0


# ── checkm2 ──────────────────────────────────────────────────────────────

def _mk_checkm2(tmp_path, name, rows):
    p = tmp_path / (name + ".tsv")
    with open(p, "w") as fh:
        fh.write("Name\tCompleteness\tContamination\n")
        for n, c, ct in rows:
            fh.write(f"{n}\t{c}\t{ct}\n")
    return str(p)


def test_merge_checkm2_reescreve_o_nome_para_o_id_namespaced(tmp_path):
    """O galah casa o relatorio pelo nome do genoma. Como o pool renomeia
    tudo, concatenar os relatorios crus deixaria o galah sem qualidade para
    NENHUM genoma -- e ele escolheria representante por outro criterio, sem
    avisar."""
    a = _mk_checkm2(tmp_path, "S1", [("binette_bin1", "95.0", "1.0")])
    b = _mk_checkm2(tmp_path, "S2", [("binette_bin1", "80.0", "3.0")])
    out = tmp_path / "merged.tsv"
    n_rows, n_src = merge_checkm2([("S1", a), ("S2", b)], str(out))

    assert (n_rows, n_src) == (2, 2)
    names = [l.split("\t")[0] for l in open(out).read().splitlines()[1:]]
    assert names == ["S1__binette_bin1", "S2__binette_bin1"]


def test_merge_checkm2_pula_fonte_ausente_ou_vazia(tmp_path):
    a = _mk_checkm2(tmp_path, "S1", [("b1", "90", "1")])
    (tmp_path / "vazio.tsv").write_text("")
    n_rows, n_src = merge_checkm2(
        [("S1", a), ("S2", str(tmp_path / "vazio.tsv")),
         ("S3", str(tmp_path / "naoexiste.tsv"))], str(tmp_path / "m.tsv"))
    assert (n_rows, n_src) == (1, 1)


def test_merge_checkm2_sem_fonte_ainda_escreve_cabecalho(tmp_path):
    out = tmp_path / "m.tsv"
    merge_checkm2([], str(out))
    assert open(out).read().startswith("Name\t")


# ── clusters e vista ─────────────────────────────────────────────────────

def test_parse_galah_converte_caminho_em_id(tmp_path):
    """O galah escreve CAMINHOS, nao ids."""
    p = tmp_path / "c.tsv"
    p.write_text("/pool/S1__bin1.fa\t/pool/S1__bin1.fa\n"
                 "/pool/S1__bin1.fa\t/pool/S2__bin7.fa\n")
    m = parse_galah_clusters(str(p))
    assert m == {"S1__bin1": "S1__bin1", "S2__bin7": "S1__bin1"}


def test_parse_galah_ignora_cabecalho_do_fallback(tmp_path):
    p = tmp_path / "c.tsv"
    p.write_text("representative\tmember\n/pool/S1__b.fa\t/pool/S1__b.fa\n")
    assert parse_galah_clusters(str(p)) == {"S1__b": "S1__b"}


def test_vista_liga_representante_a_amostra_de_origem():
    """O ponto do catalogo: analise no representante, mas sem perder de qual
    amostra cada MAG veio."""
    clusters = {"S1__b1": "S1__b1", "S2__b7": "S1__b1"}
    prov = [("S1__b1", "sample", "S1", "b1"), ("S2__b7", "sample", "S2", "b7")]
    view = representative_view(clusters, prov)
    assert ("S2", "b7", "S2__b7", "S1__b1") in view


def test_membro_sem_cluster_e_seu_proprio_representante():
    """Singleton, ou galah que falhou: nao pode sumir da vista."""
    view = representative_view({}, [("S1__b1", "sample", "S1", "b1")])
    assert view == [("S1", "b1", "S1__b1", "S1__b1")]
