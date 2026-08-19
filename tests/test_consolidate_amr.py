"""Consenso de AMR — os tres modos de falha silenciosa achados em 2026-08-19.

Todos tinham a mesma forma: a chave de juncao nao era a mesma coisa nas tres
ferramentas, nada falhava, e `n_tools >= 2` (o filtro do relatorio) ficava
inalcancavel — o painel de AMR aparecia vazio como se fosse biologia.
"""
import consolidate_amr as ca


def _w(tmp_path, name, text):
    p = tmp_path / name
    p.write_text(text)
    return str(p)


AMRF = (
    "Protein identifier\tGene symbol\tARO Accession\tDrug Class\n"
    "binette_bin1__k141_9_1\tblaIMP\tARO:3005489\tbeta-lactam\n"
)

# O RGI devolve o cabecalho INTEIRO do Prodigal em ORF_ID.
RGI = (
    "ORF_ID\tBest_Hit_ARO\tARO\tDrug Class\n"
    "binette_bin1__k141_9_1 # 63 # 314 # -1 # ID=1_1;partial=00\tIMP-84\t3005489\tbeta-lactam\n"
)

# O argNorm 1.1.0 escreve uma linha de versao ANTES do cabecalho, e no DeepARG
# a primeira coluna e o NOME do gene, nao a proteina.
DEEP = (
    "# argNorm version: 1.1.0\n"
    "#ARG\tread_id\tARO\tpredicted_ARG-class\n"
    "IMP\tbinette_bin1__k141_9_1\t3005489\tbeta-lactam\n"
)


def test_rgi_e_amrfinder_casam_apesar_da_descricao_do_prodigal(tmp_path):
    rows = ca.consolidate(_w(tmp_path, "a.tsv", AMRF),
                          _w(tmp_path, "r.tsv", RGI),
                          _w(tmp_path, "d.tsv", "#ARG\tread_id\n"))
    assert len(rows) == 1
    assert rows[0]["locus"] == "binette_bin1__k141_9_1"
    assert rows[0]["n_tools"] == 2


def test_deeparg_entra_pelo_read_id_e_nao_pelo_nome_do_gene(tmp_path):
    rows = ca.consolidate(_w(tmp_path, "a.tsv", AMRF),
                          _w(tmp_path, "r.tsv", RGI),
                          _w(tmp_path, "d.tsv", DEEP))
    assert len(rows) == 1, "um locus so — nao um locus 'IMP' inventado"
    assert rows[0]["n_tools"] == 3
    assert float(rows[0]["consensus_score"]) == 1.0
    assert "DeepARG" in rows[0]["tools_detected"]


def test_preambulo_do_argnorm_nao_engole_o_arquivo(tmp_path):
    """Sem pular '# argNorm version:', o DictReader toma essa linha como
    cabecalho e o arquivo inteiro some sem erro nenhum."""
    hits = ca._parse_deeparg_normed(_w(tmp_path, "d.tsv", DEEP))
    assert list(hits) == ["binette_bin1__k141_9_1"]


def test_cabecalho_do_deeparg_nao_e_confundido_com_comentario(tmp_path):
    """'#ARG' comeca com '#' mas E o cabecalho — um filtro generico de '#'
    apagaria o cabecalho junto com o preambulo."""
    raw = "#ARG\tread_id\tARO\nIMP\tbinette_bin1__k141_9_1\t3005489\n"
    hits = ca._parse_deeparg_normed(_w(tmp_path, "d.tsv", raw))
    assert list(hits) == ["binette_bin1__k141_9_1"]
