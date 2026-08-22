import json
import os
import pytest
from report.renderer_v2 import (
    render_html, write_report, _quebra_por_tier, _etapas,
    _funil_da_amostra, _contagem_votus_catalogo,
    build_sequencing, build_viral, _limita_explorer,
)
from report.schema import PayloadOverBudget

DADOS = {"run": {"title": "VAPOR", "samples": ["S1"]}, "overview": {"kpis": []}}


def _assets(tmp_path):
    a = tmp_path / "assets"; a.mkdir()
    (a / "report-ui.js").write_text("console.log('bundle');", encoding='utf-8')
    (a / "report-ui.css").write_text(".app{color:red}", encoding='utf-8')
    c = tmp_path / "components"; c.mkdir()
    (c / "shell_v2.html").write_text(
        "<html><head><style>{{CSS}}</style></head>"
        "<body><div id=\"vapor-root\"></div>{{DATA_JSON}}<script>{{APP_JS}}</script></body></html>",
        encoding='utf-8')
    return str(a), str(c)


def test_render_inlineia_bundle_css_e_dados(tmp_path):
    a, c = _assets(tmp_path)
    html = render_html(DADOS, a, c)
    assert "console.log('bundle');" in html
    assert ".app{color:red}" in html
    assert '"title": "VAPOR"' in html or '"title":"VAPOR"' in html
    assert "{{" not in html


def test_render_escapa_fechamento_de_script(tmp_path):
    a, c = _assets(tmp_path)
    html = render_html({"run": {"title": "</script><script>alert(1)"}}, a, c)
    assert "</script><script>alert(1)" not in html
    assert "<\\/script>" in html


def test_write_report_recusa_payload_acima_do_orcamento(tmp_path):
    a, c = _assets(tmp_path)
    gordo = {"run": {"title": "VAPOR"}, "lixo": ["x" * 300_000]}
    with pytest.raises(PayloadOverBudget):
        write_report(gordo, str(tmp_path / "r.html"), a, c, limit_mb=0.1)
    assert not os.path.exists(tmp_path / "r.html")


def test_write_report_escreve_o_arquivo(tmp_path):
    a, c = _assets(tmp_path)
    destino = str(tmp_path / "sub" / "r.html")
    assert write_report(DADOS, destino, a, c) == destino
    assert os.path.getsize(destino) > 0


def test_quebra_por_tier_separa_nunca_avaliado_de_tier_baixo(tmp_path):
    tsv = tmp_path / "viral_discarded.tsv"
    tsv.write_text(
        "contig_id\tlength\tcheckv_quality\tcheckv_completeness\tin_vrhyme_bin\tsource_id\n"
        "k141_1\t1200\tLow-quality\t12.0\tFalse\tS1\n"
        "k141_2\t900\t\t\tFalse\tS1\n"
        "k141_3\t800\tLow-quality\t8.0\tFalse\tS1\n",
        encoding='utf-8')
    quebra = {d["reason"]: d["count"] for d in _quebra_por_tier(str(tsv))}
    assert quebra == {"Low-quality": 2, "sem avaliação CheckV": 1}


def test_quebra_por_tier_sem_arquivo_e_lista_vazia(tmp_path):
    assert _quebra_por_tier(str(tmp_path / "nao_existe.tsv")) == []


def test_etapas_omite_a_etapa_cuja_fonte_nao_existe():
    etapas = _etapas({
        "contigs": 10, "candidatos virais": None,
        "sequências virais retidas": 3,
    })
    assert [e["name"] for e in etapas] == ["contigs", "sequências virais retidas"]


def test_etapas_mantem_contagem_zero_real_quando_fonte_existe():
    # fonte AUSENTE (None) some; fonte PRESENTE com zero biologico aparece com
    # valor 0 -- e a distincao da correcao 5.
    etapas = _etapas({
        "contigs": 10, "candidatos virais": 0,
        "sequências virais retidas": None,
    })
    nomes = [e["name"] for e in etapas]
    assert nomes == ["contigs", "candidatos virais"]
    valor = next(e["value"] for e in etapas if e["name"] == "candidatos virais")
    assert valor == 0


def test_etapas_inclui_a_unidade_de_cada_uma():
    etapas = _etapas({
        "contigs": 5, "candidatos virais": 5, "sequências virais retidas": 5,
    })
    unidades = {e["name"]: e["unit"] for e in etapas}
    assert unidades == {
        "contigs": "contig",
        "candidatos virais": "contig",
        "sequências virais retidas": "sequência",
    }


def test_funil_le_contagem_de_contigs_do_quast_report(tmp_path):
    # parse_quast_all devolve {rotulo: {metrica: valor}} -- a metrica fica um
    # nivel abaixo do rotulo do assembly. Este e o bug da correcao 1: sem o
    # nivel extra, n_contigs era sempre 0 e "contigs" nunca aparecia.
    sample = "S1"
    quast_dir = tmp_path / sample / "quast"
    quast_dir.mkdir(parents=True)
    (quast_dir / "report.tsv").write_text(
        "Assembly\tassembly\n# contigs\t1234\n", encoding='utf-8')
    funil = _funil_da_amostra(str(tmp_path), sample)
    etapa = next((e for e in funil["stages"] if e["name"] == "contigs"), None)
    assert etapa is not None
    assert etapa["value"] == 1234


def test_contagem_votus_catalogo_conta_ids_distintos(tmp_path):
    caminho = tmp_path / "votu_catalog" / "vOTU_clusters.tsv"
    caminho.parent.mkdir(parents=True)
    caminho.write_text(
        "votu_id\trepresentative\tmember\n"
        "vOTU_1\tk141_1\tk141_1\n"
        "vOTU_1\tk141_1\tk141_2\n"
        "vOTU_2\tk141_9\tk141_9\n",
        encoding='utf-8')
    assert _contagem_votus_catalogo(str(tmp_path)) == 2


def test_contagem_votus_catalogo_ausente_e_none(tmp_path):
    assert _contagem_votus_catalogo(str(tmp_path)) is None


def test_sequencing_ausente_omite_a_chave_em_vez_de_lista_vazia(tmp_path):
    # Fonte ausente e "trilha desligada", nao "medi zero" -- o painel precisa
    # distinguir os dois casos, entao a chave nao pode existir vazia.
    saida = build_sequencing(str(tmp_path), ["S1"])
    assert "quast" not in saida


def test_sequencing_le_quast_no_nivel_certo(tmp_path):
    d = tmp_path / "S1" / "quast"
    d.mkdir(parents=True)
    (d / "report.tsv").write_text(
        "Assembly\tassembly\n# contigs\t1234\nN50\t5678\n", encoding='utf-8')
    saida = build_sequencing(str(tmp_path), ["S1"])
    assert saida["quast"]["S1"]["# contigs"] == "1234"


def test_viral_projeta_somente_campos_declarados(tmp_path):
    d = tmp_path / "S1" / "viral" / "checkv"
    d.mkdir(parents=True)
    (d / "quality_summary.tsv").write_text(
        "contig_id\tcheckv_quality\tcompleteness\tcampo_gigante\n"
        "k141_1\tHigh-quality\t95.0\t" + "x" * 5000 + "\n", encoding='utf-8')
    saida = build_viral(str(tmp_path), ["S1"])
    tiers = saida["checkv_tiers"]["S1"]
    assert tiers == {"High-quality": 1}
    assert "campo_gigante" not in json.dumps(saida)


def test_explorer_limita_a_50_votus():
    feats = [{"votu_id": f"v{i}", "length": i, "features": []} for i in range(200)]
    assert len(_limita_explorer(feats)) == 50
    assert _limita_explorer(feats)[0]["length"] == 199   # os mais longos primeiro
