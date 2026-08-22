import os
import pytest
from report.renderer_v2 import render_html, write_report, _quebra_por_tier, _etapas
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
    etapas = _etapas({"contigs": 10, "candidatos virais": 0, "vOTUs retidos": 3})
    assert [e["name"] for e in etapas] == ["contigs", "vOTUs retidos"]
