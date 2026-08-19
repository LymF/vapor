"""Alfa-diversidade: Simpson e Chao1 sao estimadores de CONTAGEM."""
import compute_diversity as cd


def test_chao1_usa_singletons_e_doubletons():
    # 3 taxa observados, f1=1, f2=1 -> 3 + 1/(2*1) = 3.5
    assert cd.chao1([1, 2, 10]) == 3.5


def test_simpson_de_uma_comunidade_uniforme_e_alto():
    uniforme = cd.simpson([100, 100, 100, 100])
    dominada = cd.simpson([397, 1, 1, 1])
    assert uniforme > 0.7 > dominada


def test_alpha_sem_contagens_deixa_simpson_e_chao1_vazios():
    """Sem a coluna de contagem, os dois indices NAO podem ser calculados
    sobre RPKM -- sair vazio e o comportamento correto, porque um numero
    errado passaria por indice."""
    feats = ["v1", "v2"]
    data = {"S1": {"v1": 12.5, "v2": 0.4}}
    row = cd.compute_alpha(feats, data, "viral")[0]
    assert row["simpson"] == "" and row["chao1"] == ""
    assert row["richness"] == 2
    assert row["shannon"] > 0


def test_alpha_com_contagens_calcula_os_quatro():
    feats = ["v1", "v2", "v3"]
    data   = {"S1": {"v1": 12.5, "v2": 0.4, "v3": 3.0}}
    counts = {"S1": {"v1": 100.0, "v2": 1.0, "v3": 2.0}}
    row = cd.compute_alpha(feats, data, "viral", counts)[0]
    assert row["chao1"] == cd.chao1([100.0, 1.0, 2.0])
    assert row["simpson"] == cd.simpson([100.0, 1.0, 2.0])


def test_richness_e_shannon_seguem_a_metrica_configurada():
    """Riqueza e Shannon sao funcoes de proporcao e ficam certos sobre a
    metrica normalizada -- nao devem mudar por existirem contagens."""
    feats = ["v1", "v2"]
    data   = {"S1": {"v1": 10.0, "v2": 10.0}}
    counts = {"S1": {"v1": 5.0, "v2": 500.0}}
    a = cd.compute_alpha(feats, data, "viral")[0]
    b = cd.compute_alpha(feats, data, "viral", counts)[0]
    assert a["shannon"] == b["shannon"] and a["richness"] == b["richness"]


def test_loader_do_relatorio_nao_transforma_vazio_em_zero(tmp_path):
    """Um indice nao calculado nao pode virar 0.0 no grafico."""
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "scripts"))
    from report.data_loaders import load_alpha_diversity
    p = tmp_path / "alpha_diversity.tsv"
    p.write_text("sample\tdomain\trichness\tshannon\tsimpson\tchao1\n"
                 "S1\tviral\t10\t1.5\t\t\n")
    rows = load_alpha_diversity(str(p))
    idx = {r["index"] for r in rows}
    assert idx == {"observed", "shannon"}
