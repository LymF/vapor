import pytest
from report.schema import Block, project, project_strict, UndeclaredField


def test_project_mantem_apenas_campos_declarados():
    b = Block(name="mag", fields=("Bin", "Completeness"), key="Bin")
    rows = [{"Bin": "bin1", "Completeness": "95.1", "other_related_references": "x" * 5000}]
    assert project(b, rows) == [{"Bin": "bin1", "Completeness": "95.1"}]


def test_project_preenche_campo_ausente_com_vazio():
    b = Block(name="mag", fields=("Bin", "Contamination"), key="Bin")
    assert project(b, [{"Bin": "bin1"}]) == [{"Bin": "bin1", "Contamination": ""}]


def test_project_preserva_a_ordem_declarada():
    b = Block(name="mag", fields=("Bin", "Completeness"), key="Bin")
    out = project(b, [{"Completeness": "9", "Bin": "b"}])
    assert list(out[0].keys()) == ["Bin", "Completeness"]


def test_project_strict_recusa_campo_nao_declarado():
    b = Block(name="mag", fields=("Bin",), key="Bin")
    with pytest.raises(UndeclaredField) as e:
        project_strict(b, [{"Bin": "b", "surpresa": 1}])
    assert "surpresa" in str(e.value)
