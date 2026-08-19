"""Filtro viral do colapso por hospedeiro do sylph (2026-08-19)."""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                                "scripts", "reads_classify"))
from collapse_by_host import _is_viral, _leaf_rows


def test_imgvr_comeca_no_realm_nao_no_dominio():
    """O predicado antigo era contains("d__Viruses"); a taxonomia real do
    IMG/VR comeca em "r__". 1482 clados virais casavam zero vezes."""
    assert _is_viral("r__Duplodnaviria")
    assert _is_viral("r__Varidnaviria|k__Bamfordvirae|p__Preplasmiviricota")
    assert not _is_viral("d__Bacteria")
    assert not _is_viral("d__Archaea|p__Thermoproteota")


def test_outras_taxonomias_com_dominio_seguem_valendo():
    assert _is_viral("d__Viruses|p__Uroviricota")


def test_folha_em_vez_de_rank_fixo():
    """As linhagens do IMG/VR pulam especie -- a folha e "t__IMGVR_UViG_...",
    entao exigir "s__" descartava tudo."""
    clades = [
        "r__Duplodnaviria",
        "r__Duplodnaviria|k__Heunggongvirae",
        "r__Duplodnaviria|k__Heunggongvirae|t__IMGVR_UViG_1",
        "d__Bacteria",
        "d__Bacteria|p__Pseudomonadota",
    ]
    assert _leaf_rows(clades) == [False, False, True, False, True]


def test_folha_nao_confunde_prefixo_parcial():
    """"r__Ab" nao e pai de "r__Abc" -- o corte e no separador."""
    clades = ["r__Ab", "r__Abc"]
    assert _leaf_rows(clades) == [True, True]


def test_soma_nao_multiplica_por_nivel():
    """Um unico linhagem de 3 niveis contribui uma linha, nao tres."""
    clades = ["r__X", "r__X|k__Y", "r__X|k__Y|t__Z"]
    assert sum(_leaf_rows(clades)) == 1
