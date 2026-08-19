"""As duas fontes de hospedeiro do track de reads (2026-08-19)."""
import sys, os
_RC = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                   "scripts", "reads_classify")
sys.path.insert(0, _RC)
from collapse_by_host import _parse_genus, resolve_host, UNKNOWN
from build_host_map import _is_null_host, _host_column, collect


# ── host do banco: extracao do .sylphmpa ─────────────────────────────────

def test_coluna_de_host_casa_apesar_do_sufixo():
    """O sylph-tax escreve "Virus_host (if viral)", com espaco e parenteses."""
    hdr = ["clade_name", "relative_abundance", "ANI (if strain-level)",
           "Virus_host (if viral)"]
    assert _host_column(hdr) == 3
    assert _host_column(["clade_name", "relative_abundance"]) is None


def test_unknown_encadeado_nao_conta_como_host():
    """O IMG/VR emite "UNKNOWN;UNKNOWN;..." -- string longa e nao vazia, que
    passaria por qualquer teste ingenuo de "tem valor"."""
    assert _is_null_host("UNKNOWN;UNKNOWN;UNKNOWN;UNKNOWN;UNKNOWN;UNKNOWN;UNKNOWN")
    assert _is_null_host("NA")
    assert _is_null_host("")
    assert not _is_null_host("d__Bacteria;p__Proteobacteria;UNKNOWN;UNKNOWN")


def test_collect_ignora_cabecalho_de_comentario(tmp_path):
    f = tmp_path / "a.sylphmpa"
    f.write_text(
        "#SampleID\tx\tTaxonomies_used:['IMGVR_4.1']\n"
        "clade_name\trelative_abundance\tVirus_host (if viral)\n"
        "r__Duplodnaviria|t__IMGVR_1\t9.9\td__Bacteria;p__Proteobacteria;"
        "c__Gamma;o__Pseudo;f__Moraxellaceae;g__Acinetobacter\n"
        "r__Duplodnaviria|t__IMGVR_2\t1.0\tUNKNOWN;UNKNOWN;UNKNOWN;UNKNOWN;"
        "UNKNOWN;UNKNOWN\n"
    )
    hosts, n_rows, n_with = collect([str(f)])
    assert n_rows == 2 and n_with == 1
    assert set(hosts) == {"r__Duplodnaviria|t__IMGVR_1"}


# ── resolucao entre as duas fontes ───────────────────────────────────────

def test_genero_unknown_tem_uma_grafia_so():
    """O rank existe mas vem "UNKNOWN"; antes vazava cru e a mesma coluna
    ficava com "Unknown" e "UNKNOWN"."""
    assert _parse_genus("d__Bacteria;p__X;c__Y;o__Z;f__W;UNKNOWN") == UNKNOWN
    assert _parse_genus("") == UNKNOWN
    assert _parse_genus("d__Bacteria;p__X;c__Y;o__Z;f__W;g__Acinetobacter") == "Acinetobacter"


def test_host_do_banco_vem_com_fonte():
    db = {"c1": "d__B;p__;c__;o__;f__;g__Acinetobacter"}
    assert resolve_host("c1", db) == ("Acinetobacter", "db")


def test_rank_de_genero_vazio_nao_vira_atribuicao():
    """O IMG/VR preenche o rank com "UNKNOWN"; isso e ausencia, nao um genero
    chamado UNKNOWN."""
    assert resolve_host("c1", {"c1": "d__B;p__X;c__;o__;f__;UNKNOWN"}) == (UNKNOWN, "none")


def test_sem_fonte_a_origem_e_none_nao_um_palpite():
    """host_source existe para que, depois do groupby, "ninguem atribuiu" nao
    se confunda com uma atribuicao real."""
    assert resolve_host("c1", {}) == (UNKNOWN, "none")
