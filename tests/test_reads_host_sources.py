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


def test_banco_tem_precedencia_sobre_phist():
    """Uma predicao por k-mer nao pode sobrescrever atribuicao publicada."""
    db = {"c1": "d__B;p__;c__;o__;f__;g__Acinetobacter"}
    ph = {"c1": "d__B;p__;c__;o__;f__;g__Pseudomonas"}
    assert resolve_host("c1", db, ph) == ("Acinetobacter", "db")


def test_phist_preenche_onde_o_banco_cala():
    db = {"c1": "d__B;p__;c__;o__;f__;UNKNOWN"}
    ph = {"c1": "d__B;p__;c__;o__;f__;g__Pseudomonas"}
    assert resolve_host("c1", db, ph) == ("Pseudomonas", "phist")


def test_sem_nenhuma_fonte_a_origem_e_none_nao_um_palpite():
    assert resolve_host("c1", {}, {}) == (UNKNOWN, "none")
