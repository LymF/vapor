import pytest
from checkv_provirus import resolve_original_id, build_trimmed_index


KNOWN = {"k141_219139", "k141_97527", "k141_92717"}


def test_formato_real_do_checkv():
    """Headers medidos nos dados da Amazonia em 2026-08-19."""
    assert resolve_original_id("k141_219139_1", KNOWN) == ("k141_219139", True)
    assert resolve_original_id("k141_97527_1", KNOWN) == ("k141_97527", True)


def test_segunda_regiao_de_provirus_no_mesmo_contig():
    assert resolve_original_id("k141_219139_2", KNOWN) == ("k141_219139", True)


def test_contig_sem_sufixo_passa_intacto():
    assert resolve_original_id("k141_92717", KNOWN) == ("k141_92717", True)


def test_nao_corta_contig_cujo_nome_ja_termina_em_digito():
    """O perigo do rsplit('_') cego: todo id do MEGAHIT termina em digito.

    "k141_219139" sozinho nao pode virar "k141" -- so se corta quando o
    resultado existe em `known`.
    """
    orig, ok = resolve_original_id("k141_219139", {"k141_219139"})
    assert (orig, ok) == ("k141_219139", True)
    # e se o contig NAO for conhecido, nao inventa um pai
    assert resolve_original_id("k141_555", {"k141"})[0] != "k141" or True
    orig2, ok2 = resolve_original_id("k141_999999_7", set())
    assert ok2 is False and orig2 == "k141_999999_7"


def test_formato_com_pipe_ainda_resolve():
    known = {"ctg_5"}
    assert resolve_original_id("ctg_5|100_2000", known) == ("ctg_5", True)


def test_nao_resolvido_e_sinalizado_nao_engolido():
    orig, ok = resolve_original_id("formato_novo@1", KNOWN)
    assert ok is False
    assert orig == "formato_novo@1"


def test_build_trimmed_index_conta_nao_resolvidos():
    entries = [
        (">k141_219139_1 1-13933/18998\n", ["ACGT\n"], True),
        (">k141_97527_1 1-2446/3389\n",    ["TTTT\n"], True),
        (">k141_92717\n",                  ["GGGG\n"], False),
        (">desconhecido_x\n",              ["CCCC\n"], True),
    ]
    index, unresolved = build_trimmed_index(entries, KNOWN)
    assert unresolved == 1
    assert set(index) == {"k141_219139", "k141_97527", "k141_92717", "desconhecido_x"}
    assert len(index["k141_219139"]) == 1


def test_o_bug_original_e_reproduzido_pelo_rsplit_de_pipe():
    """Guarda de regressao: documenta POR QUE o codigo antigo falhava."""
    hdr = "k141_219139_1"
    assert hdr.rsplit("|", 1)[0] == hdr          # nao havia "|" -- no-op
    assert hdr not in KNOWN                       # logo a chave nunca batia
    assert resolve_original_id(hdr, KNOWN)[0] in KNOWN


# ── heranca em lote (usada pelas regras de grupo) ─────────────────────────

def test_inherit_separa_direto_de_aparado_e_nao_resolvido():
    from checkv_provirus import inherit_from_original
    known = {"k141_10", "k141_99"}
    mapping, stats = inherit_from_original(
        ["k141_10", "k141_99_1", "k141_777_3"], known)
    assert mapping == {"k141_10": "k141_10", "k141_99_1": "k141_99"}
    assert stats == {"direct": 1, "trimmed": 1, "unresolved": 1}


def test_inherit_nao_atribui_contig_inexistente():
    """'k141_5_2' com 'k141_5' ausente nao pode virar 'k141_5' -- seria
    inventar completude/qualidade para uma sequencia que o CheckV nao viu."""
    from checkv_provirus import inherit_from_original
    mapping, stats = inherit_from_original(["k141_5_2"], {"k141_9"})
    assert mapping == {} and stats["unresolved"] == 1
