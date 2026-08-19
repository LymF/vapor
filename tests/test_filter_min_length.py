"""Filtro de comprimento minimo — a trilha ONT nao tinha nenhum."""
import io
from filter_min_length import filter_fasta, iter_fasta


def _run(text, min_len):
    out = io.StringIO()
    stats = filter_fasta(io.StringIO(text), out, min_len)
    return out.getvalue(), stats


def test_descarta_abaixo_do_minimo_e_mantem_o_resto():
    fa = ">c1\nAAAAAAAAAA\n>c2\nAAAAA\n"
    out, (n_in, n_out, dropped) = _run(fa, 10)
    assert out == ">c1\nAAAAAAAAAA\n"
    assert (n_in, n_out, dropped) == (2, 1, 5)


def test_soma_linhas_quebradas_antes_de_comparar():
    """FASTA com quebra de linha e o caso normal de montador; medir so a
    primeira linha descartaria contigs longos."""
    fa = ">c1\nAAAAA\nAAAAA\nAAAAA\n"
    out, (_, n_out, _) = _run(fa, 12)
    assert n_out == 1 and out.count("\n") == 4


def test_preserva_a_descricao_do_header():
    fa = ">c1 len=99 cov=3.5\nAAAAAAAAAA\n"
    out, _ = _run(fa, 5)
    assert out.startswith(">c1 len=99 cov=3.5\n")


def test_limite_e_inclusivo():
    fa = ">c1\nAAAAA\n"
    _, (_, n_out, _) = _run(fa, 5)
    assert n_out == 1


def test_arquivo_vazio_nao_quebra():
    out, stats = _run("", 100)
    assert out == "" and stats == (0, 0, 0)


def test_iter_fasta_ignora_lixo_antes_do_primeiro_header():
    assert list(iter_fasta(io.StringIO("sujeira\n>c1\nAA\n"))) == [("c1", ["AA"])]
