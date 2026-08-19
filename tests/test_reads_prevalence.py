"""Filtro de prevalencia da trilha de reads (sylph)."""
import subprocess
import sys
import os

SCRIPT = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "scripts", "reads_classify", "filter_by_prevalence.py")


def _run(tmp_path, rows, min_prev, header="clade_name\tS1\tS2\tS3\tS4\n"):
    inp = tmp_path / "in.tsv"
    out = tmp_path / "out.tsv"
    inp.write_text(header + "".join(rows))
    subprocess.run([sys.executable, SCRIPT, str(inp), str(out), str(min_prev)],
                   check=True, capture_output=True)
    return [l.split("\t")[0] for l in out.read_text().splitlines()[1:]]


def test_taxon_exatamente_no_corte_e_mantido(tmp_path):
    """Com '>' puro, um corte de 0.25 descartava o taxon presente em
    exatamente 1 de 4 amostras -- justamente o que 'prevalencia minima de
    25%' quer manter."""
    kept = _run(tmp_path, ["r__A\t1.0\t0\t0\t0\n", "r__B\t1.0\t1.0\t0\t0\n"], 0.25)
    assert kept == ["r__A", "r__B"]


def test_taxon_abaixo_do_corte_sai(tmp_path):
    kept = _run(tmp_path, ["r__A\t1.0\t0\t0\t0\n", "r__B\t1.0\t1.0\t0\t0\n"], 0.5)
    assert kept == ["r__B"]


def test_default_zero_ainda_descarta_linhas_zeradas(tmp_path):
    """O sylph-tax merge emite linhas zeradas em toda amostra (UNKNOWN, por
    exemplo). Com min_prev 0.0 elas nao podem voltar."""
    kept = _run(tmp_path, ["UNKNOWN\t0\t0\t0\t0\n", "r__A\t0\t0.4\t0\t0\n"], 0.0)
    assert kept == ["r__A"]
