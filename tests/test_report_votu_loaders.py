import os
import pytest
from report.data_loaders import load_votu_catalog, load_votu_presence


def _make_catalog(tmp_path, presence_rows, n_pool=10):
    d = tmp_path / "votu_catalog"
    d.mkdir()
    with open(d / "vOTU_clusters.tsv", "w") as fh:
        fh.write("votu_id\trepresentative\tmember\n")
        for i in range(1, len(presence_rows) + 1):
            fh.write(f"vOTU_{i:05d}\trep{i}\trep{i}\n")
    with open(d / "provenance.tsv", "w") as fh:
        fh.write("member_id\tsource_type\tsource_id\toriginal_contig_id\n")
        for i in range(n_pool):
            fh.write(f"S1|c{i}\tsample\tS1\tc{i}\n")
    with open(d / "presence_matrix.tsv", "w") as fh:
        fh.write("votu_id\trepresentative\tS1\tS2\n")
        for i, (a, b) in enumerate(presence_rows, start=1):
            fh.write(f"vOTU_{i:05d}\trep{i}\t{a}\t{b}\n")
    return str(tmp_path)


def test_load_votu_catalog_reports_global_richness(tmp_path):
    outdir = _make_catalog(tmp_path, [("both", "absent"), ("absent", "recruited")],
                           n_pool=10)
    cat = load_votu_catalog(outdir)
    assert cat["n_votus"] == 2
    assert cat["n_pool"] == 10
    assert cat["reduction_pct"] == pytest.approx(80.0)


def test_load_votu_presence_counts_both_signals(tmp_path):
    outdir = _make_catalog(tmp_path, [
        ("both", "absent"),
        ("assembled", "recruited"),
        ("absent", "absent"),
    ])
    pres = load_votu_presence(outdir, ["S1", "S2"])
    assert pres["per_sample"]["S1"]["assembled"] == 2   # 'both' + 'assembled'
    assert pres["per_sample"]["S1"]["recruited"] == 1   # 'both'
    assert pres["per_sample"]["S1"]["total"] == 2       # presente por qualquer sinal
    assert pres["per_sample"]["S2"]["recruited"] == 1
    assert pres["per_sample"]["S2"]["total"] == 1


def test_load_votu_catalog_missing_returns_empty(tmp_path):
    cat = load_votu_catalog(str(tmp_path))
    assert cat["n_votus"] == 0
    assert cat["reduction_pct"] == 0.0


def test_load_votu_presence_missing_returns_zeros(tmp_path):
    pres = load_votu_presence(str(tmp_path), ["S1"])
    assert pres["per_sample"]["S1"]["total"] == 0
    assert pres["votus"] == []
