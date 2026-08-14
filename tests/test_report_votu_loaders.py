import os
import pytest
from report.data_loaders import (
    load_votu_catalog, load_votu_presence,
    load_tool_status, summarize_tool_status, GLOBAL_STATUS_LABEL,
)


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


def test_load_tool_status_reports_failed_global_rule(tmp_path):
    outdir = tmp_path
    catalog_dir = outdir / "votu_catalog"
    catalog_dir.mkdir()
    (catalog_dir / "done.txt").write_text("failed: empty pool\n")
    # matrices_done.txt intentionally absent -> reported as 'unknown'.

    status = load_tool_status(str(outdir), ["S1", "S2"])

    assert GLOBAL_STATUS_LABEL in status
    assert status[GLOBAL_STATUS_LABEL]["votu_catalog_reps"]["state"] == "failed"
    assert status[GLOBAL_STATUS_LABEL]["votu_catalog_reps"]["reason"] == "empty pool"
    assert status[GLOBAL_STATUS_LABEL]["votu_catalog_matrices"]["state"] == "unknown"

    rows = summarize_tool_status(status)
    global_rows = [r for r in rows if r["sample"] == GLOBAL_STATUS_LABEL]
    assert len(global_rows) == 2
    tools_reported = {r["tool"] for r in global_rows}
    assert tools_reported == {"votu_catalog_reps", "votu_catalog_matrices"}


def test_load_tool_status_global_rule_does_not_corrupt_per_sample_counts(tmp_path):
    outdir = tmp_path
    catalog_dir = outdir / "votu_catalog"
    catalog_dir.mkdir()
    (catalog_dir / "done.txt").write_text("failed: empty pool\n")
    (catalog_dir / "matrices_done.txt").write_text("ok\n")

    status = load_tool_status(str(outdir), ["S1", "S2"])

    # Per-sample entries are untouched by the global rules: only the
    # sample-scoped tools appear under each real sample key.
    assert set(status["S1"].keys()) == {
        "amrfinderplus", "rgi", "galah_derep", "gtdbtk", "vcontact3",
    }
    assert set(status["S2"].keys()) == set(status["S1"].keys())
    assert set(status.keys()) == {"S1", "S2", GLOBAL_STATUS_LABEL}

    rows = summarize_tool_status(status)
    per_sample_rows = [r for r in rows if r["sample"] != GLOBAL_STATUS_LABEL]
    # None of the sample-scoped done.txt files exist on disk -> every
    # sample-scoped tool is 'unknown', exactly 5 per sample, unaffected by
    # the global rule's separate failure.
    assert len(per_sample_rows) == 10
    assert {(r["sample"], r["state"]) for r in per_sample_rows} == {
        ("S1", "unknown"), ("S2", "unknown"),
    }


def test_load_tool_status_rejects_colliding_sample_name(tmp_path):
    with pytest.raises(ValueError):
        load_tool_status(str(tmp_path), ["S1", GLOBAL_STATUS_LABEL])
