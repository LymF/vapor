import os
import pytest
from votu_catalog import build_pool, prefixed_id
from votu_catalog import parse_skani_sparse
from votu_catalog import (
    cluster_votus, pick_representative, assign_votu_ids,
    ClusteringCollapseError, write_clusters,
)


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def test_prefixed_id_joins_with_pipe():
    assert prefixed_id("P01_RNG", "MEGAHIT_k141_10006") == "P01_RNG|MEGAHIT_k141_10006"


def test_build_pool_prefixes_and_records_provenance(tmp_path):
    a = tmp_path / "a.fasta"
    b = tmp_path / "b.fasta"
    # Same contig ID in both samples -- the collision the prefix exists to prevent.
    _write_fasta(a, [("MEGAHIT_k141_10006", "ACGT"), ("MEGAHIT_k141_2", "TTTT")])
    _write_fasta(b, [("MEGAHIT_k141_10006", "GGGG")])

    pool = tmp_path / "pool.fasta"
    prov = tmp_path / "provenance.tsv"
    stats = build_pool(
        [("sample", "S1", str(a)), ("sample", "S2", str(b))],
        str(pool), str(prov),
    )

    assert stats["n_sequences"] == 3
    assert stats["n_sources"] == 2

    names = [l[1:].strip() for l in open(pool) if l.startswith(">")]
    assert names == ["S1|MEGAHIT_k141_10006", "S1|MEGAHIT_k141_2", "S2|MEGAHIT_k141_10006"]
    assert len(set(names)) == 3          # no collision survived

    rows = [l.rstrip("\n").split("\t") for l in open(prov)]
    assert rows[0] == ["member_id", "source_type", "source_id", "original_contig_id"]
    assert rows[1] == ["S1|MEGAHIT_k141_10006", "sample", "S1", "MEGAHIT_k141_10006"]
    assert rows[3] == ["S2|MEGAHIT_k141_10006", "sample", "S2", "MEGAHIT_k141_10006"]


def test_build_pool_keeps_sequence_content(tmp_path):
    a = tmp_path / "a.fasta"
    _write_fasta(a, [("c1", "ACGTACGT")])
    pool = tmp_path / "pool.fasta"
    build_pool([("sample", "S1", str(a))], str(pool), str(tmp_path / "p.tsv"))
    assert "ACGTACGT" in open(pool).read()


def test_build_pool_skips_missing_and_empty_sources(tmp_path):
    a = tmp_path / "a.fasta"
    _write_fasta(a, [("c1", "ACGT")])
    empty = tmp_path / "empty.fasta"
    empty.write_text("")

    stats = build_pool(
        [("sample", "S1", str(a)),
         ("sample", "S2", str(empty)),
         ("group", "G1", str(tmp_path / "missing.fasta"))],
        str(tmp_path / "pool.fasta"), str(tmp_path / "p.tsv"),
    )
    assert stats["n_sequences"] == 1
    assert stats["n_sources"] == 1
    assert stats["n_skipped"] == 2


def test_build_pool_uses_first_whitespace_token_as_id(tmp_path):
    a = tmp_path / "a.fasta"
    with open(a, "w") as fh:
        fh.write(">c1 length=500 cov=3.2\nACGT\n")
    pool = tmp_path / "pool.fasta"
    build_pool([("sample", "S1", str(a))], str(pool), str(tmp_path / "p.tsv"))
    names = [l[1:].strip() for l in open(pool) if l.startswith(">")]
    assert names == ["S1|c1"]


SKANI_HEADER = ("Ref_file\tQuery_file\tANI\tAlign_fraction_ref\t"
                "Align_fraction_query\tRef_name\tQuery_name\n")


def _skani_row(ani, af_ref, af_query, ref_name, query_name):
    return (f"/path/pool.fasta\t/path/pool.fasta\t{ani}\t{af_ref}\t{af_query}\t"
            f"{ref_name}\t{query_name}\n")


def test_parse_uses_columns_five_and_six_for_names(tmp_path):
    """Names live in Ref_name/Query_name, NOT in the first two file-path columns."""
    p = tmp_path / "ani.tsv"
    p.write_text(SKANI_HEADER + _skani_row(99.0, 95.0, 95.0, "S1|a", "S1|b"))
    edges = parse_skani_sparse(str(p), 95.0, 85.0, {"S1|a", "S1|b"})
    assert edges == [("S1|a", "S1|b")]


def test_parse_applies_ani_and_af_thresholds(tmp_path):
    p = tmp_path / "ani.tsv"
    p.write_text(
        SKANI_HEADER
        + _skani_row(99.0, 95.0, 95.0, "a", "b")   # passa
        + _skani_row(94.9, 95.0, 95.0, "a", "c")   # ANI baixo
        + _skani_row(99.0, 80.0, 84.9, "a", "d")   # AF baixo nos dois lados
        + _skani_row(99.0, 84.0, 90.0, "a", "e")   # max(AF) passa
    )
    ids = {"a", "b", "c", "d", "e"}
    edges = parse_skani_sparse(str(p), 95.0, 85.0, ids)
    assert sorted(edges) == [("a", "b"), ("a", "e")]


def test_parse_drops_self_comparisons(tmp_path):
    p = tmp_path / "ani.tsv"
    p.write_text(SKANI_HEADER + _skani_row(100.0, 100.0, 100.0, "a", "a"))
    assert parse_skani_sparse(str(p), 95.0, 85.0, {"a"}) == []


def test_parse_drops_unknown_names(tmp_path):
    p = tmp_path / "ani.tsv"
    p.write_text(SKANI_HEADER + _skani_row(99.0, 95.0, 95.0, "a", "ghost"))
    assert parse_skani_sparse(str(p), 95.0, 85.0, {"a"}) == []


def test_parse_rejects_dense_matrix_format(tmp_path):
    """The dense PHYLIP matrix must yield zero edges, not garbage ones.

    This is the exact format the old parser was silently fed.
    """
    p = tmp_path / "dense.tsv"
    p.write_text("3\nS1|a\nS1|b\t0.00\nS1|c\t0.00\t0.00\n")
    edges = parse_skani_sparse(str(p), 95.0, 85.0, {"S1|a", "S1|b", "S1|c"})
    assert edges == []


def test_parse_tolerates_missing_file(tmp_path):
    assert parse_skani_sparse(str(tmp_path / "nope.tsv"), 95.0, 85.0, {"a"}) == []


def test_parse_skips_malformed_rows(tmp_path):
    p = tmp_path / "ani.tsv"
    p.write_text(SKANI_HEADER
                 + "só\tdois\n"
                 + _skani_row("NA", 95.0, 95.0, "a", "b")
                 + _skani_row(99.0, 95.0, 95.0, "a", "c"))
    assert parse_skani_sparse(str(p), 95.0, 85.0, {"a", "b", "c"}) == [("a", "c")]


def test_cluster_merges_transitively():
    ids = ["a", "b", "c", "d"]
    edges = [("a", "b"), ("b", "c")]
    clusters = cluster_votus(ids, edges, {})
    assert sorted(len(c) for c in clusters) == [1, 3]
    big = [c for c in clusters if len(c) == 3][0]
    assert set(big) == {"a", "b", "c"}


def test_cluster_keeps_singletons():
    clusters = cluster_votus(["a", "b"], [], {})
    assert sorted(clusters) == [["a"], ["b"]]


def test_pick_representative_prefers_highest_completeness():
    assert pick_representative(["a", "b", "c"], {"a": 40.0, "b": 95.0, "c": 70.0}) == "b"


def test_pick_representative_breaks_ties_by_member_order():
    assert pick_representative(["a", "b"], {"a": 50.0, "b": 50.0}) == "a"


def test_pick_representative_handles_missing_completeness():
    assert pick_representative(["a", "b"], {"b": 10.0}) == "b"


def test_clusters_sorted_by_size_then_representative():
    ids = ["z", "y", "m", "n", "o"]
    edges = [("m", "n"), ("n", "o")]          # cluster de 3
    clusters = cluster_votus(ids, edges, {})
    assert len(clusters[0]) == 3               # maior primeiro
    reps = [pick_representative(c, {}) for c in clusters[1:]]
    assert reps == sorted(reps)                # empates em ordem estavel


def test_assign_votu_ids_is_deterministic():
    ids = ["a", "b", "c", "d"]
    edges = [("a", "b")]
    rows1 = assign_votu_ids(cluster_votus(ids, edges, {}))
    rows2 = assign_votu_ids(cluster_votus(ids, edges, {}))
    assert rows1 == rows2
    votu_ids = sorted({r[0] for r in rows1})
    assert votu_ids == ["vOTU_00001", "vOTU_00002", "vOTU_00003"]


def test_assign_votu_ids_emits_one_row_per_member():
    rows = assign_votu_ids([["a", "b"], ["c"]])
    assert len(rows) == 3
    assert all(len(r) == 3 for r in rows)


def test_write_clusters_raises_when_nothing_collapsed(tmp_path):
    """N sequences -> N clusters is the signature of the historical bug."""
    clusters = [["a"], ["b"], ["c"]]
    with pytest.raises(ClusteringCollapseError, match="no clusters formed"):
        write_clusters(clusters, n_input=3, path=str(tmp_path / "c.tsv"))


def test_write_clusters_allows_single_input_sequence(tmp_path):
    """One genome cannot collapse -- that is not the bug."""
    out = tmp_path / "c.tsv"
    write_clusters([["a"]], n_input=1, path=str(out))
    assert out.read_text().startswith("votu_id\trepresentative\tmember\n")


def test_write_clusters_writes_header_and_rows(tmp_path):
    out = tmp_path / "c.tsv"
    write_clusters([["a", "b"], ["c"]], n_input=3, path=str(out),
                   completeness={"a": 10.0, "b": 90.0})
    rows = [l.rstrip("\n").split("\t") for l in open(out)]
    assert rows[0] == ["votu_id", "representative", "member"]
    assert rows[1] == ["vOTU_00001", "b", "a"]
    assert rows[2] == ["vOTU_00001", "b", "b"]
    assert rows[3] == ["vOTU_00002", "c", "c"]
