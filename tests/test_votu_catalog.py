import os
import pytest
from votu_catalog import build_pool, prefixed_id


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
