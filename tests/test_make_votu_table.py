"""
Unit tests for scripts/make_votu_table.py's per-sample filtering of the
global vOTU catalog.

vOTU_clusters.tsv (produced by rules/votu_catalog.smk) namespaces every
member ID as "{source_id}|{contig_id}" because clustering happens once over
the pooled viral sets of every sample and co-assembly group. This module's
`load_votu_clusters(clusters_tsv, sample)` must:
  - keep votu_id / representative namespaced (global identity — the
    representative of a vOTU may belong to a different sample)
  - drop members that don't belong to the requested sample
  - strip the prefix off of surviving members only
  - drop a vOTU entirely when it has no member in the requested sample
"""

import csv

from make_votu_table import load_votu_clusters


def _write_clusters_tsv(path, rows):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["votu_id", "representative", "member"])
        for votu_id, rep, mem in rows:
            writer.writerow([votu_id, rep, mem])


def test_member_from_other_sample_is_excluded(tmp_path):
    clusters_tsv = tmp_path / "vOTU_clusters.tsv"
    # vOTU_1 has members from S1 and S2; only S1's member should survive
    # when filtering for sample=S1.
    _write_clusters_tsv(clusters_tsv, [
        ("vOTU_1", "S1|contigA", "S1|contigA"),
        ("vOTU_1", "S1|contigA", "S2|contigB"),
    ])

    votu_order, votu_rep, votu_members, votu_total_size = load_votu_clusters(
        str(clusters_tsv), "S1")

    assert votu_order == ["vOTU_1"]
    assert votu_members["vOTU_1"] == ["contigA"]
    # The S2 member must not leak into S1's member list under any key.
    all_members = [m for members in votu_members.values() for m in members]
    assert "S2|contigB" not in all_members
    assert "contigB" not in all_members


def test_surviving_member_id_has_no_prefix(tmp_path):
    clusters_tsv = tmp_path / "vOTU_clusters.tsv"
    _write_clusters_tsv(clusters_tsv, [
        ("vOTU_1", "S1|contigA", "S1|contigA"),
    ])

    _, _, votu_members, _ = load_votu_clusters(str(clusters_tsv), "S1")

    assert votu_members["vOTU_1"] == ["contigA"]
    for members in votu_members.values():
        for m in members:
            assert "|" not in m


def test_representative_keeps_prefix_even_from_another_sample(tmp_path):
    clusters_tsv = tmp_path / "vOTU_clusters.tsv"
    # The vOTU's representative is S2's contig, but S1 still has a member
    # in the same cluster -- this is the scenario the naive prefix-strip
    # would have destroyed the link for.
    _write_clusters_tsv(clusters_tsv, [
        ("vOTU_1", "S2|contigX", "S2|contigX"),
        ("vOTU_1", "S2|contigX", "S1|contigA"),
    ])

    votu_order, votu_rep, votu_members, votu_total_size = load_votu_clusters(
        str(clusters_tsv), "S1")

    assert votu_order == ["vOTU_1"]
    # Representative stays fully namespaced -- NOT stripped, NOT replaced
    # with a local ID -- even though it belongs to a different sample.
    assert votu_rep["vOTU_1"] == "S2|contigX"
    # Only S1's (bare) member is present as a row-producing member.
    assert votu_members["vOTU_1"] == ["contigA"]
    # Global cluster size still reflects both samples' members.
    assert votu_total_size["vOTU_1"] == 2


def test_votu_with_no_member_in_sample_produces_no_row(tmp_path):
    clusters_tsv = tmp_path / "vOTU_clusters.tsv"
    _write_clusters_tsv(clusters_tsv, [
        ("vOTU_1", "S2|contigX", "S2|contigX"),
        ("vOTU_2", "S1|contigA", "S1|contigA"),
    ])

    votu_order, votu_rep, votu_members, votu_total_size = load_votu_clusters(
        str(clusters_tsv), "S1")

    assert votu_order == ["vOTU_2"]
    assert "vOTU_1" not in votu_rep
    assert "vOTU_1" not in votu_members


# ── PHIST: o prefixo do nome do arquivo e a heranca do representante ──────

def _phist_csv(tmp_path, rows):
    p = tmp_path / "phist_results.csv"
    with open(p, "w", newline="") as fh:
        fh.write("phage,host,#common-kmers,pvalue,adj-pvalue\n")
        for r in rows:
            fh.write(",".join(r) + "\n")
    return str(p)


def test_phist_ignora_o_prefixo_contig_do_nome_do_arquivo(tmp_path):
    """split_viral_fastas.py grava 'contig_{id}.fasta'. Sem tirar o prefixo,
    a chave vira 'contig_S1|k141_10' e nada casa -- host_bin sai vazio em
    toda amostra, sem erro."""
    from make_votu_table import load_phist
    p = _phist_csv(tmp_path, [
        ["contig_S1|k141_10.fasta", "binette_bin2.fa", "12", "1e-12", "1e-10"],
    ])
    hosts = load_phist(p)
    assert list(hosts) == ["S1|k141_10"]
    assert hosts["S1|k141_10"]["host_bin"] == "binette_bin2"


def test_phist_mantem_o_melhor_hit_por_virus(tmp_path):
    from make_votu_table import load_phist
    p = _phist_csv(tmp_path, [
        ["contig_S1|k141_10.fasta", "binette_bin2.fa", "12", "1e-3", "1e-2"],
        ["contig_S1|k141_10.fasta", "binette_bin9.fa", "40", "1e-30", "1e-28"],
    ])
    assert load_phist(p)["S1|k141_10"]["host_bin"] == "binette_bin9"


def test_phist_nao_perde_bins_do_vrhyme(tmp_path):
    """Os bins do vRhyme entram no PHIST com o proprio nome de arquivo, sem
    o prefixo 'contig_' -- nao podem ser mutilados pela remocao."""
    from make_votu_table import load_phist
    p = _phist_csv(tmp_path, [
        ["vRhyme_bin_3.fasta", "binette_bin2.fa", "12", "1e-12", "1e-10"],
    ])
    assert list(load_phist(p)) == ["vRhyme_bin_3"]
