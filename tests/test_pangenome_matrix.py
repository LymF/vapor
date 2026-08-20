# tests/test_pangenome_matrix.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from pangenome_matrix import build_matrix, summarize_clusters

MEMBERS = {"R1": ["m1", "m2", "m3", "m4"]}
COMPLETE = {"m1": 95.0, "m2": 88.0, "m3": 91.0, "m4": 42.0}   # m4 e ruim
HITS = {"m1": {"RM_Type_II", "Gabija"}, "m2": {"RM_Type_II"},
        "m3": {"RM_Type_II", "Gabija"}, "m4": set()}


def _rows():
    return {r["gene"]: r for r in build_matrix(["R1"], MEMBERS, HITS, COMPLETE)}


class TestBuildMatrix:
    def test_incomplete_member_is_unassessable_not_absent(self):
        # O ponto cientifico central: m4 esta 42% completo. Um "." ali
        # afirmaria que o organismo nao tem o gene, quando o que houve foi
        # a regiao nao ter montado.
        assert _rows()["Gabija"]["states"]["m4"] == "?"
        assert _rows()["RM_Type_II"]["states"]["m4"] == "?"

    def test_unassessable_member_is_out_of_the_denominator(self):
        # Gabija em 2 de 3 AVALIAVEIS, nao 2 de 4.
        row = _rows()["Gabija"]
        assert row["n_present"] == 2
        assert row["n_evaluable"] == 3
        assert row["freq"] == "2/3"

    def test_present_and_absent_are_distinguished(self):
        row = _rows()["Gabija"]
        assert row["states"]["m1"] == "x"
        assert row["states"]["m2"] == "."

    def test_gene_in_every_evaluable_member(self):
        row = _rows()["RM_Type_II"]
        assert row["n_present"] == 3
        assert row["freq"] == "3/3"


class TestSummarizeClusters:
    def test_core_uses_90_percent_not_99(self):
        # RM_Type_II esta em 3/3 avaliaveis -> core. Gabija em 2/3 (67%)
        # -> variavel. O limiar de 99% zeraria o core com MAG.
        rows = build_matrix(["R1"], MEMBERS, HITS, COMPLETE)
        summary = {s["representative_id"]: s
                   for s in summarize_clusters(rows, MEMBERS, COMPLETE)}["R1"]
        assert summary["n_genes_core"] == 1
        assert summary["n_genes_variaveis"] == 1

    def test_evaluable_member_count_excludes_the_incomplete_one(self):
        rows = build_matrix(["R1"], MEMBERS, HITS, COMPLETE)
        summary = summarize_clusters(rows, MEMBERS, COMPLETE)[0]
        assert summary["n_members"] == 4
        assert summary["n_members_avaliaveis"] == 3

    def test_singleton_gene_is_counted(self):
        hits = dict(HITS, m1=HITS["m1"] | {"CBASS"})
        rows = build_matrix(["R1"], MEMBERS, hits, COMPLETE)
        summary = summarize_clusters(rows, MEMBERS, COMPLETE)[0]
        assert summary["n_genes_singleton"] == 1

    def test_cluster_with_no_evaluable_member_does_not_divide_by_zero(self):
        rows = build_matrix(["R1"], MEMBERS, HITS,
                            {m: 10.0 for m in MEMBERS["R1"]})
        summary = summarize_clusters(rows, MEMBERS,
                                     {m: 10.0 for m in MEMBERS["R1"]})[0]
        assert summary["n_members_avaliaveis"] == 0
        assert summary["n_genes_core"] == 0
