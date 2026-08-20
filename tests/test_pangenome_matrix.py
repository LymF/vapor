# tests/test_pangenome_matrix.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from pangenome_matrix import CORE_FRACTION, build_matrix, summarize_clusters

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
        # Com apenas 3-6 avaliaveis (regime real da fase 1), 0.90*n_eval
        # arredonda para n_eval e o limiar de 99% teria o MESMO efeito --
        # o teste original (3 avaliaveis) nao discriminava nada. Com 10
        # avaliaveis os dois limiares DIVERGEM de verdade: 9/10 e core sob
        # 0.90 (9 >= 9.0) e seria variavel sob 0.99 (9 < 9.9).
        members = [f"m{i}" for i in range(1, 11)]
        members_by_rep = {"R1": members}
        completeness = {m: 95.0 for m in members}
        hits = {m: {"RM_Type_II"} for m in members[:9]}
        hits[members[9]] = set()

        rows = build_matrix(["R1"], members_by_rep, hits, completeness)
        summary = summarize_clusters(rows, members_by_rep, completeness)[0]

        # sob o limiar real (0.90) o gene em 9/10 e core.
        assert summary["n_genes_core"] == 1
        assert summary["n_genes_variaveis"] == 0
        # prova de que 0.90 e 0.99 realmente divergem neste tamanho --
        # sob 0.99 o mesmo 9/10 cairia em variavel, nao core.
        assert 9 >= CORE_FRACTION * 10
        assert 9 < 0.99 * 10

    def test_core_at_small_cluster_size_is_operationally_all_evaluable(self):
        # Caso pequeno (regime real): RM_Type_II em 3/3 avaliaveis -> core.
        # Gabija em 2/3 (67%) -> variavel. Documenta o efeito descrito no
        # docstring do modulo: com poucos avaliaveis "core" e "presente em
        # todos", sem tolerancia de fato.
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

    def test_completeness_precedes_gene_presence(self):
        # A propriedade cientifica central da fase 1: hoje a garantia vem
        # so da ordem do if/elif em build_matrix, sem teste proprio. m4
        # esta abaixo do piso de 70% (42.0) mas TEM um hit real de Gabija
        # em gene_hits -- se a checagem de completude nao viesse primeiro,
        # isso viraria "x" e um MAG incompleto que por acaso montou o gene
        # seria lido como evidencia de presenca, sem o denominador dizer
        # que ele era avaliavel.
        hits = dict(HITS, m4={"Gabija"})
        row = _rows_with(hits)["Gabija"]
        assert row["states"]["m4"] == "?"
        assert row["states"]["m4"] != "x"
        # e nao entra no numerador nem no denominador da frequencia
        assert row["n_present"] == 2
        assert row["n_evaluable"] == 3


def _rows_with(hits):
    return {r["gene"]: r for r in build_matrix(["R1"], MEMBERS, hits, COMPLETE)}


class TestMultipleClusters:
    MEMBERS2 = {
        "R1": ["a1", "a2", "a3"],
        "R2": ["b1", "b2", "b3", "b4"],
    }
    COMPLETE2 = {"a1": 95.0, "a2": 90.0, "a3": 92.0,
                "b1": 88.0, "b2": 85.0, "b3": 91.0, "b4": 20.0}
    HITS2 = {
        "a1": {"RM_Type_II"}, "a2": {"RM_Type_II"}, "a3": set(),
        "b1": {"Gabija", "CBASS"}, "b2": {"Gabija"},
        "b3": {"Gabija", "CBASS"}, "b4": {"Gabija"},
    }

    def test_matrix_and_summary_keep_clusters_separate(self):
        # Um teste so, cobrindo matriz e sumario: dois clusters com
        # representative_id distintos, cada um com seus proprios membros,
        # nao podem se misturar em nenhuma das duas saidas.
        rows = build_matrix(["R1", "R2"], self.MEMBERS2, self.HITS2,
                            self.COMPLETE2)
        r1_rows = [r for r in rows if r["representative_id"] == "R1"]
        r2_rows = [r for r in rows if r["representative_id"] == "R2"]

        # genes de R2 (Gabija/CBASS) nao aparecem nas linhas de R1
        assert {r["gene"] for r in r1_rows} == {"RM_Type_II"}
        assert {r["gene"] for r in r2_rows} == {"Gabija", "CBASS"}

        # os estados de R1 so citam membros de R1
        r1_gene = [r for r in r1_rows if r["gene"] == "RM_Type_II"][0]
        assert set(r1_gene["states"]) == {"a1", "a2", "a3"}

        summary = {s["representative_id"]: s
                   for s in summarize_clusters(rows, self.MEMBERS2,
                                               self.COMPLETE2)}

        # R1: 3 membros, todos avaliaveis
        assert summary["R1"]["n_members"] == 3
        assert summary["R1"]["n_members_avaliaveis"] == 3

        # R2: 4 membros, b4 abaixo do piso -> 3 avaliaveis. O
        # denominador de R2 nao pode vazar para R1 nem vice-versa.
        assert summary["R2"]["n_members"] == 4
        assert summary["R2"]["n_members_avaliaveis"] == 3
