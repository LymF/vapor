# tests/test_pangenome_matrix.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from pangenome_matrix import (CORE_FRACTION, build_matrix, summarize_clusters,
                              load_defensefinder_summary)

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


class TestAnnotatedGatesAllTypes:
    # O caso "membro fora do manifesto" (prodigal falhou ou o genoma nunca
    # chegou ao pool): '?' em TODAS as linhas do membro, defesa e amr.
    MEMBERS3 = {"R1": ["m1", "m2", "m3"]}
    COMPLETE3 = {"m1": 95.0, "m2": 90.0, "m3": 92.0}
    HITS3 = {
        "m1": {("defesa", "Gabija"), ("amr", "blaTEM")},
        "m2": {("defesa", "Gabija")},
        "m3": {("amr", "blaTEM")},
    }

    def test_member_outside_annotated_is_unassessable_in_every_row(self):
        # m3 nunca aparece no manifesto do pangenoma (prodigal falhou ou o
        # genoma ficou fora do pool) -- annotated exclui m3.
        annotated = {"m1", "m2"}
        rows = build_matrix(["R1"], self.MEMBERS3, self.HITS3, self.COMPLETE3,
                            annotated=annotated)
        by_gene_tipo = {(r["tipo"], r["gene"]): r for r in rows}

        assert by_gene_tipo[("defesa", "Gabija")]["states"]["m3"] == "?"
        assert by_gene_tipo[("amr", "blaTEM")]["states"]["m3"] == "?"

    def test_member_outside_annotated_leaves_both_denominators(self):
        annotated = {"m1", "m2"}
        rows = build_matrix(["R1"], self.MEMBERS3, self.HITS3, self.COMPLETE3,
                            annotated=annotated)
        by_gene_tipo = {(r["tipo"], r["gene"]): r for r in rows}

        # annotated={m1,m2} e o mesmo denominador base para AMBOS os tipos
        # (m3 esta fora do manifesto, entao sai de toda linha, defesa e
        # amr) -- e exatamente o ponto do teste: nao ha diferenca de
        # denominador entre tipos aqui, so quando defense_failed entra em
        # jogo (ver TestDefenseFailedOnlyGatesDefenseRows).
        assert by_gene_tipo[("defesa", "Gabija")]["n_evaluable"] == 2
        assert by_gene_tipo[("amr", "blaTEM")]["n_evaluable"] == 2

    def test_this_case_is_unaffected_by_defense_failed(self):
        # annotated e defense_failed sao independentes: um membro fora do
        # manifesto continua '?' em amr mesmo que defense_failed esteja
        # vazio -- a causa aqui e prodigal, nao DefenseFinder.
        annotated = {"m1", "m2"}
        rows = build_matrix(["R1"], self.MEMBERS3, self.HITS3, self.COMPLETE3,
                            annotated=annotated, defense_failed=set())
        row = {(r["tipo"], r["gene"]): r for r in rows}[("amr", "blaTEM")]
        assert row["states"]["m3"] == "?"


class TestDefenseFailedOnlyGatesDefenseRows:
    # m2 tem prodigal OK (esta em `annotated`) mas o DefenseFinder falhou
    # naquele genoma especifico -- so as linhas tipo='defesa' de m2 podem
    # virar '?'; AMR vem de outra ferramenta e continua valido.
    MEMBERS = {"R1": ["m1", "m2", "m3"]}
    COMPLETE = {"m1": 95.0, "m2": 90.0, "m3": 92.0}
    ANNOTATED = {"m1", "m2", "m3"}
    HITS = {
        "m1": {("defesa", "Gabija"), ("amr", "blaTEM")},
        "m2": {("defesa", "Gabija"), ("amr", "blaTEM")},
        "m3": {("amr", "blaTEM")},
    }
    DEFENSE_FAILED = {"m2"}

    def _rows(self):
        return {(r["tipo"], r["gene"]): r
                for r in build_matrix(["R1"], self.MEMBERS, self.HITS,
                                      self.COMPLETE, annotated=self.ANNOTATED,
                                      defense_failed=self.DEFENSE_FAILED)}

    def test_defense_row_marks_failed_member_unassessable(self):
        row = self._rows()[("defesa", "Gabija")]
        assert row["states"]["m2"] == "?"

    def test_amr_row_keeps_failed_member_real_state(self):
        # m2 tem o ARG de verdade: precisa continuar 'x', nao virar '?'
        # so porque o DefenseFinder quebrou nele.
        row = self._rows()[("amr", "blaTEM")]
        assert row["states"]["m2"] == "x"
        assert row["states"]["m3"] == "x"

    def test_defense_denominator_excludes_failed_member(self):
        row = self._rows()[("defesa", "Gabija")]
        assert row["n_evaluable"] == 2   # m1, m3 (m3 nao tem hit mas e avaliavel)
        assert row["n_present"] == 1     # so m1
        assert row["freq"] == "1/2"

    def test_amr_denominator_keeps_failed_member(self):
        row = self._rows()[("amr", "blaTEM")]
        assert row["n_evaluable"] == 3
        assert row["n_present"] == 3
        assert row["freq"] == "3/3"

    def test_summarize_agrees_with_matrix_on_defense_row(self):
        # A propriedade central do item A.2: build_matrix e
        # summarize_clusters NUNCA podem discordar sobre quantos membros
        # contam para uma linha de defesa com falha do DefenseFinder.
        rows = build_matrix(["R1"], self.MEMBERS, self.HITS, self.COMPLETE,
                            annotated=self.ANNOTATED,
                            defense_failed=self.DEFENSE_FAILED)
        summary = summarize_clusters(rows, self.MEMBERS, self.COMPLETE,
                                     annotated=self.ANNOTATED)[0]

        gabija = {(r["tipo"], r["gene"]): r for r in rows}[("defesa", "Gabija")]
        # Gabija em 1/2 avaliaveis (50%) fica abaixo de CORE_FRACTION -> variavel.
        assert gabija["n_present"] / gabija["n_evaluable"] < CORE_FRACTION
        assert summary["n_genes_variaveis"] >= 1

    def test_summarize_amr_denominator_unaffected_by_defense_failure(self):
        rows = build_matrix(["R1"], self.MEMBERS, self.HITS, self.COMPLETE,
                            annotated=self.ANNOTATED,
                            defense_failed=self.DEFENSE_FAILED)
        summary = summarize_clusters(rows, self.MEMBERS, self.COMPLETE,
                                     annotated=self.ANNOTATED)[0]
        # n_members_avaliaveis e o denominador "largo" (prodigal+completude),
        # o mesmo que as linhas amr usam -- nao encolhe por causa do
        # DefenseFinder.
        assert summary["n_members_avaliaveis"] == 3

    def test_defense_failed_member_still_present_in_manifest_gate(self):
        # m2 esta em `annotated` (prodigal OK): a falha e so do
        # DefenseFinder, entao m2 NAO pode virar '?' nas linhas amr por
        # causa de defense_failed sozinho.
        row = self._rows()[("amr", "blaTEM")]
        assert "m2" in row["states"]
        assert row["states"]["m2"] != "?"


class TestLoadDefenseFinderSummary:
    def test_parses_ok_and_failed(self, tmp_path):
        p = tmp_path / "defensefinder_summary.tsv"
        p.write_text("genome\tstatus\ng1\tok\ng2\tfailed\ng3\tok\n")
        assert load_defensefinder_summary(str(p)) == {"g2"}

    def test_unknown_status_counts_as_failed(self, tmp_path):
        # Um status que este parser nunca viu (nem 'ok' nem 'failed') nao
        # pode ser lido como sucesso silencioso.
        p = tmp_path / "defensefinder_summary.tsv"
        p.write_text("genome\tstatus\ng1\tok\ng2\ttimeout\n")
        assert load_defensefinder_summary(str(p)) == {"g2"}

    def test_header_only_file_returns_empty_set(self, tmp_path):
        # Caminho de skip da regra (desabilitada ou manifesto vazio):
        # arquivo so com cabecalho, nenhum genoma marcado como falho.
        p = tmp_path / "defensefinder_summary.tsv"
        p.write_text("genome\tstatus\n")
        assert load_defensefinder_summary(str(p)) == set()

    def test_all_ok_returns_empty_set(self, tmp_path):
        p = tmp_path / "defensefinder_summary.tsv"
        p.write_text("genome\tstatus\ng1\tok\ng2\tok\n")
        assert load_defensefinder_summary(str(p)) == set()
