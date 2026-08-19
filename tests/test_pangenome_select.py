import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from pangenome_select import select_clusters

MEMB = {
    "S1__binette_bin1": [{"member_id": f"S{i}__binette_bin1"} for i in range(1, 5)],
    "S9__binette_bin2": [{"member_id": "S9__binette_bin2"}, {"member_id": "S8__binette_bin7"}],
}
NO_EVIDENCE = {"n_islands": 0, "n_systems": 0, "n_args": 0, "n_plasmid": 0}


def _sel(evidence, membership=MEMB, **kw):
    return {c["representative_id"]: c
            for c in select_clusters(membership, evidence, **kw)}


class TestSelectClusters:
    def test_island_with_enough_members_is_eligible(self):
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_islands=1)}
        out = _sel(ev)
        assert out["S1__binette_bin1"]["eligible"] is True
        assert out["S1__binette_bin1"]["criterio"] == "ilha"

    def test_two_members_is_never_eligible_however_strong_the_evidence(self):
        # O piso de 3 e absoluto: com 2 nao ha frequencia que signifique nada.
        ev = {"S9__binette_bin2": dict(NO_EVIDENCE, n_islands=5, n_args=9)}
        assert _sel(ev)["S9__binette_bin2"]["eligible"] is False

    def test_three_defense_systems_without_island_is_eligible(self):
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_systems=3)}
        out = _sel(ev)["S1__binette_bin1"]
        assert out["eligible"] is True
        assert out["criterio"] == "sistemas"

    def test_two_defense_systems_without_island_is_not(self):
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_systems=2)}
        assert _sel(ev)["S1__binette_bin1"]["eligible"] is False

    def test_consensus_arg_alone_is_eligible(self):
        # O AMR entra por direito proprio: um cluster sem ilha e com ARG e
        # exatamente um caso que um portao so-de-defesa perderia.
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_args=1)}
        out = _sel(ev)["S1__binette_bin1"]
        assert out["eligible"] is True
        assert out["criterio"] == "amr"

    def test_plasmid_alone_is_NOT_eligible(self):
        # Plasmidio sem defesa nem ARG e sinal de mobilidade, nao motivo.
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_plasmid=4)}
        assert _sel(ev)["S1__binette_bin1"]["eligible"] is False

    def test_plasmid_reinforces_but_island_names_the_criterion(self):
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_islands=1, n_plasmid=2)}
        out = _sel(ev)["S1__binette_bin1"]
        assert out["eligible"] is True
        assert out["criterio"] == "ilha"

    def test_cluster_without_evidence_row_is_reported_not_dropped(self):
        # Toda linha aparece no candidates.tsv: a selecao tem de ser
        # auditavel, nao um numero que apareceu.
        out = _sel({})
        assert set(out) == {"S1__binette_bin1", "S9__binette_bin2"}
        assert all(c["eligible"] is False for c in out.values())

    def test_n_members_counts_the_membership_rows(self):
        assert _sel({})["S1__binette_bin1"]["n_members"] == 4


import pytest

from pangenome_select import load_completeness, load_membership


class TestLoadMembership:
    def test_groups_members_by_representative(self, tmp_path):
        # nomes realistas: um bin Binette e um inteiro nu do VAMB, que e o
        # que o mag_catalog_pool produz para grupos de co-assembly.
        path = tmp_path / "mag_membership.tsv"
        path.write_text(
            "source_id\toriginal_bin_id\tmember_id\trepresentative_id\n"
            "S1\tbinette_bin1\tS1__binette_bin1\tS1__binette_bin1\n"
            "S2\tbinette_bin3\tS2__binette_bin3\tS1__binette_bin1\n"
            "G1\t7\tG1__7\tS1__binette_bin1\n"
        )
        out = load_membership(str(path))
        assert set(out) == {"S1__binette_bin1"}
        members = out["S1__binette_bin1"]
        assert len(members) == 3
        assert members[0] == {
            "source_id": "S1",
            "original_bin_id": "binette_bin1",
            "member_id": "S1__binette_bin1",
        }
        assert {m["member_id"] for m in members} == {
            "S1__binette_bin1", "S2__binette_bin3", "G1__7",
        }

    def test_skips_row_with_empty_representative_or_member(self, tmp_path):
        path = tmp_path / "mag_membership.tsv"
        path.write_text(
            "source_id\toriginal_bin_id\tmember_id\trepresentative_id\n"
            "S1\tbinette_bin1\t\tS1__binette_bin1\n"
            "S2\tbinette_bin2\tS2__binette_bin2\t\n"
            "S3\tbinette_bin3\tS3__binette_bin3\tS3__binette_bin3\n"
        )
        out = load_membership(str(path))
        assert set(out) == {"S3__binette_bin3"}
        assert len(out["S3__binette_bin3"]) == 1

    def test_accumulates_multiple_members_same_representative(self, tmp_path):
        path = tmp_path / "mag_membership.tsv"
        path.write_text(
            "source_id\toriginal_bin_id\tmember_id\trepresentative_id\n"
            "S1\tbinette_bin1\tS1__binette_bin1\tS1__binette_bin1\n"
            "S2\tbinette_bin1\tS2__binette_bin1\tS1__binette_bin1\n"
            "S3\tbinette_bin1\tS3__binette_bin1\tS1__binette_bin1\n"
            "S4\tbinette_bin1\tS4__binette_bin1\tS1__binette_bin1\n"
        )
        out = load_membership(str(path))
        assert len(out["S1__binette_bin1"]) == 4


class TestLoadCompleteness:
    def test_reads_name_and_completeness(self, tmp_path):
        path = tmp_path / "checkm2_quality_report.tsv"
        path.write_text(
            "Name\tCompleteness\tContamination\n"
            "S1__binette_bin1\t98.5\t1.2\n"
            "S2__binette_bin3\t76.0\t0.5\n"
        )
        out = load_completeness(str(path))
        assert out == {"S1__binette_bin1": 98.5, "S2__binette_bin3": 76.0}

    def test_skips_row_with_unparseable_completeness_keeps_others(self, tmp_path):
        path = tmp_path / "checkm2_quality_report.tsv"
        path.write_text(
            "Name\tCompleteness\n"
            "S1__binette_bin1\tNA\n"
            "S2__binette_bin3\t76.0\n"
        )
        out = load_completeness(str(path))
        assert out == {"S2__binette_bin3": 76.0}

    def test_missing_file_propagates_error(self, tmp_path):
        # Input declarado do Snakemake: faltar e bug de pipeline, nao deve
        # ser mascarado como matriz de completude vazia.
        missing = tmp_path / "does_not_exist.tsv"
        with pytest.raises(OSError):
            load_completeness(str(missing))
