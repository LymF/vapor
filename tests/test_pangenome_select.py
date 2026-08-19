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
