import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from defense_islands import find_islands


def _genes(n, contig="k141_1"):
    """n genes consecutivos com coordenadas de 1000 bp cada."""
    return {contig: [{"Protein": f"{contig}_{i+1}", "Start": i * 1000 + 1,
                      "End": (i + 1) * 1000, "Strand": 1} for i in range(n)]}


class TestFindIslands:
    def test_run_of_defense_genes_is_an_island(self):
        genes = _genes(10)
        prot_to_sys = {f"k141_1_{i}": ("bin1", f"Sys{i}", f"id{i}")
                       for i in range(1, 6)}
        islands = find_islands(genes, prot_to_sys, min_genes=5, min_systems=3)
        assert len(islands) == 1
        assert islands[0]["n_genes"] == 5
        assert islands[0]["n_systems"] == 5
        assert islands[0]["Contig"] == "k141_1"

    def test_too_few_systems_is_not_an_island(self):
        # 5 genes, mas todos do MESMO sistema: nao e ilha.
        genes = _genes(10)
        prot_to_sys = {f"k141_1_{i}": ("bin1", "RM", "id1") for i in range(1, 6)}
        assert find_islands(genes, prot_to_sys, min_genes=5, min_systems=3) == []

    def test_genes_further_apart_than_window_split_into_two_clusters(self):
        genes = _genes(40)
        # 3 genes no inicio, 3 no fim: nenhum grupo alcanca min_genes=5
        hits = [1, 2, 3, 30, 31, 32]
        prot_to_sys = {f"k141_1_{i}": ("bin1", f"Sys{i}", f"id{i}") for i in hits}
        assert find_islands(genes, prot_to_sys, min_genes=5, min_systems=3) == []

    def test_island_carries_genomic_extent(self):
        genes = _genes(10)
        prot_to_sys = {f"k141_1_{i}": ("bin1", f"Sys{i}", f"id{i}")
                       for i in range(1, 6)}
        island = find_islands(genes, prot_to_sys, min_genes=5, min_systems=3)[0]
        assert island["Start_bp"] == 1
        assert island["End_bp"] == 5000

    def test_no_defense_genes_yields_nothing(self):
        assert find_islands(_genes(10), {}, min_genes=5, min_systems=3) == []
