import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from defense_islands import find_islands, genes_by_contig


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


class TestGenesByContig:
    def test_parses_coordinates_and_strand(self, tmp_path):
        # header real do Prodigal: >{seqid}_{n} # {start} # {end} # {strand} # ID=...
        faa = tmp_path / "genome.faa"
        faa.write_text(
            ">k141_1_1 # 3 # 1400 # 1 # ID=1_1;partial=00;start_type=ATG\n"
            "MSEQ\n"
            ">k141_1_2 # 1500 # 2600 # -1 # ID=1_2;partial=00;start_type=ATG\n"
            "MSEQ\n"
        )
        by_contig = genes_by_contig(str(faa))
        genes = by_contig["k141_1"]
        assert genes[0] == {"Protein": "k141_1_1", "Start": 3, "End": 1400,
                            "Strand": 1}
        assert genes[1] == {"Protein": "k141_1_2", "Start": 1500, "End": 2600,
                            "Strand": -1}

    def test_groups_by_contig_in_file_order(self, tmp_path):
        faa = tmp_path / "genome.faa"
        faa.write_text(
            ">k141_1_1 # 1 # 100 # 1 # ID=1_1\n"
            "MSEQ\n"
            ">k141_2_1 # 1 # 100 # 1 # ID=2_1\n"
            "MSEQ\n"
            ">k141_1_2 # 200 # 300 # 1 # ID=1_2\n"
            "MSEQ\n"
        )
        by_contig = genes_by_contig(str(faa))
        assert list(by_contig) == ["k141_1", "k141_2"]
        assert [g["Protein"] for g in by_contig["k141_1"]] == [
            "k141_1_1", "k141_1_2"]
        assert [g["Protein"] for g in by_contig["k141_2"]] == ["k141_2_1"]

    def test_header_without_coordinate_fields_keeps_gene_with_none(self, tmp_path):
        # header sem os campos de coordenada: o gene tem de ficar, com
        # Start/End/Strand em None, nao ser descartado.
        faa = tmp_path / "genome.faa"
        faa.write_text(">k141_1_1 some free-text description\nMSEQ\n")
        by_contig = genes_by_contig(str(faa))
        genes = by_contig["k141_1"]
        assert len(genes) == 1
        assert genes[0] == {"Protein": "k141_1_1", "Start": None,
                            "End": None, "Strand": None}

    def test_missing_file_returns_empty(self, tmp_path):
        missing = tmp_path / "does_not_exist.faa"
        assert genes_by_contig(str(missing)) == {}
