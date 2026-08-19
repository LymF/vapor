"""Testes das funcoes puras de scripts/annotation_tables.py.

Cobrem os dois formatos verificados contra arquivos reais em 2026-08-19:

  - a coluna CAZy do eggNOG-mapper (coluna 19 de
    `eggnog_annotations.emapper.annotations`, valores como `CBM48,GH13`);
  - o input do `give_completeness -i`, que e largo
    (`nome<TAB>KO<TAB>KO...`) enquanto o `ko_per_mag.tsv` e longo.
"""

import os
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from annotation_tables import cazy_class, ko_long_to_wide, parse_cazy_field


class TestParseCazyField:
    def test_single_family(self):
        assert parse_cazy_field("GH13") == ["GH13"]

    def test_multiple_families_are_comma_separated(self):
        # Valor real do ERR4682430: 34 proteinas com exatamente este campo.
        assert parse_cazy_field("CBM48,GH13") == ["CBM48", "GH13"]

    def test_empty_and_dash_yield_nothing(self):
        # O emapper escreve "-" para ausencia, nao vazio.
        assert parse_cazy_field("-") == []
        assert parse_cazy_field("") == []
        assert parse_cazy_field("   ") == []

    def test_whitespace_is_stripped(self):
        assert parse_cazy_field(" GT2 , GH3 ") == ["GT2", "GH3"]

    def test_subfamily_suffix_is_preserved(self):
        # GH13_20 e uma subfamilia; jogar fora o sufixo perderia resolucao.
        assert parse_cazy_field("GH13_20") == ["GH13_20"]


class TestCazyClass:
    @pytest.mark.parametrize(
        "family,expected",
        [
            ("GH13", "GH"),
            ("GT51", "GT"),
            ("PL1", "PL"),
            ("CE4", "CE"),
            ("AA10", "AA"),
            ("CBM48", "CBM"),
            ("GH13_20", "GH"),
        ],
    )
    def test_known_classes(self, family, expected):
        assert cazy_class(family) == expected

    def test_cbm_wins_over_ce(self):
        # "CBM48" comeca com "C"; um prefixo mais curto casado primeiro
        # devolveria "CE" para nada e "CBM" para nada. Ordem importa.
        assert cazy_class("CBM48") == "CBM"

    def test_unknown_family_is_other(self):
        assert cazy_class("XYZ9") == "Other"
        assert cazy_class("") == "Other"


class TestKoLongToWide:
    def test_groups_kos_by_mag(self):
        rows = [("S1__binette_bin1", "K00844"), ("S1__binette_bin1", "K01810"),
                ("G2__3", "K00845")]
        assert ko_long_to_wide(rows) == {
            "S1__binette_bin1": ["K00844", "K01810"],
            "G2__3": ["K00845"],
        }

    def test_deduplicates_and_sorts(self):
        # Varias proteinas do mesmo MAG carregam o mesmo KO; o
        # give_completeness nao precisa da multiplicidade.
        rows = [("m", "K00002"), ("m", "K00001"), ("m", "K00002")]
        assert ko_long_to_wide(rows) == {"m": ["K00001", "K00002"]}

    def test_mag_names_with_double_underscore_survive_intact(self):
        # O nome do catalogo e {source}__{bin} e o do MAG membro pode ser
        # um inteiro do VAMB. Nada aqui pode cortar em separador.
        rows = [("S1__binette_bin1", "K00844")]
        assert list(ko_long_to_wide(rows)) == ["S1__binette_bin1"]

    def test_empty_input(self):
        assert ko_long_to_wide([]) == {}
