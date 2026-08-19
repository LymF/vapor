#!/usr/bin/env python3
"""Converte ko_per_mag.tsv (longo) para o input do `give_completeness -i`.

    entrada: mag<TAB>ko          (uma linha por par, com cabecalho)
    saida:   mag<TAB>KO<TAB>KO   (uma linha por MAG, sem cabecalho)

Uso: ko_to_wide.py <ko_per_mag.tsv> <saida.tsv>

O nome do MAG e escrito verbatim. O `give_completeness` v1.4.4 devolve esse
nome intacto na coluna `contig` -- verificado com `S1__binette_bin1`, que
sobrevive ao round-trip sem ser cortado em separador nenhum.
"""

import sys
import os

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from annotation_tables import ko_long_to_wide


def main(argv):
    if len(argv) != 3:
        sys.exit(__doc__)
    src, dest = argv[1], argv[2]

    rows = []
    with open(src) as f:
        for i, line in enumerate(f):
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 2:
                continue
            if i == 0 and parts[0] == "mag":   # cabecalho
                continue
            rows.append((parts[0], parts[1]))

    wide = ko_long_to_wide(rows)
    with open(dest, "w") as f:
        for mag in sorted(wide):
            f.write("\t".join([mag] + wide[mag]) + "\n")

    print(f"[ko_to_wide] {len(rows)} pares -> {len(wide)} MAGs")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
