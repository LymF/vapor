"""vrhyme_bins.py — localizar os bins do vRhyme e recuperar o ID do contig.

Existe porque o layout real do vRhyme nao e o que o codigo desta pipeline
assumia. Verificado contra saida real (7 grupos de co-assembly, 2870 headers
de bin, 2026-08-18):

  diretorio : {vrhyme_dir}/vRhyme_best_bins_fasta/
  arquivo   : vRhyme_bin_100.fasta        (junto de .faa e .ffn, que NAO sao bins)
  header    : >vRhyme_100__k141_246312    (prefixo vRhyme_<bin>__ no nome do contig)

O codigo antigo fazia glob("{vrhyme_dir}/vRhyme_best_bins.*.fasta") -- diretorio
errado E padrao de nome errado. Casava ZERO arquivos em 100% dos casos. Os
`vRhyme_best_bins.19.*` que existem na raiz sao .membership.tsv e .summary.tsv,
nunca .fasta, o que fazia o padrao parecer plausivel numa leitura rapida.

Consequencias do bug, todas silenciosas: nenhum vMAG entrava no conjunto
nao-redundante; o contador de bins do log reportava 0; e o dedup do
split_viral_fastas.py (PHIST) nunca disparava, duplicando linha por sequencia.
"""

import glob
import os
import re

_BIN_HEADER = re.compile(r"^vRhyme_\d+__(.+)$")


def bin_fastas(vrhyme_dir):
    """FASTAs de bin dentro de <vrhyme_dir>, ordenados. Lista vazia se nao ha bins."""
    return sorted(glob.glob(
        os.path.join(vrhyme_dir, "vRhyme_best_bins_fasta", "*.fasta")))


def contig_from_bin_header(header):
    """'vRhyme_100__k141_246312' -> 'k141_246312'.

    Aceita o header com ou sem '>' e com ou sem descricao depois do ID.
    Um nome sem o prefixo volta inalterado -- versoes do vRhyme podem nao
    prefixar, e nesse caso o nome ja e o ID do contig.
    """
    name = header[1:] if header.startswith(">") else header
    name = name.split()[0] if name.split() else name
    m = _BIN_HEADER.match(name)
    return m.group(1) if m else name


def binned_contigs(vrhyme_dir):
    """Conjunto dos IDs de contig que estao em algum bin, ja sem o prefixo."""
    out = set()
    for path in bin_fastas(vrhyme_dir):
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    out.add(contig_from_bin_header(line))
    return out
