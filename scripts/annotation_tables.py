"""Funcoes puras das tabelas funcionais derivadas do eggNOG-mapper.

Duas tabelas saem do MESMO arquivo que o `mag_eggnog_prok` ja grava
(`eggnog_annotations.tsv`), sem ferramenta nem banco novo:

  - `ko_per_mag.tsv`   (coluna 12, `KEGG_ko`)   -- ja existia
  - `cazy_per_mag.tsv` (coluna 19, `CAZy`)      -- entrou em 2026-08-19

Os indices sao 1-based porque foram conferidos contra o cabecalho real de
`eggnog_annotations.emapper.annotations` em `~/global/results`, nao contra a
documentacao do emapper. Em 46.694 proteinas daquele MAG, 641 tem CAZy
nao-vazio -- o emapper escreve "-" para ausencia, nunca campo vazio.

`ko_long_to_wide` existe porque o `give_completeness -i` quer uma linha por
genoma (`nome<TAB>KO<TAB>KO...`) e o `ko_per_mag.tsv` e longo (um par por
linha). Verificado com a ferramenta de verdade (v1.4.4): ela devolve o nome
na coluna `contig` **intacto**, entao nomes como `S1__binette_bin1` passam
sem o corte em separador que assombra o resto desta pipeline.
"""

CAZY_COL = 19   # 1-based, conferido no cabecalho real
KEGG_KO_COL = 12

# Ordem importa: a checagem e por prefixo e "CBM" tem de ser testada antes de
# "CE"/"C*" qualquer. Sao as seis classes do CAZy.
_CAZY_CLASSES = ("CBM", "GH", "GT", "PL", "CE", "AA")


def parse_cazy_field(value):
    """Familias CAZy de um campo do emapper, em ordem de aparicao.

    Multi-familia vem separada por virgula (`CBM48,GH13`, 34 proteinas reais
    em ERR4682430). Ausencia e "-". O sufixo de subfamilia (`GH13_20`) fica:
    joga-lo fora perderia a resolucao que motiva anotar CAZy.
    """
    if not value:
        return []
    text = value.strip()
    if not text or text == "-":
        return []
    return [fam.strip() for fam in text.split(",") if fam.strip() and fam.strip() != "-"]


def cazy_class(family):
    """Classe CAZy de uma familia: GH, GT, PL, CE, AA, CBM ou 'Other'."""
    if not family:
        return "Other"
    fam = family.strip().upper()
    for cls in _CAZY_CLASSES:
        if fam.startswith(cls):
            return cls
    return "Other"


def ko_long_to_wide(rows):
    """`[(mag, ko), ...]` -> `{mag: [ko ordenado e sem repeticao]}`.

    Formato de entrada do `give_completeness -i`. A multiplicidade nao
    interessa: o que a completude de modulo pergunta e se o KO esta
    presente no genoma, nao em quantas copias.
    """
    wide = {}
    for mag, ko in rows:
        wide.setdefault(mag, set()).add(ko)
    return {mag: sorted(kos) for mag, kos in wide.items()}
