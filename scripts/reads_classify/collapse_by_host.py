"""
Colapsar abundancia viral (sylph) por genero do hospedeiro predito.

FONTE DO HOSPEDEIRO: a anotacao do BANCO de referencia -- coluna
"Virus_host (if viral)" dos .sylphmpa por amostra, recuperada por
build_host_map.py, ja que o `sylph-tax merge` a descarta. Existe para IMG/VR e
UHGV, e o IMG/VR deixa a maioria como UNKNOWN (medido: 8 de 1407 taxa virais
com genero real nos dados da Amazonia).

`host_source` acompanha `host_genus` de proposito: sem ele, "o banco atribuiu
Acinetobacter" e "ninguem atribuiu nada" ficam indistinguiveis depois da
agregacao. Um genero sem proveniencia e um numero sem procedencia.

Predicao de hospedeiro por PHIST NAO existe neste track -- ver o bloco de
comentario em rules/reads_classify.smk para a decisao e o motivo.

Uso:
    python collapse_by_host.py <merged.tsv> <saida.tsv> [--host-map F]
                               [--assignments F]

Colunas de entrada:  clade_name | sample1 | sample2 | ...
Colunas de saida:    host_genus | host_source | n_viral_taxa | sample1 | ...
"""
import argparse
import csv
import sys

import pandas as pd


_NULL = {"", "na", "nan", "unknown", "none", "-"}
UNKNOWN = "Unknown"


def _is_viral(clade: str) -> bool:
    """Uma linha de clado e viral?

    Era `contains("d__Viruses")`. A taxonomia do IMG/VR no sylph-tax nao usa
    dominio para virus: ela comeca no REALM, "r__Duplodnaviria". Medido nos
    dados da Amazonia em 2026-08-19 -- 1482 clados virais na tabela mesclada,
    ZERO com "d__Viruses", e o viral_abundance_by_host.tsv resultante tinha
    exatamente uma linha, o cabecalho.
    """
    return clade.startswith("r__") or "d__Viruses" in clade


def _leaf_rows(clades):
    """Quais linhas sao folha da hierarquia.

    A tabela do sylph-tax traz uma linha por nivel ("r__X", "r__X|k__Y", ...),
    entao somar todas multiplicaria a abundancia. O filtro antigo pegava um
    nivel fixo exigindo "s__" -- mas as linhagens do IMG/VR pulam especie: a
    folha e "t__IMGVR_UViG_...". Nenhum dos 1482 clados virais tinha "s__",
    o que sozinho ja zerava a saida.

    Escolher a folha em vez de um rank fixo funciona para as duas taxonomias
    e nao exige saber qual banco gerou a tabela.
    """
    clades = list(clades)
    parents = {c.rsplit("|", 1)[0] for c in clades if "|" in c}
    return [c not in parents for c in clades]


def _parse_genus(host_str) -> str:
    """Genero a partir de uma linhagem estilo GTDB separada por ';'.

    Devolve UNKNOWN (grafia unica) quando o rank de genero existe mas esta
    preenchido com "UNKNOWN" -- o IMG/VR faz isso com frequencia, e antes o
    valor cru vazava para a saida, criando duas grafias de desconhecido na
    mesma coluna ("Unknown" e "UNKNOWN").
    """
    host_str = (host_str or "").strip()
    if host_str.lower() in _NULL:
        return UNKNOWN
    parts = host_str.split(";")
    # genero e o indice 5 em d__;p__;c__;o__;f__;g__
    genus_raw = parts[5] if len(parts) > 5 else parts[-1]
    genus = genus_raw.replace("g__", "").strip()
    return UNKNOWN if genus.lower() in _NULL else genus


def _load_map(path):
    """TSV de duas colunas (clade_name, linhagem) -> dict. Vazio se sem path."""
    if not path:
        return {}
    out = {}
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 2 or row[0].startswith("#") or row[0] == "clade_name":
                continue
            if row[0].strip():
                out[row[0].strip()] = row[1].strip()
    return out


def resolve_host(clade, db_map):
    """(genero, fonte) para um clado. Fonte e "db" ou "none".

    Devolver a fonte junto, e nao so o genero, e o ponto: depois do groupby um
    "Unknown" sem procedencia se confunde com um hospedeiro que o banco de fato
    nao conhece.
    """
    g_db = _parse_genus(db_map.get(clade, ""))
    return (g_db, "db") if g_db != UNKNOWN else (UNKNOWN, "none")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("input")
    ap.add_argument("output")
    ap.add_argument("--host-map", default="",
                    help="TSV clade_name/host_db (build_host_map.py)")
    ap.add_argument("--assignments", default="",
                    help="TSV auditavel por taxon: as duas fontes lado a lado")
    args = ap.parse_args()

    db_map = _load_map(args.host_map)
    sys.stderr.write("[collapse_by_host] host_db: %d clados com hospedeiro\n"
                     % len(db_map))

    df = pd.read_csv(args.input, sep="\t", comment="#")
    taxon_col   = df.columns[0] if len(df.columns) else "clade_name"
    sample_cols = [c for c in df.columns[1:]]

    out_cols = ["host_genus", "host_source", "n_viral_taxa"] + list(sample_cols)
    if df.empty:
        pd.DataFrame(columns=out_cols).to_csv(args.output, sep="\t", index=False)
        return

    is_viral = df[taxon_col].fillna("").map(_is_viral)
    is_leaf  = pd.Series(_leaf_rows(df[taxon_col].fillna("")), index=df.index)
    viral = df[is_viral & is_leaf].copy()
    sys.stderr.write("[collapse_by_host] %d linhas, %d virais, %d folhas virais\n"
                     % (len(df), int(is_viral.sum()), len(viral)))

    if viral.empty:
        sys.stderr.write(
            "[collapse_by_host] AVISO: nenhuma linha viral. Se a tabela tem "
            "clados virais, o predicado de _is_viral nao reconhece esta "
            "taxonomia -- nao trate a saida vazia como ausencia de virus.\n")
        pd.DataFrame(columns=out_cols).to_csv(args.output, sep="\t", index=False)
        return

    resolved = [resolve_host(c, db_map) for c in viral[taxon_col]]
    viral["host_genus"]  = [r[0] for r in resolved]
    viral["host_source"] = [r[1] for r in resolved]

    if args.assignments:
        aud = pd.DataFrame({
            "clade_name":  viral[taxon_col].values,
            "host_db":     [db_map.get(c, "") for c in viral[taxon_col]],
            "host_genus":  viral["host_genus"].values,
            "host_source": viral["host_source"].values,
        })
        aud.to_csv(args.assignments, sep="\t", index=False)

    # Agrupa por GENERO so. Agrupar por (genero, fonte) partiria o mesmo
    # genero em duas linhas quando parte dos seus virus vem do banco e parte
    # do PHIST -- no grafico isso apareceria como dois taxa distintos com o
    # mesmo nome. A proveniencia vira um campo do grupo ("db", "phist" ou
    # "db,phist"); o detalhe por taxon fica no sidecar de assignments.
    grp = viral.groupby("host_genus")
    agg = grp[sample_cols].sum()
    agg.insert(0, "n_viral_taxa", grp[taxon_col].count())
    sources = grp["host_source"].apply(
        lambda col: ",".join(sorted(set(col) - {"none"}) or ["none"]))  # noqa: E501
    agg.insert(0, "host_source", sources)
    agg = agg.reset_index()
    if sample_cols:
        agg = agg.sort_values(sample_cols[0], ascending=False)
    agg[out_cols].to_csv(args.output, sep="\t", index=False)

    by_src = viral["host_source"].value_counts().to_dict()
    sys.stderr.write("[collapse_by_host] %d generos; taxa por fonte: %s -> %s\n"
                     % (len(agg), by_src, args.output))


if __name__ == "__main__":
    main()
