"""
Collapse viral relative abundances by predicted bacterial host genus.

Works with the sylph-tax .sylphmpa merged table when the -a flag was used
(host column present), OR falls back to species-level output without grouping.

Usage:
    python collapse_by_host.py <merged_abundance.tsv> <output.tsv>

Input columns:  clade_name | host (optional) | sample1 | sample2 | ...
Output columns: host_genus | n_viral_taxa | sample1 | sample2 | ...

NOTA sobre a coluna `host`: o `sylph-tax merge` (rules/reads_classify.smk
sylph_merge) NAO propaga a coluna "Virus_host (if viral)" dos .sylphmpa por
amostra para a tabela mesclada. Entao, por este caminho, `has_host` e sempre
False e o colapso por hospedeiro nunca acontece de fato -- escreve-se a saida
em nivel de taxon. Verificado nos dados da Amazonia em 2026-08-19: os
.sylphmpa por amostra TEM a coluna (107 linhas com valor), a mesclada nao.
Levar o host ate aqui exige mudar o merge, o que e decisao de layout e nao
foi feito junto.
"""
import pandas as pd
import sys


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
    """Indices das linhas que sao folha da hierarquia.

    A tabela do sylph-tax traz uma linha por nivel ("r__X", "r__X|k__Y", ...),
    entao somar todas multiplicaria a abundancia. O filtro antigo pegava um
    nivel fixo exigindo "s__" -- mas as linhagens do IMG/VR pulam especie: a
    folha e "t__IMGVR_UViG_...". Nenhum dos 1482 clados virais tinha "s__",
    o que sozinho ja zerava a saida.

    Escolher a folha em vez de um rank fixo funciona para as duas taxonomias
    e nao precisa saber qual banco gerou a tabela.
    """
    clades = list(clades)
    parents = {c.rsplit("|", 1)[0] for c in clades if "|" in c}
    return [c not in parents for c in clades]


def _parse_genus(host_str: str) -> str:
    """Extract genus from a semicolon-delimited GTDB-style lineage."""
    if not host_str or host_str in ("Unknown", "NA", "nan"):
        return "Unknown"
    parts = host_str.split(";")
    # genus is index 5 in d__;p__;c__;o__;f__;g__
    genus_raw = parts[5] if len(parts) > 5 else parts[-1]
    return genus_raw.replace("g__", "").strip() or "Unknown"


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: collapse_by_host.py input_merged.tsv output.tsv")

    df = pd.read_csv(sys.argv[1], sep="\t", comment="#")

    if df.empty:
        df.to_csv(sys.argv[2], sep="\t", index=False)
        return

    taxon_col   = df.columns[0]
    has_host    = "host" in df.columns
    sample_cols = [c for c in df.columns if c not in (taxon_col, "host")]

    # Uma linha por taxon viral: folhas da hierarquia, nao um rank fixo.
    is_viral = df[taxon_col].fillna("").map(_is_viral)
    is_leaf  = _leaf_rows(df[taxon_col].fillna(""))
    viral = df[is_viral & pd.Series(is_leaf, index=df.index)].copy()
    sys.stderr.write(
        "[collapse_by_host] %d linhas, %d virais, %d folhas virais mantidas\n"
        % (len(df), int(is_viral.sum()), len(viral))
    )

    if viral.empty:
        sys.stderr.write(
            "[collapse_by_host] AVISO: nenhuma linha viral. Se a tabela tem "
            "clados virais, o predicado de _is_viral nao reconhece esta "
            "taxonomia -- nao trate a saida vazia como ausencia de virus.\n"
        )
        pd.DataFrame(columns=["host_genus", "n_viral_taxa"] + list(sample_cols)).to_csv(
            sys.argv[2], sep="\t", index=False
        )
        return

    if not has_host:
        sys.stderr.write(
            "[collapse_by_host] no 'host' column found "
            "(re-run sylph-tax with -a for pre-built viral DBs); "
            "writing species-level output unchanged\n"
        )
        viral[[taxon_col] + sample_cols].to_csv(sys.argv[2], sep="\t", index=False)
        return

    viral["host_genus"] = viral["host"].fillna("Unknown").apply(_parse_genus)

    # Sum abundances per host genus, add taxon count
    agg = viral.groupby("host_genus")[sample_cols].sum()
    agg.insert(0, "n_viral_taxa", viral.groupby("host_genus")[taxon_col].count())
    agg = agg.reset_index().sort_values(sample_cols[0], ascending=False)
    agg.to_csv(sys.argv[2], sep="\t", index=False)
    print(f"[collapse_by_host] {len(agg)} host genera written to {sys.argv[2]}")


if __name__ == "__main__":
    main()
