"""
Filter a sylph-tax merged abundance table by minimum prevalence.

Usage:
    python filter_by_prevalence.py <input.tsv> <output.tsv> <min_prevalence>

min_prevalence: fraction 0–1 (e.g. 0.1 = taxon must be present in >=10% of samples).
Input TSV: first column is clade_name, remaining columns are samples.

O corte e >= min_prevalence E prevalencia > 0. Os dois testes sao precisos:
o `>=` e o que a documentacao (e a literatura) chama de prevalencia minima --
com `>` puro, um corte de 0.1 descartava justamente os taxa em exatamente 10%
das amostras; e o `> 0` mantem o comportamento util do default 0.0, que e
jogar fora as linhas zeradas em toda amostra (o sylph-tax merge as emite).
"""
import pandas as pd
import sys


def main():
    if len(sys.argv) != 4:
        sys.exit("Usage: filter_by_prevalence.py input.tsv output.tsv min_prevalence")

    path_in, path_out, min_prev = sys.argv[1], sys.argv[2], float(sys.argv[3])

    df = pd.read_csv(path_in, sep="\t", comment="#")
    if df.empty:
        df.to_csv(path_out, sep="\t", index=False)
        return

    taxon_col  = df.columns[0]
    sample_cols = df.columns[1:]

    df["prevalence"] = (df[sample_cols] > 0).sum(axis=1) / len(sample_cols)
    keep = (df["prevalence"] >= min_prev) & (df["prevalence"] > 0)
    df = df[keep].drop(columns="prevalence")
    df.to_csv(path_out, sep="\t", index=False)
    print(f"[filter_by_prevalence] kept {len(df)} taxa "
          f"(prevalence >= {min_prev} and > 0)")


if __name__ == "__main__":
    main()
