"""
Filter a sylph-tax merged abundance table by minimum prevalence.

Usage:
    python filter_by_prevalence.py <input.tsv> <output.tsv> <min_prevalence>

min_prevalence: fraction 0–1 (e.g. 0.1 = taxon must be present in ≥10% of samples).
Input TSV: first column is clade_name, remaining columns are samples.
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
    df = df[df["prevalence"] > min_prev].drop(columns="prevalence")
    df.to_csv(path_out, sep="\t", index=False)
    print(f"[filter_by_prevalence] kept {len(df)} taxa (prevalence > {min_prev})")


if __name__ == "__main__":
    main()
