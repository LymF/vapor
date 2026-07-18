"""
Reformat a sylph-tax merged table to standard OTU-table format.

Usage:
    python make_otu.py <input.tsv> <output.tsv>

Renames the first column (clade_name) to #OTU_ID so the table is
compatible with QIIME2, phyloseq, and other downstream tools.
"""
import pandas as pd
import sys


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: make_otu.py input.tsv output.tsv")

    df = pd.read_csv(sys.argv[1], sep="\t", comment="#")
    df = df.rename(columns={df.columns[0]: "#OTU_ID"})
    df.to_csv(sys.argv[2], sep="\t", index=False)
    print(f"[make_otu] wrote {len(df)} OTUs to {sys.argv[2]}")


if __name__ == "__main__":
    main()
