"""
Collapse viral relative abundances by predicted bacterial host genus.

Works with the sylph-tax .sylphmpa merged table when the -a flag was used
(host column present), OR falls back to species-level output without grouping.

Usage:
    python collapse_by_host.py <merged_abundance.tsv> <output.tsv>

Input columns:  clade_name | host (optional) | sample1 | sample2 | ...
Output columns: host_genus | n_viral_taxa | sample1 | sample2 | ...
"""
import pandas as pd
import sys


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

    # Keep only viral species rows
    viral = df[
        df[taxon_col].str.contains("d__Viruses", na=False) &
        df[taxon_col].str.contains("s__", na=False)
    ].copy()

    if viral.empty:
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
