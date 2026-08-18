"""
Convert IMG/VR metadata TSV to sylph-tax custom taxonomy format.

Usage:
    python build_imgvr_taxonomy.py <IMG_VR_meta.tsv> <output_taxonomy.tsv>

Input:  IMG/VR _meta.tsv  (columns: accession, phylum, class, order, family, genus, organism, domain)
Output: two-column TSV for sylph-tax taxprof -t  (genome_id <TAB> taxonomy_string)

Taxonomy string format (sylph-tax custom):
    d__Viruses;p__Uroviricota;c__Caudoviricetes;o__Caudovirales;f__Straboviridae;g__Tequatrovirus;s__organism_name
"""
import sys
import pandas as pd


def _clean(val: str) -> str:
    return str(val).strip() if pd.notna(val) and str(val).strip() not in ("", "nan") else ""


def build_lineage(row) -> str:
    domain  = _clean(row.get("domain",  "Viruses"))
    phylum  = _clean(row.get("phylum",  ""))
    cls     = _clean(row.get("class",   ""))
    order   = _clean(row.get("order",   ""))
    family  = _clean(row.get("family",  ""))
    genus   = _clean(row.get("genus",   ""))
    species = _clean(row.get("organism",""))

    parts = [f"d__{domain or 'Viruses'}"]
    parts.append(f"p__{phylum}" if phylum else "p__")
    parts.append(f"c__{cls}"    if cls    else "c__")
    parts.append(f"o__{order}"  if order  else "o__")
    parts.append(f"f__{family}" if family else "f__")
    parts.append(f"g__{genus}"  if genus  else "g__")
    parts.append(f"s__{species}" if species else "s__")

    return ";".join(parts)


def main():
    if len(sys.argv) != 3:
        sys.exit("Usage: build_imgvr_taxonomy.py <meta.tsv> <output_taxonomy.tsv>")

    meta_path, out_path = sys.argv[1], sys.argv[2]

    df = pd.read_csv(meta_path, sep="\t", dtype=str, low_memory=False)

    if "accession" not in df.columns:
        sys.exit(f"[build_imgvr_taxonomy] 'accession' column not found in {meta_path}")

    df["taxonomy_string"] = df.apply(build_lineage, axis=1)

    out = df[["accession", "taxonomy_string"]].copy()
    out.columns = ["genome_id", "taxonomy"]
    out.to_csv(out_path, sep="\t", index=False, header=False)

    print(f"[build_imgvr_taxonomy] {len(out):,} entries written to {out_path}")
    print(f"  Example: {out.iloc[0].tolist()}")


if __name__ == "__main__":
    main()
