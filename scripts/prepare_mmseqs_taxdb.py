#!/usr/bin/env python3
"""
prepare_mmseqs_taxdb.py — Build an MMseqs2 seqTaxDB (custom taxonomy,
LCA-capable) from the IMG NR protein FASTA + taxonOId2Taxonomy.tsv lineage
file already used by prepare_diamond_db.py.

Why: diamond_custom_prok's best-hit + majority-vote approach is prone to
"spurious specificity" -- a single highly-conserved or HGT-acquired region
can give a confident species-level best hit even when the organism is
actually novel/divergent (documented for best-hit-only methods, von
Meijenfeldt et al. 2019, the CAT/BAT paper). This matters more, not less,
for IMG_NR specifically: it's a custom database in the first place because
standard reference databases under-represent the environmental/divergent
organisms it's meant to classify. `mmseqs taxonomy` computes a real
lowest-common-ancestor per query across all its hits instead, resolving to
the correct (coarser) rank when hits disagree, rather than committing to
one bin's species call by majority vote of single best-hits.

nodes.dmp/names.dmp use the minimal custom-taxonomy format MMseqs2 itself
documents for non-NCBI taxonomies (confirmed against a real working
example, GTDB-style conversion: github.com/soedinglab/MMseqs2/issues/976):
    nodes.dmp: "{taxid}\t|\t{parent_taxid}\t|\t{rank}\t|\t-\t|"
    names.dmp: "{taxid}\t|\t{name}\t|\t-\t|\tscientific name\t|"

IMG lineage strings are semicolon-separated, always in this rank order
(confirmed against a real taxonOId2Taxonomy.tsv on litrp4):
    Domain;Phylum;Class;Order;Family;Genus;Species;Strain
Shorter lineages (fewer than 8 levels) are handled -- whatever's missing
just isn't assigned a node, the leaf taxid is whatever rank is deepest.

Two-pass design, memory bounded by the number of GENOMES (taxonOId2Taxonomy
.tsv, tens-hundreds of thousands of rows), not the number of PROTEINS
(the raw FASTA can have hundreds of millions of headers) -- pass 1 builds
the small taxon_oid -> leaf_taxid dict from the lineage file; pass 2
streams the FASTA headers and writes the mapping file line by line,
never holding all protein IDs in memory at once.

Usage:
  python3 prepare_mmseqs_taxdb.py \\
      --faa img_unrestricted_isolates_nr.faa \\
      --img-tax taxonOId2Taxonomy.tsv \\
      --out /path/to/img_nr_mmseqs \\
      --threads 32
"""
import argparse, os, subprocess, sys

RANKS = ['domain', 'phylum', 'class', 'order', 'family', 'genus', 'species', 'strain']


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--faa', required=True,
                   help='Raw IMG protein FASTA (img_unrestricted_isolates_nr.faa) -- '
                        'headers must be ">HASH length TAXON_OID:GENE_OID"')
    p.add_argument('--img-tax', required=True, help='IMG taxonOId2Taxonomy.tsv')
    p.add_argument('--out', required=True, help='Output dir for the seqTaxDB + taxdump files')
    p.add_argument('--threads', type=int, default=16)
    return p.parse_args()


def build_taxdump(img_tax_path, taxdump_dir):
    """Parse taxonOId2Taxonomy.tsv lineages into nodes.dmp/names.dmp.
    Returns {taxon_oid: leaf_taxid}."""
    os.makedirs(taxdump_dir, exist_ok=True)
    nodes_path  = os.path.join(taxdump_dir, 'nodes.dmp')
    names_path  = os.path.join(taxdump_dir, 'names.dmp')
    merged_path = os.path.join(taxdump_dir, 'merged.dmp')
    # mmseqs createtaxdb's createbintaxonomy step requires merged.dmp to
    # exist even for a from-scratch custom taxonomy with no merged taxon
    # IDs to report (confirmed live on litrp4: "Input .../merged.dmp does
    # not exist" -- the GTDB-style conversion example that documents the
    # minimal nodes.dmp/names.dmp format doesn't mention needing this file,
    # but a real run against mmseqs2 18.8cc5c does). Empty file, NCBI
    # format ("old_taxid\t|\tnew_taxid\t|" per line) -- zero lines means
    # no merges, which is exactly true here.
    open(merged_path, 'w').close()

    node_id = {(): 1}   # path-tuple (ancestor names + this name) -> synthetic taxid
    next_id = 2
    toid_to_leaf = {}

    with open(nodes_path, 'w') as nf, open(names_path, 'w') as mf, \
         open(img_tax_path) as f:
        nf.write("1\t|\t1\t|\tno rank\t|\t-\t|\n")
        mf.write("1\t|\troot\t|\t-\t|\tscientific name\t|\n")

        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 2:
                continue
            taxon_oid = parts[0].strip()
            lineage = [x.strip() for x in parts[1].split(';') if x.strip()]
            if not lineage:
                continue

            path = ()
            parent_id = 1
            leaf_id = 1
            for depth, name in enumerate(lineage):
                path = path + (name,)
                existing = node_id.get(path)
                if existing is not None:
                    parent_id = existing
                    leaf_id = existing
                    continue
                this_id = next_id
                next_id += 1
                node_id[path] = this_id
                rank = RANKS[depth] if depth < len(RANKS) else 'no rank'
                nf.write(f"{this_id}\t|\t{parent_id}\t|\t{rank}\t|\t-\t|\n")
                mf.write(f"{this_id}\t|\t{name}\t|\t-\t|\tscientific name\t|\n")
                parent_id = this_id
                leaf_id = this_id
            toid_to_leaf[taxon_oid] = leaf_id

    print(f"[prepare_mmseqs_taxdb] {next_id - 1:,} taxonomy nodes written, "
          f"{len(toid_to_leaf):,} genomes (taxonOIds) mapped")
    return toid_to_leaf


def write_mapping(faa_path, toid_to_leaf, mapping_path):
    """Stream the raw FASTA headers (>HASH length TAXON_OID:GENE_OID),
    write '{hash}\t{leaf_taxid}' per protein -- never holds all headers in
    memory (the FASTA can have hundreds of millions of them)."""
    n_written = n_skipped = 0
    with open(faa_path) as f, open(mapping_path, 'w') as out:
        for line in f:
            if not line.startswith('>'):
                continue
            parts = line[1:].split()
            if len(parts) < 3:
                n_skipped += 1
                continue
            h = parts[0]
            taxon_oid = parts[2].split(':')[0]
            taxid = toid_to_leaf.get(taxon_oid)
            if taxid is None:
                n_skipped += 1
                continue
            out.write(f"{h}\t{taxid}\n")
            n_written += 1
    print(f"[prepare_mmseqs_taxdb] mapping written: {n_written:,} proteins mapped, "
          f"{n_skipped:,} skipped (no taxonomy entry for their taxonOId)")
    return n_written


def main():
    args = parse_args()
    os.makedirs(args.out, exist_ok=True)
    taxdump_dir  = os.path.join(args.out, 'taxdump')
    mapping_path = os.path.join(args.out, 'mapping.tsv')
    seqdb_path   = os.path.join(args.out, 'seqTaxDB')
    tmp_dir      = os.path.join(args.out, 'tmp')

    toid_to_leaf = build_taxdump(args.img_tax, taxdump_dir)
    n_mapped = write_mapping(args.faa, toid_to_leaf, mapping_path)
    if n_mapped == 0:
        print("[prepare_mmseqs_taxdb] ERROR: no proteins mapped to a taxonomy node -- "
              "check --faa header format and --img-tax content", file=sys.stderr)
        sys.exit(1)

    if os.path.exists(seqdb_path + '.dbtype'):
        print(f"[prepare_mmseqs_taxdb] {seqdb_path} already exists -- "
              "skipping mmseqs createdb (delete it to force a rebuild)")
    else:
        print("[prepare_mmseqs_taxdb] Building MMseqs2 sequence DB (mmseqs createdb)...")
        subprocess.run(['mmseqs', 'createdb', args.faa, seqdb_path], check=True)

    print("[prepare_mmseqs_taxdb] Running mmseqs createtaxdb "
          "(can take a while for large protein sets)...")
    subprocess.run(['mmseqs', 'createtaxdb', seqdb_path, tmp_dir,
                     '--ncbi-tax-dump', taxdump_dir,
                     '--tax-mapping-file', mapping_path,
                     '--threads', str(args.threads)], check=True)

    print(f"\n[prepare_mmseqs_taxdb] Done!")
    print(f"  seqTaxDB: {seqdb_path}")
    print(f"\nAdd to config.yaml:")
    print(f"  custom_prok_mmseqs_db: \"{seqdb_path}\"")


if __name__ == '__main__':
    main()
