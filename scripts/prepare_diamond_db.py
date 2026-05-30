#!/usr/bin/env python3
"""
prepare_diamond_db.py — Prepare a Diamond DB + standardised metadata TSV
for use with the custom taxonomy tiers in the metagenomics pipeline.

Supports input formats:
  --format img    : IMG NR (hash headers + taxonOId2Taxonomy.tsv)
  --format ncbi   : NCBI RefSeq/GenBank (>accession.version description [organism])
  --format simple : Already has header format >accession [family=X] [genus=Y]
  --format tsv    : Custom FASTA + separate accession→taxonomy TSV

Output:
  <out>.dmnd       Diamond database
  <out>_meta.tsv   Standard metadata: accession\tphylum\tclass\torder\tfamily\tgenus\torganism

Usage examples:
  # IMG NR (viral + prokaryote)
  python3 prepare_diamond_db.py \\
      --faa img_unrestricted_isolates_nr.faa \\
      --format img \\
      --img-tax taxonOId2Taxonomy.tsv \\
      --out /path/to/imgnr \\
      --threads 32

  # NCBI RefSeq viral
  python3 prepare_diamond_db.py \\
      --faa viral.1.protein.faa \\
      --format ncbi \\
      --ncbi-tax taxonomy.tsv \\
      --out /path/to/refseq_viral \\
      --threads 32

  # Filter to only viral entries (for CUSTOM_VIRAL_DMND)
  python3 prepare_diamond_db.py \\
      --faa img_unrestricted_isolates_nr.faa \\
      --format img \\
      --img-tax taxonOId2Taxonomy.tsv \\
      --filter-kingdom Viruses \\
      --out /path/to/imgnr_viral \\
      --threads 32
"""

import argparse, csv, os, sys, subprocess, tempfile
from collections import defaultdict

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument('--faa',       required=True, help='Input protein FASTA')
    p.add_argument('--format',    required=True,
                   choices=['img','ncbi','simple','tsv'],
                   help='Input format (see usage examples)')
    p.add_argument('--img-tax',   help='IMG taxonOId2Taxonomy.tsv (required for --format img)')
    p.add_argument('--ncbi-tax',  help='NCBI accession2taxid or lineage TSV (optional for --format ncbi)')
    p.add_argument('--custom-tsv',help='Custom accession→taxonomy TSV (required for --format tsv)')
    p.add_argument('--filter-kingdom', default='',
                   help='Only keep entries matching this kingdom (e.g. Viruses, Bacteria)')
    p.add_argument('--out',       required=True,
                   help='Output prefix (e.g. /path/to/imgnr_viral). '
                        'Creates <out>.dmnd and <out>_meta.tsv')
    p.add_argument('--threads',   type=int, default=16)
    p.add_argument('--min-len',   type=int, default=0,
                   help='Minimum protein length to include (default: no filter)')
    return p.parse_args()


def parse_lineage(lineage_str):
    """Parse semicolon-separated lineage into standard fields."""
    parts = [p.strip() for p in lineage_str.split(';') if p.strip()]
    # IMG/standard format: Domain;Phylum;Class;Order;Family;Genus;Species;Strain
    def get(idx, default=''):
        return parts[idx] if idx < len(parts) else default
    return {
        'domain':   get(0),
        'phylum':   get(1),
        'class':    get(2),
        'order':    get(3),
        'family':   get(4),
        'genus':    get(5),
        'organism': get(6) or get(5),
    }


def load_img_taxonomy(tax_path):
    """Load IMG taxonOId2Taxonomy.tsv → {taxon_oid: lineage_dict}"""
    print(f"[prepare_diamond_db] Loading IMG taxonomy from {tax_path}...")
    tax = {}
    with open(tax_path) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                tax[parts[0]] = parse_lineage(parts[1])
    print(f"[prepare_diamond_db] Loaded {len(tax):,} taxonomy entries")
    return tax


def load_ncbi_taxonomy(tax_path):
    """Load NCBI lineage TSV → {accession: lineage_dict}"""
    print(f"[prepare_diamond_db] Loading NCBI taxonomy from {tax_path}...")
    tax = {}
    with open(tax_path) as f:
        rdr = csv.DictReader(f, delimiter='\t')
        for row in rdr:
            acc = (row.get('accession','') or row.get('Accession','')).strip()
            if not acc: continue
            tax[acc] = {
                'domain':   row.get('domain',  row.get('Domain','')),
                'phylum':   row.get('phylum',  row.get('Phylum','')),
                'class':    row.get('class',   row.get('Class','')),
                'order':    row.get('order',   row.get('Order','')),
                'family':   row.get('family',  row.get('Family','')),
                'genus':    row.get('genus',   row.get('Genus','')),
                'organism': row.get('organism',row.get('Organism',
                            row.get('species', row.get('Species','')))),
            }
    print(f"[prepare_diamond_db] Loaded {len(tax):,} taxonomy entries")
    return tax


def process_img(faa_path, img_tax, filter_kingdom, min_len):
    """
    Parse IMG FAA headers: >HASH length TAXON_OID:GENE_OID
    Build hash→taxon_oid map, join with taxonomy.
    Returns: (filtered_faa_path, meta_dict {hash: lineage_dict})
    """
    print("[prepare_diamond_db] Parsing IMG FAA headers...")
    hash2taxon = {}
    with open(faa_path) as f:
        for line in f:
            if line.startswith('>'):
                parts = line[1:].split()
                if len(parts) >= 3:
                    h = parts[0]
                    taxon_oid = parts[2].split(':')[0]
                    hash2taxon[h] = taxon_oid

    print(f"[prepare_diamond_db] {len(hash2taxon):,} protein headers parsed")

    meta = {}
    kept = set()
    filter_lower = filter_kingdom.lower()

    for h, toid in hash2taxon.items():
        lin = img_tax.get(toid, {})
        if filter_lower and filter_lower not in lin.get('domain','').lower():
            continue
        meta[h] = lin
        kept.add(h)

    print(f"[prepare_diamond_db] {len(kept):,} proteins after kingdom filter "
          f"('{filter_kingdom or 'none'}')")
    return meta, kept


def process_ncbi(faa_path, ncbi_tax, filter_kingdom, min_len):
    """Parse NCBI FAA: >accession.version description"""
    print("[prepare_diamond_db] Parsing NCBI FAA headers...")
    meta = {}
    kept = set()
    filter_lower = filter_kingdom.lower()

    with open(faa_path) as f:
        for line in f:
            if line.startswith('>'):
                acc = line[1:].split()[0]
                # Try versioned and unversioned accession
                lin = (ncbi_tax.get(acc) or
                       ncbi_tax.get(acc.split('.')[0]) or {})
                if filter_lower and filter_lower not in lin.get('domain','').lower():
                    continue
                meta[acc] = lin
                kept.add(acc)

    print(f"[prepare_diamond_db] {len(kept):,} proteins")
    return meta, kept


def write_filtered_faa(faa_path, kept_ids, out_faa, min_len=0):
    """Write filtered FASTA keeping only kept_ids."""
    print(f"[prepare_diamond_db] Writing filtered FAA → {out_faa}")
    written = 0; current_id = None; current_seq = []

    def flush():
        nonlocal written
        if current_id and current_id in kept_ids:
            seq = ''.join(current_seq)
            if len(seq) >= min_len:
                out.write('>' + current_id + '\n' + seq + '\n')
                written += 1

    with open(faa_path) as f, open(out_faa, 'w') as out:
        for line in f:
            if line.startswith('>'):
                flush()
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line.strip())
        flush()

    print(f"[prepare_diamond_db] {written:,} sequences written")
    return written


def write_meta_tsv(meta, out_path):
    """Write standardised metadata TSV."""
    print(f"[prepare_diamond_db] Writing metadata → {out_path}")
    fields = ['accession','phylum','class','order','family','genus','organism','domain']
    with open(out_path, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        w.writeheader()
        for acc, lin in meta.items():
            w.writerow({'accession': acc,
                        'phylum':   lin.get('phylum',''),
                        'class':    lin.get('class',''),
                        'order':    lin.get('order',''),
                        'family':   lin.get('family',''),
                        'genus':    lin.get('genus',''),
                        'organism': lin.get('organism',''),
                        'domain':   lin.get('domain',''),
                       })
    print(f"[prepare_diamond_db] {len(meta):,} entries in metadata TSV")


def build_diamond(faa_path, dmnd_path, threads):
    """Build Diamond database."""
    print(f"[prepare_diamond_db] Building Diamond DB → {dmnd_path}")
    cmd = ['diamond', 'makedb',
           '--in', faa_path,
           '--db', dmnd_path,
           '--threads', str(threads)]
    ret = subprocess.run(cmd)
    if ret.returncode != 0:
        print("[prepare_diamond_db] ERROR: diamond makedb failed", file=sys.stderr)
        sys.exit(1)
    print(f"[prepare_diamond_db] Diamond DB built: {dmnd_path}.dmnd")


def main():
    args = parse_args()

    out_dmnd = args.out + '.dmnd' if not args.out.endswith('.dmnd') else args.out
    out_meta = args.out + '_meta.tsv'
    out_faa  = args.out + '_filtered.faa'

    # ── Parse format ─────────────────────────────────────────────────
    if args.format == 'img':
        if not args.img_tax:
            print("ERROR: --img-tax required for --format img", file=sys.stderr)
            sys.exit(1)
        img_tax = load_img_taxonomy(args.img_tax)
        meta, kept = process_img(args.faa, img_tax,
                                 args.filter_kingdom, args.min_len)

    elif args.format == 'ncbi':
        ncbi_tax = load_ncbi_taxonomy(args.ncbi_tax) if args.ncbi_tax else {}
        meta, kept = process_ncbi(args.faa, ncbi_tax,
                                  args.filter_kingdom, args.min_len)

    elif args.format == 'tsv':
        if not args.custom_tsv:
            print("ERROR: --custom-tsv required for --format tsv", file=sys.stderr)
            sys.exit(1)
        meta = {}; kept = set()
        with open(args.custom_tsv) as f:
            for row in csv.DictReader(f, delimiter='\t'):
                acc = row.get('accession','').strip()
                if not acc: continue
                fl = args.filter_kingdom.lower()
                if fl and fl not in row.get('domain','').lower(): continue
                meta[acc] = row; kept.add(acc)

    elif args.format == 'simple':
        # Headers already have accession as first word, no filtering needed
        meta = {}; kept = set()
        with open(args.faa) as f:
            for line in f:
                if line.startswith('>'):
                    acc = line[1:].split()[0]
                    meta[acc] = {}; kept.add(acc)

    # ── Write filtered FAA ────────────────────────────────────────────
    n = write_filtered_faa(args.faa, kept, out_faa, args.min_len)
    if n == 0:
        print("[prepare_diamond_db] ERROR: No sequences passed filters", file=sys.stderr)
        sys.exit(1)

    # ── Write metadata ────────────────────────────────────────────────
    write_meta_tsv(meta, out_meta)

    # ── Build Diamond DB ──────────────────────────────────────────────
    build_diamond(out_faa, out_dmnd.replace('.dmnd',''), args.threads)

    print(f"\n[prepare_diamond_db] Done!")
    print(f"  Diamond DB : {out_dmnd}")
    print(f"  Metadata   : {out_meta}")
    print(f"\nAdd to Snakefile config:")
    print(f"  CUSTOM_VIRAL_DMND = \"{out_dmnd}\"")
    print(f"  CUSTOM_VIRAL_META = \"{out_meta}\"")
    print(f"  # or for prokaryotes:")
    print(f"  CUSTOM_PROK_DMND  = \"{out_dmnd}\"")
    print(f"  CUSTOM_PROK_META  = \"{out_meta}\"")


if __name__ == '__main__':
    main()
