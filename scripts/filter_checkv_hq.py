#!/usr/bin/env python3
"""
filter_checkv_hq.py — Filter viral FASTA for vConTACT3 input.

Strategy (two tiers to maximise coverage):
  Tier 1 — CheckV MQ+: Complete / High-quality / Medium-quality OR completeness >= 50%
  Tier 2 — Any contig with a length >= MIN_LEN (default 5000 bp) regardless of CheckV
            These short/fragmented contigs still contribute to protein-sharing network

Usage: python3 filter_checkv_hq.py <checkv_tsv> <viral_fasta> <output_fasta> [min_len]
"""
import csv, sys

checkv_p = sys.argv[1]
viral_p  = sys.argv[2]
out_p    = sys.argv[3]
min_len  = int(sys.argv[4]) if len(sys.argv) > 4 else 5000

# Tier 1: CheckV quality filter
keep_quality = set()
with open(checkv_p) as f:
    for row in csv.DictReader(f, delimiter='\t'):
        q = row.get('checkv_quality', '')
        try:
            comp = float(row.get('completeness', '0') or 0)
        except:
            comp = 0.0
        if q in ('Complete', 'High-quality', 'Medium-quality') or comp >= 50:
            keep_quality.add(row.get('contig_id', ''))

# Tier 2: length-based (long enough to have informative protein content)
written = 0; tier1 = 0; tier2 = 0
with open(viral_p) as fin, open(out_p, 'w') as fout:
    hdr = None; seq = []
    def flush(hdr, seq):
        global written, tier1, tier2
        if not hdr: return
        length = sum(len(s.strip()) for s in seq)
        in_t1 = hdr in keep_quality
        in_t2 = length >= min_len
        if in_t1 or in_t2:
            fout.write('>' + hdr + '\n' + ''.join(seq))
            written += 1
            if in_t1: tier1 += 1
            elif in_t2: tier2 += 1

    for line in fin:
        if line.startswith('>'):
            flush(hdr, seq)
            hdr = line[1:].strip().split()[0]; seq = []
        else:
            seq.append(line)
    flush(hdr, seq)

print(f'[filter_checkv_hq] Written: {written} genomes '
      f'(tier1 CheckV MQ+: {tier1}, tier2 len>={min_len}bp: {tier2})')
print(f'[filter_checkv_hq] CheckV keep list: {len(keep_quality)} genomes')
