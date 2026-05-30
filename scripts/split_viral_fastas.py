#!/usr/bin/env python3
"""
split_viral_fastas.py — Split viral consensus FASTA into per-genome FASTA files.
Bins-first strategy: vRhyme bins are copied as-is (each bin = one viral genome).
Individual FASTA files are written ONLY for contigs NOT already in a vRhyme bin,
eliminating the kmer-db redundancy that caused duplicate PHIST rows.

Usage: python3 split_viral_fastas.py <viral_fasta> <vrhyme_dir> <output_dir>
"""
import os, sys, shutil, glob

viral_fa   = sys.argv[1]
vrhyme_dir = sys.argv[2]
out_dir    = sys.argv[3]

os.makedirs(out_dir, exist_ok=True)

# 1. vRhyme bins — copy bin FASTAs and record which contigs are already binned
binned_contigs = set()
vbins = glob.glob(os.path.join(vrhyme_dir, 'vRhyme_best_bins.*.fasta'))
bins_copied = 0
for vbin in vbins:
    bname = os.path.basename(vbin)
    shutil.copy(vbin, os.path.join(out_dir, bname))
    bins_copied += 1
    with open(vbin) as f:
        for line in f:
            if line.startswith('>'):
                binned_contigs.add(line[1:].split()[0])
print(f'[split_viral_fastas] vRhyme bins copied: {bins_copied} ({len(binned_contigs)} binned contigs)')

# 2. Individual FASTAs for unbinned contigs only (bins-first, no redundancy)
contigs_written = 0
skipped = 0
if os.path.exists(viral_fa):
    seq_name = None
    seq_lines = []
    with open(viral_fa) as f:
        for line in f:
            if line.startswith('>'):
                if seq_name and seq_lines:
                    if seq_name not in binned_contigs:
                        safe = seq_name.replace('/', '_').replace(' ', '_')
                        out = os.path.join(out_dir, 'contig_' + safe + '.fasta')
                        with open(out, 'w') as fout:
                            fout.write('>' + seq_name + '\n' + ''.join(seq_lines))
                        contigs_written += 1
                    else:
                        skipped += 1
                seq_name = line[1:].split()[0]
                seq_lines = []
            else:
                seq_lines.append(line)
        if seq_name and seq_lines:
            if seq_name not in binned_contigs:
                safe = seq_name.replace('/', '_').replace(' ', '_')
                out = os.path.join(out_dir, 'contig_' + safe + '.fasta')
                with open(out, 'w') as fout:
                    fout.write('>' + seq_name + '\n' + ''.join(seq_lines))
                contigs_written += 1
            else:
                skipped += 1

total = len(os.listdir(out_dir))
print(f'[split_viral_fastas] Unbinned contigs written: {contigs_written} (skipped {skipped} already binned)')
print(f'[split_viral_fastas] Total fastas in output dir: {total}')
