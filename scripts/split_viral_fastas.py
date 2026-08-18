#!/usr/bin/env python3
"""
split_viral_fastas.py — Split viral consensus FASTA into per-genome FASTA files.
Bins-first strategy: vRhyme bins are copied as-is (each bin = one viral genome).
Individual FASTA files are written ONLY for contigs NOT already in a vRhyme bin,
eliminating the kmer-db redundancy that caused duplicate PHIST rows.

NAMESPACE: <viral_fasta> is the GLOBAL vOTU catalog, whose headers are
namespaced "{source_id}|{contig_id}". The vRhyme bins are built from the
per-source viral set, so their headers are BARE "{contig_id}". Comparing the
two directly makes the "already binned" test never fire, every binned contig
gets written a second time as its own FASTA, and PHIST emits duplicate rows
for the same sequence -- exactly the redundancy this script exists to remove
(downstream, build_host_collapse then double-counts n_viruses and total_rpkm).
That is why <source_id> is required: bare ids are only unique WITHIN a source,
so the bin names are lifted into the namespaced form instead of the catalog
names being stripped down -- stripping would let a contig from another source
collide with this source's bins.

Usage: python3 split_viral_fastas.py <viral_fasta> <vrhyme_dir> <output_dir> <source_id>
"""
import os, sys, shutil, glob

viral_fa   = sys.argv[1]
vrhyme_dir = sys.argv[2]
out_dir    = sys.argv[3]
source_id  = sys.argv[4] if len(sys.argv) > 4 else ''

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
                name = line[1:].split()[0]
                binned_contigs.add(name)
                if source_id:
                    binned_contigs.add(f'{source_id}|{name}')
print(f'[split_viral_fastas] vRhyme bins copied: {bins_copied} '
      f'({len(binned_contigs)} binned contig keys, source_id={source_id or "<none>"})')

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
# As duas pontas do join, explicitas: "0 pulados" com bins presentes e o
# sintoma de descasamento de namespace, nao de um conjunto viral sem bins.
if bins_copied and skipped == 0:
    print('[split_viral_fastas] WARNING: %d vRhyme bins carregados mas NENHUM '
          'contig casou como ja binado -- suspeitar de descasamento de '
          'namespace entre os headers dos bins e os do FASTA viral' % bins_copied)
print(f'[split_viral_fastas] Total fastas in output dir: {total}')
