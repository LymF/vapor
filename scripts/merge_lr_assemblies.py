#!/usr/bin/env python3
"""
merge_lr_assemblies.py — Merge long-read assembly FASTAs.
Labels contigs with assembler prefix for traceability.

Usage (3 assemblers):
  python3 merge_lr_assemblies.py <flye.fa> <hifiasm.fa> <metaMDBG.fa> <output.fa> <sample>

Usage (2 assemblers — legacy):
  python3 merge_lr_assemblies.py <flye.fa> <hifiasm.fa> <output.fa> <sample>
"""
import sys, os

args = sys.argv[1:]
if len(args) == 5:
    flye_fa, hifiasm_fa, mdbg_fa, out_fa, sample = args
elif len(args) == 4:
    flye_fa, hifiasm_fa, out_fa, sample = args
    mdbg_fa = ""
else:
    print("ERROR: expected 4 or 5 positional arguments", file=sys.stderr)
    sys.exit(1)

sources = [
    (flye_fa,    "FLYE"),
    (hifiasm_fa, "HIFIASM"),
    (mdbg_fa,    "MDBG"),
]

written = 0
counts = {}
with open(out_fa, 'w') as fout:
    for fa, prefix in sources:
        if not fa or not os.path.exists(fa) or not os.path.getsize(fa):
            counts[prefix] = 0
            continue
        n = 0
        with open(fa) as fin:
            for line in fin:
                if line.startswith('>'):
                    name = line[1:].split()[0]
                    fout.write(f'>{prefix}_{sample}_{name}\n')
                    n += 1
                    written += 1
                else:
                    fout.write(line)
        counts[prefix] = n

summary = " + ".join(f"{counts[p]} {p}" for _, p in sources if counts.get(p, 0) > 0)
print(f'[merge_lr] {written} total contigs merged: {summary}')
