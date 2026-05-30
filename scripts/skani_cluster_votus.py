#!/usr/bin/env python3
"""
Greedy single-linkage vOTU clustering from a skani triangle TSV.

Usage:
    skani_cluster_votus.py <skani_triangle.tsv> <viral_fasta.fa> <ani_min> <af_min> <out_clusters.tsv>

Two genomes are placed in the same cluster when both
    ani  >= ani_min   AND   max(af_q, af_r) >= af_min
This is the ICTV / Roux 2019 vOTU definition (defaults 95.0 / 85.0).

The cluster representative is the first member encountered in FASTA order
(the upstream pipeline writes vRhyme bins first, then unbinned contigs;
within each block FASTA order roughly follows length).
"""
import sys


def main():
    if len(sys.argv) != 6:
        sys.stderr.write(__doc__)
        sys.exit(2)
    ani_tsv, fasta, ani_min_s, af_min_s, out_tsv = sys.argv[1:6]
    ani_min = float(ani_min_s)
    af_min = float(af_min_s)

    ids = []
    with open(fasta) as f:
        for line in f:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])

    neigh = {i: set() for i in ids}
    try:
        with open(ani_tsv) as f:
            f.readline()  # header
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                q, r = parts[0], parts[1]
                try:
                    ani = float(parts[2])
                    afq = float(parts[3])
                    afr = float(parts[4])
                except ValueError:
                    continue
                if ani >= ani_min and max(afq, afr) >= af_min:
                    if q in neigh and r in neigh:
                        neigh[q].add(r)
                        neigh[r].add(q)
    except FileNotFoundError:
        pass

    seen = set()
    clusters = []
    for n in ids:
        if n in seen:
            continue
        comp = []
        stack = [n]
        while stack:
            x = stack.pop()
            if x in seen:
                continue
            seen.add(x)
            comp.append(x)
            for y in neigh.get(x, ()):
                if y not in seen:
                    stack.append(y)
        clusters.append(comp)

    with open(out_tsv, "w") as f:
        f.write("representative\tmember\n")
        for comp in clusters:
            rep = comp[0]
            for m in comp:
                f.write(f"{rep}\t{m}\n")

    sys.stderr.write(
        f"[skani_cluster_votus] genomes={len(ids)} clusters={len(clusters)} "
        f"ani>={ani_min} af>={af_min}\n"
    )


if __name__ == "__main__":
    main()
