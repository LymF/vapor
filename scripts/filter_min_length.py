#!/usr/bin/env python3
"""
filter_min_length.py — descarta contigs abaixo de um comprimento minimo.

Existe por causa de uma assimetria entre os montadores: o MEGAHIT
(`--min-contig-len`) e o metaMDBG (`--min-contig-length`) aplicam o
`min_contig` do config sozinhos, mas o **Flye nao tem flag equivalente** --
entao, ate 2026-08-19, uma corrida ONT levava contigs de qualquer tamanho
para deteccao viral, mapeamento e binning, enquanto short-read e HiFi nao
levavam. O `min_contig` do config, documentado como "comprimento minimo
pos-montagem", simplesmente nao valia para uma das tres trilhas.

Uso:
    python3 filter_min_length.py <entrada.fasta> <saida.fasta> <min_bp>

Escreve no stderr quantos contigs entraram, sairam e quantas bases foram
descartadas -- um filtro silencioso e indistinguivel de uma montagem ruim.
"""
import sys


def iter_fasta(handle):
    """(header_sem_'>', [linhas_de_sequencia]) por registro."""
    header, seq = None, []
    for line in handle:
        if line.startswith(">"):
            if header is not None:
                yield header, seq
            header, seq = line[1:].rstrip("\n"), []
        elif header is not None:
            seq.append(line.rstrip("\n"))
    if header is not None:
        yield header, seq


def filter_fasta(fin, fout, min_len):
    """Devolve (n_entrada, n_saida, bases_descartadas)."""
    n_in = n_out = dropped = 0
    for header, seq in iter_fasta(fin):
        length = sum(len(s) for s in seq)
        n_in += 1
        if length < min_len:
            dropped += length
            continue
        n_out += 1
        fout.write(">" + header + "\n")
        for chunk in seq:
            fout.write(chunk + "\n")
    return n_in, n_out, dropped


def main():
    if len(sys.argv) != 4:
        sys.exit("Uso: filter_min_length.py <entrada.fasta> <saida.fasta> <min_bp>")
    src, dst, min_len = sys.argv[1], sys.argv[2], int(sys.argv[3])
    with open(src) as fin, open(dst, "w") as fout:
        n_in, n_out, dropped = filter_fasta(fin, fout, min_len)
    sys.stderr.write(
        "[filter_min_length] %d contigs -> %d (>= %d bp); %d descartados, "
        "%d bases\n" % (n_in, n_out, min_len, n_in - n_out, dropped))


if __name__ == "__main__":
    main()
