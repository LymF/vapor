# ══════════════════════════════════════════════════════════════════════
# rules/merge_dedup.smk — BLOCK 3: Merge, Filter and Deduplication
#
# merge_contigs : takes the MEGAHIT contigs, adds the MEGAHIT_ prefix and
#                 filters < MIN_CONTIG bp
# mmseqs2       : deduplication at 95% identity → {sample}_rep_seq.fasta
#                 (used as input by ALL downstream rules)
#
# For long reads the single assembler's output (Flye+Medaka for ONT,
# metaMDBG for HiFi) is passed directly to mmseqs2 (merge_contigs is SR-only).
# ══════════════════════════════════════════════════════════════════════


rule merge_contigs:
    """
    Prefix and length-filter the MEGAHIT contigs (short reads only).
    - Adds the MEGAHIT_ prefix to each header (kept so downstream contig IDs
      stay stable now that metaSPAdes/metaviralSPAdes are gone)
    - Filters contigs shorter than MIN_CONTIG bp
    """
    input:
        megahit   = rules.megahit.output.contigs,
    output:
        merged = f"{OUTDIR}/{{sample}}/mmseqs/{{sample}}_merged_filtered.fasta",
    log:
        f"{OUTDIR}/{{sample}}/logs/merge_contigs.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/merge_contigs.tsv"
    run:
        def iter_fasta(path, prefix):
            with open(path) as f:
                header, seq = "", []
                for line in f:
                    line = line.rstrip()
                    if line.startswith(">"):
                        if header:
                            yield header, "".join(seq)
                        header = line.replace(">", f">{prefix}_", 1)
                        seq = []
                    else:
                        seq.append(line)
                if header:
                    yield header, "".join(seq)

        kept = total = 0
        with open(output.merged, "w") as out, open(log[0], "w") as lf:
            sources = [
                (input.megahit,   "MEGAHIT"),
            ]
            for src, pfx in sources:
                for header, seq in iter_fasta(src, pfx):
                    total += 1
                    if len(seq) >= MIN_CONTIG:
                        out.write(header + "\n" + seq + "\n")
                        kept += 1
            msg = f"Kept {kept}/{total} contigs >= {MIN_CONTIG} bp\n"
            lf.write(msg)
            print(msg)


rule mmseqs2:
    """
    Deduplicate merged contigs with MMseqs2 easy-linclust.
    Output: {sample}_rep_seq.fasta — non-redundant contig set used by ALL
    downstream rules (viral detection, mapping, binning, taxonomy).
    Input:
      Short reads       → merge_contigs (MEGAHIT)
      Long reads (ONT)  → medaka_lr     (Flye, polished)
      Long reads (HiFi) → metaMDBG_lr   (metaMDBG, no polishing)
    """
    input:
        contigs = (
            (
                f"{OUTDIR}/{{sample}}/assembly/lr/flye_polished/assembly.fasta"
                if LR_TECH == "ont"
                else f"{OUTDIR}/{{sample}}/assembly/lr/metaMDBG/assembly.fasta"
            )
            if LONG_READS
            else rules.merge_contigs.output.merged
        ),
    output:
        rep     = f"{OUTDIR}/{{sample}}/mmseqs/{{sample}}_rep_seq.fasta",
        cluster = f"{OUTDIR}/{{sample}}/mmseqs/{{sample}}_cluster.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/mmseqs2.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/mmseqs2.tsv"
    conda:      "../envs/env_assembly.yaml"
    container:  CONTAINERS.get("mmseqs2")
    threads: THREADS
    params:
        prefix = f"{OUTDIR}/{{sample}}/mmseqs/{{sample}}",
        tmp    = f"{OUTDIR}/{{sample}}/mmseqs/tmp",
    shell:
        """
        mkdir -p {params.tmp}
        mmseqs easy-linclust \
            {input.contigs} \
            {params.prefix} \
            {params.tmp} \
            --min-seq-id {MIN_SEQ_ID} \
            -c 0.9 --cov-mode 2 \
            --threads {threads} \
            > {log} 2>&1
        """
