# ══════════════════════════════════════════════════════════════════════
# rules/merge_dedup.smk — BLOCK 3: Merge, Filter and Deduplication
#
# merge_contigs : combines MEGAHIT + metaSPAdes + metaviralSPAdes,
#                 adds MEGAHIT_/SPADES_/METAVIRAL_ prefix, filters < MIN_CONTIG bp
# mmseqs2       : deduplication at 95% identity → {sample}_rep_seq.fasta
#                 (used as input by ALL downstream rules)
#
# For long reads the merged FASTA comes from rules.merge_lr.output.merged
# and is passed directly to mmseqs2 (merge_contigs is SR-only).
# ══════════════════════════════════════════════════════════════════════


rule merge_contigs:
    """
    Merge MEGAHIT + metaSPAdes + metaviralSPAdes contigs (short reads only).
    - Adds MEGAHIT_ / SPADES_ / METAVIRAL_ prefix to each header to prevent collisions
    - Filters contigs shorter than MIN_CONTIG bp
    - metaviralSPAdes may produce an empty file (no viral contigs) — handled gracefully
    """
    input:
        megahit   = rules.megahit.output.contigs,
        spades    = rules.metaspades.output.contigs,
        metaviral = rules.metaviral_spades.output.contigs,
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
                (input.spades,    "SPADES"),
                (input.metaviral, "METAVIRAL"),
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
      Short reads → merge_contigs (MEGAHIT + metaSPAdes + metaviralSPAdes)
      Long reads  → merge_lr     (Flye + hifiasm)
    """
    input:
        contigs = (
            f"{OUTDIR}/{{sample}}/assembly/lr/merged/{{sample}}_lr_merged.fasta"
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
    conda: "envs/env_assembly.yaml"
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
