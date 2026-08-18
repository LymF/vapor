# ══════════════════════════════════════════════════════════════════════
# rules/quast.smk — BLOCK 4: Assembly Quality Assessment
#
# Short reads: avalia 3 conjuntos (MEGAHIT, merged+filtrado, dedup)
# Long reads : avalia apenas o conjunto final deduplicated
# ══════════════════════════════════════════════════════════════════════


rule quast:
    """
    QUAST assembly quality assessment.
    Short reads: 3 contig sets (MEGAHIT, merged+filtered, deduplicated).
    Long reads : only the deduplicated final contigs (no MEGAHIT to compare).
    Key metrics: N50, L50, total length, # contigs, GC content.
    """
    input:
        megahit   = (rules.megahit.output.contigs           if not LONG_READS else []),
        merged    = (rules.merge_contigs.output.merged       if not LONG_READS else []),
        rep       = rules.mmseqs2.output.rep,
    output:
        report = f"{OUTDIR}/{{sample}}/quast/report.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/quast.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/quast.tsv"
    conda:      "../envs/env_qc.yaml"
    container:  CONTAINERS.get("quast")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/quast"
    shell:
        """
        if [ "{LONG_READS}" = "True" ]; then
            quast.py \
                {input.rep} \
                --labels "lr_deduplicated" \
                --min-contig {MIN_CONTIG} \
                --threads {threads} \
                -o {params.outdir} \
                > {log} 2>&1
        else
            quast.py \
                {input.megahit} \
                {input.merged} \
                {input.rep} \
                --labels "MEGAHIT,merged_filtered,deduplicated" \
                --min-contig {MIN_CONTIG} \
                --threads {threads} \
                -o {params.outdir} \
                > {log} 2>&1
        fi
        """
