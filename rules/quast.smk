# ══════════════════════════════════════════════════════════════════════
# rules/quast.smk — BLOCK 4: Assembly Quality Assessment
#
# Since (d): a single contig set per track (no merge, no dedup), so QUAST
# evaluates just the assembly itself -- short reads (MEGAHIT) and long
# reads (Flye+Medaka for ONT, metaMDBG for HiFi) alike.
# ══════════════════════════════════════════════════════════════════════


rule quast:
    """
    QUAST assembly quality assessment on the sample's assembly (the hub
    since item (d) of docs/ROADMAP_SIMPLIFICACAO.md).
    Key metrics: N50, L50, total length, # contigs, GC content.
    """
    input:
        assembly  = _sample_contigs,
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
        quast.py \
            {input.assembly} \
            --labels "assembly" \
            --min-contig {MIN_CONTIG} \
            --threads {threads} \
            -o {params.outdir} \
            > {log} 2>&1
        """
