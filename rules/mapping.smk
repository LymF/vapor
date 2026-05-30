# ══════════════════════════════════════════════════════════════════════
# rules/mapping.smk — BLOCK 6: Read Mapping
#
# Short reads : bwa_index → bwa_mem → calc_depth
# Long reads  : minimap2_lr → calc_depth
#
# calc_depth output (jgi_summarize_bam_contig_depths) is consumed by
# MetaBAT2, VAMB, SemiBin2, vRhyme.
# ══════════════════════════════════════════════════════════════════════


rule bwa_index:
    """Index deduplicated contigs for BWA-MEM2."""
    input:
        contigs = rules.mmseqs2.output.rep,
    output:
        idx = f"{OUTDIR}/{{sample}}/mapping/contigs_index.bwt.2bit.64",
    log:
        f"{OUTDIR}/{{sample}}/logs/bwa_index.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/bwa_index.tsv"
    conda:      "../envs/env_mapping.yaml"
    container:  CONTAINERS.get("bwa_mem2")
    params:
        prefix = f"{OUTDIR}/{{sample}}/mapping/contigs_index"
    shell:
        """
        mkdir -p {OUTDIR}/{wildcards.sample}/mapping
        bwa-mem2 index -p {params.prefix} {input.contigs} > {log} 2>&1
        """


rule bwa_mem:
    """
    Short read mapping with BWA-MEM2. Skipped when LONG_READS=True
    (minimap2_lr takes precedence via ruleorder in Snakefile).
    Outputs done.txt sentinel to decouple from minimap2_lr BAM output.
    """
    input:
        tr1 = rules.fastp.output.tr1,
        tr2 = rules.fastp.output.tr2,
        idx = rules.bwa_index.output.idx,
    output:
        bam   = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bai   = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam.bai",
        stats = f"{OUTDIR}/{{sample}}/mapping/flagstat.txt",
        done  = f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/bwa_mem.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/bwa_mem.tsv"
    conda:      "../envs/env_mapping.yaml"
    container:  CONTAINERS.get("bwa_mem2")
    threads: THREADS
    params:
        prefix = f"{OUTDIR}/{{sample}}/mapping/contigs_index",
    shell:
        """
        if [ "{LONG_READS}" = "True" ]; then
            mkdir -p $(dirname {output.bam})
            touch {output.bam} {output.bai} {output.stats} {output.done}; exit 0
        fi
        mkdir -p $(dirname {output.bam})
        bwa-mem2 mem \
            -t {threads} \
            {params.prefix} \
            {input.tr1} {input.tr2} \
            2> {log} \
        | samtools sort -@ {threads} -o {output.bam}
        samtools index {output.bam}
        samtools flagstat {output.bam} > {output.stats} 2>> {log}
        touch {output.done}
        """


if LONG_READS:

    rule minimap2_lr:
        """
        Long read mapping with minimap2.
        Replaces BWA-MEM2 for long reads (ruleorder in Snakefile).
        ONT: -ax map-ont; HiFi: -ax map-hifi.
        Output: sorted BAM identical in format to BWA-MEM2 output.
        """
        input:
            contigs = rules.mmseqs2.output.rep,
            reads   = f"{OUTDIR}/{{sample}}/lr_filtered/{{sample}}_filtered.fastq.gz",
        output:
            bam   = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
            bai   = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam.bai",
            stats = f"{OUTDIR}/{{sample}}/mapping/flagstat.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/minimap2_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/minimap2_lr.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("minimap2")
        threads: THREADS
        shell:
            """
            mkdir -p $(dirname {output.bam})
            if [ "{LR_TECH}" = "hifi" ]; then
                PRESET="map-hifi"
            else
                PRESET="map-ont"
            fi
            minimap2 -ax $PRESET \
                -t {threads} \
                {input.contigs} {input.reads} | \
            samtools sort -@ {threads} -o {output.bam} -
            samtools index {output.bam}
            samtools flagstat {output.bam} > {output.stats} 2>> {log}
            """


rule calc_depth:
    """
    Per-contig coverage depth from BAM (jgi_summarize_bam_contig_depths).
    Output columns: contigName, contigLen, totalAvgDepth, bam_avg, bam_var.
    Used by: MetaBAT2, VAMB, SemiBin2, vRhyme.
    """
    input:
        bam      = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bwa_done = (f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt"
                    if not LONG_READS else []),
    output:
        depth = f"{OUTDIR}/{{sample}}/mapping/{{sample}}_depth.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/calc_depth.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/calc_depth.tsv"
    conda:      "../envs/env_mapping.yaml"
    container:  CONTAINERS.get("metabat2")
    shell:
        """
        jgi_summarize_bam_contig_depths \
            --outputDepth {output.depth} \
            {input.bam} \
            > {log} 2>&1
        """
