# ══════════════════════════════════════════════════════════════════════
# rules/mapping.smk — BLOCK 6: Read Mapping
#
# Short reads : bwa_index → bwa_mem → samtools_sort → calc_depth
# Long reads  : minimap2_lr → samtools_sort_lr → calc_depth
#
# Each rule uses exactly one container (bwa-mem2 / minimap2 / samtools).
# Piped bwa-mem2|samtools was split to avoid needing a combined image.
# Intermediate SAM files are marked temp() and auto-deleted by Snakemake.
#
# calc_depth output (jgi_summarize_bam_contig_depths) is consumed by
# MetaBAT2, SemiBin2, vRhyme.
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
    Short read alignment with BWA-MEM2. Outputs an unsorted SAM (temp).
    Sorted by the downstream samtools_sort rule.
    SE mode: maps R1 only. PE mode: maps R1 + R2.
    """
    input:
        tr1 = _clean_r1,
        tr2 = _clean_r2,
        idx = rules.bwa_index.output.idx,
    output:
        sam = temp(f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sam"),
    log:
        f"{OUTDIR}/{{sample}}/logs/bwa_mem.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/bwa_mem.tsv"
    conda:      "../envs/env_mapping.yaml"
    container:  CONTAINERS.get("bwa_mem2")
    threads: THREADS
    params:
        prefix     = f"{OUTDIR}/{{sample}}/mapping/contigs_index",
        single_end = SINGLE_END,
    shell:
        """
        mkdir -p $(dirname {output.sam})
        if [ "{params.single_end}" = "True" ]; then
            bwa-mem2 mem \
                -t {threads} \
                {params.prefix} \
                {input.tr1} \
                > {output.sam} 2> {log}
        else
            bwa-mem2 mem \
                -t {threads} \
                {params.prefix} \
                {input.tr1} {input.tr2} \
                > {output.sam} 2> {log}
        fi
        """


rule samtools_sort:
    """
    Sort BWA-MEM2 SAM → BAM, index, and compute flagstat.
    Receives the temp SAM from bwa_mem; deletes it on completion.
    """
    input:
        sam = rules.bwa_mem.output.sam,
    output:
        bam   = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bai   = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam.bai",
        stats = f"{OUTDIR}/{{sample}}/mapping/flagstat.txt",
        done  = f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/samtools_sort.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/samtools_sort.tsv"
    conda:     "../envs/env_mapping.yaml"
    container: CONTAINERS.get("samtools")
    threads: THREADS
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
        samtools index {output.bam} 2>> {log}
        samtools flagstat {output.bam} > {output.stats} 2>> {log}
        touch {output.done}
        """


if LONG_READS:

    rule minimap2_lr:
        """
        Long read alignment with minimap2. Outputs an unsorted SAM (temp).
        Sorted by the downstream samtools_sort_lr rule.
        ONT: -ax map-ont; HiFi: -ax map-hifi.
        """
        input:
            contigs = rules.mmseqs2.output.rep,
            reads   = _clean_lr,
        output:
            sam = temp(f"{OUTDIR}/{{sample}}/mapping/{{sample}}_lr.sam"),
        log:   f"{OUTDIR}/{{sample}}/logs/minimap2_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/minimap2_lr.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("minimap2")
        threads: THREADS
        shell:
            """
            mkdir -p $(dirname {output.sam})
            if [ "{LR_TECH}" = "hifi" ]; then
                PRESET="map-hifi"
            else
                PRESET="map-ont"
            fi
            minimap2 -ax $PRESET \
                -t {threads} \
                {input.contigs} {input.reads} \
                > {output.sam} 2> {log}
            """

    rule samtools_sort_lr:
        """
        Sort minimap2 SAM → BAM, index, and compute flagstat (long reads).
        Receives the temp SAM from minimap2_lr.
        """
        input:
            sam = rules.minimap2_lr.output.sam,
        output:
            bam   = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
            bai   = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam.bai",
            stats = f"{OUTDIR}/{{sample}}/mapping/flagstat.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/samtools_sort_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/samtools_sort_lr.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("samtools")
        threads: THREADS
        shell:
            """
            samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
            samtools index {output.bam} 2>> {log}
            samtools flagstat {output.bam} > {output.stats} 2>> {log}
            """


rule calc_depth:
    """
    Per-contig coverage depth from BAM (jgi_summarize_bam_contig_depths).
    Output columns: contigName, contigLen, totalAvgDepth, bam_avg, bam_var.
    Used by: MetaBAT2, SemiBin2, vRhyme.
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
