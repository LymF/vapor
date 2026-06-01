# ══════════════════════════════════════════════════════════════════════
# rules/host_removal.smk — BLOCK 1.5: Host Genome Removal (optional)
#
# Only active when USE_HOST_REMOVAL=True (host_genome set in config).
# SR: bwa-mem2 → keep pairs where BOTH mates are unmapped (samtools -f 12 -F 256)
# LR: minimap2 → keep reads that are unmapped (samtools -f 4)
#
# Both rules reuse env_mapping (bwa-mem2, minimap2, samtools).
# ══════════════════════════════════════════════════════════════════════

if USE_HOST_REMOVAL:

    if not HOST_INDEX:

        rule build_host_index:
            """
            Build bwa-mem2 index for the host reference genome.
            Runs once and is shared across all samples.
            Skipped when host_index is pre-set in config.
            """
            input:
                genome = HOST_GENOME,
            output:
                idx = f"{OUTDIR}/host_index/host.bwt.2bit.64",
            log:
                f"{OUTDIR}/host_index/bwa_index.log"
            benchmark:
                f"{OUTDIR}/host_index/bwa_index.tsv"
            conda:      "../envs/env_mapping.yaml"
            container:  CONTAINERS.get("bwa_mem2")
            threads: min(THREADS, 8)
            params:
                prefix = f"{OUTDIR}/host_index/host",
            shell:
                """
                mkdir -p {OUTDIR}/host_index
                bwa-mem2 index -p {params.prefix} {input.genome} > {log} 2>&1
                """

    rule remove_host_sr:
        """
        Remove host contamination from short reads.
        PE: maps R1+R2; keeps pairs where BOTH mates are unmapped (-f 12 -F 256).
        SE: maps R1 only; keeps reads that are unmapped (-f 4).
        clean_r2 is an empty sentinel file in SE mode.
        """
        input:
            tr1 = f"{OUTDIR}/{{sample}}/trimmed/{{sample}}_R1_fastp.fq.gz",
            tr2 = (f"{OUTDIR}/{{sample}}/trimmed/{{sample}}_R2_fastp.fq.gz"
                   if not SINGLE_END else []),
            idx = (f"{HOST_INDEX}.bwt.2bit.64"
                   if HOST_INDEX
                   else f"{OUTDIR}/host_index/host.bwt.2bit.64"),
        output:
            clean_r1 = f"{OUTDIR}/{{sample}}/host_removed/{{sample}}_R1_clean.fq.gz",
            clean_r2 = f"{OUTDIR}/{{sample}}/host_removed/{{sample}}_R2_clean.fq.gz",
            stats    = f"{OUTDIR}/{{sample}}/host_removed/host_stats.txt",
        log:
            f"{OUTDIR}/{{sample}}/logs/remove_host_sr.log"
        benchmark:
            f"{OUTDIR}/{{sample}}/benchmarks/remove_host_sr.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("bwa_mem2")
        threads: THREADS
        params:
            idx_prefix = (HOST_INDEX if HOST_INDEX
                          else f"{OUTDIR}/host_index/host"),
            tmpbam     = f"{OUTDIR}/{{sample}}/host_removed/tmp_namesort.bam",
            single_end = SINGLE_END,
        shell:
            """
            mkdir -p $(dirname {output.clean_r1})
            if [ "{params.single_end}" = "True" ]; then
                bwa-mem2 mem -t {threads} {params.idx_prefix} \
                    {input.tr1} 2>> {log} \
                | samtools sort -@ {threads} -o {params.tmpbam} - 2>> {log}
                samtools flagstat {params.tmpbam} > {output.stats} 2>> {log}
                samtools view -@ {threads} -f 4 -b {params.tmpbam} \
                | samtools fastq - 2>> {log} \
                | gzip > {output.clean_r1}
                touch {output.clean_r2}
            else
                bwa-mem2 mem -t {threads} {params.idx_prefix} \
                    {input.tr1} {input.tr2} 2>> {log} \
                | samtools sort -@ {threads} -n \
                    -o {params.tmpbam} - 2>> {log}
                samtools flagstat {params.tmpbam} > {output.stats} 2>> {log}
                samtools view -@ {threads} -f 12 -F 256 -b {params.tmpbam} \
                | samtools fastq \
                    -1 {output.clean_r1} \
                    -2 {output.clean_r2} \
                    -0 /dev/null -s /dev/null \
                    2>> {log}
            fi
            rm -f {params.tmpbam}
            """

    if LONG_READS:

        rule remove_host_lr:
            """
            Remove host contamination from long reads.
            Maps reads to the host reference with minimap2;
            keeps only reads that are unmapped (samtools -f 4).
            Outputs host-removed FASTQ + flagstat for reporting.
            """
            input:
                reads  = f"{OUTDIR}/{{sample}}/lr_filtered/{{sample}}_filtered.fastq.gz",
                genome = HOST_GENOME,
            output:
                clean = f"{OUTDIR}/{{sample}}/host_removed/{{sample}}_lr_clean.fastq.gz",
                stats = f"{OUTDIR}/{{sample}}/host_removed/host_lr_stats.txt",
            log:
                f"{OUTDIR}/{{sample}}/logs/remove_host_lr.log"
            benchmark:
                f"{OUTDIR}/{{sample}}/benchmarks/remove_host_lr.tsv"
            conda:      "../envs/env_mapping.yaml"
            container:  CONTAINERS.get("minimap2")
            threads: THREADS
            params:
                tmpbam = f"{OUTDIR}/{{sample}}/host_removed/tmp_lr_host.bam",
            shell:
                """
                mkdir -p $(dirname {output.clean})
                if [ "{LR_TECH}" = "hifi" ]; then
                    PRESET="map-hifi"
                else
                    PRESET="map-ont"
                fi
                minimap2 -ax $PRESET -t {threads} {input.genome} {input.reads} \
                    2>> {log} \
                | samtools sort -@ {threads} -o {params.tmpbam} - 2>> {log}
                samtools flagstat {params.tmpbam} > {output.stats} 2>> {log}
                samtools view -@ {threads} -f 4 -b {params.tmpbam} \
                | samtools fastq - 2>> {log} \
                | gzip > {output.clean}
                rm -f {params.tmpbam}
                """
