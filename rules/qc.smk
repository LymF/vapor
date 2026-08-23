# ══════════════════════════════════════════════════════════════════════
# rules/qc.smk — BLOCK 1: Quality Control
#
# Short reads : fastp (adapter trimming + QC report in one step)
# Long reads  : nanoplot_lr → porechop_lr → filtlong_lr
#
# All variables (OUTDIR, THREADS, LONG_READS, …) come from the main
# Snakefile namespace — no imports needed.
# ══════════════════════════════════════════════════════════════════════


# ── Short read QC ─────────────────────────────────────────────────────

rule fastp:
    """
    Adapter trimming + quality filtering + QC report with fastp.
    PE mode: paired -i/-I -o/-O with --detect_adapter_for_pe.
    SE mode: single -i -o only; tr2 is an empty sentinel file.
    """
    input:
        r1 = lambda wc: SAMPLES[wc.sample]["R1"],
        r2 = lambda wc: SAMPLES[wc.sample].get("R2", []),
    output:
        tr1  = f"{OUTDIR}/{{sample}}/trimmed/{{sample}}_R1_fastp.fq.gz",
        tr2  = f"{OUTDIR}/{{sample}}/trimmed/{{sample}}_R2_fastp.fq.gz",
        json = f"{OUTDIR}/{{sample}}/qc_raw/{{sample}}_fastp.json",
        html = f"{OUTDIR}/{{sample}}/qc_raw/{{sample}}_fastp.html",
        done = f"{OUTDIR}/{{sample}}/qc_raw/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/fastp.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/fastp.tsv"
    conda:      "../envs/env_qc.yaml"
    container:  CONTAINERS.get("fastp")
    threads: min(THREADS, 16)
    params:
        single_end = SINGLE_END,
    shell:
        """
        mkdir -p {OUTDIR}/{wildcards.sample}/trimmed \
                 {OUTDIR}/{wildcards.sample}/qc_raw
        if [ "{params.single_end}" = "True" ]; then
            fastp \
                -i {input.r1} \
                -o {output.tr1} \
                --json {output.json} \
                --html {output.html} \
                --thread {threads} \
                --qualified_quality_phred 20 \
                --length_required 50 \
                --low_complexity_filter \
                --complexity_threshold 30 \
                2> {log}
            touch {output.tr2}
        else
            fastp \
                -i {input.r1} -I {input.r2} \
                -o {output.tr1} -O {output.tr2} \
                --json {output.json} \
                --html {output.html} \
                --thread {threads} \
                --detect_adapter_for_pe \
                --qualified_quality_phred 20 \
                --length_required 50 \
                --low_complexity_filter \
                --complexity_threshold 30 \
                2> {log}
        fi
        # 'ok' explicito, nao arquivo vazio: _read_status_file() le
        # done.txt vazio como 'unknown', entao um fastp que terminou
        # bem aparecia no report como lacuna.
        echo ok > {output.done}
        """


# ── QC alternativo por amostra: cutadapt ──────────────────────────────
# TEMPORARIO. Existe porque algumas libs paired-end grandes nao terminam no
# fastp. So e definida quando `qc_cutadapt_samples` tem alguem, e o
# wildcard_constraints limita a regra a essas amostras -- as demais continuam
# no fastp, sem ambiguidade de DAG.
#
# NAO e um fastp equivalente, e a diferenca importa na hora de comparar
# numeros entre amostras:
#   - `-q 20` do cutadapt apara PONTAS de baixa qualidade; o
#     `--qualified_quality_phred 20` do fastp e o limiar do filtro de
#     percentual de bases nao qualificadas. Nao sao a mesma operacao.
#   - o cutadapt nao tem equivalente de `--low_complexity_filter`, entao
#     leituras de baixa complexidade que o fastp descartaria passam aqui.
#   - o adaptador vai fixo (TruSeq), que foi o que o proprio fastp detectou
#     nestas libs; nao ha auto-deteccao equivalente ao
#     `--detect_adapter_for_pe`.
# Por isso o relatorio NAO recebe um `{sample}_fastp.json` falso: a aba de QC
# mostra lacuna para estas amostras, em vez de exibir numero de cutadapt sob
# o nome do fastp.
if QC_CUTADAPT_SAMPLES:

    ruleorder: cutadapt > fastp

    rule cutadapt:
        """Adapter/quality trimming com cutadapt, para as amostras listadas
        em `qc_cutadapt_samples`. Escreve os MESMOS FASTQ aparados que o
        fastp escreveria, entao tudo a jusante (`_clean_r1`/`_clean_r2`)
        segue sem saber a diferenca."""
        wildcard_constraints:
            sample = "|".join(re.escape(s) for s in QC_CUTADAPT_SAMPLES) or "^$",
        input:
            r1 = lambda wc: SAMPLES[wc.sample]["R1"],
            r2 = lambda wc: SAMPLES[wc.sample].get("R2", []),
        output:
            tr1  = f"{OUTDIR}/{{sample}}/trimmed/{{sample}}_R1_fastp.fq.gz",
            tr2  = f"{OUTDIR}/{{sample}}/trimmed/{{sample}}_R2_fastp.fq.gz",
            json = f"{OUTDIR}/{{sample}}/qc_raw/{{sample}}_cutadapt.json",
            done = f"{OUTDIR}/{{sample}}/qc_raw/done.txt",
        log:
            f"{OUTDIR}/{{sample}}/logs/cutadapt.log"
        benchmark:
            f"{OUTDIR}/{{sample}}/benchmarks/cutadapt.tsv"
        conda:      "../envs/env_cutadapt.yaml"
        container:  CONTAINERS.get("cutadapt")
        threads: min(THREADS, 16)
        params:
            single_end = SINGLE_END,
            # TruSeq, os mesmos que o --detect_adapter_for_pe do fastp
            # identificou no log destas libs.
            a1 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
            a2 = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
        shell:
            """
            mkdir -p {OUTDIR}/{wildcards.sample}/trimmed \
                     {OUTDIR}/{wildcards.sample}/qc_raw
            if [ "{params.single_end}" = "True" ]; then
                cutadapt \
                    -a {params.a1} \
                    -q 20 -m 50 \
                    -j {threads} \
                    --json {output.json} \
                    -o {output.tr1} \
                    {input.r1} \
                    > {log} 2>&1
                touch {output.tr2}
            else
                cutadapt \
                    -a {params.a1} -A {params.a2} \
                    -q 20 -m 50 \
                    -j {threads} \
                    --json {output.json} \
                    -o {output.tr1} -p {output.tr2} \
                    {input.r1} {input.r2} \
                    > {log} 2>&1
            fi
            echo "ok: cutadapt (substituto temporario do fastp)" > {output.done}
            """


# ── Long read QC ──────────────────────────────────────────────────────
# Rules below are only active when LONG_READS=True.

if LONG_READS:

    rule nanoplot_lr:
        """Long read QC with NanoPlot. Replaces FastQC for long reads."""
        input:
            reads = (
                lambda wc: SAMPLES[wc.sample]["LR"]
                if LONG_READS else []
            ),
        output:
            done = f"{OUTDIR}/{{sample}}/qc_lr/nanoplot_done.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/nanoplot_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/nanoplot_lr.tsv"
        conda:      "../envs/env_lr_utils.yaml"
        container:  CONTAINERS.get("nanoplot")
        threads: min(THREADS, 8)
        params:
            outdir = f"{OUTDIR}/{{sample}}/qc_lr",
        shell:
            """
            if [ "{LONG_READS}" != "True" ]; then touch {output.done}; exit 0; fi
            mkdir -p {params.outdir}
            # NanoPlot accepts space-separated list of fastq files
            READS=$(echo {input.reads} | tr '\n' ' ')
            NanoPlot \
                --fastq $READS \
                --outdir {params.outdir} \
                --threads {threads} \
                --N50 \
                --loglength \
                --tsv_stats \
                --plots dot \
                --drop_outliers \
                > {log} 2>&1
            touch {output.done}
            """

    rule porechop_lr:
        """ONT adapter removal with Porechop_ABI. Skipped for HiFi."""
        input:
            reads = (
                lambda wc: SAMPLES[wc.sample]["LR"]
                if LONG_READS else []
            ),
        output:
            trimmed = f"{OUTDIR}/{{sample}}/lr_trimmed/{{sample}}_porechop.fastq.gz",
        log:   f"{OUTDIR}/{{sample}}/logs/porechop_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/porechop_lr.tsv"
        conda:      "../envs/env_lr_utils.yaml"
        container:  CONTAINERS.get("porechop_abi")
        threads: THREADS
        shell:
            """
            if [ "{LONG_READS}" != "True" ] || [ "{LR_TECH}" != "ont" ]; then
                mkdir -p $(dirname {output.trimmed})
                cat {input.reads} > {output.trimmed}; exit 0
            fi
            mkdir -p $(dirname {output.trimmed})
            # --no_split: keep chimeric reads intact (split in Filtlong instead)
            porechop_abi \
                --input {input.reads} \
                --output {output.trimmed} \
                --threads {threads} \
                --no_split \
                > {log} 2>&1
            """

    rule filtlong_lr:
        """
        Quality/length filtering with Filtlong. Applies to both ONT and HiFi.
        Input is always the porechop_lr output: for HiFi that rule is a
        pass-through (it only cats the raw reads), so both technologies
        share one upstream path. It used to point at a
        {sample}_placeholder.fastq.gz for HiFi -- a file no rule produces,
        which killed every HiFi run with MissingInputException.
        """
        input:
            reads = f"{OUTDIR}/{{sample}}/lr_trimmed/{{sample}}_porechop.fastq.gz",
        output:
            filtered = f"{OUTDIR}/{{sample}}/lr_filtered/{{sample}}_filtered.fastq.gz",
        log:   f"{OUTDIR}/{{sample}}/logs/filtlong_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/filtlong_lr.tsv"
        conda:      "../envs/env_lr_utils.yaml"
        container:  CONTAINERS.get("filtlong")
        threads: 2
        params:
            min_len      = LR_MIN_LEN,
            min_mean_q   = LR_MIN_MEAN_Q,
            min_window_q = 5,
            keep_pct     = 90,
        shell:
            """
            if [ "{LONG_READS}" != "True" ]; then
                mkdir -p $(dirname {output.filtered})
                cp {input.reads} {output.filtered}; exit 0
            fi
            mkdir -p $(dirname {output.filtered})
            if [ "{LR_TECH}" = "ont" ]; then
                # ONT: filter by length + mean quality + window quality
                filtlong \
                    --min_length {params.min_len} \
                    --min_mean_q {params.min_mean_q} \
                    --min_window_q {params.min_window_q} \
                    --keep_percent {params.keep_pct} \
                    {input.reads} 2> {log} | gzip > {output.filtered}
            else
                # HiFi: already Q20+, only filter by length
                filtlong \
                    --min_length {params.min_len} \
                    {input.reads} 2> {log} | gzip > {output.filtered}
            fi
            """
