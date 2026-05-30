# ══════════════════════════════════════════════════════════════════════
# rules/report.smk — BLOCK 12/13: HTML Report + MultiQC
#
# generate_report — relatório Plotly standalone (todos os samples)
#   Tabs: Overview | Read QC | Assembly QC | Viral Analysis |
#         Taxonomy | Bin Quality | Abundance | About
#   Script: scripts/generate_report.py
#   Requer: pip install plotly (no env do snakemake)
#
# multiqc — consolida FastQC + Trim Galore + QUAST (SR only)
# ══════════════════════════════════════════════════════════════════════


rule generate_report:
    """
    Relatório HTML standalone com Plotly — abre em qualquer browser.
    Agrega todos os samples numa única página multi-tab.
    """
    input:
        quast             = expand(f"{OUTDIR}/{{sample}}/quast/report.tsv",                            sample=SAMPLES),
        checkv            = expand(f"{OUTDIR}/{{sample}}/viral/checkv/quality_summary.tsv",            sample=SAMPLES),
        checkv_vrhyme     = expand(f"{OUTDIR}/{{sample}}/viral/checkv_vrhyme/quality_summary.tsv",     sample=SAMPLES),
        checkm2           = expand(f"{OUTDIR}/{{sample}}/bins/checkm2/quality_report.tsv",             sample=SAMPLES),
        support           = expand(f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_tool_support.tsv", sample=SAMPLES),
        depth             = expand(f"{OUTDIR}/{{sample}}/mapping/{{sample}}_depth.txt",                sample=SAMPLES),
        binette           = expand(f"{OUTDIR}/{{sample}}/bins/binette/binette_results.tsv",            sample=SAMPLES),
        taxonomy          = expand(f"{OUTDIR}/{{sample}}/viral/taxonomy/viral_taxonomy_merged.tsv",    sample=SAMPLES),
        vcontact3         = expand(f"{OUTDIR}/{{sample}}/viral/vcontact3/genome_clusters.tsv",         sample=SAMPLES),
        custom_prok       = expand(f"{OUTDIR}/{{sample}}/bins/diamond_custom_prok/diamond_vs_custom.tsv", sample=SAMPLES),
        gtdbtk_bac        = expand(f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.bac120.summary.tsv", sample=SAMPLES),
        gtdbtk_arc        = expand(f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.ar53.summary.tsv",   sample=SAMPLES),
        phist             = expand(f"{OUTDIR}/{{sample}}/viral/phist/phist_results.csv",               sample=SAMPLES),
        org               = expand(f"{OUTDIR}/{{sample}}/final/done.txt",                              sample=SAMPLES),
        benchmark_summary = f"{OUTDIR}/benchmarks/pipeline_timing_summary.tsv",
    output:
        html = f"{OUTDIR}/report.html",
    params:
        samples         = list(SAMPLES.keys()),
        outdir          = OUTDIR,
        threads         = THREADS,
        spades_mem      = SPADES_MEM,
        megahit_mem     = int(MEGAHIT_MEM // 1_000_000_000),
        min_contig      = MIN_CONTIG,
        min_seq_id      = MIN_SEQ_ID,
        semibin_env     = SEMIBIN_ENV,
        min_viral_tools = MIN_VIRAL_TOOLS,
        checkv_db       = CHECKV_DB,
        genomad_db      = GENOMAD_DB,
        inphared_db     = INPHARED_DB,
        gtdbtk_db        = GTDBTK_DB,
        custom_prok_meta = CUSTOM_PROK_META,
    benchmark:
        f"{OUTDIR}/benchmarks/generate_report.tsv"
    script:
        "../scripts/generate_report.py"


rule multiqc:
    """MultiQC consolida fastp + QUAST (short reads only)."""
    input:
        expand(f"{OUTDIR}/{{sample}}/qc_raw/{{sample}}_fastp.json", sample=SAMPLES),
        expand(f"{OUTDIR}/{{sample}}/quast/report.tsv",             sample=SAMPLES),
    output:
        report = f"{OUTDIR}/multiqc_report/multiqc_report.html",
    log:
        f"{OUTDIR}/logs/multiqc.log"
    benchmark:
        f"{OUTDIR}/benchmarks/multiqc.tsv"
    conda: "envs/env_qc.yaml"
    container:  CONTAINERS.get("multiqc")
    shell:
        """
        multiqc {OUTDIR} \
            -o {OUTDIR}/multiqc_report \
            --force \
            > {log} 2>&1
        """
