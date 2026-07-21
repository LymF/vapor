# ══════════════════════════════════════════════════════════════════════
# rules/report.smk — BLOCK 12/13: HTML Report + MultiQC
#
# generate_report — VAPOR standalone HTML report (ECharts + D3)
#   Tabs: Overview | Sequencing | Viral | Prokaryotic | Diversity | Annotation | About
#   Entry: scripts/generate_report.py → scripts/report/renderer.py
#
# multiqc — consolida fastp + QUAST (SR only)
# ══════════════════════════════════════════════════════════════════════


rule generate_report:
    """
    VAPOR HTML report — standalone, no CDN required (ECharts + D3 inlined).
    Aggregates all samples in one multi-tab page.
    """
    input:
        # Viral track
        **({
            "checkv":             expand(f"{OUTDIR}/{{sample}}/viral/checkv/quality_summary.tsv",                sample=SAMPLES),
            "checkv_vrhyme":      expand(f"{OUTDIR}/{{sample}}/viral/checkv_vrhyme/quality_summary.tsv",         sample=SAMPLES),
            "support":            expand(f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_tool_support.tsv",     sample=SAMPLES),
            "taxonomy":           expand(f"{OUTDIR}/{{sample}}/viral/taxonomy/viral_taxonomy_merged.tsv",        sample=SAMPLES),
            "vcontact3":          expand(f"{OUTDIR}/{{sample}}/viral/vcontact3/genome_clusters.tsv",             sample=SAMPLES),
            "antidefense_viral":  expand(f"{OUTDIR}/{{sample}}/viral/defensefinder/viral_antidefense_systems.tsv", sample=SAMPLES),
            "dbapis_viral":       expand(f"{OUTDIR}/{{sample}}/viral/dbapis/dbapis_hits.tsv",                     sample=SAMPLES),
        } if TRACK_VIRAL else {}),
        # Prokaryotic track
        **({
            "checkm2":                expand(f"{OUTDIR}/{{sample}}/bins/checkm2/quality_report.tsv",                 sample=SAMPLES),
            "binette":                expand(f"{OUTDIR}/{{sample}}/bins/binette/binette_results.tsv",                sample=SAMPLES),
            "mmseqs_prok":            expand(f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/taxonomy.tsv",          sample=SAMPLES),
            "gtdbtk_bac":             expand(f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.bac120.summary.tsv", sample=SAMPLES),
            "gtdbtk_arc":             expand(f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.ar53.summary.tsv",   sample=SAMPLES),
            "defensefinder":          expand(f"{OUTDIR}/{{sample}}/bins/defensefinder/defensefinder_systems.tsv",     sample=SAMPLES),
            "antidefensefinder":      expand(f"{OUTDIR}/{{sample}}/bins/defensefinder/antidefensefinder_systems.tsv",  sample=SAMPLES),
            "amr_consensus":          expand(f"{OUTDIR}/{{sample}}/bins/amr_consensus/amr_consensus.tsv",             sample=SAMPLES),
            "prok_protein_manifest":  expand(f"{OUTDIR}/{{sample}}/bins/proteins/manifest.txt",                   sample=SAMPLES),
        } if TRACK_PROK else {}),
        # Host prediction / integration
        **({
            "phist": expand(f"{OUTDIR}/{{sample}}/viral/phist/phist_results.csv", sample=SAMPLES),
        } if INTEGRATION_ENABLED else {}),
        # Assembly QC + mapping depth + diversity + finalize (either track)
        **({
            "quast":             expand(f"{OUTDIR}/{{sample}}/quast/report.tsv",             sample=SAMPLES),
            "depth":             expand(f"{OUTDIR}/{{sample}}/mapping/{{sample}}_depth.txt",  sample=SAMPLES),
            "alpha_div":         f"{OUTDIR}/diversity/alpha_diversity.tsv",
            "pcoa_viral":        f"{OUTDIR}/diversity/beta_pcoord_viral.tsv",
            "pcoa_prok":         f"{OUTDIR}/diversity/beta_pcoord_prok.tsv",
            "pcoa_combined":     f"{OUTDIR}/diversity/beta_pcoord_combined.tsv",
            "org":               expand(f"{OUTDIR}/{{sample}}/final/done.txt",                sample=SAMPLES),
            "benchmark_summary": f"{OUTDIR}/benchmarks/pipeline_timing_summary.tsv",
        } if (TRACK_VIRAL or TRACK_PROK) else {}),
        # reads_classify module (optional — only wired when enabled)
        **({
            "reads_classify_abundance": f"{OUTDIR}/reads_classify/merged_relative_abundance_filtered.tsv",
            "reads_classify_host":      f"{OUTDIR}/reads_classify/viral_abundance_by_host.tsv",
        } if READS_CLASSIFY_ENABLED else {}),
        # coassembly module (optional — group MAGs report tab)
        **({
            "coassembly_sentinel": expand(
                f"{OUTDIR}/coassembly/{{group}}/gtdbtk/done.txt", group=list(GROUPS.keys()))
        } if (COASSEMBLY_ENABLED and COASSEMBLY_BINNING and not LONG_READS and GROUPS) else {}),
    output:
        html = f"{OUTDIR}/report.html",
    params:
        samples          = list(SAMPLES.keys()),
        outdir           = OUTDIR,
        threads          = THREADS,
        spades_mem       = SPADES_MEM,
        megahit_mem      = int(MEGAHIT_MEM // 1_000_000_000),
        min_contig       = MIN_CONTIG,
        min_seq_id       = MIN_SEQ_ID,
        semibin_env      = SEMIBIN_ENV,
        min_viral_tools  = MIN_VIRAL_TOOLS,
        checkv_db        = CHECKV_DB,
        genomad_db       = GENOMAD_DB,
        inphared_db      = INPHARED_DB,
        gtdbtk_db        = GTDBTK_DB,
        # Same fallback as rule dbapis_viral's params.apis_dir -- the
        # seed_and_familyrep_all_infor.tsv family->gene/defense-system
        # mapping lives wherever that rule actually downloaded it.
        apis_db_dir      = APIS_DB or f"{OUTDIR}/dbapis_db",
        low_depth_mode   = LOW_DEPTH_MODE,
        coassembly_groups = list(GROUPS.keys()),
        tracks = {
            "reads":       bool(TRACK_READS or READS_CLASSIFY_ENABLED),
            "viral":       bool(TRACK_VIRAL),
            "prok":        bool(TRACK_PROK),
            "integration": bool(INTEGRATION_ENABLED),
            "coassembly":  bool(COASSEMBLY_ENABLED),
        },
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
    conda: "../envs/env_qc.yaml"
    container:  CONTAINERS.get("multiqc")
    shell:
        """
        multiqc {OUTDIR} \
            -o {OUTDIR}/multiqc_report \
            --force \
            > {log} 2>&1
        """
