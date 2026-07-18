# ══════════════════════════════════════════════════════════════════════
# rules/reads_classify.smk — BLOCK 14: Reads-only classification
#
# Ultrafast taxonomic profiling directly from raw reads using Sylph.
# Runs independently of the assembly pipeline — no contigs required.
#
# Supported databases (configured in config.yaml):
#   imgvr   — IMG/VR 4.1  (2.9M viral genomes,  ~2 GB)
#   uhgv    — UHGV        (171k gut vOTUs,      ~0.4 GB)
#   gtdb    — GTDB r226+  (prokaryotes,          ~4 GB)
#   custom  — any user-built .syldb (INPHARED or other)
#
# Short reads (PE/SE) and long reads (ONT/HiFi) use the same rule.
# ══════════════════════════════════════════════════════════════════════

import os as _os

_RC_SCRIPTS = _os.path.join(SCRIPTS_DIR, "reads_classify")

# ── Derive list of active DB paths and their sylph-tax identifiers ─────
_RC_DBS_CFG    = config.get("reads_classify_dbs", {})
_RC_KEY_ORDER  = ["imgvr", "uhgv", "gtdb", "custom"]
_GTDB_VER      = _RC_DBS_CFG.get("gtdb_version", "r226")
_RC_TAXID_MAP  = {
    "imgvr":  "IMGVR_4.1",
    "uhgv":   "UHGV",
    "gtdb":   f"GTDB_{_GTDB_VER}",
    "custom": _RC_DBS_CFG.get("custom_taxonomy", ""),
}

# Active: key → abs path (only non-empty values)
_RC_ACTIVE_DBS = {
    k: _expand(v)
    for k in _RC_KEY_ORDER
    if (v := _RC_DBS_CFG.get(k, "")) and k not in ("gtdb_version", "custom_taxonomy")
}

# Ordered list of paths for sylph profile command
_RC_DB_PATHS = list(_RC_ACTIVE_DBS.values())

# Ordered list of tax identifiers for sylph-tax taxprof -t
_RC_TAX_IDS = [
    _RC_TAXID_MAP[k]
    for k in _RC_KEY_ORDER
    if k in _RC_ACTIVE_DBS and _RC_TAXID_MAP.get(k)
]

# Whether to pass -a (host annotation) — only works for pre-built viral DBs
_RC_HAS_PREBUILT_VIRAL = any(k in _RC_ACTIVE_DBS for k in ("imgvr", "uhgv"))

# Taxonomy download dir
_RC_TAX_DIR = _expand(config.get("reads_classify_tax_dir", "")) \
    or _os.path.join(OUTDIR, "reads_classify", "taxonomy")

# Optional genome FASTA for BACPHLIP lifestyle prediction
_RC_GENOME_FASTA = _expand(config.get("reads_classify_genome_fasta", "")) \
    if config.get("reads_classify_genome_fasta", "") else ""

_RC_VIR_THRESHOLD    = config.get("reads_classify_virulence_threshold", 0.5)
_RC_MIN_PREVALENCE   = config.get("reads_classify_min_prevalence", 0.0)


# ── One-time taxonomy file download ────────────────────────────────────

rule sylph_tax_download:
    """Download sylph-tax taxonomy files (~50 MB). Runs once per installation."""
    output:
        done = _os.path.join(_RC_TAX_DIR, ".tax_downloaded"),
    log:
        _os.path.join(OUTDIR, "logs", "sylph_tax_download.log"),
    conda:      "../envs/env_reads_classify.yaml"
    container:  CONTAINERS.get("sylph_tax")
    params:
        tax_dir = _RC_TAX_DIR,
    shell:
        """
        mkdir -p {params.tax_dir}
        export SYLPH_TAXONOMY_CONFIG={params.tax_dir}/config.json
        sylph-tax download --download-to {params.tax_dir} 2>&1 | tee {log}
        touch {output.done}
        """


# ── Per-sample: profile reads against all active DBs ──────────────────

rule sylph_profile:
    """
    Profile raw reads against configured sylph databases.
    Handles PE short reads, SE short reads, and long reads (ONT/HiFi).
    """
    input:
        r1   = lambda wc: (
            SAMPLES[wc.sample].get("LR") or SAMPLES[wc.sample]["R1"]
        ),
        r2   = lambda wc: SAMPLES[wc.sample].get("R2", []),
        dbs  = _RC_DB_PATHS,
    output:
        tsv  = f"{OUTDIR}/{{sample}}/reads_classify/sylph_results.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/sylph_profile.log",
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/sylph_profile.tsv",
    conda:      "../envs/env_reads_classify.yaml"
    container:  CONTAINERS.get("sylph")
    threads: THREADS
    params:
        long_reads  = LONG_READS,
        single_end  = SINGLE_END,
        min_kmers   = config.get("reads_classify_min_kmers", 20),
    shell:
        """
        mkdir -p $(dirname {output.tsv})

        if [ "{params.long_reads}" = "True" ] || [ "{params.single_end}" = "True" ]; then
            reads_arg="-r {input.r1}"
        else
            reads_arg="-1 {input.r1} -2 {input.r2}"
        fi

        sylph profile \
            {input.dbs} \
            $reads_arg \
            --min-number-kmers {params.min_kmers} \
            -t {threads} \
            -o {output.tsv} \
            2> {log}
        """


# ── Per-sample: annotate with taxonomy ────────────────────────────────

rule sylph_tax:
    """
    Annotate sylph profile output with full taxonomy using sylph-tax.
    Produces one .sylphmpa file per configured database.
    """
    input:
        tsv      = f"{OUTDIR}/{{sample}}/reads_classify/sylph_results.tsv",
        tax_done = _os.path.join(_RC_TAX_DIR, ".tax_downloaded"),
    output:
        done = f"{OUTDIR}/{{sample}}/reads_classify/taxprof_done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/sylph_tax.log",
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/sylph_tax.tsv",
    conda:      "../envs/env_reads_classify.yaml"
    container:  CONTAINERS.get("sylph_tax")
    params:
        tax_ids  = _RC_TAX_IDS,
        prefix   = lambda wc: f"{OUTDIR}/{wc.sample}/reads_classify/{wc.sample}_",
        ann_flag = "-a" if _RC_HAS_PREBUILT_VIRAL else "",
        tax_dir  = _RC_TAX_DIR,
    shell:
        """
        export SYLPH_TAXONOMY_CONFIG={params.tax_dir}/config.json

        sylph-tax taxprof {input.tsv} \
            -t {params.tax_ids} \
            {params.ann_flag} \
            -o {params.prefix} \
            2> {log}

        touch {output.done}
        """


# ── Cross-sample merge ─────────────────────────────────────────────────

rule sylph_merge:
    """
    Merge per-sample sylphmpa files into cross-sample abundance tables.
    Produces two tables: relative taxonomic abundance and sequence abundance.
    """
    input:
        done = expand(
            f"{OUTDIR}/{{sample}}/reads_classify/taxprof_done.txt",
            sample=SAMPLES,
        ),
    output:
        abundance = f"{OUTDIR}/reads_classify/merged_relative_abundance.tsv",
        seqabund  = f"{OUTDIR}/reads_classify/merged_sequence_abundance.tsv",
    log:
        f"{OUTDIR}/logs/sylph_merge.log",
    conda:      "../envs/env_reads_classify.yaml"
    container:  CONTAINERS.get("sylph_tax")
    shell:
        """
        mkdir -p {OUTDIR}/reads_classify

        mapfile -t sylphmpa_files < <(
            find {OUTDIR}/*/reads_classify/ -name "*.sylphmpa" 2>/dev/null | sort
        )

        if [ "${{#sylphmpa_files[@]}}" -gt 0 ]; then
            sylph-tax merge "${{sylphmpa_files[@]}}" \
                --column relative_abundance \
                -o {output.abundance} 2>> {log}

            sylph-tax merge "${{sylphmpa_files[@]}}" \
                --column sequence_abundance \
                -o {output.seqabund} 2>> {log}
        else
            echo "No .sylphmpa files found" | tee -a {log}
            touch {output.abundance} {output.seqabund}
        fi
        """


# ── Filter by prevalence (optional) ───────────────────────────────────

rule reads_filter_prevalence:
    """Remove taxa present in fewer than min_prevalence fraction of samples."""
    input:
        f"{OUTDIR}/reads_classify/merged_relative_abundance.tsv",
    output:
        f"{OUTDIR}/reads_classify/merged_relative_abundance_filtered.tsv",
    log:
        f"{OUTDIR}/logs/reads_filter_prevalence.log",
    conda:      "../envs/env_reads_classify.yaml"
    container:  CONTAINERS.get("sylph_tax")
    params:
        min_prev = _RC_MIN_PREVALENCE,
        script   = _os.path.join(_RC_SCRIPTS, "filter_by_prevalence.py"),
    shell:
        """
        python {params.script} \
            {input} \
            {output} \
            {params.min_prev} \
            2> {log}
        """


# ── Export OTU table ───────────────────────────────────────────────────

rule reads_make_otu:
    """
    Reformat merged table to standard OTU format (#OTU_ID as first column).
    Compatible with QIIME2, phyloseq, and other diversity analysis tools.
    """
    input:
        f"{OUTDIR}/reads_classify/merged_relative_abundance_filtered.tsv",
    output:
        f"{OUTDIR}/reads_classify/otu_table.tsv",
    log:
        f"{OUTDIR}/logs/reads_make_otu.log",
    conda:      "../envs/env_reads_classify.yaml"
    container:  CONTAINERS.get("sylph_tax")
    params:
        script = _os.path.join(_RC_SCRIPTS, "make_otu.py"),
    shell:
        "python {params.script} {input} {output} 2> {log}"


# ── Collapse viral abundance by host genus ─────────────────────────────

rule reads_collapse_host:
    """
    Aggregate viral relative abundance by predicted bacterial host genus.
    Requires sylph-tax -a host annotation (available for IMGVR/UHGV databases).
    """
    input:
        f"{OUTDIR}/reads_classify/merged_relative_abundance_filtered.tsv",
    output:
        f"{OUTDIR}/reads_classify/viral_abundance_by_host.tsv",
    log:
        f"{OUTDIR}/logs/reads_collapse_host.log",
    conda:      "../envs/env_reads_classify.yaml"
    container:  CONTAINERS.get("sylph_tax")
    params:
        script = _os.path.join(_RC_SCRIPTS, "collapse_by_host.py"),
    shell:
        "python {params.script} {input} {output} 2> {log}"


# ── BACPHLIP lifestyle prediction (optional) ───────────────────────────
# Only active when reads_classify_genome_fasta is set in config.

if _RC_GENOME_FASTA:

    rule reads_bacphlip:
        """
        Predict viral lifestyle (virulent/temperate) for sylph-detected genomes.
        Requires reads_classify_genome_fasta: path to reference genome FASTA.
        """
        input:
            tsv   = f"{OUTDIR}/{{sample}}/reads_classify/sylph_results.tsv",
            fasta = _RC_GENOME_FASTA,
        output:
            done  = f"{OUTDIR}/{{sample}}/reads_classify/bacphlip/done.txt",
        log:
            f"{OUTDIR}/{{sample}}/logs/reads_bacphlip.log",
        benchmark:
            f"{OUTDIR}/{{sample}}/benchmarks/reads_bacphlip.tsv",
        conda:      "../envs/env_reads_classify.yaml"
        container:  CONTAINERS.get("sylph_tax")
        params:
            out_dir   = f"{OUTDIR}/{{sample}}/reads_classify/bacphlip",
            vir_thresh = _RC_VIR_THRESHOLD,
            script    = _os.path.join(_RC_SCRIPTS, "bacphlip_lifestyle.py"),
        shell:
            """
            python {params.script} \
                {input.tsv} \
                {input.fasta} \
                {params.out_dir} \
                {params.vir_thresh} \
                2> {log}
            touch {output.done}
            """


# ── Final sentinel ─────────────────────────────────────────────────────

localrules: reads_classify_done

rule reads_classify_done:
    """Sentinel marking successful completion of the reads-classify module."""
    input:
        otu     = f"{OUTDIR}/reads_classify/otu_table.tsv",
        host    = f"{OUTDIR}/reads_classify/viral_abundance_by_host.tsv",
        *(expand(
            f"{OUTDIR}/{{sample}}/reads_classify/bacphlip/done.txt",
            sample=SAMPLES,
        ) if _RC_GENOME_FASTA else []),
    output:
        f"{OUTDIR}/reads_classify/reads_classify_done.txt",
    shell:
        "touch {output}"
