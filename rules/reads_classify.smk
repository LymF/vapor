# ══════════════════════════════════════════════════════════════════════
# rules/reads_classify.smk — BLOCK 14: Reads-only classification
#
# Ultrafast taxonomic profiling from reads using Sylph — sem montagem, mas
# NAO sem pre-processamento: o perfil roda sobre os reads APARADOS (fastp /
# filtlong) e, quando `host_genome` esta configurado, sobre os reads SEM
# HOSPEDEIRO. E a mesma entrada que a montagem usa, pelos mesmos helpers
# (`_clean_r1`/`_clean_r2`/`_clean_lr`, Snakefile).
#
# Ate 2026-08-19 esta trilha lia os FASTQ CRUS: adaptadores, caudas de baixa
# qualidade e reads de hospedeiro entravam no perfil, e o `host_genome`
# configurado nao tinha efeito nenhum aqui. Consequencia pratica: rodar so a
# trilha de reads agora dispara o trimming (e o mapeamento contra o hospedeiro,
# se configurado) -- e assim que tem de ser.
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
    Profile reads against configured sylph databases.
    Handles PE short reads, SE short reads, and long reads (ONT/HiFi).

    Entrada: reads aparados e, se `host_genome` estiver configurado, sem
    hospedeiro -- os mesmos `_clean_*` que alimentam a montagem. Ver a nota
    no topo deste arquivo.
    """
    input:
        r1   = lambda wc: _clean_lr(wc) if LONG_READS else _clean_r1(wc),
        r2   = lambda wc: [] if LONG_READS else _clean_r2(wc),
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
        # NOTE: 2026-08-26 — the -a flag (host annotation) was causing sylph-tax
        # taxprof to fail silently with exit code 1. Temporarily disabled to
        # unblock the pipeline; investigate sylph-tax 1.9.1 compatibility.
        # ann_flag = "-a" if _RC_HAS_PREBUILT_VIRAL else "",
        ann_flag = "",
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
    params:
        # Diretorios das amostras DESTA corrida, nao "{OUTDIR}/*/". O glob
        # varria o disco: um sample renomeado ou removido do config deixa seu
        # .sylphmpa antigo no OUTDIR e ele entrava na mesclagem em silencio,
        # como se fosse uma amostra a mais.
        dirs = [f"{OUTDIR}/{s}/reads_classify" for s in SAMPLES],
    shell:
        """
        mkdir -p {OUTDIR}/reads_classify

        mapfile -t sylphmpa_files < <(
            find {params.dirs} -maxdepth 1 -name "*.sylphmpa" 2>/dev/null | sort
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

rule reads_host_map:
    """
    Recuperar o hospedeiro anotado pelo BANCO dos .sylphmpa por amostra.

    O `sylph_merge` roda `sylph-tax merge --column relative_abundance`, entao
    a coluna "Virus_host (if viral)" nao chega a tabela mesclada -- e era so
    de la que o `reads_collapse_host` tentava ler o hospedeiro, o que fazia o
    agrupamento nunca acontecer. Esta regra le a fonte original.
    """
    input:
        done = expand(
            f"{OUTDIR}/{{sample}}/reads_classify/taxprof_done.txt",
            sample=SAMPLES,
        ),
    output:
        f"{OUTDIR}/reads_classify/host_map_db.tsv",
    log:
        f"{OUTDIR}/logs/reads_host_map.log",
    conda:      "../envs/env_reads_classify.yaml"
    container:  CONTAINERS.get("sylph_tax")
    params:
        script = _os.path.join(_RC_SCRIPTS, "build_host_map.py"),
        # mesma razao do `params.dirs` em sylph_merge: nada de varrer o disco.
        dirs   = [f"{OUTDIR}/{s}/reads_classify" for s in SAMPLES],
    shell:
        """
        mkdir -p {OUTDIR}/reads_classify
        mapfile -t mpa < <(find {params.dirs} -maxdepth 1 -name "*.sylphmpa" 2>/dev/null | sort)
        if [ "${{#mpa[@]}}" -eq 0 ]; then
            printf "clade_name\thost_db\n" > {output}
            echo "[reads_host_map] nenhum .sylphmpa encontrado" > {log}
        else
            python {params.script} {output} "${{mpa[@]}}" 2> {log}
        fi
        """


rule reads_collapse_host:
    """
    Agregar abundancia viral por genero do hospedeiro anotado pelo BANCO.

    Fonte unica: a coluna "Virus_host (if viral)" dos .sylphmpa, recuperada por
    `reads_host_map` -- o `sylph-tax merge` a descarta. O PHIST sobre genomas de
    referencia foi descartado neste track; ver o bloco de comentario abaixo.

    O sidecar `viral_host_assignments.tsv` traz a atribuicao por taxon, para a
    proveniencia do numero agregado ser auditavel: `host_source` distingue "o
    banco atribuiu" de "ninguem atribuiu".
    """
    input:
        table    = f"{OUTDIR}/reads_classify/merged_relative_abundance_filtered.tsv",
        host_map = rules.reads_host_map.output,
    output:
        table       = f"{OUTDIR}/reads_classify/viral_abundance_by_host.tsv",
        assignments = f"{OUTDIR}/reads_classify/viral_host_assignments.tsv",
    log:
        f"{OUTDIR}/logs/reads_collapse_host.log",
    conda:      "../envs/env_reads_classify.yaml"
    container:  CONTAINERS.get("sylph_tax")
    params:
        script    = _os.path.join(_RC_SCRIPTS, "collapse_by_host.py"),
    shell:
        """
        python {params.script} {input.table} {output.table} \
            --host-map {input.host_map} \
            --assignments {output.assignments} 2> {log}
        """


# -- BACPHLIP e PHIST NAO existem neste track, por decisao (2026-08-19) ------
#
# Ambos precisam das SEQUENCIAS dos genomas de referencia, e o que ha em disco
# sao apenas os sketches .syldb do sylph -- um .syldb e um esboco de k-mers, do
# qual nao se recupera sequencia. Habilita-los exigiria baixar o IMG/VR em
# FASTA e, para o PHIST, os representantes do GTDB r232: dezenas de GB numa
# pipeline que ja pede ~500 GB de bancos.
#
# E o custo nao compraria nada novo. As duas perguntas ja tem resposta sobre os
# virus DO USUARIO, e nao sobre genomas de referencia publicos:
#   lifestyle  -> `bacphlip_votu`  (votu_catalog.smk, sobre all_fasta)
#   hospedeiro -> `rule phist`     (host_prediction.smk, mq_fasta vs os MAGs)
# A versao deste track responderia as mesmas perguntas sobre o IMG/VR, que nao
# e o objeto de estudo. O hospedeiro anotado pelo BANCO continua disponivel,
# via `reads_host_map` -> coluna host_db.
#
# Se um dia voltar: o PHIST exige um FASTA POR GENOMA dos dois lados (o kmer-db
# emite uma linha por ARQUIVO, nao por sequencia) mais dois arquivos de lista
# com um caminho por linha -- ver split_viral_fastas.py. E atencao: os
# resultados do `rule phist` NAO servem aqui, porque descrevem contigs montados
# ("k141_...") e nao os genomas de referencia do IMG/VR ("t__IMGVR_UViG_...");
# a juncao sairia vazia, mais um caso da familia de bugs de namespace.


# ── Final sentinel ─────────────────────────────────────────────────────

localrules: reads_classify_done

rule reads_classify_done:
    """Sentinel marking successful completion of the reads-classify module."""
    input:
        otu     = f"{OUTDIR}/reads_classify/otu_table.tsv",
        host    = f"{OUTDIR}/reads_classify/viral_abundance_by_host.tsv",
    output:
        f"{OUTDIR}/reads_classify/reads_classify_done.txt",
    shell:
        "touch {output}"
