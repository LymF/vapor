# ══════════════════════════════════════════════════════════════════════
# rules/host_prediction.smk — BLOCK 10: Phage-Host Prediction
#
# phist — fast k-mer similarity vs sample MAGs (local hosts)
#
# phist uses split_viral_fastas.py to create per-genome FASTA files
# because kmer-db assigns one result row per input file.
#
# Host-side defense/antidefense/AMR annotation lives in defense_amr.smk
# (BLOCK 10.5) and is cross-linked to these PHIST results in the report.
# ══════════════════════════════════════════════════════════════════════


rule phist:
    """
    PHIST: fast phage-host prediction using k-mer similarity.
    Uses sample Binette MAGs as candidate hosts.
    Per-genome mode: each viral contig/vRhyme bin = individual FASTA
    so kmer-db assigns one result row per genome.
    Script: scripts/split_viral_fastas.py
    """
    input:
        viral  = rules.viral_votu_reps.output.mq_fasta,
        gtdbtk = rules.gtdbtk.output.done,
    output:
        done    = f"{OUTDIR}/{{sample}}/viral/phist/done.txt",
        results = f"{OUTDIR}/{{sample}}/viral/phist/phist_results.csv",
    log:
        f"{OUTDIR}/{{sample}}/logs/phist.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/phist.tsv"
    conda: "../envs/env_phist.yaml"
    container:  CONTAINERS.get("phist")
    threads: THREADS
    params:
        bins_dir    = lambda wc: f"{OUTDIR}/{wc.sample}/bins/binette/final_bins",
        vrhyme_dir  = lambda wc: f"{OUTDIR}/{wc.sample}/bins/vrhyme",
        outdir      = f"{OUTDIR}/{{sample}}/viral/phist",
        scripts_dir = SCRIPTS_DIR,
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        N_BINS=$(find {params.bins_dir} -maxdepth 1 -name "*.fa" 2>/dev/null | wc -l)
        if [ "$N_BINS" -eq 0 ] || [ ! -s {input.viral} ]; then
            echo "[phist] No MAGs or no viral sequences — skipping" | tee {log}
            printf "phage,host,#common-kmers,pvalue,adj-pvalue\n" > {output.results}
            touch {output.done}; exit 0
        fi

        # Split viral sequences into individual FASTA files (per-genome kmer-db mode)
        VFASTA_DIR={params.outdir}/viral_fastas
        mkdir -p "$VFASTA_DIR"
        rm -f "$VFASTA_DIR"/*.fasta
        python3 {params.scripts_dir}/split_viral_fastas.py \
            {input.viral} {params.vrhyme_dir} "$VFASTA_DIR" \
            >> {log} 2>&1
        echo "[phist] Total viral genomes: $(find $VFASTA_DIR -maxdepth 1 -name "*.fasta" 2>/dev/null | wc -l)" \
            | tee -a {log}

        # Build k-mer DB — file list mode (one file per genome = one row in output)
        ls "$VFASTA_DIR"/*.fasta > {params.outdir}/phage.list 2>/dev/null
        if [ ! -s {params.outdir}/phage.list ]; then
            echo "[phist] No viral fastas found" | tee -a {log}
            printf "phage,host,#common-kmers,pvalue,adj-pvalue\n" > {output.results}
            touch {output.done}; exit 0
        fi

        kmer-db build -k 25 -t {threads} \
            {params.outdir}/phage.list \
            {params.outdir}/phages.db \
            >> {log} 2>&1

        ls {params.bins_dir}/*.fa > {params.outdir}/bacteria.list 2>/dev/null

        kmer-db new2all -sparse -t {threads} \
            {params.outdir}/phages.db \
            {params.outdir}/bacteria.list \
            {params.outdir}/kmers.csv \
            >> {log} 2>&1

        if [ ! -s {params.outdir}/kmers.csv ] || \
           [ "$(wc -l < {params.outdir}/kmers.csv)" -lt 2 ]; then
            echo "[phist] No shared k-mers found" | tee -a {log}
            printf "phage,host,#common-kmers,pvalue,adj-pvalue\n" > {output.results}
        else
            phist {params.outdir}/kmers.csv {output.results} >> {log} 2>&1 || \
                printf "phage,host,#common-kmers,pvalue,adj-pvalue\n" > {output.results}
        fi
        touch {output.done}
        """
