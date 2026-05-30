# ══════════════════════════════════════════════════════════════════════
# rules/prok_binning.smk — BLOCK 8: Prokaryote Binning
#
# Binners:
#   metabat2          — tetranucleotide + coverage
#   vamb              — Variational Autoencoder (v5 interface)
#   semibin2          — semi-supervised with environment priors
#
# Integration:
#   prepare_scaffold2bin — convert binner outputs to scaffold2bin.tsv
#   binette              — best-bin selection (successor to DAS Tool)
#
# QC + taxonomy:
#   checkm2           — completeness / contamination (MIMAG standards)
#   gtdbtk            — GTDB taxonomy for MAGs
# ══════════════════════════════════════════════════════════════════════


rule metabat2:
    """
    MetaBAT2 — tetranucleotide + coverage binning.
    --minContig enforced minimum is 1500 (tool limit) — we use max(MIN_CONTIG, 1500).
    """
    input:
        contigs = rules.mmseqs2.output.rep,
        depth   = rules.calc_depth.output.depth,
    output:
        done = f"{OUTDIR}/{{sample}}/bins/metabat2/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/metabat2.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/metabat2.tsv"
    conda: "envs/env_binning.yaml"
    container:  CONTAINERS.get("metabat2")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/bins/metabat2",
    shell:
        """
        mkdir -p {params.outdir}
        METABAT_MIN=$(( {MIN_CONTIG} > 1500 ? {MIN_CONTIG} : 1500 ))
        metabat2 \
            -i {input.contigs} \
            -a {input.depth} \
            -o {params.outdir}/bin \
            -t {threads} \
            --minContig $METABAT_MIN \
            --unbinned \
            > {log} 2>&1
        touch {output.done}
        """


rule vamb:
    """
    VAMB v5 — Variational Autoencoder binning.
    IMPORTANT v5 changes from v4:
      - Subcommand: vamb bin default
      - No --bamfiles; requires --abundance_tsv
      - TSV header must be: contigname <TAB> sample_name
    NOTE: VAMB creates its output dir itself — rm -rf before run.
    """
    input:
        contigs = rules.mmseqs2.output.rep,
        depth   = rules.calc_depth.output.depth,
    output:
        done = f"{OUTDIR}/{{sample}}/bins/vamb/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/vamb.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/vamb.tsv"
    conda: "envs/env_binning.yaml"
    container:  CONTAINERS.get("vamb")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/bins/vamb",
    shell:
        """
        rm -rf {params.outdir}
        echo -e "contigname\\t{wildcards.sample}" > {params.outdir}_abundance.tsv
        awk 'NR>1 {{print $1"\\t"$3}}' {input.depth} \
            >> {params.outdir}_abundance.tsv
        CUDA_FLAG=""
        if [ "{USE_GPU}" = "True" ]; then CUDA_FLAG="--cuda"; fi
        vamb bin default \
            --outdir {params.outdir} \
            --fasta {input.contigs} \
            --abundance_tsv {params.outdir}_abundance.tsv \
            --minfasta 200000 \
            -p {threads} \
            $CUDA_FLAG \
            > {log} 2>&1
        touch {output.done}
        """


rule semibin2:
    """
    SemiBin2 — semi-supervised binning with environment-specific priors.
    Set SEMIBIN_ENV to match your sample type:
    soil / ocean / gut / wastewater / global (use global if unsure).
    """
    input:
        contigs  = rules.mmseqs2.output.rep,
        bam      = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bwa_done = (f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt"
                    if not LONG_READS else []),
    output:
        done = f"{OUTDIR}/{{sample}}/bins/semibin2/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/semibin2.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/semibin2.tsv"
    conda: "envs/env_binning.yaml"
    container:  CONTAINERS.get("semibin")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/bins/semibin2",
    shell:
        """
        mkdir -p {params.outdir}
        ENGINE="cpu"
        if [ "{USE_GPU}" = "True" ]; then ENGINE="gpu"; fi
        SemiBin2 single_easy_bin \
            -i {input.contigs} \
            -b {input.bam} \
            -o {params.outdir} \
            --environment {SEMIBIN_ENV} \
            --engine $ENGINE \
            -t {threads} \
            > {log} 2>&1
        touch {output.done}
        """


rule comebin:
    """
    COMEBin — transformer embedding + contrastive learning binner.
    Rank 1 in 2025 metagenome binning benchmark (Nature Communications).
    Added as 4th binner alongside MetaBAT2 + VAMB + SemiBin2.
    Binette uses all 4 scaffold2bin files for refined bin selection.
    Set comebin_enabled: false in config.yaml to skip.
    """
    input:
        contigs  = rules.mmseqs2.output.rep,
        bam      = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bai      = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam.bai",
        bwa_done = (f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt"
                    if not LONG_READS else []),
    output:
        done = f"{OUTDIR}/{{sample}}/bins/comebin/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/comebin.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/comebin.tsv"
    conda: "envs/env_comebin.yaml"
    container:  CONTAINERS.get("comebin")
    threads: THREADS
    params:
        outdir  = f"{OUTDIR}/{{sample}}/bins/comebin",
        bam_dir = f"{OUTDIR}/{{sample}}/mapping",
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{COMEBIN_ENABLED}" != "True" ]; then
            echo "[COMEBin] Disabled via config (comebin_enabled: false)" | tee {log}
            touch {output.done}; exit 0
        fi
        # COMEBin: -p is the directory containing BAM files, not the BAM file itself
        #          -t is threads; -n is number of views (contrastive learning, keep default 6)
        #          GPU: set CUDA_VISIBLE_DEVICES=0 before run_comebin.sh
        if [ "{USE_GPU}" = "True" ]; then
            export CUDA_VISIBLE_DEVICES=0
        else
            export CUDA_VISIBLE_DEVICES=""
        fi
        run_comebin.sh \
            -a {input.contigs} \
            -o {params.outdir} \
            -p {params.bam_dir} \
            -t {threads} \
            >> {log} 2>&1
        touch {output.done}
        N=$(ls {params.outdir}/bins/*.fa 2>/dev/null | wc -l)
        echo "[COMEBin] $N bins produced" | tee -a {log}
        """


rule prepare_scaffold2bin:
    """
    Convert each binner output to scaffold2bin.tsv (contig <TAB> bin).
    Each binner has a different format:
      MetaBAT2 : bin.X.fa files
      VAMB v5  : clusters.tsv — REVERSED columns (bin_id TAB contig_id)
      SemiBin2 : contig_bins.tsv — contig TAB bin_number
      COMEBin  : bins/*.fa files (if comebin_enabled)
    """
    input:
        mb2     = rules.metabat2.output.done,
        vamb    = rules.vamb.output.done,
        sb2     = rules.semibin2.output.done,
        comebin = rules.comebin.output.done,
    output:
        done = f"{OUTDIR}/{{sample}}/bins/scaffold2bin/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/scaffold2bin.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/prepare_scaffold2bin.tsv"
    params:
        s = f"{OUTDIR}/{{sample}}",
    run:
        import os, glob

        outdir = f"{params.s}/bins/scaffold2bin"
        os.makedirs(outdir, exist_ok=True)

        def write_s2b(bin_files, outfile, strip_ext):
            with open(outfile, "w") as out:
                for bf in bin_files:
                    bin_name = os.path.basename(bf).replace(strip_ext, "")
                    with open(bf) as f:
                        for line in f:
                            if line.startswith(">"):
                                contig = line[1:].strip().split()[0]
                                out.write(f"{contig}\t{bin_name}\n")

        mb_bins = glob.glob(f"{params.s}/bins/metabat2/bin.*.fa")
        if mb_bins:
            write_s2b(mb_bins, f"{outdir}/metabat2_s2b.tsv", ".fa")

        vamb_clusters = f"{params.s}/bins/vamb/clusters.tsv"
        if os.path.exists(vamb_clusters):
            with open(vamb_clusters) as fin, open(f"{outdir}/vamb_s2b.tsv", "w") as fout:
                for line in fin:
                    parts = line.strip().split("\t")
                    if len(parts) == 2 and not parts[0].startswith("clusterid"):
                        fout.write(f"{parts[1]}\t{parts[0]}\n")

        # SemiBin2: use contig_bins.tsv directly (bins are SemiBin_N.fa.gz)
        sb_tsv = f"{params.s}/bins/semibin2/contig_bins.tsv"
        if os.path.exists(sb_tsv):
            with open(sb_tsv) as fin, open(f"{outdir}/semibin2_s2b.tsv", "w") as fout:
                header = True
                for line in fin:
                    if header:
                        header = False
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        contig, bin_num = parts
                        fout.write(f"{contig}\tSemiBin_{bin_num}\n")

        # COMEBin: bins/*.fa files (same format as MetaBAT2)
        cb_bins = glob.glob(f"{params.s}/bins/comebin/bins/*.fa")
        if cb_bins:
            write_s2b(cb_bins, f"{outdir}/comebin_s2b.tsv", ".fa")

        with open(output.done, "w") as f:
            f.write("ok\n")


rule binette:
    """
    Binette — best-bin selection across multiple binners (2023).
    Successor to DAS Tool: faster, better scoring, same scaffold2bin input format.
    Outputs: final_bins/ (FASTA) + binette_results.tsv (quality table).
    Filters scaffold2bin to only contigs present in the assembly (safety check).
    """
    input:
        contigs = rules.mmseqs2.output.rep,
        s2b     = rules.prepare_scaffold2bin.output.done,
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/binette/done.txt",
        summary = f"{OUTDIR}/{{sample}}/bins/binette/binette_results.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/binette.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/binette.tsv"
    conda: "envs/env_binette.yaml"
    container:  CONTAINERS.get("binette")
    threads: THREADS
    params:
        s2b_dir = f"{OUTDIR}/{{sample}}/bins/scaffold2bin",
        outdir  = f"{OUTDIR}/{{sample}}/bins/binette",
    shell:
        """
        mkdir -p {params.outdir}

        VALID_CONTIGS=$(mktemp)
        grep "^>" {input.contigs} | sed 's/^>//' | awk '{{print $1}}' > "$VALID_CONTIGS"
        echo "[binette] $(wc -l < $VALID_CONTIGS) assembly contigs" | tee {log}

        FILTERED_DIR={params.s2b_dir}/filtered
        mkdir -p "$FILTERED_DIR"
        rm -f "$FILTERED_DIR"/*.tsv
        for S2B in {params.s2b_dir}/*_s2b.tsv; do
            [ -f "$S2B" ] || continue
            BNAME=$(basename "$S2B")
            awk 'NR==FNR{{valid[$1]=1;next}} ($1 in valid){{print}}' \
                "$VALID_CONTIGS" "$S2B" > "$FILTERED_DIR/$BNAME"
            [ -s "$FILTERED_DIR/$BNAME" ] || rm -f "$FILTERED_DIR/$BNAME"
        done
        rm -f "$VALID_CONTIGS"

        S2B_FILES=$(ls "$FILTERED_DIR"/*_s2b.tsv 2>/dev/null | tr '\n' ' ')
        if [ -z "$(echo $S2B_FILES | tr -d ' ')" ]; then
            echo "[WARNING] No valid s2b files — skipping Binette" | tee -a {log}
            mkdir -p {params.outdir}/final_bins
            printf "Name\tCompleteness\tContamination\tScore\n" > {output.summary}
            touch {output.done}; exit 0
        fi

        binette -b $S2B_FILES -c {input.contigs} -o {params.outdir} -t {threads} --checkm2_db {CHECKM2_DB} >> {log} 2>&1

        # Locate summary (Binette may name it differently across versions)
        for candidate in {params.outdir}/binette_results.tsv \
                          {params.outdir}/final_bins_quality.tsv \
                          {params.outdir}/results.tsv; do
            [ -f "$candidate" ] && cp "$candidate" {output.summary} && break
        done
        [ -f {output.summary} ] || \
            printf "Name\tCompleteness\tContamination\tScore\n" > {output.summary}
        touch {output.done}
        """


rule checkm2:
    """
    CheckM2 — completeness and contamination of Binette final bins.
    Uses universal ML model — no lineage-specific marker database needed.
    MIMAG thresholds:
      High-quality   : completeness >= 90%, contamination <= 5%
      Medium-quality : completeness >= 50%, contamination <= 10%
    """
    input:
        done = rules.binette.output.done,
    output:
        report = f"{OUTDIR}/{{sample}}/bins/checkm2/quality_report.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/checkm2.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/checkm2.tsv"
    conda: "envs/env_checkm2.yaml"
    container:  CONTAINERS.get("checkm2")
    threads: THREADS
    params:
        bins_dir = lambda wc: f"{OUTDIR}/{wc.sample}/bins/binette/final_bins",
        outdir   = f"{OUTDIR}/{{sample}}/bins/checkm2",
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        N_BINS=$(ls {params.bins_dir}/*.fa 2>/dev/null | wc -l)
        if [ "$N_BINS" -eq 0 ]; then
            echo "[checkm2] No bins found — skipping" | tee {log}
            printf "Name\tCompleteness\tContamination\tGenome_Size\n" > {output.report}
            exit 0
        fi

        checkm2 predict \
            --threads {threads} \
            --input {params.bins_dir} \
            --output-directory {params.outdir} \
            -x fa \
            --database_path {CHECKM2_DB} \
            > {log} 2>&1 || echo "[checkm2] WARNING: predict failed" | tee -a {log}

        # Ensure output exists
        [ -f {output.report} ] || \
            printf "Name\tCompleteness\tContamination\tGenome_Size\n" > {output.report}
        """


rule gtdbtk:
    """
    GTDB-Tk classify_wf: assign GTDB taxonomy to Binette MAGs.
    Replaces NCBI taxonomy for metagenome-assembled genomes.
    Outputs: bac120.summary.tsv + ar53.summary.tsv
    Creates empty outputs if no bins are available (safe fallback).
    """
    input:
        done = rules.checkm2.output.report,
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/gtdbtk/done.txt",
        bac_tsv = f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.bac120.summary.tsv",
        ar_tsv  = f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/gtdbtk.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/gtdbtk.tsv"
    conda: "envs/env_gtdbtk.yaml"
    container:  CONTAINERS.get("gtdbtk")
    threads: THREADS
    params:
        bins_dir = lambda wc: f"{OUTDIR}/{wc.sample}/bins/binette/final_bins",
        outdir   = f"{OUTDIR}/{{sample}}/bins/gtdbtk",
    shell:
        """
        mkdir -p {params.outdir}

        # If no bins, create empty outputs
        N_BINS=$(ls {params.bins_dir}/*.fa 2>/dev/null | wc -l)
        if [ "$N_BINS" -eq 0 ]; then
            echo "[gtdbtk] No bins found — skipping" | tee {log}
            mkdir -p {params.outdir}/classify
            printf "user_genome\tclassification\n" > {output.bac_tsv}
            printf "user_genome\tclassification\n" > {output.ar_tsv}
            touch {output.done}; exit 0
        fi

        export GTDBTK_DATA_PATH={GTDBTK_DB}
        gtdbtk classify_wf \
            --genome_dir {params.bins_dir} \
            --out_dir    {params.outdir} \
            --cpus       {threads} \
            --extension  fa \
            >> {log} 2>&1 || echo "[gtdbtk] WARNING: classify_wf failed — creating empty outputs" | tee -a {log}

        # Always ensure output files exist — gtdbtk may fail if bins are low quality
        mkdir -p {params.outdir}/classify
        [ -f {output.bac_tsv} ] || printf "user_genome\tclassification\n" > {output.bac_tsv}
        [ -f {output.ar_tsv}  ] || printf "user_genome\tclassification\n" > {output.ar_tsv}
        touch {output.done}
        """
