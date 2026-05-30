# ══════════════════════════════════════════════════════════════════════
# rules/assembly.smk — BLOCK 2: Assembly
#
# Short reads : megahit + metaspades (parallel, both outputs feed merge_dedup)
# Long reads  : flye_lr + hifiasm_lr → medaka_lr (ONT only) → merge_lr
# ══════════════════════════════════════════════════════════════════════


# ── Short read assembly ───────────────────────────────────────────────

rule megahit:
    """
    MEGAHIT assembly. Fast, memory-efficient.
    MEGAHIT_PRESET controls k-mer strategy:
      "default"        — MEGAHIT auto k-mers (21,29,39,59,79,99,119,141)
      "meta-sensitive" — --presets meta-sensitive (same k-list, more sensitive)
      "meta-large"     — --presets meta-large (larger step, big/complex metagenomes)
      "custom"         — pass MEGAHIT_CUSTOM_PARAMS verbatim (e.g. --k-list 21,55,99)
    NOTE: fails if output dir already exists — removed before run.
    NOTE: -m is in bytes (MEGAHIT_MEM is defined in bytes in the config).
    """
    input:
        tr1 = rules.fastp.output.tr1,
        tr2 = rules.fastp.output.tr2,
    output:
        contigs = f"{OUTDIR}/{{sample}}/assembly/megahit/final.contigs.fa",
    log:
        f"{OUTDIR}/{{sample}}/logs/megahit.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/megahit.tsv"
    conda:      "envs/env_assembly.yaml"
    container:  CONTAINERS.get("megahit")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/assembly/megahit",
        preset = (
            f"--presets {MEGAHIT_PRESET}" if MEGAHIT_PRESET in ("meta-sensitive", "meta-large")
            else (MEGAHIT_CUSTOM_PARAMS   if MEGAHIT_PRESET == "custom" else "")
        ),
    shell:
        """
        rm -rf {params.outdir}
        megahit \
            -1 {input.tr1} -2 {input.tr2} \
            -t {threads} \
            -m {MEGAHIT_MEM} \
            --min-contig-len {MIN_CONTIG} \
            {params.preset} \
            -o {params.outdir} \
            > {log} 2>&1
        """

rule metaspades:
    """
    metaSPAdes assembly. Better for low-abundance organisms; more memory.
    Also outputs assembly graph (.gfa) — required by GraphMB.
    SPADES_KMERS controls k-mer selection:
      "default" — SPAdes auto-selects based on read length
      "custom"  — uses SPADES_KMER_LIST (e.g. "21,33,55,77,99,127")
    """
    input:
        tr1 = rules.fastp.output.tr1,
        tr2 = rules.fastp.output.tr2,
    output:
        contigs = f"{OUTDIR}/{{sample}}/assembly/metaspades/contigs.fasta",
        graph   = f"{OUTDIR}/{{sample}}/assembly/metaspades/assembly_graph_with_scaffolds.gfa",
    log:
        f"{OUTDIR}/{{sample}}/logs/metaspades.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/metaspades.tsv"
    conda:      "envs/env_assembly.yaml"
    container:  CONTAINERS.get("spades")
    threads: THREADS
    params:
        kmers = f"-k {SPADES_KMER_LIST}" if SPADES_KMERS == "custom" else "",
    shell:
        """
        spades.py --meta \
            -1 {input.tr1} -2 {input.tr2} \
            -t {threads} -m {SPADES_MEM} \
            {params.kmers} \
            -o {OUTDIR}/{wildcards.sample}/assembly/metaspades \
            > {log} 2>&1
        """


rule metaviral_spades:
    """
    metaviralSPAdes: viral-optimised assembly that reuses the metaSPAdes GFA graph
    (--assembly-graph), skipping expensive k-mer/graph construction (~70% less CPU).
    Extracts circular genomes and linear paths with terminal repeats (viral signatures).
    Output is merged alongside MEGAHIT + metaSPAdes in merge_contigs → MMseqs2.
    Falls back to an empty FASTA if no viral contigs are produced.
    """
    input:
        tr1   = rules.fastp.output.tr1,
        tr2   = rules.fastp.output.tr2,
        graph = rules.metaspades.output.graph,
    output:
        contigs = f"{OUTDIR}/{{sample}}/assembly/metaviral/contigs.fasta",
    log:
        f"{OUTDIR}/{{sample}}/logs/metaviral_spades.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/metaviral_spades.tsv"
    conda:      "envs/env_assembly.yaml"
    container:  CONTAINERS.get("spades")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/assembly/metaviral",
    shell:
        """
        mkdir -p {params.outdir}
        spades.py --metaviral \
            --assembly-graph {input.graph} \
            -1 {input.tr1} -2 {input.tr2} \
            -t {threads} -m {SPADES_MEM} \
            -o {params.outdir} \
            > {log} 2>&1 || true
        if [ -f {params.outdir}/scaffolds.fasta ] && [ -s {params.outdir}/scaffolds.fasta ]; then
            cp {params.outdir}/scaffolds.fasta {output.contigs}
        elif [ -f {params.outdir}/contigs.fasta ] && [ -s {params.outdir}/contigs.fasta ]; then
            cp {params.outdir}/contigs.fasta {output.contigs}
        else
            echo "[metaviral_spades] WARNING: no contigs produced — writing empty FASTA" \
                | tee -a {log}
            touch {output.contigs}
        fi
        """


# ── Long read assembly ────────────────────────────────────────────────
# Rules below are only active when LONG_READS=True.

if LONG_READS:

    rule flye_lr:
        """
        metaFlye assembly for long reads (ONT and HiFi).
        --meta: metagenome-specific graph simplification.
        ONT: --nano-raw or --nano-hq depending on read quality (LR_ONT_CHEM).
        HiFi: --pacbio-hifi for CCS/HiFi reads.
        """
        input:
            reads = f"{OUTDIR}/{{sample}}/lr_filtered/{{sample}}_filtered.fastq.gz",
        output:
            fasta = f"{OUTDIR}/{{sample}}/assembly/lr/flye/assembly.fasta",
            done  = f"{OUTDIR}/{{sample}}/assembly/lr/flye/done.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/flye_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/flye_lr.tsv"
        conda:      "envs/env_flye.yaml"
        container:  CONTAINERS.get("flye")
        threads: THREADS
        params:
            outdir  = f"{OUTDIR}/{{sample}}/assembly/lr/flye",
            overlap = LR_FLYE_OVERLAP,
        shell:
            """
            if [ "{LONG_READS}" != "True" ]; then
                mkdir -p {params.outdir}
                touch {output.fasta} {output.done}; exit 0
            fi
            mkdir -p {params.outdir}
            if [ "{LR_TECH}" = "hifi" ]; then
                READ_FLAG="--pacbio-hifi"
            elif [ "{LR_ONT_CHEM}" = "hq" ]; then
                READ_FLAG="--nano-hq"    # R10 / Q20+ reads
            else
                READ_FLAG="--nano-raw"   # R9 / older chemistry
            fi
            flye $READ_FLAG {input.reads} \
                --out-dir  {params.outdir} \
                --meta \
                --min-overlap {params.overlap} \
                --threads {threads} \
                --scaffold \
                --iterations 2 \
                >> {log} 2>&1 || true
            [ -f {params.outdir}/assembly.fasta ] && \
                cp {params.outdir}/assembly.fasta {output.fasta} || \
                touch {output.fasta}
            touch {output.done}
            """

    rule hifiasm_lr:
        """
        hifiasm-meta: metagenome-specific assembler for long reads.
        HiFi: hifiasm_meta (dedicated binary).
        ONT: regular hifiasm --ont (requires >= v0.21).
        Produces GFA graph — converted to FASTA via awk.
        """
        input:
            reads = f"{OUTDIR}/{{sample}}/lr_filtered/{{sample}}_filtered.fastq.gz",
        output:
            fasta = f"{OUTDIR}/{{sample}}/assembly/lr/hifiasm/assembly.fasta",
            done  = f"{OUTDIR}/{{sample}}/assembly/lr/hifiasm/done.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/hifiasm_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/hifiasm_lr.tsv"
        conda:      "envs/env_flye.yaml"
        container:  CONTAINERS.get("hifiasm")
        threads: THREADS
        params:
            outdir = f"{OUTDIR}/{{sample}}/assembly/lr/hifiasm",
            prefix = f"{OUTDIR}/{{sample}}/assembly/lr/hifiasm/{{sample}}",
        shell:
            """
            if [ "{LONG_READS}" != "True" ]; then
                mkdir -p {params.outdir}
                touch {output.fasta} {output.done}; exit 0
            fi
            mkdir -p {params.outdir}
            if [ "{LR_TECH}" = "hifi" ]; then
                hifiasm_meta \
                    -t {threads} \
                    -o {params.prefix} \
                    {input.reads} \
                    > {log} 2>&1 || true
            else
                hifiasm \
                    --ont \
                    -t {threads} \
                    -o {params.prefix} \
                    --sc-cut 10 \
                    {input.reads} \
                    > {log} 2>&1 || true
            fi
            # Convert primary contigs GFA → FASTA
            GFA=""
            for pattern in \
                "{params.outdir}/*.bp.p_ctg.gfa" \
                "{params.outdir}/*.p_ctg.gfa"; do
                for f in $pattern; do
                    [ -f "$f" ] && [ -s "$f" ] && GFA="$f" && break 2
                done
            done
            if [ -n "$GFA" ]; then
                awk '/^S/{{print ">"$2; print $3}}' "$GFA" > {output.fasta}
                echo "[hifiasm_lr] $(grep -c '^>' {output.fasta}) contigs from $GFA" | tee -a {log}
            else
                echo "[hifiasm_lr] WARNING: no GFA found" | tee -a {log}
                touch {output.fasta}
            fi
            touch {output.done}
            """

    rule medaka_lr:
        """
        ONT consensus polishing with Medaka v2+. Skipped for HiFi.
        LR_MEDAKA_MODEL = "auto" → --bacteria flag (recommended for metagenomics).
        Or set an explicit model: "r1041_e82_400bps_sup_v5.2.0"
        """
        input:
            flye_fa    = f"{OUTDIR}/{{sample}}/assembly/lr/flye/assembly.fasta",
            hifiasm_fa = f"{OUTDIR}/{{sample}}/assembly/lr/hifiasm/assembly.fasta",
            reads      = f"{OUTDIR}/{{sample}}/lr_filtered/{{sample}}_filtered.fastq.gz",
        output:
            flye_pol    = f"{OUTDIR}/{{sample}}/assembly/lr/flye_polished/assembly.fasta",
            hifiasm_pol = f"{OUTDIR}/{{sample}}/assembly/lr/hifiasm_polished/assembly.fasta",
            done        = f"{OUTDIR}/{{sample}}/assembly/lr/medaka_done.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/medaka_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/medaka_lr.tsv"
        conda:      "envs/env_medaka.yaml"
        container:  CONTAINERS.get("medaka")
        threads: THREADS
        params:
            model       = LR_MEDAKA_MODEL,
            flye_out    = f"{OUTDIR}/{{sample}}/assembly/lr/flye_polished",
            hifiasm_out = f"{OUTDIR}/{{sample}}/assembly/lr/hifiasm_polished",
        shell:
            """
            if [ "{LONG_READS}" != "True" ] || [ "{LR_TECH}" != "ont" ]; then
                mkdir -p {params.flye_out} {params.hifiasm_out}
                cp {input.flye_fa}    {output.flye_pol}    2>/dev/null || touch {output.flye_pol}
                cp {input.hifiasm_fa} {output.hifiasm_pol} 2>/dev/null || touch {output.hifiasm_pol}
                touch {output.done}; exit 0
            fi
            mkdir -p {params.flye_out} {params.hifiasm_out}
            if [ "{params.model}" = "auto" ]; then
                MODEL_FLAG="--bacteria"
            else
                MODEL_FLAG="-m {params.model}"
            fi
            polish_assembly() {{
                local draft=$1 outdir=$2 outfile=$3
                [ -s "$draft" ] || {{ touch "$outfile"; return; }}
                mkdir -p "$outdir"
                medaka_consensus \
                    -i {input.reads} -d "$draft" -o "$outdir" \
                    $MODEL_FLAG -t {threads} -f >> {log} 2>&1 || true
                if [ -f "$outdir/consensus.fasta" ]; then
                    cp "$outdir/consensus.fasta" "$outfile"
                else
                    echo "[medaka] polishing failed — using unpolished" | tee -a {log}
                    cp "$draft" "$outfile"
                fi
            }}
            polish_assembly {input.flye_fa}    {params.flye_out}    {output.flye_pol}
            polish_assembly {input.hifiasm_fa} {params.hifiasm_out} {output.hifiasm_pol}
            touch {output.done}
            """

    rule metaMDBG_lr:
        """
        metaMDBG — 3rd long-read assembler (HiFi or ONT R10+).
        2× more circular MAGs, 1.5–12× faster, 1/10 RAM vs Flye+hifiasm.
        Nature Biotechnology 2024.

        Activation:
          LR_TECH="hifi"                    → --preset hifi
          LR_TECH="ont" + LR_ONT_CHEM="hq" → --preset ont  (R10+ / Q20+)
          LR_TECH="ont" + LR_ONT_CHEM="raw" → SKIP (R9 not supported)
        Set lr_metaMDBG: false in config.yaml to disable entirely.
        """
        input:
            reads = f"{OUTDIR}/{{sample}}/lr_filtered/{{sample}}_filtered.fastq.gz",
        output:
            fasta = f"{OUTDIR}/{{sample}}/assembly/lr/metaMDBG/assembly.fasta",
            done  = f"{OUTDIR}/{{sample}}/assembly/lr/metaMDBG/done.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/metaMDBG_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/metaMDBG_lr.tsv"
        conda:      "envs/env_assembly.yaml"
        container:  CONTAINERS.get("metamdbg")
        threads: THREADS
        params:
            outdir = f"{OUTDIR}/{{sample}}/assembly/lr/metaMDBG",
        shell:
            """
            mkdir -p {params.outdir}
            if [ "{USE_METAMDBG}" != "True" ]; then
                echo "[metaMDBG] Skipped (USE_METAMDBG=False)" | tee {log}
                touch {output.fasta} {output.done}; exit 0
            fi
            if [ "{LR_TECH}" = "hifi" ]; then
                metaMDBG asm \
                    --in-hifi {input.reads} \
                    --out-dir {params.outdir} \
                    --threads {threads} \
                    --min-contig-length {MIN_CONTIG} \
                    > {log} 2>&1 || true
            else
                metaMDBG asm \
                    --in-ont {input.reads} \
                    --out-dir {params.outdir} \
                    --threads {threads} \
                    --min-contig-length {MIN_CONTIG} \
                    > {log} 2>&1 || true
            fi
            # metaMDBG outputs contigs.fasta or assembly.fasta depending on version
            for CANDIDATE in \
                {params.outdir}/contigs.fasta \
                {params.outdir}/assembly.fasta \
                {params.outdir}/*.fasta; do
                [ -f "$CANDIDATE" ] && [ -s "$CANDIDATE" ] && \
                    cp "$CANDIDATE" {output.fasta} && break
            done
            [ -f {output.fasta} ] || touch {output.fasta}
            touch {output.done}
            N=$(grep -c '^>' {output.fasta} 2>/dev/null || echo 0)
            echo "[metaMDBG] $N contigs assembled" | tee -a {log}
            """

    rule merge_lr:
        """
        Merge Flye + hifiasm (+ metaMDBG if USE_METAMDBG) assemblies.
        Labels contigs with assembler prefix (FLYE_ / HIFIASM_ / MDBG_) for traceability.
        Equivalent to merge_contigs for short reads — both feed MMseqs2.
        """
        input:
            flye_fa = lambda wc: (
                f"{OUTDIR}/{wc.sample}/assembly/lr/flye_polished/assembly.fasta"
                if LR_TECH == "ont"
                else f"{OUTDIR}/{wc.sample}/assembly/lr/flye/assembly.fasta"
            ),
            hifiasm_fa = lambda wc: (
                f"{OUTDIR}/{wc.sample}/assembly/lr/hifiasm_polished/assembly.fasta"
                if LR_TECH == "ont"
                else f"{OUTDIR}/{wc.sample}/assembly/lr/hifiasm/assembly.fasta"
            ),
            mdbg_fa = f"{OUTDIR}/{{sample}}/assembly/lr/metaMDBG/assembly.fasta",
            mdbg_done = f"{OUTDIR}/{{sample}}/assembly/lr/metaMDBG/done.txt",
        output:
            merged = f"{OUTDIR}/{{sample}}/assembly/lr/merged/{{sample}}_lr_merged.fasta",
            done   = f"{OUTDIR}/{{sample}}/assembly/lr/merged/done.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/merge_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/merge_lr.tsv"
        conda:      "envs/env_assembly.yaml"
        container:  CONTAINERS.get("mmseqs2")
        params:
            scripts_dir = SCRIPTS_DIR,
        shell:
            """
            if [ "{LONG_READS}" != "True" ]; then
                mkdir -p $(dirname {output.merged})
                touch {output.merged} {output.done}; exit 0
            fi
            mkdir -p $(dirname {output.merged})
            python3 {params.scripts_dir}/merge_lr_assemblies.py \
                {input.flye_fa} {input.hifiasm_fa} {input.mdbg_fa} \
                {output.merged} {wildcards.sample} \
                >> {log} 2>&1
            touch {output.done}
            echo "[merge_lr] Done" | tee -a {log}
            """
