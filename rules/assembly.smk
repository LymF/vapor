# ══════════════════════════════════════════════════════════════════════
# rules/assembly.smk — BLOCK 2: Assembly
#
# Short reads : megahit (single assembler — metaSPAdes/metaviralSPAdes removed)
# Long reads  : single assembler chosen by lr_tech —
#               ONT  → flye_lr → medaka_lr (polishing)
#               HiFi → metaMDBG_lr (no polishing)
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
    PE: -1/-2; SE: -r (single reads).
    NOTE: fails if output dir already exists — removed before run.
    NOTE: -m is in bytes (MEGAHIT_MEM is defined in bytes in the config).
    """
    input:
        tr1 = _clean_r1,
        tr2 = _clean_r2,
    output:
        contigs = f"{OUTDIR}/{{sample}}/assembly/megahit/final.contigs.fa",
    log:
        f"{OUTDIR}/{{sample}}/logs/megahit.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/megahit.tsv"
    conda:      "../envs/env_assembly.yaml"
    container:  CONTAINERS.get("megahit")
    threads: THREADS
    params:
        outdir     = f"{OUTDIR}/{{sample}}/assembly/megahit",
        preset     = (
            f"--presets {MEGAHIT_PRESET}" if MEGAHIT_PRESET in ("meta-sensitive", "meta-large")
            else (MEGAHIT_CUSTOM_PARAMS   if MEGAHIT_PRESET == "custom" else "")
        ),
        single_end = SINGLE_END,
    shell:
        """
        rm -rf {params.outdir}
        if [ "{params.single_end}" = "True" ]; then
            megahit \
                -r {input.tr1} \
                -t {threads} \
                -m {MEGAHIT_MEM} \
                --min-contig-len {MIN_CONTIG} \
                {params.preset} \
                -o {params.outdir} \
                > {log} 2>&1
        else
            megahit \
                -1 {input.tr1} -2 {input.tr2} \
                -t {threads} \
                -m {MEGAHIT_MEM} \
                --min-contig-len {MIN_CONTIG} \
                {params.preset} \
                -o {params.outdir} \
                > {log} 2>&1
        fi
        """

# ── Long read assembly ────────────────────────────────────────────────
# Rules below are only active when LONG_READS=True.

if LONG_READS:

  if LR_TECH == "ont":

    rule flye_lr:
        """
        metaFlye assembly — sole long-read assembler for ONT.
        --meta: metagenome-specific graph simplification.
        --nano-raw or --nano-hq depending on read chemistry (LR_ONT_CHEM).
        HiFi does not reach this rule: it is assembled by metaMDBG_lr.
        """
        input:
            reads = _clean_lr,
        output:
            fasta = f"{OUTDIR}/{{sample}}/assembly/lr/flye/assembly.fasta",
            done  = f"{OUTDIR}/{{sample}}/assembly/lr/flye/done.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/flye_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/flye_lr.tsv"
        conda:      "../envs/env_flye.yaml"
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
            if [ "{LR_ONT_CHEM}" = "hq" ]; then
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
                >> {log} 2>&1 && RC=0 || RC=$?
            [ -f {params.outdir}/assembly.fasta ] && \
                cp {params.outdir}/assembly.fasta {output.fasta} || \
                touch {output.fasta}
            if [ "$RC" -ne 0 ]; then
                echo "failed: flye exit $RC" > {output.done}
            else
                echo "ok" > {output.done}
            fi
            """

    rule medaka_lr:
        """
        ONT consensus polishing with Medaka v2+. Only instantiated for ONT
        (HiFi goes through metaMDBG_lr, which needs no polishing).
        LR_MEDAKA_MODEL = "auto" → --bacteria flag (recommended for metagenomics).
        Or set an explicit model: "r1041_e82_400bps_sup_v5.2.0"
        """
        input:
            flye_fa    = f"{OUTDIR}/{{sample}}/assembly/lr/flye/assembly.fasta",
            reads      = _clean_lr,
        output:
            flye_pol    = f"{OUTDIR}/{{sample}}/assembly/lr/flye_polished/assembly.fasta",
            done        = f"{OUTDIR}/{{sample}}/assembly/lr/medaka_done.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/medaka_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/medaka_lr.tsv"
        conda:      "../envs/env_medaka.yaml"
        container:  CONTAINERS.get("medaka")
        threads: THREADS
        params:
            model       = LR_MEDAKA_MODEL,
            flye_out    = f"{OUTDIR}/{{sample}}/assembly/lr/flye_polished",
        shell:
            """
            if [ "{LONG_READS}" != "True" ] || [ "{LR_TECH}" != "ont" ]; then
                mkdir -p {params.flye_out}
                cp {input.flye_fa}    {output.flye_pol}    2>/dev/null || touch {output.flye_pol}
                touch {output.done}; exit 0
            fi
            mkdir -p {params.flye_out}
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
            touch {output.done}
            """

  else:

    rule metaMDBG_lr:
        """
        metaMDBG — sole long-read assembler for HiFi.
        2× more circular MAGs, 1.5–12× faster, 1/10 RAM vs Flye.
        Nature Biotechnology 2024. No polishing step (HiFi reads are
        already near Q30+ accuracy).
        """
        input:
            reads = _clean_lr,
        output:
            fasta = f"{OUTDIR}/{{sample}}/assembly/lr/metaMDBG/assembly.fasta",
            done  = f"{OUTDIR}/{{sample}}/assembly/lr/metaMDBG/done.txt",
        log:   f"{OUTDIR}/{{sample}}/logs/metaMDBG_lr.log"
        benchmark: f"{OUTDIR}/{{sample}}/benchmarks/metaMDBG_lr.tsv"
        conda:      "../envs/env_assembly.yaml"
        container:  CONTAINERS.get("metamdbg")
        threads: THREADS
        params:
            outdir = f"{OUTDIR}/{{sample}}/assembly/lr/metaMDBG",
        shell:
            """
            mkdir -p {params.outdir}
            metaMDBG asm \
                --in-hifi {input.reads} \
                --out-dir {params.outdir} \
                --threads {threads} \
                --min-contig-length {MIN_CONTIG} \
                > {log} 2>&1 || true
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
