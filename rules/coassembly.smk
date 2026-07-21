# co-assembly / co-binning rules (Plan 2)
#
# Canonical per-group assembly output:  {OUTDIR}/coassembly/{group}/contigs.fa
#   short reads (PE/SE) → megahit_coassembly   (if not LONG_READS)
#   long reads   (LR)   → flye_coassembly      (if LONG_READS)
# Downstream rules (map-back, co-binning) consume the canonical path, so they
# are mode-agnostic.

# Group → list of each sample's cleaned reads (post host-removal / trimming).
def _group_r1(wc):
    return [_clean_r1(type("W", (), {"sample": s})) for s in GROUPS[wc.group]]

def _group_r2(wc):
    if SINGLE_END:
        return []
    return [_clean_r2(type("W", (), {"sample": s})) for s in GROUPS[wc.group]]

def _group_lr(wc):
    return [_clean_lr(type("W", (), {"sample": s})) for s in GROUPS[wc.group]]


_group_samples = sorted({s for members in GROUPS.values() for s in members}) if GROUPS else []

# Multi-split co-binning (Task 5b) maps EVERY sample against one global
# catalog, regardless of co-assembly grouping — so the `sample` wildcard
# constraint must also cover samples that a metadata-based grouping might
# not have assigned to any group.
_multisplit_samples = sorted(SAMPLES.keys()) if COBINNING_MULTISPLIT and not LONG_READS else []
_sample_constraint_set = sorted(set(_group_samples) | set(_multisplit_samples))

if GROUPS:
    _grp_re = "|".join(re.escape(g) for g in GROUPS)
    wildcard_constraints:
        group = _grp_re

if _sample_constraint_set:
    _smp_re = "|".join(re.escape(s) for s in _sample_constraint_set)
    wildcard_constraints:
        sample = _smp_re


def _coassembly_prok_contigs(wc):
    """FASTA used as input to co-binning (VAMB co-binning).
    Returns the free-living-virus-filtered FASTA (mirrors per-sample
    `_prok_input_contigs` / `filter_viral_for_prok`) when COASSEMBLY_VIRAL is
    enabled; otherwise falls back to the raw group co-assembly contigs (no
    viral detection ran, so no filter is possible).
    NOTE: built as a literal path (not `rules.coassembly_filter_viral_for_prok.output`)
    because that rule is defined later in this file, inside the
    `if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:` block — Snakemake resolves
    input functions at DAG-build time, so the literal string is safe even
    though the rule object isn't in scope yet at this point in the file.
    """
    if COASSEMBLY_VIRAL:
        return f"{OUTDIR}/coassembly/{wc.group}/prok_input/{wc.group}_contigs_nonviral.fasta"
    return f"{OUTDIR}/coassembly/{wc.group}/contigs.fa"


if not LONG_READS:

    rule megahit_coassembly:
        """
        MEGAHIT co-assembly of all samples in a group (short reads).
        Mirrors `rule megahit` (rules/assembly.smk): same env/container/flags,
        but inputs are the comma-joined reads of every sample in the group.
        PE: -1/-2 (comma-joined per-sample files); SE: -r (comma-joined).
        Writes the canonical {group}/contigs.fa.
        NOTE: -m is in bytes (MEGAHIT_MEM is defined in bytes in the config).
        """
        input:
            r1 = _group_r1,
            r2 = _group_r2,
        output:
            contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/megahit.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/megahit.tsv"
        conda:      "../envs/env_assembly.yaml"
        container:  CONTAINERS.get("megahit")
        threads: THREADS
        params:
            outdir     = f"{OUTDIR}/coassembly/{{group}}/megahit",
            preset     = (
                f"--presets {MEGAHIT_PRESET}" if MEGAHIT_PRESET in ("meta-sensitive", "meta-large")
                else (MEGAHIT_CUSTOM_PARAMS   if MEGAHIT_PRESET == "custom" else "")
            ),
            single_end = SINGLE_END,
        shell:
            """
            mkdir -p $(dirname {log})
            rm -rf {params.outdir}
            if [ "{params.single_end}" = "True" ]; then
                R1=$(echo {input.r1} | tr ' ' ',')
                megahit \
                    -r "$R1" \
                    -t {threads} \
                    -m {MEGAHIT_MEM} \
                    --min-contig-len {MIN_CONTIG} \
                    {params.preset} \
                    -o {params.outdir} \
                    > {log} 2>&1
            else
                R1=$(echo {input.r1} | tr ' ' ',')
                R2=$(echo {input.r2} | tr ' ' ',')
                megahit \
                    -1 "$R1" -2 "$R2" \
                    -t {threads} \
                    -m {MEGAHIT_MEM} \
                    --min-contig-len {MIN_CONTIG} \
                    {params.preset} \
                    -o {params.outdir} \
                    > {log} 2>&1
            fi
            cp {params.outdir}/final.contigs.fa {output.contigs}
            """


    rule coassembly_index:
        """
        Index the group's co-assembly contigs for BWA-MEM2.
        Mirrors `rule bwa_index` (rules/mapping.smk): same env/container/flags.
        """
        input:
            contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
        output:
            idx = f"{OUTDIR}/coassembly/{{group}}/mapping/contigs_index.bwt.2bit.64",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/bwa_index.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/bwa_index.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("bwa_mem2")
        params:
            prefix = f"{OUTDIR}/coassembly/{{group}}/mapping/contigs_index",
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p {OUTDIR}/coassembly/{wildcards.group}/mapping
            bwa-mem2 index -p {params.prefix} {input.contigs} > {log} 2>&1
            """


    rule coassembly_map:
        """
        Map one sample's cleaned short reads back to its group's co-assembly.
        Outputs an unsorted SAM (temp). Mirrors `rule bwa_mem`
        (rules/mapping.smk): same env/container/flags. SE maps R1 only; PE
        maps R1 + R2.
        """
        input:
            idx     = rules.coassembly_index.output.idx,
            tr1     = lambda wc: _clean_r1(type("W", (), {"sample": wc.sample})),
            tr2     = lambda wc: ([] if SINGLE_END
                                   else _clean_r2(type("W", (), {"sample": wc.sample}))),
        output:
            sam = temp(f"{OUTDIR}/coassembly/{{group}}/mapping/{{sample}}.sam"),
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/{{sample}}_bwa_mem.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/{{sample}}_bwa_mem.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("bwa_mem2")
        threads: THREADS
        params:
            prefix     = f"{OUTDIR}/coassembly/{{group}}/mapping/contigs_index",
            single_end = SINGLE_END,
        shell:
            """
            mkdir -p $(dirname {output.sam})
            mkdir -p $(dirname {log})
            if [ "{params.single_end}" = "True" ]; then
                bwa-mem2 mem \
                    -t {threads} \
                    {params.prefix} \
                    {input.tr1} \
                    > {output.sam} 2> {log}
            else
                bwa-mem2 mem \
                    -t {threads} \
                    {params.prefix} \
                    {input.tr1} {input.tr2} \
                    > {output.sam} 2> {log}
            fi
            """


    rule coassembly_sort:
        """
        Sort a sample's co-assembly SAM → BAM and index it.
        Mirrors `rule samtools_sort` (rules/mapping.smk): same env/container/flags.
        """
        input:
            sam = rules.coassembly_map.output.sam,
        output:
            bam = f"{OUTDIR}/coassembly/{{group}}/mapping/{{sample}}.sorted.bam",
            bai = f"{OUTDIR}/coassembly/{{group}}/mapping/{{sample}}.sorted.bam.bai",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/{{sample}}_samtools_sort.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/{{sample}}_samtools_sort.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("samtools")
        threads: THREADS
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p $(dirname {output.bam})
            samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
            samtools index {output.bam} 2>> {log}
            """


    rule coassembly_mapback:
        """
        Per-contig coverage depth for one sample against its group co-assembly
        (jgi_summarize_bam_contig_depths). Mirrors `rule calc_depth`
        (rules/mapping.smk): same env/container/flags. Consumed, per group, by
        `coassembly_abundance` to build the multi-sample VAMB abundance matrix.
        """
        input:
            bam = rules.coassembly_sort.output.bam,
        output:
            depth = f"{OUTDIR}/coassembly/{{group}}/mapping/{{sample}}_depth.txt",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/{{sample}}_calc_depth.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/{{sample}}_calc_depth.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("metabat2")
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p $(dirname {output.depth})
            jgi_summarize_bam_contig_depths \
                --outputDepth {output.depth} \
                {input.bam} \
                > {log} 2>&1
            """


    rule coassembly_abundance:
        """
        Join every group sample's per-contig depth into a single VAMB v5
        abundance matrix keyed by contigName. VAMB v5 --abundance_tsv expects
        `contigname\\t<s1>\\t<s2>...`; values are each sample's totalAvgDepth
        column (index 2) from jgi_summarize_bam_contig_depths output.
        """
        input:
            depths = lambda wc: expand(
                f"{OUTDIR}/coassembly/{wc.group}/mapping/{{sample}}_depth.txt",
                sample=GROUPS[wc.group]),
        output:
            matrix = f"{OUTDIR}/coassembly/{{group}}/vamb/abundance.tsv",
        run:
            import csv, os
            samples = GROUPS[wildcards.group]
            # jgi depth files: columns contigName, contigLen, totalAvgDepth, <bam>-var...
            cov = {}
            contig_order = []
            for s, path in zip(samples, input.depths):
                with open(path) as fh:
                    r = csv.reader(fh, delimiter="\t")
                    next(r)  # header
                    for row in r:
                        c = row[0]
                        if c not in cov:
                            cov[c] = {}
                            contig_order.append(c)
                        cov[c][s] = row[2]
            os.makedirs(os.path.dirname(output.matrix), exist_ok=True)
            with open(output.matrix, "w") as out:
                out.write("contigname\t" + "\t".join(samples) + "\n")
                for c in contig_order:
                    out.write(c + "\t" + "\t".join(cov[c].get(s, "0") for s in samples) + "\n")


    rule vamb_cobinning:
        """
        VAMB v5 co-binning on the group co-assembly (short reads).
        Mirrors the (now-removed) per-sample `rule vamb` (see git history,
        commit ad1f599): same env/container/flags, `vamb bin default`,
        --minfasta 200000, GPU flag from USE_GPU. KEY DIFFERENCE: the
        abundance TSV is the multi-sample matrix already built by
        `coassembly_abundance` (one column per group member) — used directly,
        no per-sample TSV is rebuilt here.
        VAMB creates its own output dir, so it is written under a `run/`
        subdir (not directly under `{group}/vamb/`) to avoid clobbering the
        sibling `{group}/vamb/abundance.tsv` matrix with `rm -rf`.
        Bin FASTAs consumed by downstream CheckM2 (Task 6) live at:
            {OUTDIR}/coassembly/{group}/vamb/run/bins/

        VIRAL FILTERING (Plan 4, Task 1): input contigs are the free-living-
        virus-filtered FASTA (`_coassembly_prok_contigs`, produced by
        `coassembly_filter_viral_for_prok`) when COASSEMBLY_VIRAL is enabled
        — matching per-sample behavior (`filter_viral_for_prok`). Because
        `coassembly_abundance`'s matrix spans ALL contigs (viral + non-viral)
        but VAMB requires --abundance_tsv rows to exactly match --fasta
        contigs, the abundance matrix is filtered down to the fasta's contig
        set inside the shell block before being passed to VAMB (mirrors the
        old per-sample `rule vamb`'s grep/awk filtering trick, see git show
        275340e:rules/prok_binning.smk). When COASSEMBLY_VIRAL is off, the
        fasta is the unfiltered contigs.fa (== abundance matrix contigs), so
        no filtering is needed.
        """
        input:
            contigs   = _coassembly_prok_contigs,
            abundance = rules.coassembly_abundance.output.matrix,
        output:
            done = f"{OUTDIR}/coassembly/{{group}}/vamb/done.txt",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/vamb_cobinning.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/vamb_cobinning.tsv"
        conda:      "../envs/env_binning.yaml"
        container:  CONTAINERS.get("vamb")
        threads: THREADS
        params:
            outdir       = f"{OUTDIR}/coassembly/{{group}}/vamb/run",
            low_depth    = LOW_DEPTH_MODE,
            viral_filter = COASSEMBLY_VIRAL,
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p $(dirname {params.outdir})
            if [ "{params.low_depth}" = "True" ]; then
                echo "[VAMB co-binning] Skipped -- low_depth_mode enabled" | tee {log}
                touch {output.done}; exit 0
            fi
            rm -rf {params.outdir}

            ABUNDANCE_TSV={input.abundance}
            if [ "{params.viral_filter}" = "True" ]; then
                # Filter the multi-sample abundance matrix down to the contigs
                # present in the (viral-filtered) FASTA — VAMB requires an
                # exact match between --fasta and --abundance_tsv contigs.
                grep "^>" {input.contigs} | sed 's/^>//' | awk '{{print $1}}' \
                    > {params.outdir}_contig_names.txt
                head -n1 {input.abundance} > {params.outdir}_abundance_filtered.tsv
                awk 'NR==FNR{{keep[$1]=1; next}} FNR>1 && ($1 in keep)' \
                    {params.outdir}_contig_names.txt {input.abundance} \
                    >> {params.outdir}_abundance_filtered.tsv
                rm -f {params.outdir}_contig_names.txt
                ABUNDANCE_TSV={params.outdir}_abundance_filtered.tsv
            fi

            CUDA_FLAG=""
            if [ "{USE_GPU}" = "True" ]; then CUDA_FLAG="--cuda"; fi
            vamb bin default \
                --outdir {params.outdir} \
                --fasta {input.contigs} \
                --abundance_tsv $ABUNDANCE_TSV \
                --minfasta 200000 \
                -p {threads} \
                $CUDA_FLAG \
                > {log} 2>&1
            touch {output.done}
            """


    rule checkm2_group:
        """
        CheckM2 — completeness and contamination of VAMB co-binning MAGs
        (group-level, short reads). Mirrors `rule checkm2` (rules/prok_binning.smk):
        same env/container/flags, universal ML model, MIMAG thresholds.
        NOTE: unlike the per-sample rule (Binette bins, `.fa`), VAMB v5 writes
        bin FASTAs as `.fna` — extension differs accordingly (-x fna). Confirm
        the VAMB output extension on a real run.
        """
        input:
            done = rules.vamb_cobinning.output.done,
        output:
            report = f"{OUTDIR}/coassembly/{{group}}/checkm2/quality_report.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/checkm2.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/checkm2.tsv"
        conda: "../envs/env_checkm2.yaml"
        container:  CONTAINERS.get("checkm2")
        threads: THREADS
        params:
            bins_dir = f"{OUTDIR}/coassembly/{{group}}/vamb/run/bins",
            outdir   = f"{OUTDIR}/coassembly/{{group}}/checkm2",
        shell:
            """
            mkdir -p $(dirname {log})
            rm -rf {params.outdir}
            mkdir -p {params.outdir}

            N_BINS=$(find {params.bins_dir} -maxdepth 1 -name "*.fna" 2>/dev/null | wc -l)
            if [ "$N_BINS" -eq 0 ]; then
                echo "[checkm2_group] No bins found — skipping" | tee {log}
                printf "Name\tCompleteness\tContamination\tGenome_Size\n" > {output.report}
                exit 0
            fi

            checkm2 predict \
                --threads {threads} \
                --input {params.bins_dir} \
                --output-directory {params.outdir} \
                -x fna \
                --database_path {CHECKM2_DB} \
                > {log} 2>&1 || echo "[checkm2_group] WARNING: predict failed" | tee -a {log}

            # Ensure output exists
            [ -f {output.report} ] || \
                printf "Name\tCompleteness\tContamination\tGenome_Size\n" > {output.report}
            """


    rule gtdbtk_group:
        """
        GTDB-Tk classify_wf — assign GTDB taxonomy to VAMB co-binning MAGs
        (group-level, short reads). Mirrors `rule gtdbtk` (rules/prok_binning.smk):
        same env/container/flags. No dereplication step for co-binning MAGs
        (galah_derep is per-sample only), so bins_dir is the raw VAMB bins dir.
        NOTE: extension is `.fna` (VAMB v5 output), unlike the per-sample rule's
        `.fa` (Binette/galah output). Confirm on a real run.
        Creates empty outputs if no bins are available (safe fallback).
        """
        input:
            report = rules.checkm2_group.output.report,
        output:
            done    = f"{OUTDIR}/coassembly/{{group}}/gtdbtk/done.txt",
            bac_tsv = f"{OUTDIR}/coassembly/{{group}}/gtdbtk/classify/gtdbtk.bac120.summary.tsv",
            ar_tsv  = f"{OUTDIR}/coassembly/{{group}}/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/gtdbtk.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/gtdbtk.tsv"
        conda: "../envs/env_gtdbtk.yaml"
        container:  CONTAINERS.get("gtdbtk")
        threads: THREADS
        params:
            bins_dir = f"{OUTDIR}/coassembly/{{group}}/vamb/run/bins",
            outdir   = f"{OUTDIR}/coassembly/{{group}}/gtdbtk",
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p {params.outdir}

            # If no bins, create empty outputs
            N_BINS=$(find {params.bins_dir} -maxdepth 1 -name "*.fna" 2>/dev/null | wc -l)
            if [ "$N_BINS" -eq 0 ]; then
                echo "[gtdbtk_group] No bins found — skipping" | tee {log}
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
                --extension  fna \
                >> {log} 2>&1 || echo "[gtdbtk_group] WARNING: classify_wf failed — creating empty outputs" | tee -a {log}

            # Always ensure output files exist — gtdbtk may fail if bins are low quality
            mkdir -p {params.outdir}/classify
            [ -f {output.bac_tsv} ] || printf "user_genome\tclassification\n" > {output.bac_tsv}
            [ -f {output.ar_tsv}  ] || printf "user_genome\tclassification\n" > {output.ar_tsv}
            touch {output.done}
            """


if LONG_READS:

    rule flye_coassembly:
        """
        metaFlye co-assembly of all samples in a group (long reads).
        Mirrors `rule flye_lr` (rules/assembly.smk): --meta graph simplification,
        read-type flag from LR_TECH / LR_ONT_CHEM. All group LR reads are passed
        as positional inputs (Flye accepts multiple read files). Writes the
        canonical {group}/contigs.fa.
        """
        input:
            reads = _group_lr,
        output:
            contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/flye.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/flye.tsv"
        conda:      "../envs/env_flye.yaml"
        container:  CONTAINERS.get("flye")
        threads: THREADS
        params:
            outdir  = f"{OUTDIR}/coassembly/{{group}}/flye",
            overlap = LR_FLYE_OVERLAP,
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p {params.outdir}
            if [ "{LR_TECH}" = "hifi" ]; then
                READ_FLAG="--pacbio-hifi"
            elif [ "{LR_ONT_CHEM}" = "hq" ]; then
                READ_FLAG="--nano-hq"
            else
                READ_FLAG="--nano-raw"
            fi
            flye $READ_FLAG {input.reads} \
                --out-dir {params.outdir} \
                --meta \
                --min-overlap {params.overlap} \
                --threads {threads} \
                --scaffold \
                --iterations 2 \
                >> {log} 2>&1 || true
            [ -f {params.outdir}/assembly.fasta ] && \
                cp {params.outdir}/assembly.fasta {output.contigs} || \
                touch {output.contigs}
            """


# ══════════════════════════════════════════════════════════════════════
#  VAMB multi-split co-binning (Task 5b, Plan 2)
#
#  Independent of co-assembly: concatenates every sample's per-sample
#  dereplicated assembly (rep_seq.fasta) into ONE global catalog (contigs
#  renamed with a sample prefix so bins stay traceable per-sample), maps
#  every sample's reads to that catalog, and runs VAMB on the resulting
#  multi-sample abundance matrix. Short reads only.
#
#  Everything lives under a fixed "multisplit" pseudo-group path so the
#  report tab (rules/report.smk `load_coassembly`) can reuse it exactly
#  like any other co-assembly group, with no loader changes required:
#      {OUTDIR}/coassembly/multisplit/...
# ══════════════════════════════════════════════════════════════════════

if COBINNING_MULTISPLIT and not LONG_READS:

    rule multisplit_catalog:
        """
        Concatenate every sample's dereplicated assembly (mmseqs rep_seq)
        into one global catalog for pooled VAMB co-binning. Headers are
        renamed `>X` -> `>S{sample}C{X}` so bins stay traceable back to the
        originating sample/contig. Pooled VAMB co-binning over the
        concatenated per-sample assembly catalog. (True per-sample
        bin-splitting via VAMB's --binsplit_separator is a possible
        refinement — verify VAMB v5 binsplit behavior on a real run before
        enabling.)
        """
        input:
            rep_seqs = expand(f"{OUTDIR}/{{sample}}/mmseqs/{{sample}}_rep_seq.fasta",
                               sample=list(SAMPLES.keys())),
        output:
            catalog = f"{OUTDIR}/coassembly/multisplit/catalog.fasta",
        log:
            f"{OUTDIR}/coassembly/multisplit/logs/catalog.log"
        run:
            import os
            os.makedirs(os.path.dirname(output.catalog), exist_ok=True)
            os.makedirs(os.path.dirname(log[0]), exist_ok=True)
            samples = list(SAMPLES.keys())
            with open(output.catalog, "w") as out, open(log[0], "w") as lg:
                for s, path in zip(samples, input.rep_seqs):
                    n = 0
                    with open(path) as fh:
                        for line in fh:
                            if line.startswith(">"):
                                header = line[1:].strip().split()[0]
                                out.write(f">S{s}C{header}\n")
                                n += 1
                            else:
                                out.write(line)
                    lg.write(f"{s}: {n} contigs\n")


    rule multisplit_index:
        """
        Index the multi-split catalog for BWA-MEM2.
        Mirrors `rule coassembly_index`: same env/container/flags.
        """
        input:
            catalog = rules.multisplit_catalog.output.catalog,
        output:
            idx = f"{OUTDIR}/coassembly/multisplit/mapping/contigs_index.bwt.2bit.64",
        log:
            f"{OUTDIR}/coassembly/multisplit/logs/bwa_index.log"
        benchmark:
            f"{OUTDIR}/coassembly/multisplit/benchmarks/bwa_index.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("bwa_mem2")
        params:
            prefix = f"{OUTDIR}/coassembly/multisplit/mapping/contigs_index",
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p {OUTDIR}/coassembly/multisplit/mapping
            bwa-mem2 index -p {params.prefix} {input.catalog} > {log} 2>&1
            """


    rule multisplit_map:
        """
        Map one sample's cleaned short reads to the multi-split catalog.
        Outputs an unsorted SAM (temp). Mirrors `rule coassembly_map`:
        same env/container/flags. SE maps R1 only; PE maps R1 + R2.
        """
        input:
            idx = rules.multisplit_index.output.idx,
            tr1 = lambda wc: _clean_r1(type("W", (), {"sample": wc.sample})),
            tr2 = lambda wc: ([] if SINGLE_END
                               else _clean_r2(type("W", (), {"sample": wc.sample}))),
        output:
            sam = temp(f"{OUTDIR}/coassembly/multisplit/mapping/{{sample}}.sam"),
        log:
            f"{OUTDIR}/coassembly/multisplit/logs/{{sample}}_bwa_mem.log"
        benchmark:
            f"{OUTDIR}/coassembly/multisplit/benchmarks/{{sample}}_bwa_mem.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("bwa_mem2")
        threads: THREADS
        params:
            prefix     = f"{OUTDIR}/coassembly/multisplit/mapping/contigs_index",
            single_end = SINGLE_END,
        shell:
            """
            mkdir -p $(dirname {output.sam})
            mkdir -p $(dirname {log})
            if [ "{params.single_end}" = "True" ]; then
                bwa-mem2 mem \
                    -t {threads} \
                    {params.prefix} \
                    {input.tr1} \
                    > {output.sam} 2> {log}
            else
                bwa-mem2 mem \
                    -t {threads} \
                    {params.prefix} \
                    {input.tr1} {input.tr2} \
                    > {output.sam} 2> {log}
            fi
            """


    rule multisplit_sort:
        """
        Sort a sample's multi-split SAM -> BAM and index it.
        Mirrors `rule coassembly_sort`: same env/container/flags.
        """
        input:
            sam = rules.multisplit_map.output.sam,
        output:
            bam = f"{OUTDIR}/coassembly/multisplit/mapping/{{sample}}.sorted.bam",
            bai = f"{OUTDIR}/coassembly/multisplit/mapping/{{sample}}.sorted.bam.bai",
        log:
            f"{OUTDIR}/coassembly/multisplit/logs/{{sample}}_samtools_sort.log"
        benchmark:
            f"{OUTDIR}/coassembly/multisplit/benchmarks/{{sample}}_samtools_sort.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("samtools")
        threads: THREADS
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p $(dirname {output.bam})
            samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
            samtools index {output.bam} 2>> {log}
            """


    rule multisplit_depth:
        """
        Per-contig coverage depth for one sample against the multi-split
        catalog (jgi_summarize_bam_contig_depths). Mirrors
        `rule coassembly_mapback`: same env/container/flags. Consumed by
        `multisplit_abundance` to build the multi-sample VAMB abundance
        matrix.
        """
        input:
            bam = rules.multisplit_sort.output.bam,
        output:
            depth = f"{OUTDIR}/coassembly/multisplit/mapping/{{sample}}_depth.txt",
        log:
            f"{OUTDIR}/coassembly/multisplit/logs/{{sample}}_calc_depth.log"
        benchmark:
            f"{OUTDIR}/coassembly/multisplit/benchmarks/{{sample}}_calc_depth.tsv"
        conda:      "../envs/env_mapping.yaml"
        container:  CONTAINERS.get("metabat2")
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p $(dirname {output.depth})
            jgi_summarize_bam_contig_depths \
                --outputDepth {output.depth} \
                {input.bam} \
                > {log} 2>&1
            """


    rule multisplit_abundance:
        """
        Join ALL samples' per-contig depth into a single VAMB v5 abundance
        matrix keyed by contigName. Mirrors `rule coassembly_abundance`,
        but spans every sample in the run rather than one group's members.
        """
        input:
            depths = expand(f"{OUTDIR}/coassembly/multisplit/mapping/{{sample}}_depth.txt",
                             sample=list(SAMPLES.keys())),
        output:
            matrix = f"{OUTDIR}/coassembly/multisplit/vamb/abundance.tsv",
        run:
            import csv, os
            samples = list(SAMPLES.keys())
            # jgi depth files: columns contigName, contigLen, totalAvgDepth, <bam>-var...
            cov = {}
            contig_order = []
            for s, path in zip(samples, input.depths):
                with open(path) as fh:
                    r = csv.reader(fh, delimiter="\t")
                    next(r)  # header
                    for row in r:
                        c = row[0]
                        if c not in cov:
                            cov[c] = {}
                            contig_order.append(c)
                        cov[c][s] = row[2]
            os.makedirs(os.path.dirname(output.matrix), exist_ok=True)
            with open(output.matrix, "w") as out:
                out.write("contigname\t" + "\t".join(samples) + "\n")
                for c in contig_order:
                    out.write(c + "\t" + "\t".join(cov[c].get(s, "0") for s in samples) + "\n")


    rule multisplit_vamb:
        """
        Pooled VAMB co-binning over the concatenated per-sample assembly
        catalog. (True per-sample bin-splitting via VAMB's
        --binsplit_separator is a possible refinement — verify VAMB v5
        binsplit behavior on a real run before enabling.)
        Mirrors `rule vamb_cobinning`: same env/container/flags, `vamb bin
        default`, --minfasta 200000, GPU flag from USE_GPU. The abundance
        TSV is the multi-sample matrix built by `multisplit_abundance`
        (one column per SAMPLE in the whole run). This is NOT true
        per-sample bin-splitting — the `S{sample}C` header prefix is for
        traceability only; bins can still mix contigs from multiple
        samples.
        Bin FASTAs consumed by downstream CheckM2/GTDB-Tk live at:
            {OUTDIR}/coassembly/multisplit/vamb/run/bins/
        """
        input:
            catalog   = rules.multisplit_catalog.output.catalog,
            abundance = rules.multisplit_abundance.output.matrix,
        output:
            done = f"{OUTDIR}/coassembly/multisplit/vamb/done.txt",
        log:
            f"{OUTDIR}/coassembly/multisplit/logs/vamb_multisplit.log"
        benchmark:
            f"{OUTDIR}/coassembly/multisplit/benchmarks/vamb_multisplit.tsv"
        conda:      "../envs/env_binning.yaml"
        container:  CONTAINERS.get("vamb")
        threads: THREADS
        params:
            outdir    = f"{OUTDIR}/coassembly/multisplit/vamb/run",
            low_depth = LOW_DEPTH_MODE,
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p $(dirname {params.outdir})
            if [ "{params.low_depth}" = "True" ]; then
                echo "[VAMB multi-split] Skipped -- low_depth_mode enabled" | tee {log}
                touch {output.done}; exit 0
            fi
            rm -rf {params.outdir}
            CUDA_FLAG=""
            if [ "{USE_GPU}" = "True" ]; then CUDA_FLAG="--cuda"; fi
            vamb bin default \
                --outdir {params.outdir} \
                --fasta {input.catalog} \
                --abundance_tsv {input.abundance} \
                --minfasta 200000 \
                -p {threads} \
                $CUDA_FLAG \
                > {log} 2>&1
            touch {output.done}
            """


    rule multisplit_checkm2:
        """
        CheckM2 — completeness and contamination of multi-split VAMB MAGs.
        Mirrors `rule checkm2_group`: same env/container/flags. VAMB v5
        writes bin FASTAs as `.fna`.
        """
        input:
            done = rules.multisplit_vamb.output.done,
        output:
            report = f"{OUTDIR}/coassembly/multisplit/checkm2/quality_report.tsv",
        log:
            f"{OUTDIR}/coassembly/multisplit/logs/checkm2.log"
        benchmark:
            f"{OUTDIR}/coassembly/multisplit/benchmarks/checkm2.tsv"
        conda: "../envs/env_checkm2.yaml"
        container:  CONTAINERS.get("checkm2")
        threads: THREADS
        params:
            bins_dir = f"{OUTDIR}/coassembly/multisplit/vamb/run/bins",
            outdir   = f"{OUTDIR}/coassembly/multisplit/checkm2",
        shell:
            """
            mkdir -p $(dirname {log})
            rm -rf {params.outdir}
            mkdir -p {params.outdir}

            N_BINS=$(find {params.bins_dir} -maxdepth 1 -name "*.fna" 2>/dev/null | wc -l)
            if [ "$N_BINS" -eq 0 ]; then
                echo "[multisplit_checkm2] No bins found — skipping" | tee {log}
                printf "Name\tCompleteness\tContamination\tGenome_Size\n" > {output.report}
                exit 0
            fi

            checkm2 predict \
                --threads {threads} \
                --input {params.bins_dir} \
                --output-directory {params.outdir} \
                -x fna \
                --database_path {CHECKM2_DB} \
                > {log} 2>&1 || echo "[multisplit_checkm2] WARNING: predict failed" | tee -a {log}

            # Ensure output exists
            [ -f {output.report} ] || \
                printf "Name\tCompleteness\tContamination\tGenome_Size\n" > {output.report}
            """


    rule multisplit_gtdbtk:
        """
        GTDB-Tk classify_wf — assign GTDB taxonomy to multi-split VAMB MAGs.
        Mirrors `rule gtdbtk_group`: same env/container/flags. No
        dereplication step (galah_derep is per-sample only), so bins_dir is
        the raw VAMB bins dir. Extension is `.fna` (VAMB v5 output).
        Creates empty outputs if no bins are available (safe fallback).
        """
        input:
            report = rules.multisplit_checkm2.output.report,
        output:
            done    = f"{OUTDIR}/coassembly/multisplit/gtdbtk/done.txt",
            bac_tsv = f"{OUTDIR}/coassembly/multisplit/gtdbtk/classify/gtdbtk.bac120.summary.tsv",
            ar_tsv  = f"{OUTDIR}/coassembly/multisplit/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
        log:
            f"{OUTDIR}/coassembly/multisplit/logs/gtdbtk.log"
        benchmark:
            f"{OUTDIR}/coassembly/multisplit/benchmarks/gtdbtk.tsv"
        conda: "../envs/env_gtdbtk.yaml"
        container:  CONTAINERS.get("gtdbtk")
        threads: THREADS
        params:
            bins_dir = f"{OUTDIR}/coassembly/multisplit/vamb/run/bins",
            outdir   = f"{OUTDIR}/coassembly/multisplit/gtdbtk",
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p {params.outdir}

            # If no bins, create empty outputs
            N_BINS=$(find {params.bins_dir} -maxdepth 1 -name "*.fna" 2>/dev/null | wc -l)
            if [ "$N_BINS" -eq 0 ]; then
                echo "[multisplit_gtdbtk] No bins found — skipping" | tee {log}
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
                --extension  fna \
                >> {log} 2>&1 || echo "[multisplit_gtdbtk] WARNING: classify_wf failed — creating empty outputs" | tee -a {log}

            # Always ensure output files exist — gtdbtk may fail if bins are low quality
            mkdir -p {params.outdir}/classify
            [ -f {output.bac_tsv} ] || printf "user_genome\tclassification\n" > {output.bac_tsv}
            [ -f {output.ar_tsv}  ] || printf "user_genome\tclassification\n" > {output.ar_tsv}
            touch {output.done}
            """


    # The "multisplit" pseudo-group path (coassembly/multisplit/...) shares its
    # output filename patterns with the per-group co-assembly rules above
    # (mapping/, vamb/, checkm2/, gtdbtk/ subpaths). `parse_groups` rejects a
    # real group literally named "multisplit" (see coassembly_groups.py), so
    # this collision can only happen when the `group` wildcard is otherwise
    # unconstrained (e.g. co-assembly is off, so GROUPS is empty and no
    # `group` wildcard_constraints is emitted). Disambiguate explicitly so
    # Snakemake doesn't raise AmbiguousRuleException in that case.
    ruleorder: multisplit_index > coassembly_index
    ruleorder: multisplit_map > coassembly_map
    ruleorder: multisplit_sort > coassembly_sort
    ruleorder: multisplit_depth > coassembly_mapback
    ruleorder: multisplit_abundance > coassembly_abundance
    ruleorder: multisplit_vamb > vamb_cobinning
    ruleorder: multisplit_checkm2 > checkm2_group
    ruleorder: multisplit_gtdbtk > gtdbtk_group


# ══════════════════════════════════════════════════════════════════════
#  Co-assembly viral detection + consensus (Plan 3, Task 1)
#
#  Mode-agnostic (contig-based): mirrors rules/viral_detection.smk exactly,
#  but runs on the group co-assembly contigs (megahit for SR, flye for LR)
#  instead of the per-sample mmseqs2 rep_seq.fasta. NOT gated by
#  `if not LONG_READS:` — the viral consumer runs for either mode, since
#  `{OUTDIR}/coassembly/{group}/contigs.fa` is produced by whichever
#  assembler is active.
# ══════════════════════════════════════════════════════════════════════

if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:

    rule coassembly_virsorter2:
        """
        VirSorter2 viral detection on the group co-assembly.
        Mirrors `rule virsorter2` (rules/viral_detection.smk): same
        env/container/flags (--include-groups, --min-score,
        --hallmark-required-on-short).
        """
        input:
            contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
        output:
            viral = f"{OUTDIR}/coassembly/{{group}}/viral/virsorter2/final-viral-combined.fa",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/virsorter2.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/virsorter2.tsv"
        conda: "../envs/env_viral.yaml"
        container:  CONTAINERS.get("virsorter")
        threads: THREADS
        shell:
            """
            virsorter run \
                -i {input.contigs} \
                -w {OUTDIR}/coassembly/{wildcards.group}/viral/virsorter2 \
                --db-dir {VS2_DB} \
                --include-groups dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae \
                --min-score 0.5 \
                --hallmark-required-on-short \
                -j {threads} all \
                > {log} 2>&1
            """


    rule coassembly_genomad:
        """
        GeNomad viral (and plasmid) detection on the group co-assembly.
        Mirrors `rule genomad` (rules/viral_detection.smk): same
        env/container/flags. Output: *_summary/*_virus_summary.tsv —
        original contig names, no renaming.
        """
        input:
            contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
        output:
            done = f"{OUTDIR}/coassembly/{{group}}/viral/genomad/done.txt",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/genomad.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/genomad.tsv"
        conda: "../envs/env_genomad.yaml"
        container:  CONTAINERS.get("genomad")
        threads: THREADS
        params:
            outdir = f"{OUTDIR}/coassembly/{{group}}/viral/genomad",
        shell:
            """
            rm -rf {params.outdir}
            # --disable-nn-classification bypasses the PyTorch NN (faster on CPU,
            # avoids execstack kernel warnings). With use_gpu: true the NN is enabled
            # and PyTorch auto-detects CUDA for a moderate speedup.
            NN_FLAG="--disable-nn-classification"
            if [ "{USE_GPU}" = "True" ]; then NN_FLAG=""; fi
            genomad end-to-end \
                {input.contigs} \
                {params.outdir} \
                {GENOMAD_DB} \
                --threads {threads} \
                --splits 8 \
                --min-score 0.7 \
                --enable-score-calibration \
                $NN_FLAG \
                > {log} 2>&1
            touch {output.done}
            """


    rule coassembly_vibrant:
        """
        VIBRANT viral detection on the group co-assembly. Mirrors
        `rule vibrant` (rules/viral_detection.smk): same env/container/flags.
        VIBRANT outputs to CWD — cd to outdir first, using absolute paths.
        """
        input:
            contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
        output:
            done = f"{OUTDIR}/coassembly/{{group}}/viral/vibrant/done.txt",
        log:   f"{OUTDIR}/coassembly/{{group}}/logs/vibrant.log"
        benchmark: f"{OUTDIR}/coassembly/{{group}}/benchmarks/vibrant.tsv"
        conda: "../envs/phage_vibrant.yaml"
        container:  CONTAINERS.get("vibrant")
        threads: THREADS
        params:
            outdir    = f"{OUTDIR}/coassembly/{{group}}/viral/vibrant",
            minlen    = MIN_CONTIG,
            db_dir    = f"{_VIBRANT_BASE}/databases",
            files_dir = f"{_VIBRANT_BASE}/files",
        shell:
            """
            mkdir -p {params.outdir}
            cp {input.contigs} {params.outdir}/input.fasta
            cd {params.outdir}
            VIBRANT_run.py \
                -i input.fasta \
                -f nucl \
                -t {threads} \
                -l {params.minlen} \
                -no_plot \
                -d {params.db_dir} \
                -m {params.files_dir} \
                > {log} 2>&1 || true
            touch {output.done}
            """


    rule coassembly_viral_consensus:
        """
        Integrate results from all 3 viral detection tools on the group
        co-assembly. Mirrors `rule viral_consensus` (rules/viral_detection.smk)
        exactly, operating on `{OUTDIR}/coassembly/{group}/contigs.fa` instead
        of the per-sample mmseqs2 rep_seq.fasta.

        CONTIG NAME NORMALIZATION:
        - VirSorter2 : appends ||full or ||lt0.5 — stripped with split("||")[0]
        - GeNomad    : uses original names directly; provirus "contig|prov_X_Y" → "contig"
        - VIBRANT    : uses original names directly

        CONSENSUS STRATEGIES (VIRAL_CONSENSUS_MODE):
        - "count"  : keep contigs called by >= MIN_VIRAL_TOOLS tools (default: 2 of 3)
        - "score"  : keep contigs with any tool score >= threshold
        - "hybrid" : count OR single high-confidence tool (score mode)

        OUTPUT:
          *_viral_consensus.fasta — confirmed viral contigs
          *_tool_support.tsv      — all contigs × tools matrix
        """
        input:
            contigs      = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
            vs2_done     = rules.coassembly_virsorter2.output.viral,
            genomad_done = rules.coassembly_genomad.output.done,
            vibrant_done = rules.coassembly_vibrant.output.done,
        output:
            fasta   = f"{OUTDIR}/coassembly/{{group}}/viral/consensus/{{group}}_viral_consensus.fasta",
            support = f"{OUTDIR}/coassembly/{{group}}/viral/consensus/{{group}}_tool_support.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/viral_consensus.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/viral_consensus.tsv"
        params:
            vs2_dir     = f"{OUTDIR}/coassembly/{{group}}/viral/virsorter2",
            genomad_dir = f"{OUTDIR}/coassembly/{{group}}/viral/genomad",
            vibrant_dir = f"{OUTDIR}/coassembly/{{group}}/viral/vibrant",
            outdir      = f"{OUTDIR}/coassembly/{{group}}/viral/consensus",
        run:
            import os, glob, csv
            from collections import defaultdict

            os.makedirs(params.outdir, exist_ok=True)

            def iter_fasta_names(path):
                names = set()
                if not os.path.exists(path):
                    return names
                with open(path) as f:
                    for line in f:
                        if line.startswith(">"):
                            names.add(line[1:].strip().split()[0])
                return names

            # ── VirSorter2 ────────────────────────────────────────────────
            vs2_names = {
                n.split("||")[0]
                for n in iter_fasta_names(
                    os.path.join(params.vs2_dir, "final-viral-combined.fa")
                )
            }

            # ── GeNomad ───────────────────────────────────────────────────
            # Summary TSV: {outdir}/{prefix}_summary/{prefix}_virus_summary.tsv
            # Column 'seq_name' = original contig name (no renaming)
            # Provirus names: "contig|provirus_X_Y" → normalize to "contig"
            genomad_names = set()
            genomad_prefix = os.path.basename(input.contigs).replace(".fasta","").replace(".fa","")
            genomad_summary = os.path.join(
                params.genomad_dir,
                f"{genomad_prefix}_summary",
                f"{genomad_prefix}_virus_summary.tsv"
            )
            if os.path.exists(genomad_summary):
                with open(genomad_summary) as f:
                    header = None
                    seq_col = 0
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        parts = line.split("\t")
                        if header is None:
                            header = [h.lower() for h in parts]
                            for c in ["seq_name", "sequence", "contig"]:
                                if c in header:
                                    seq_col = header.index(c)
                                    break
                            continue
                        if seq_col < len(parts):
                            raw = parts[seq_col].strip()
                            name = raw.split("|")[0] if "|" in raw else raw
                            genomad_names.add(name)

            # ── VIBRANT ───────────────────────────────────────────────────
            vibrant_names = set()
            for fa in (glob.glob(os.path.join(str(params.vibrant_dir), "**", "*phages_combined*"), recursive=True) +
                       glob.glob(os.path.join(str(params.vibrant_dir), "**", "VIBRANT_phages_*", "*.fna"), recursive=True)):
                try:
                    for l in open(fa):
                        if l.startswith(">"): vibrant_names.add(l[1:].strip().split()[0])
                except: pass

            # ── Count tool support per contig ─────────────────────────────
            tool_hits = defaultdict(list)
            for n in vs2_names:     tool_hits[n].append("VirSorter2")
            for n in genomad_names: tool_hits[n].append("GeNomad")
            for n in vibrant_names: tool_hits[n].append("VIBRANT")

            # ── Consensus strategy ────────────────────────────────────────
            def _safe_float(v, d=0.0):
                try: return float(v)
                except: return d

            if VIRAL_CONSENSUS_MODE in ("score", "hybrid"):
                high_conf = set()

                # VirSorter2: final-viral-score.tsv — column max_score (0-1)
                vs2_sf = os.path.join(params.vs2_dir, "final-viral-score.tsv")
                if os.path.exists(vs2_sf):
                    with open(vs2_sf) as _f:
                        for _r in csv.DictReader(_f, delimiter="\t"):
                            if _safe_float(_r.get("max_score",0)) >= SCORE_VS2_MIN:
                                high_conf.add(_r.get("seqname","").strip())
                for _gf in glob.glob(os.path.join(params.genomad_dir,"**","*_virus_summary.tsv"),recursive=True):
                    with open(_gf) as _f:
                        for _r in csv.DictReader(_f, delimiter="\t"):
                            if _safe_float(_r.get("virus_score",0)) >= SCORE_GENOMAD_MIN:
                                high_conf.add(_r.get("seq_name", "").strip())
                for _vf in glob.glob(os.path.join(params.vibrant_dir, "**", "VIBRANT_genome_quality_*.tsv"), recursive=True):
                    try:
                        with open(_vf) as _f:
                            for _r in csv.DictReader(_f, delimiter="\t"):
                                if _r.get("type", "").lower() in ("lytic", "lysogenic"):
                                    high_conf.add(_r.get("scaffold", "").strip())
                    except Exception: pass
            else:
                high_conf = set()

            if VIRAL_CONSENSUS_MODE == "count":
                consensus = {n: t for n, t in tool_hits.items() if len(t) >= MIN_VIRAL_TOOLS}
            elif VIRAL_CONSENSUS_MODE == "score":
                consensus = {n: t for n, t in tool_hits.items() if n in high_conf}
            elif VIRAL_CONSENSUS_MODE == "hybrid":
                consensus = {n: t for n, t in tool_hits.items()
                             if len(t) >= MIN_VIRAL_TOOLS or n in high_conf}
            else:
                consensus = {n: t for n, t in tool_hits.items() if len(t) >= MIN_VIRAL_TOOLS}

            # ── Write tool support table ──────────────────────────────────
            with open(output.support, "w") as tsv:
                tsv.write("contig\tn_tools\ttools\n")
                for n, t in sorted(tool_hits.items(), key=lambda x: -len(x[1])):
                    tsv.write(f"{n}\t{len(t)}\t{','.join(sorted(t))}\n")

            # ── Extract FASTA sequences ───────────────────────────────────
            kept = 0
            with open(input.contigs) as fin, open(output.fasta, "w") as fout:
                header, seq = "", []
                for line in fin:
                    line = line.rstrip()
                    if line.startswith(">"):
                        if header and header in consensus:
                            fout.write(">" + header + "\n" + "".join(seq) + "\n")
                            kept += 1
                        header = line[1:].strip().split()[0]
                        seq = []
                    else:
                        seq.append(line)
                if header and header in consensus:
                    fout.write(">" + header + "\n" + "".join(seq) + "\n")
                    kept += 1

            # ── Summary log ───────────────────────────────────────────────
            with open(log[0], "a") as lf:
                lf.write(f"\nVirSorter2 : {len(vs2_names)}\n")
                lf.write(f"GeNomad    : {len(genomad_names)}\n")
                lf.write(f"VIBRANT    : {len(vibrant_names)}\n")
                lf.write(f"Union total: {len(tool_hits)}\n")
                lf.write(f"Consensus mode={VIRAL_CONSENSUS_MODE} (min_tools={MIN_VIRAL_TOOLS}): {len(consensus)}\n")
                lf.write(f"FASTA output: {kept} → {output.fasta}\n")

            print(f"[coassembly_viral_consensus] {wildcards.group}: {kept} consensus viral contigs")


    rule coassembly_checkv:
        """
        CheckV — quality assessment of consensus viral contigs on the group
        co-assembly. Mirrors `rule checkv` (rules/viral_binning.smk), operating
        on the group viral consensus fasta instead of the per-sample
        (COBRA-extended | viral_consensus) input.
        Classifies: Complete / High-quality / Medium-quality / Low-quality / Not-determined.
        NOTE: always removes output dir before running — CheckV skips gene calling
        if it finds existing files, causing KeyError when contig names changed.
        Outputs viruses.fna (complete viruses) and proviruses.fna (trimmed provirus
        regions — host DNA flanks removed) for use in viral_nonredundant.
        """
        input:
            viral = rules.coassembly_viral_consensus.output.fasta,
        output:
            summary   = f"{OUTDIR}/coassembly/{{group}}/viral/checkv/quality_summary.tsv",
            viruses   = f"{OUTDIR}/coassembly/{{group}}/viral/checkv/viruses.fna",
            proviruses = f"{OUTDIR}/coassembly/{{group}}/viral/checkv/proviruses.fna",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/checkv.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/checkv.tsv"
        conda: "../envs/env_viral.yaml"
        container:  CONTAINERS.get("checkv")
        threads: THREADS
        shell:
            """
            rm -rf {OUTDIR}/coassembly/{wildcards.group}/viral/checkv
            checkv end_to_end \
                {input.viral} \
                {OUTDIR}/coassembly/{wildcards.group}/viral/checkv \
                -d {CHECKV_DB} \
                -t {threads} \
                > {log} 2>&1
            touch {output.viruses} {output.proviruses}
            """


    rule coassembly_viral_trimmed:
        """
        CheckV-trim proviruses in the group co-assembly viral consensus set.
        Mirrors the CheckV-trim portion of `rule viral_nonredundant`
        (rules/viral_binning.smk) — NOT the multi-assembler dedup / bins-first
        / MMseqs parts, which don't apply at the co-assembly level (there is
        no per-sample vRhyme binning input here; the group vOTU chain does its
        own skani-based clustering downstream).

        For every contig in the consensus fasta: if CheckV trimmed it
        (provirus), the trimmed sequence from viruses.fna / proviruses.fna
        replaces the original so that host DNA flanking proviruses is always
        removed before vOTU clustering. Sequences not processed by CheckV
        (rare, e.g. below CheckV internal length threshold) fall back to the
        original. One contig may yield >1 trimmed entry when CheckV found
        multiple provirus regions within it.
        """
        input:
            consensus         = rules.coassembly_viral_consensus.output.fasta,
            checkv_viruses    = rules.coassembly_checkv.output.viruses,
            checkv_proviruses = rules.coassembly_checkv.output.proviruses,
        output:
            fasta = f"{OUTDIR}/coassembly/{{group}}/viral/checkv/{{group}}_viral_trimmed.fasta",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/viral_trimmed.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/viral_trimmed.tsv"
        run:
            import os
            from collections import defaultdict

            # Build CheckV trimmed dict: orig_contig_id -> list of (header_line, seq_lines)
            # viruses.fna: header = original ID (no trimming needed, but use this file to
            #   be consistent — avoids re-reading the consensus for non-provirus sequences)
            # proviruses.fna: header = "orig_id|start_end" (trimmed region only)
            trimmed = defaultdict(list)
            for fna_path, is_provirus in [
                (str(input.checkv_viruses), False),
                (str(input.checkv_proviruses), True),
            ]:
                if not os.path.exists(fna_path) or os.path.getsize(fna_path) == 0:
                    continue
                curr_hdr, curr_seq = None, []
                with open(fna_path) as fh:
                    for line in fh:
                        if line.startswith('>'):
                            if curr_hdr:
                                hdr_id = curr_hdr[1:].split()[0]
                                orig = hdr_id.rsplit('|', 1)[0] if is_provirus else hdr_id
                                trimmed[orig].append((curr_hdr, curr_seq))
                            curr_hdr = line; curr_seq = []
                        else:
                            curr_seq.append(line)
                    if curr_hdr:
                        hdr_id = curr_hdr[1:].split()[0]
                        orig = hdr_id.rsplit('|', 1)[0] if is_provirus else hdr_id
                        trimmed[orig].append((curr_hdr, curr_seq))

            # Fallback: original sequences for anything CheckV didn't process
            orig_seqs = {}
            with open(str(input.consensus)) as fh:
                curr_hdr, curr_seq = None, []
                for line in fh:
                    if line.startswith('>'):
                        if curr_hdr:
                            orig_seqs[curr_hdr[1:].split()[0]] = (curr_hdr, curr_seq)
                        curr_hdr = line; curr_seq = []
                    else:
                        curr_seq.append(line)
                if curr_hdr:
                    orig_seqs[curr_hdr[1:].split()[0]] = (curr_hdr, curr_seq)

            def emit(name, out_lines):
                if name in trimmed:
                    for hdr, seq in trimmed[name]:
                        out_lines.append(hdr)
                        out_lines.extend(seq)
                elif name in orig_seqs:
                    hdr, seq = orig_seqs[name]
                    out_lines.append(hdr)
                    out_lines.extend(seq)

            out_lines = []
            for name in orig_seqs:
                emit(name, out_lines)

            with open(str(output.fasta), 'w') as fh:
                fh.writelines(out_lines)

            n_total = sum(1 for l in out_lines if l.startswith('>'))
            n_trimmed = sum(1 for n in orig_seqs if n in trimmed)
            with open(str(log[0]), 'w') as lf:
                lf.write(f'CheckV-trimmed sequences: {n_trimmed}\n')
                lf.write(f'Total trimmed set: {n_total}\n')


    rule coassembly_filter_viral_for_prok:
        """
        Remove free-living viral contigs from the co-binning (VAMB) input —
        group-level mirror of the per-sample `filter_viral_for_prok`
        (rules/prok_binning.smk:40). Fixes a Plan 2 gap: `vamb_cobinning`
        previously binned the FULL group co-assembly (`contigs.fa`),
        including free-living viral contigs, so group MAGs could be
        virus-contaminated.

        Strategy (identical to the per-sample rule):
          1. coassembly_viral_consensus fasta        → set of viral contigs
          2. CheckV `provirus=Yes` + GeNomad `|provirus_` suffix → provirus set
          3. remove = viral_consensus MINUS provirus
          4. {group}/contigs.fa MINUS remove → {group}_contigs_nonviral.fasta

        Provirus-bearing contigs always stay in the prok input (the provirus
        region is integrated within a bacterial host contig — removing it
        would drop the host genome with the prophage). Free-living viruses
        (no host chromosome context) are removed to clean up MAGs.

        Only defined when COASSEMBLY_VIRAL is enabled — it requires the
        group's viral detection outputs (consensus/CheckV/GeNomad). When
        COASSEMBLY_VIRAL is off, `_coassembly_prok_contigs` falls back to the
        unfiltered `contigs.fa` directly (no filter possible without viral
        detection).
        """
        input:
            contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
            viral   = rules.coassembly_viral_consensus.output.fasta,
            checkv  = rules.coassembly_checkv.output.summary,
            genomad = rules.coassembly_genomad.output.done,
        output:
            filtered = f"{OUTDIR}/coassembly/{{group}}/prok_input/{{group}}_contigs_nonviral.fasta",
            stats    = f"{OUTDIR}/coassembly/{{group}}/prok_input/filter_stats.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/filter_viral_for_prok.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/filter_viral_for_prok.tsv"
        params:
            genomad_dir = lambda wc: f"{OUTDIR}/coassembly/{wc.group}/viral/genomad",
        run:
            import os, glob, csv

            viral_set = set()
            with open(input.viral) as f:
                for line in f:
                    if line.startswith(">"):
                        viral_set.add(line[1:].strip().split()[0])

            # Always keep provirus-bearing contigs in prok binning
            provirus_set = set()
            try:
                with open(input.checkv) as f:
                    rdr = csv.DictReader(f, delimiter="\t")
                    for row in rdr:
                        if (row.get("provirus", "") or "").strip().lower() == "yes":
                            cid = (row.get("contig_id", "") or "").strip()
                            if cid:
                                provirus_set.add(cid)
            except Exception:
                pass

            for gf in glob.glob(os.path.join(str(params.genomad_dir),
                                             "**", "*_virus_summary.tsv"),
                                recursive=True):
                try:
                    with open(gf) as f:
                        rdr = csv.DictReader(f, delimiter="\t")
                        for row in rdr:
                            seq = (row.get("seq_name", "") or "").strip()
                            if "|provirus_" in seq:
                                provirus_set.add(seq.split("|")[0])
                except Exception:
                    pass

            # Contigs we will remove from the prok input
            remove_set = viral_set - provirus_set

            os.makedirs(os.path.dirname(output.filtered), exist_ok=True)
            kept = removed = total = 0
            with open(input.contigs) as fin, open(output.filtered, "w") as fout:
                write = False
                header = None
                for line in fin:
                    if line.startswith(">"):
                        if header is not None:
                            # finalize previous record (already written if write=True)
                            pass
                        header = line[1:].strip().split()[0]
                        total += 1
                        if header in remove_set:
                            write = False
                            removed += 1
                        else:
                            write = True
                            kept += 1
                            fout.write(line)
                    else:
                        if write:
                            fout.write(line)

            with open(output.stats, "w") as s:
                s.write("metric\tcount\n")
                s.write(f"total_contigs\t{total}\n")
                s.write(f"viral_total\t{len(viral_set)}\n")
                s.write(f"provirus_kept\t{len(provirus_set)}\n")
                s.write(f"removed\t{removed}\n")
                s.write(f"kept_for_prok\t{kept}\n")

            with open(log[0], "w") as lf:
                lf.write(f"Total contigs        : {total}\n")
                lf.write(f"Viral consensus      : {len(viral_set)}\n")
                lf.write(f"Provirus (kept)      : {len(provirus_set)}\n")
                lf.write(f"Removed (viral-only) : {removed}\n")
                lf.write(f"Kept for prok binning: {kept}\n")


    rule coassembly_skani_votu:
        """
        skani triangle — pairwise ANI matrix for group co-assembly viral genomes.
        Mirrors `rule skani_votu` (rules/viral_binning.smk), operating on the
        group consensus fasta (rules.coassembly_viral_consensus) instead of the
        per-sample viral_nonredundant dedup — vOTU clustering itself dereplicates
        by ANI, so the consensus + checkv set is a consistent contig set.
        Clustering is done by the downstream coassembly_skani_cluster rule (pure
        Python, runs in Snakemake's own interpreter — no container needed).
        """
        input:
            fasta = rules.coassembly_viral_trimmed.output.fasta,
        output:
            ani = f"{OUTDIR}/coassembly/{{group}}/viral/votu/skani_ani.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/skani_votu.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/skani_votu.tsv"
        conda: "../envs/env_derep.yaml"
        container: CONTAINERS.get("skani")
        threads: THREADS
        params:
            outdir  = f"{OUTDIR}/coassembly/{{group}}/viral/votu",
            enabled = VOTU_CLUSTERING_ENABLED,
        shell:
            """
            mkdir -p {params.outdir}
            if [ "{params.enabled}" != "True" ]; then
                echo "[skani_votu] Disabled via config" | tee {log}
                printf "qname\trname\tani\taf_q\taf_r\n" > {output.ani}
                exit 0
            fi
            N_SEQ=$(grep -c '^>' {input.fasta} 2>/dev/null || echo 0)
            if [ "$N_SEQ" -eq 0 ]; then
                echo "[skani_votu] Empty viral set" | tee {log}
                printf "qname\trname\tani\taf_q\taf_r\n" > {output.ani}
                exit 0
            fi
            skani triangle \
                -i {input.fasta} \
                -o {output.ani} \
                -t {threads} \
                --slow \
                >> {log} 2>&1 || echo "[skani_votu] WARNING: triangle failed" | tee -a {log}
            """


    rule coassembly_skani_cluster:
        """
        Greedy single-linkage vOTU clustering from the skani ANI matrix, for the
        group co-assembly viral set. Mirrors `rule skani_cluster`
        (rules/viral_binning.smk). Pure Python (stdlib only) — runs in
        Snakemake's interpreter, no container needed.
        ICTV / Roux 2019 definition: ANI >= VOTU_ANI AND max(af_q, af_r) >= VOTU_AF.
        The cluster representative is the member with the highest CheckV
        completeness (ties broken by FASTA order) — not simply the first contig
        encountered, since that's an arbitrary assembly-order artifact and the
        representative's sequence/length is what downstream genome maps and
        vOTU abundance use.
        """
        input:
            ani    = rules.coassembly_skani_votu.output.ani,
            fasta  = rules.coassembly_viral_trimmed.output.fasta,
            checkv = rules.coassembly_checkv.output.summary,
        output:
            clusters = f"{OUTDIR}/coassembly/{{group}}/viral/votu/vOTU_clusters.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/skani_cluster.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/skani_cluster.tsv"
        params:
            ani_min = VOTU_ANI,
            af_min  = VOTU_AF,
            enabled = VOTU_CLUSTERING_ENABLED,
        run:
            import csv
            import sys
            with open(log[0], "w") as _lf:
                if not params.enabled:
                    _lf.write("[skani_cluster] Disabled via config\n")
                    with open(output.clusters, "w") as f:
                        f.write("representative\tmember\n")
                else:
                    ids = []
                    with open(input.fasta) as f:
                        for line in f:
                            if line.startswith(">"):
                                ids.append(line[1:].strip().split()[0])

                    completeness = {}
                    with open(input.checkv) as f:
                        for row in csv.DictReader(f, delimiter="\t"):
                            try:
                                completeness[row["contig_id"]] = float(row.get("completeness", "0") or 0)
                            except (KeyError, ValueError):
                                continue

                    neigh = {i: set() for i in ids}
                    try:
                        with open(input.ani) as f:
                            f.readline()
                            for line in f:
                                parts = line.rstrip("\n").split("\t")
                                if len(parts) < 5:
                                    continue
                                q, r = parts[0], parts[1]
                                try:
                                    ani = float(parts[2])
                                    afq = float(parts[3])
                                    afr = float(parts[4])
                                except ValueError:
                                    continue
                                if ani >= params.ani_min and max(afq, afr) >= params.af_min:
                                    if q in neigh and r in neigh:
                                        neigh[q].add(r)
                                        neigh[r].add(q)
                    except FileNotFoundError:
                        pass

                    seen = set()
                    clusters = []
                    for n in ids:
                        if n in seen:
                            continue
                        comp = []
                        stack = [n]
                        while stack:
                            x = stack.pop()
                            if x in seen:
                                continue
                            seen.add(x)
                            comp.append(x)
                            for y in neigh.get(x, ()):
                                if y not in seen:
                                    stack.append(y)
                        clusters.append(comp)

                    with open(output.clusters, "w") as f:
                        f.write("representative\tmember\n")
                        for comp in clusters:
                            # Representative = highest CheckV completeness (ties -> FASTA order)
                            rep = max(comp, key=lambda m: (completeness.get(m, 0.0), -comp.index(m)))
                            for m in comp:
                                f.write(f"{rep}\t{m}\n")

                    msg = (f"[skani_cluster] genomes={len(ids)} clusters={len(clusters)} "
                           f"ani>={params.ani_min} af>={params.af_min}\n")
                    _lf.write(msg)
                    print(msg, end="")


    rule coassembly_viral_votu_reps:
        """
        Extract vOTU representative sequences from the group co-assembly viral
        consensus fasta and produce quality-filtered subsets used by downstream
        analyses. Mirrors `rule viral_votu_reps` (rules/viral_binning.smk).

        Outputs:
          all_fasta      — all representatives (one per vOTU cluster), for reporting
          mq_fasta       — MQ+ (Complete/HQ/MQ or completeness>=50%) representatives,
                           for prodigal_viral, PHIST, Pharokka, genome maps
          hq_10kb_fasta  — HQ+/Complete AND >= 10 kb representatives, for vConTACT3

        Why only representatives?
          skani_cluster picks the highest-completeness member per cluster as the
          canonical genome. Running taxonomy/PHIST/etc. on all members would inflate
          compute by ~20–40% and produce duplicate annotations for near-identical
          sequences, with no scientific gain.
        """
        input:
            fasta    = rules.coassembly_viral_trimmed.output.fasta,
            clusters = rules.coassembly_skani_cluster.output.clusters,
            checkv   = rules.coassembly_checkv.output.summary,
        output:
            all_fasta     = f"{OUTDIR}/coassembly/{{group}}/viral/votu/votu_all_reps.fasta",
            mq_fasta      = f"{OUTDIR}/coassembly/{{group}}/viral/votu/votu_mq_reps.fasta",
            hq_10kb_fasta = f"{OUTDIR}/coassembly/{{group}}/viral/votu/votu_hq_10kb_reps.fasta",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/viral_votu_reps.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/viral_votu_reps.tsv"
        params:
            hq_min_len = 10000,
        run:
            import csv, os

            # Collect unique representatives from cluster TSV
            reps = set()
            with open(str(input.clusters)) as fh:
                for row in csv.DictReader(fh, delimiter='\t'):
                    reps.add(row['representative'])

            # CheckV quality sets for filtering
            mq_plus = set()   # Complete / HQ / MQ or completeness >= 50%
            hq_plus = set()   # Complete / HQ only
            with open(str(input.checkv)) as fh:
                for row in csv.DictReader(fh, delimiter='\t'):
                    cid = row.get('contig_id', '').strip()
                    q   = row.get('checkv_quality', '').strip()
                    try:
                        comp = float(row.get('completeness', '0') or 0)
                    except (ValueError, TypeError):
                        comp = 0.0
                    if q in ('Complete', 'High-quality', 'Medium-quality') or comp >= 50:
                        mq_plus.add(cid)
                    if q in ('Complete', 'High-quality'):
                        hq_plus.add(cid)

            # Parse viral consensus fasta into dict
            seqs = {}
            with open(str(input.fasta)) as fh:
                curr_hdr, curr_seq = None, []
                for line in fh:
                    if line.startswith('>'):
                        if curr_hdr:
                            seqs[curr_hdr[1:].split()[0]] = (curr_hdr, curr_seq)
                        curr_hdr = line; curr_seq = []
                    else:
                        curr_seq.append(line)
                if curr_hdr:
                    seqs[curr_hdr[1:].split()[0]] = (curr_hdr, curr_seq)

            def write_filtered(names, path, length_min=0):
                written = 0
                with open(path, 'w') as out:
                    for name in names:
                        if name not in seqs:
                            continue
                        hdr, seq_lines = seqs[name]
                        if length_min > 0:
                            length = sum(len(s.rstrip()) for s in seq_lines)
                            if length < length_min:
                                continue
                        out.write(hdr)
                        out.writelines(seq_lines)
                        written += 1
                return written

            n_all      = write_filtered(reps, str(output.all_fasta))
            mq_reps    = reps & mq_plus
            n_mq       = write_filtered(mq_reps, str(output.mq_fasta))
            hq_reps    = reps & hq_plus
            n_hq_10kb  = write_filtered(hq_reps, str(output.hq_10kb_fasta),
                                        length_min=params.hq_min_len)

            with open(str(log[0]), 'w') as lf:
                lf.write(f'[viral_votu_reps] Total vOTU reps: {n_all}\n')
                lf.write(f'[viral_votu_reps] MQ+ reps (taxonomy/PHIST/Pharokka): {n_mq}\n')
                lf.write(f'[viral_votu_reps] HQ+/>=10kb reps (vConTACT3): {n_hq_10kb}\n')


    rule coassembly_prodigal_viral:
        """
        Predict ORFs from the group vOTU MQ+ representatives for taxonomy
        searches. Mirrors `rule prodigal_viral` (rules/taxonomy.smk).
        """
        input:
            viral = rules.coassembly_viral_votu_reps.output.mq_fasta,
        output:
            faa  = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/viral_proteins.faa",
            done = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/prodigal_done.txt",
        log:   f"{OUTDIR}/coassembly/{{group}}/logs/prodigal_viral.log"
        benchmark: f"{OUTDIR}/coassembly/{{group}}/benchmarks/prodigal_viral.tsv"
        conda: "../envs/env_viral.yaml"
        container:  CONTAINERS.get("prodigal")
        threads: 1
        shell:
            """
            mkdir -p $(dirname {output.faa})
            if [ ! -s {input.viral} ]; then
                touch {output.faa} {output.done}; exit 0
            fi
            prodigal -i {input.viral} -a {output.faa} -p meta -f gff > {log} 2>&1
            touch {output.done}
            """


    rule coassembly_mmseqs_taxonomy_viral:
        """
        Real per-query LCA against INPHARED-derived seqTaxDB on the group
        vOTU proteins. Mirrors `rule mmseqs_taxonomy_viral` (rules/taxonomy.smk)
        exactly -- same seqTaxDB path, same skip guards. vConTACT3 and
        custom-MMseqs sources are out of scope for the co-assembly viral
        taxonomy core; only GeNomad + this rule feed coassembly_viral_taxonomy.
        """
        input:
            faa  = rules.coassembly_prodigal_viral.output.faa,
            done = rules.coassembly_prodigal_viral.output.done,
        output:
            hits = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/mmseqs_vs_inphared.tsv",
            done = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/mmseqs_inphared_done.txt",
        log:   f"{OUTDIR}/coassembly/{{group}}/logs/mmseqs_taxonomy_viral.log"
        benchmark: f"{OUTDIR}/coassembly/{{group}}/benchmarks/mmseqs_taxonomy_viral.tsv"
        conda: "../envs/env_assembly.yaml"
        container:  CONTAINERS.get("mmseqs2")
        threads: THREADS
        params:
            seqtaxdb = f"{INPHARED_DB}/inphared_mmseqs_taxdb/seqTaxDB",
            outdir   = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/mmseqs_inphared",
            querydb  = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/mmseqs_inphared/queryDB",
            result   = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/mmseqs_inphared/result",
            tmp      = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/mmseqs_inphared/tmp",
        run:
            import os
            from pathlib import Path

            os.makedirs(params.outdir, exist_ok=True)
            header = "qseqid\ttaxid\trank\tname\tlineage\n"

            def write_empty(msg):
                with open(str(log[0]), "a") as lf:
                    lf.write(msg + "\n")
                Path(str(output.hits)).write_text(header)
                Path(str(output.done)).touch()

            if not os.path.exists(str(params.seqtaxdb) + ".dbtype"):
                write_empty(
                    "[coassembly_mmseqs_taxonomy_viral] No seqTaxDB at " + str(params.seqtaxdb) +
                    " -- run scripts/prepare_mmseqs_taxdb.py --format inphared once first (see INSTALL.md). Skipping."
                )
                return

            if not os.path.exists(str(input.faa)) or os.path.getsize(str(input.faa)) == 0:
                write_empty("[coassembly_mmseqs_taxonomy_viral] No viral proteins -- skipping")
                return

            # mmseqs taxonomy refuses to run if its output DB already exists
            # ("result.dbtype exists already!") -- same fix as mmseqs_taxonomy_viral.
            shell("rm -rf {params.tmp} {params.result}*; mkdir -p {params.tmp}")
            shell("mmseqs createdb {input.faa} {params.querydb} >> {log} 2>&1")
            shell(
                "mmseqs taxonomy {params.querydb} {params.seqtaxdb} {params.result} {params.tmp} "
                "--threads {threads} --tax-lineage 1 >> {log} 2>&1"
            )
            shell(
                "mmseqs createtsv {params.querydb} {params.result} {output.hits}.raw >> {log} 2>&1"
            )
            Path(str(output.hits)).write_text(header)
            if os.path.exists(str(output.hits) + ".raw"):
                with open(str(output.hits) + ".raw") as f, open(str(output.hits), "a") as out:
                    out.writelines(f)
                os.remove(str(output.hits) + ".raw")
            Path(str(output.done)).touch()


    rule coassembly_viral_taxonomy:
        """
        Merge taxonomy from GeNomad + MMseqs2/INPHARED into one table per
        contig, on the group co-assembly vOTU representatives. Mirrors
        `rule viral_taxonomy` (rules/taxonomy.smk), but scoped to only the
        two sources that are in-scope for co-assembly viral taxonomy --
        vConTACT3 and custom-MMseqs are NOT wired here. With those two
        sources absent, the effective priority becomes
        mmseqs_inphared > genomad (still resolved by deepest-rank-wins,
        same as the per-sample rule; ties just have fewer contenders here).
        """
        input:
            genomad_done = rules.coassembly_genomad.output.done,
            mmseqs_hits  = rules.coassembly_mmseqs_taxonomy_viral.output.hits,
            mmseqs_done  = rules.coassembly_mmseqs_taxonomy_viral.output.done,
            viral        = rules.coassembly_viral_votu_reps.output.mq_fasta,
        output:
            tsv  = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/viral_taxonomy_merged.tsv",
            done = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/taxonomy_done.txt",
        log:   f"{OUTDIR}/coassembly/{{group}}/logs/viral_taxonomy.log"
        benchmark: f"{OUTDIR}/coassembly/{{group}}/benchmarks/viral_taxonomy.tsv"
        conda: "../envs/env_viral.yaml"
        container:  CONTAINERS.get("diamond")
        threads: 1
        run:
            import csv, os, glob, collections
            from pathlib import Path

            lf = open(str(log[0]), "w")

            # All viral contigs
            contigs = []
            if os.path.exists(str(input.viral)):
                with open(str(input.viral)) as f:
                    for line in f:
                        if line.startswith(">"): contigs.append(line[1:].split()[0])
            lf.write(f"Total viral contigs: {len(contigs)}\n")

            # vConTACT3 and custom-MMseqs are out of scope for co-assembly
            # viral taxonomy -- these sources are always empty here.
            vc3_tax = {}
            custom_tax = {}

            # ── MMseqs2/INPHARED (real per-query LCA) ───────────────────
            _RANKS = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus']
            mmseqs_tax = _mmseqs_lca_rollup(str(input.mmseqs_hits), _RANKS)
            lf.write(f"MMseqs2/INPHARED: {len(mmseqs_tax)} contigs\n")

            # ── GeNomad ───────────────────────────────────────────────────
            genomad_tax = {}
            gdir   = os.path.dirname(str(input.genomad_done))
            gfiles = glob.glob(os.path.join(gdir, "**", "*_virus_summary.tsv"), recursive=True)
            gpath  = gfiles[0] if gfiles else ""
            if gpath and os.path.getsize(gpath) > 0:
                with open(gpath) as f:
                    for row in csv.DictReader(f, delimiter="\t"):
                        name  = row.get("seq_name","")
                        tax   = row.get("taxonomy","")
                        score = row.get("virus_score","0")
                        if not name or not tax: continue
                        parts = [p.strip() for p in tax.split(";")
                                 if p.strip() and p.strip() not in ("Viruses", "")]
                        family=""; genus=""; order=""; cls=""; best=""
                        # High-rank names (phylum/kingdom) must not be used as fallback
                        _high = set()
                        for p in parts:
                            if p.endswith("viridae") or p.endswith("virnae"):
                                family = p
                            elif p.endswith("virales"):
                                order  = p
                            elif p.endswith("viricetes"):
                                cls    = p
                            elif p.endswith("virus") or p.endswith("phage"):
                                genus  = p
                            elif any(p.endswith(s) for s in
                                     ("viricota", "virae", "viria", "virites")):
                                _high.add(p)  # phylum/kingdom — skip as taxonomy fallback
                        _low = [p for p in parts if p not in _high]
                        best = genus or family or order or cls or (
                            _low[-1] if _low else (parts[-1] if parts else ""))
                        genomad_tax[name] = {
                            "family": family, "genus": genus, "order": order,
                            "class": cls, "best": best, "lineage": tax,
                            "score": float(score or 0),
                        }
            lf.write(f"GeNomad: {len(genomad_tax)} contigs\n")

            # ── Build final table ─────────────────────────────────────────
            # Priority order retained for parity with the per-sample rule,
            # but vcontact3/mmseqs_custom never appear as candidates here.
            _PRIORITY = {"vcontact3": 0, "mmseqs_inphared": 1,
                         "mmseqs_custom": 2, "genomad": 3}

            def _depth(genus, family, order):
                if genus:  return 3
                if family: return 2
                if order:  return 1
                return 0

            rows = []; stats = collections.Counter()
            for contig in contigs:
                vc3 = vc3_tax.get(contig, {})
                mms = mmseqs_tax.get(contig, {})
                gmd = genomad_tax.get(contig, {})
                cms = custom_tax.get(contig, {})

                candidates = []  # (source, ff, fg, fo, lin, conf, best)

                if mms and (mms.get("family") or mms.get("genus") or mms.get("order")):
                    ff, fg, fo = mms.get("family",""), mms.get("genus",""), mms.get("order","")
                    candidates.append(("mmseqs_inphared", ff, fg, fo, mms.get("lineage",""),
                                        f"{mms.get('rank','')} ({mms.get('n_proteins',0)} proteins)",
                                        fg or ff or fo))

                if gmd and gmd.get("family"):  # only true family; no class/order via this candidate
                    ff, fg, fo = gmd.get("family",""), gmd.get("genus",""), gmd.get("order","")
                    candidates.append(("genomad", ff, fg, fo, gmd.get("lineage",""),
                                        f"{gmd['score']:.3f}", gmd.get("best","")))

                if candidates:
                    candidates.sort(key=lambda c: (-_depth(c[2], c[1], c[3]), _PRIORITY[c[0]]))
                    source, ff, fg, fo, lin, conf, best = candidates[0]
                elif gmd:
                    # GeNomad's only signal is class/order/higher -- still better than nothing
                    source = "genomad"
                    ff, fg, fo = "", "", gmd.get("order","")
                    best   = gmd.get("best","")
                    conf   = f"{gmd['score']:.3f}"
                    lin    = gmd.get("lineage","")
                else:
                    source = "unclassified"
                    ff = fg = fo = conf = lin = best = ""

                stats[source] += 1
                rows.append({
                    "seq_name":      contig,
                    "final_family":  ff,
                    "final_genus":   fg,
                    "final_order":   fo,
                    "best_taxonomy": best,
                    "source":        source,
                    "confidence":    conf,
                    "lineage":       lin,
                    "vc3_status":    vc3.get("status",""),
                    "vc3_novel_anchor": vc3.get("novel_anchor",""),
                    "genomad_best":  gmd.get("best",""),
                    "genomad_class": gmd.get("class",""),
                    "genomad_score": gmd.get("score",""),
                    "mmseqs_rank":         mms.get("rank",""),
                    "mmseqs_lineage":      mms.get("lineage",""),
                    "mmseqs_n_proteins":   mms.get("n_proteins",""),
                    "custom_rank":         cms.get("rank",""),
                    "custom_lineage":      cms.get("lineage",""),
                    "custom_n_proteins":   cms.get("n_proteins",""),
                })

            lf.write("\nSummary:\n")
            for k, v in stats.most_common():
                lf.write(f"  {k}: {v}\n")
            total = len(rows)
            unclass = stats.get("unclassified", 0)
            lf.write(f"  Novel (unclassified): {unclass}/{total} = {100*unclass/total:.1f}%\n" if total else "")
            lf.close()

            fields = ["seq_name","final_family","final_genus","final_order","best_taxonomy",
                      "source","confidence","lineage",
                      "vc3_status","vc3_novel_anchor",
                      "genomad_best","genomad_class","genomad_score",
                      "mmseqs_rank","mmseqs_lineage","mmseqs_n_proteins",
                      "custom_rank","custom_lineage","custom_n_proteins"]
            with open(str(output.tsv), "w", newline="") as f:
                w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
                w.writeheader(); w.writerows(rows)
            Path(str(output.done)).write_text("ok\n")


# ── Group vRhyme (viral vMAGs) — short reads only (needs coverage) ──────────────
# Mirrors rule vrhyme / checkv_vrhyme (rules/viral_binning.smk) but bins the group's
# CheckV-trimmed viral genomes using MULTI-SAMPLE differential coverage (all the
# group's per-sample BAMs mapped to the co-assembly by the Plan 2 map-back).
if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL and not LONG_READS:

    rule coassembly_vrhyme:
        """
        vRhyme — group viral contigs into vMAGs using coverage + protein homology.
        Uses ALL of the group's per-sample BAMs (differential coverage), analogous
        to VAMB co-binning on the prokaryotic side.
        NOTE: vRhyme creates its output dir itself — rm -rf before run.
        """
        input:
            viral = rules.coassembly_viral_trimmed.output.fasta,
            bams  = lambda wc: expand(
                f"{OUTDIR}/coassembly/{wc.group}/mapping/{{sample}}.sorted.bam",
                sample=GROUPS[wc.group]),
        output:
            done = f"{OUTDIR}/coassembly/{{group}}/bins/vrhyme/done.txt",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/vrhyme.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/vrhyme.tsv"
        conda: "../envs/env_vrhyme.yaml"
        container:  CONTAINERS.get("vrhyme")
        threads: THREADS
        params:
            outdir = f"{OUTDIR}/coassembly/{{group}}/bins/vrhyme",
        shell:
            """
            rm -rf {params.outdir}
            vRhyme \
                -i {input.viral} \
                -b {input.bams} \
                -o {params.outdir} \
                -t {threads} \
                -l {MIN_CONTIG} \
                > {log} 2>&1 || true
            mkdir -p {params.outdir}
            touch {output.done}
            """


    rule coassembly_checkv_vrhyme:
        """CheckV on the group vRhyme vMAGs (empty summary if no bins)."""
        input:
            done = rules.coassembly_vrhyme.output.done,
        output:
            summary = f"{OUTDIR}/coassembly/{{group}}/viral/checkv_vrhyme/quality_summary.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/checkv_vrhyme.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/checkv_vrhyme.tsv"
        conda: "../envs/env_viral.yaml"
        container:  CONTAINERS.get("checkv")
        threads: THREADS
        params:
            bin_dir  = f"{OUTDIR}/coassembly/{{group}}/bins/vrhyme/vRhyme_best_bins_fasta",
            out_dir  = f"{OUTDIR}/coassembly/{{group}}/viral/checkv_vrhyme",
            combined = f"{OUTDIR}/coassembly/{{group}}/viral/checkv_vrhyme/vrhyme_combined.fasta",
        shell:
            """
            rm -rf {params.out_dir}
            mkdir -p {params.out_dir}
            shopt -s nullglob
            fastas=({params.bin_dir}/*.fasta)
            if [ ${{#fastas[@]}} -gt 0 ]; then
                cat "${{fastas[@]}}" > {params.combined}
                checkv end_to_end \
                    {params.combined} {params.out_dir} \
                    -d {CHECKV_DB} -t {threads} >> {log} 2>&1
            else
                echo "No vRhyme bins — skipping CheckV" > {log}
                echo -e "contig_id\tcheckv_quality\tcompleteness\tcontig_length" > {output.summary}
            fi
            """


# ── Group prok functional foundation: protein prediction (Plan 5) ──────────────
# Prodigal per group MAG — feeds group AMR/defense/annotation. VAMB bins are *.fna
# (not *.fa like per-sample Binette). Reuses _read_manifest/_concat_proteins from
# rules/prok_binning.smk (included before this file).
if COASSEMBLY_ENABLED and COASSEMBLY_BINNING and not LONG_READS:

    rule coassembly_prok_bin_proteins:
        """Per-genome Prodigal on the group VAMB MAGs (mirrors prok_bin_proteins)."""
        input:
            done = rules.vamb_cobinning.output.done,
        output:
            manifest = f"{OUTDIR}/coassembly/{{group}}/bins/proteins/manifest.txt",
            done     = f"{OUTDIR}/coassembly/{{group}}/bins/proteins/done.txt",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/prok_bin_proteins.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/prok_bin_proteins.tsv"
        conda: "../envs/env_viral.yaml"
        container:  CONTAINERS.get("prodigal")
        threads: 1
        params:
            bins_dir = f"{OUTDIR}/coassembly/{{group}}/vamb/run/bins",
            outdir   = f"{OUTDIR}/coassembly/{{group}}/bins/proteins",
            enabled  = DEFENSE_AMR_ENABLED,
        run:
            import glob, os
            from pathlib import Path
            os.makedirs(params.outdir, exist_ok=True)
            manifest_rows = []
            with open(str(log[0]), "w") as lf:
                if not params.enabled:
                    lf.write("[coassembly_prok_bin_proteins] defense_amr disabled -- skipping\n")
                else:
                    bins = sorted(glob.glob(os.path.join(params.bins_dir, "*.fna")))
                    if bins:
                        lf.write(f"[coassembly_prok_bin_proteins] {len(bins)} MAGs -- per-genome prodigal\n")
                        for bin_fa in bins:
                            name = os.path.splitext(os.path.basename(bin_fa))[0]
                            faa  = os.path.join(params.outdir, f"{name}.faa")
                            gff  = os.path.join(params.outdir, f"{name}.gff")
                            shell("prodigal -i {bin_fa} -a {faa} -f gff -o {gff} -p single -q >> {log} 2>&1 || true")
                            if os.path.exists(faa) and os.path.getsize(faa) > 0:
                                manifest_rows.append((name, "bins", bin_fa, faa, gff))
                    else:
                        lf.write("[coassembly_prok_bin_proteins] No MAGs -- skipping\n")
                with open(str(output.manifest), "w") as mf:
                    for name, mode, fna, faa, gff in manifest_rows:
                        mf.write(f"{name}\t{mode}\t{fna}\t{faa}\t{gff}\n")
                lf.write(f"[coassembly_prok_bin_proteins] {len(manifest_rows)} genome unit(s)\n")
            Path(str(output.done)).touch()
