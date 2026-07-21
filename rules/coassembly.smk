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
        """
        input:
            contigs   = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
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
            outdir    = f"{OUTDIR}/coassembly/{{group}}/vamb/run",
            low_depth = LOW_DEPTH_MODE,
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p $(dirname {params.outdir})
            if [ "{params.low_depth}" = "True" ]; then
                echo "[VAMB co-binning] Skipped -- low_depth_mode enabled" | tee {log}
                touch {output.done}; exit 0
            fi
            rm -rf {params.outdir}
            CUDA_FLAG=""
            if [ "{USE_GPU}" = "True" ]; then CUDA_FLAG="--cuda"; fi
            vamb bin default \
                --outdir {params.outdir} \
                --fasta {input.contigs} \
                --abundance_tsv {input.abundance} \
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
