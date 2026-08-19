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
            # Fail here rather than let an empty assembly surface as geNomad's
            # cryptic "empty or duplicate identifiers" three rules later.
            if [ ! -s {params.outdir}/final.contigs.fa ]; then
                echo "[megahit_coassembly] ERROR: group {wildcards.group} produced no contigs >= {MIN_CONTIG} bp." >&2
                echo "  Lower min_contig, or drop this group from the metadata TSV." >&2
                exit 1
            fi
            # tmp + atomic mv: an interrupted `cp` leaves a 0-byte destination
            # that Snakemake accepts as a finished output.
            cp {params.outdir}/final.contigs.fa {output.contigs}.tmp
            mv {output.contigs}.tmp {output.contigs}
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
            viral_filter = COASSEMBLY_VIRAL,
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p $(dirname {params.outdir})
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


    # `coassembly_bakta`, `coassembly_eggnog_prok` e
    # `coassembly_extract_kegg_kos` foram APAGADAS em 2026-08-19: bakta,
    # eggNOG e a extracao de KO rodam uma vez nas representantes do catalogo
    # (rules/annotation.smk) e `mag_views_group` (rules/defense_amr.smk)
    # escreve o sumario do Bakta no caminho de sempre.

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
            filter_script = os.path.join(SCRIPTS_DIR, "filter_min_length.py"),
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
            # Same guard as megahit_coassembly: never promote an empty assembly.
            if [ ! -s {params.outdir}/assembly.fasta ]; then
                echo "[flye_coassembly] ERROR: group {wildcards.group} produced no contigs (see {log})." >&2
                exit 1
            fi
            # MIN_CONTIG: o Flye nao tem flag para isso, ao contrario do
            # MEGAHIT e do metaMDBG (ver scripts/filter_min_length.py).
            python3 {params.filter_script} \
                {params.outdir}/assembly.fasta {output.contigs}.tmp {MIN_CONTIG} \
                2>> {log}
            mv {output.contigs}.tmp {output.contigs}
            """


# ══════════════════════════════════════════════════════════════════════
#  VAMB multi-split co-binning (Task 5b, Plan 2)
#
#  Independent of co-assembly: concatenates every sample's assembly
#  (MEGAHIT -- the hub since item (d), no merge, no dedup) into ONE global
#  catalog (contigs renamed with a sample prefix so bins stay traceable
#  per-sample), maps
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
        Concatenate every sample's assembly (MEGAHIT -- the hub since item
        (d), no merge, no dedup) into one global catalog for pooled VAMB
        co-binning. Headers are renamed `>X` -> `>S{sample}C{X}` so bins
        stay traceable back to the originating sample/contig. Pooled VAMB
        co-binning over the concatenated per-sample assembly catalog. (True
        per-sample bin-splitting via VAMB's --binsplit_separator is a
        possible refinement — verify VAMB v5 binsplit behavior on a real
        run before enabling.)
        """
        input:
            # _contigs_path (Snakefile): a arvore de decisao SR/ONT/HiFi mora
            # num lugar so. Este bloco e SR-only, mas chamar o helper evita a
            # quarta copia do caminho literal.
            rep_seqs = [_contigs_path(s) for s in SAMPLES],
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
        shell:
            """
            mkdir -p $(dirname {log})
            mkdir -p $(dirname {params.outdir})
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
        A trilha `multisplit` NAO entra no catalogo global de MAGs
        (rules/mag_catalog.smk) -- seus bins nao estao no pool -- entao segue
        computando o proprio GTDB-Tk sobre os bins crus do VAMB. E coerente:
        e um experimento de co-binning alternativo, nao uma fonte de MAGs
        finais. Extensao `.fna` (saida do VAMB v5).
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
#  instead of the per-sample assembly (_sample_contigs). NOT gated by
#  `if not LONG_READS:` — the viral consumer runs for either mode, since
#  `{OUTDIR}/coassembly/{group}/contigs.fa` is produced by whichever
#  assembler is active.
# ══════════════════════════════════════════════════════════════════════

if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:

    # VirSorter2 on the group co-assembly.
    # Inherits the ENTIRE body (shell, conda, container, threads, params) from
    # `rule virsorter2` in rules/viral_detection.smk — only the paths are
    # overridden. There is no second copy of the command to keep in sync.
    use rule virsorter2 as coassembly_virsorter2 with:
        input:
            contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
        output:
            viral = f"{OUTDIR}/coassembly/{{group}}/viral/virsorter2/final-viral-combined.fa",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/virsorter2.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/virsorter2.tsv"


    # geNomad no co-assembly do grupo. Herda corpo/conda/container/params
    # de `rule genomad` (rules/viral_detection.smk); so os caminhos mudam.
    use rule genomad as coassembly_genomad with:
        input:
            contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
        output:
            done = f"{OUTDIR}/coassembly/{{group}}/viral/genomad/done.txt",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/genomad.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/genomad.tsv"


    # Consenso das ferramentas no co-assembly. Herda `rule viral_consensus`
    # inteiro, incluindo o bloco run: de ~180 linhas que antes existia em
    # duplicata aqui.
    use rule viral_consensus as coassembly_viral_consensus with:
        input:
            contigs      = f"{OUTDIR}/coassembly/{{group}}/contigs.fa",
            vs2_done     = rules.coassembly_virsorter2.output.viral,
            genomad_done = rules.coassembly_genomad.output.done,
        output:
            fasta   = f"{OUTDIR}/coassembly/{{group}}/viral/consensus/{{group}}_viral_consensus.fasta",
            support = f"{OUTDIR}/coassembly/{{group}}/viral/consensus/{{group}}_tool_support.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/viral_consensus.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/viral_consensus.tsv"


    # CheckV no consenso viral do grupo. Herda `rule checkv`
    # (rules/viral_binning.smk); so o input e os caminhos mudam — o
    # per-sample entra pelo FASTA COBRA-estendido, o grupo pelo consenso.
    use rule checkv as coassembly_checkv with:
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

            import sys as _sys
            _sys.path.insert(0, SCRIPTS_DIR)
            from checkv_provirus import build_trimmed_index

            # Consenso PRIMEIRO: seus ids sao o `known` contra o qual o header
            # de provirus e resolvido (ver scripts/checkv_provirus.py). Gemeo
            # exato de `viral_trimmed` em rules/viral_binning.smk.
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

            # viruses.fna    : header = id original.
            # proviruses.fna : header = "{contig}_{n}" (so a regiao aparada).
            def _iter_fna(path, is_provirus):
                if not os.path.exists(path) or os.path.getsize(path) == 0:
                    return
                curr_hdr, curr_seq = None, []
                with open(path) as fh:
                    for line in fh:
                        if line.startswith('>'):
                            if curr_hdr:
                                yield (curr_hdr, curr_seq, is_provirus)
                            curr_hdr = line; curr_seq = []
                        else:
                            curr_seq.append(line)
                    if curr_hdr:
                        yield (curr_hdr, curr_seq, is_provirus)

            _entries = list(_iter_fna(str(input.checkv_viruses), False)) + \
                       list(_iter_fna(str(input.checkv_proviruses), True))
            trimmed, _unresolved = build_trimmed_index(_entries, set(orig_seqs))
            _n_prov = sum(1 for e in _entries if e[2])
            if _n_prov and _unresolved == _n_prov:
                raise RuntimeError(
                    "nenhum dos %d headers de provirus do CheckV foi resolvido "
                    "para um contig do consenso -- o formato do proviruses.fna "
                    "mudou. Ver scripts/checkv_provirus.py." % _n_prov)

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


    # Remove contigs virais antes do co-binning. Herda
    # `rule filter_viral_for_prok` (rules/prok_binning.smk): os corpos eram
    # identicos. O grupo entra pelos contigs do co-assembly (nao ha etapa de
    # dereplicacao MMseqs2), dai o nome de output `_contigs_nonviral`.
    use rule filter_viral_for_prok as coassembly_filter_viral_for_prok with:
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
        NOTE: vRhyme refuses -l below 2000 ("minimum scaffold length cannot be
              set below 2000 for binning"), so MIN_CONTIG is clamped to that
              floor — same pattern as rule vrhyme (viral_binning.smk) and
              MetaBAT2's 1500 bp clamp (prok_binning.smk). Below 2026-08-18
              this rule passed MIN_CONTIG unclamped and swallowed the failure
              with `|| true` + an unconditional `touch done.txt`, so once
              min_contig dropped to 1000 the group viral track would have
              produced zero bins in silence.
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
            VRHYME_MIN=$(( {MIN_CONTIG} > 2000 ? {MIN_CONTIG} : 2000 ))
            vRhyme \
                -i {input.viral} \
                -b {input.bams} \
                -o {params.outdir} \
                -t {threads} \
                -l $VRHYME_MIN \
                > {log} 2>&1 && RC=0 || RC=$?
            mkdir -p {params.outdir}
            if [ "$RC" -ne 0 ]; then
                echo "failed: vRhyme exit $RC" > {output.done}
            else
                echo "ok" > {output.done}
            fi
            """


    # CheckV nos vMAGs do vRhyme do grupo. Herda `rule checkv_vrhyme`
    # (rules/viral_binning.smk): corpos identicos.
    use rule checkv_vrhyme as coassembly_checkv_vrhyme with:
        input:
            done = rules.coassembly_vrhyme.output.done,
        output:
            summary = f"{OUTDIR}/coassembly/{{group}}/viral/checkv_vrhyme/quality_summary.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/checkv_vrhyme.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/checkv_vrhyme.tsv"


    rule coassembly_viral_nonredundant:
        """
        Group equivalent of `rule viral_nonredundant` (rules/viral_binning.smk):
        bins-first vRhyme + the composite length/quality/bin gate (item (e),
        docs/ROADMAP_SIMPLIFICACAO.md, 2026-08-18) applied to unbinned
        sequences — restoring sample/group parity (docs/AUDITORIA_COASSEMBLY_PARES.md).
        Before this rule, `_catalog_sources()` (rules/votu_catalog.smk) pointed
        the group source at `coassembly_viral_trimmed`'s PRE-binning output, so
        `coassembly_vrhyme`'s group vMAGs were computed but never reached the
        global vOTU catalog. This rule is now that source instead.

        NOT inherited from `rule viral_nonredundant` via `use rule ... as ...`:
        the per-sample rule does its OWN CheckV-trim inline (its `viral` input
        is the untrimmed consensus, trimmed via checkv_viruses/checkv_proviruses
        read fresh here). The group track already CheckV-trims upstream in
        `rule coassembly_viral_trimmed`, and that trimmed output is what feeds
        `rule coassembly_vrhyme` — so vRhyme's group bin headers already live in
        the TRIMMED namespace, not the untrimmed one `rule viral_nonredundant`
        expects. Re-running the per-sample rule's trim logic here would trim
        twice, and worse, look up vRhyme bin membership against the wrong
        (untrimmed) id namespace. A `use rule` here would need to suppress half
        the parent rule's logic; a small sibling run: block is clearer than
        forcing one rule body to serve two different upstream shapes.

        Quality lookup: `coassembly_checkv`'s quality_summary.tsv is keyed by
        the PRE-trim (consensus) contig_id. A trimmed provirus header in the
        already-trimmed input fasta carries the "{orig_id}|{start}_{end}"
        suffix `rule coassembly_viral_trimmed` writes (same convention as
        `rule viral_nonredundant`'s proviruses.fna handling) — stripped here
        the same way (rsplit on the last "|") to recover the lookup key.

        Discard audit sidecar (same schema as `rule viral_nonredundant`,
        rules/viral_binning.smk — kept identical on purpose, so the two
        tracks' TSVs concatenate without reconciling columns): every dropped
        sequence goes to `discarded_fasta` (bare contig_id header, no reason
        encoded in it -- so nothing downstream keying off the header breaks)
        and one row of `discarded_tsv`, columns
        `viral_length_gate.DISCARD_TSV_COLUMNS` = contig_id, length,
        checkv_quality, checkv_completeness, in_vrhyme_bin (always False
        here), source_id (this group's name).
        """
        input:
            trimmed        = rules.coassembly_viral_trimmed.output.fasta,
            done           = rules.coassembly_vrhyme.output.done,
            checkv_summary = rules.coassembly_checkv.output.summary,
        output:
            fasta           = f"{OUTDIR}/coassembly/{{group}}/viral/consensus/{{group}}_viral_nonredundant.fasta",
            discarded_fasta = f"{OUTDIR}/coassembly/{{group}}/viral/consensus/{{group}}_viral_discarded.fasta",
            discarded_tsv   = f"{OUTDIR}/coassembly/{{group}}/viral/consensus/{{group}}_viral_discarded.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/viral_nonredundant.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/viral_nonredundant.tsv"
        params:
            vrhyme_dir = lambda wc, input: os.path.dirname(str(input.done)),
            min_length = VIRAL_MIN_CONTIG,
        run:
            import csv, glob, os
            import sys as _sys
            _sys.path.insert(0, SCRIPTS_DIR)
            from vrhyme_bins import bin_fastas, contig_from_bin_header
            from viral_length_gate import (
                passes_gate, summarize, format_discard_row, DISCARD_TSV_COLUMNS,
            )

            # vRhyme falhou? Entao o braco "esta num bin" e INDISPONIVEL, nao
            # falso -- e o portao degradaria para comprimento-puro, descartando
            # em silencio contigs curtos que teriam sido resgatados por um bin.
            # Perda silenciosa e exatamente a classe de falha que este repo ja
            # pagou caro (docs/AUDITORIA_COASSEMBLY_PARES.md §0), entao aborta.
            # "ok" com zero bins e legitimo (comunidade diversa demais para
            # binar) e segue adiante -- so o crash aborta.
            _vr_status = ""
            if os.path.exists(str(input.done)):
                with open(str(input.done)) as _fh:
                    _vr_status = _fh.read().strip()
            if _vr_status.startswith("failed:"):
                raise RuntimeError(
                    "vRhyme falhou (%s) -- o braco de bin do portao viral ficaria "
                    "indisponivel e contigs curtos binnaveis seriam descartados em "
                    "silencio. Corrija o vRhyme antes de seguir." % _vr_status)


            # CheckV quality/completeness, keyed by the PRE-trim contig_id.
            quality = {}
            completeness = {}
            with open(str(input.checkv_summary)) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    cid = (row.get("contig_id") or "").strip()
                    if not cid:
                        continue
                    quality[cid] = (row.get("checkv_quality") or "").strip()
                    try:
                        completeness[cid] = float(row.get("completeness") or 0)
                    except ValueError:
                        completeness[cid] = 0.0

            import sys as _sys
            _sys.path.insert(0, SCRIPTS_DIR)
            from checkv_provirus import resolve_original_id

            def lookup_key(name):
                """Chave pre-trim do quality_summary para um header ja aparado.

                Usava `rsplit('|', 1)`, que nunca casava: o CheckV sufixa
                "_{n}", nao "|start_end". Efeito -- todo provirus caia para
                completude 0.0 e reprovava o braco de qualidade do portao.
                """
                orig, _ok = resolve_original_id(name, set(quality))
                return orig

            # Read the already-trimmed fasta: name -> (header, seq_lines, length)
            seqs = {}
            curr_hdr, curr_seq = None, []
            def _flush():
                if curr_hdr:
                    nm = curr_hdr[1:].split()[0]
                    ln = sum(len(l.strip()) for l in curr_seq)
                    seqs[nm] = (curr_hdr, curr_seq, ln)
            with open(str(input.trimmed)) as fh:
                for line in fh:
                    if line.startswith('>'):
                        _flush()
                        curr_hdr, curr_seq = line, []
                    else:
                        curr_seq.append(line)
                _flush()

            # vRhyme bin membership — headers are in the trimmed namespace
            # (coassembly_vrhyme's -i is the already-trimmed fasta).
            vbins_dir = params.vrhyme_dir
            binned = set()
            for vbin in bin_fastas(vbins_dir):
                with open(vbin) as fh:
                    for line in fh:
                        if line.startswith('>'):
                            binned.add(line[1:].split()[0])

            out_lines = []
            discard_out_lines = []
            discard_rows = []
            gate_results = []
            for name, (hdr, seq, length) in seqs.items():
                is_binned = name in binned
                key = lookup_key(name)
                res = passes_gate(
                    binned=is_binned,
                    quality_tier=quality.get(key, ""),
                    completeness=completeness.get(key, 0.0),
                    length=length,
                    min_length=params.min_length)
                gate_results.append(res)
                if res.kept:
                    out_lines.append(hdr)
                    out_lines.extend(seq)
                else:
                    frag_id = hdr[1:].split()[0]
                    discard_out_lines.append(f'>{frag_id}\n')
                    discard_out_lines.extend(seq)
                    discard_rows.append(format_discard_row(
                        contig_id=frag_id,
                        length=length,
                        quality_tier=quality.get(key) if key in quality else None,
                        completeness=completeness.get(key) if key in quality else None,
                        source_id=wildcards.group,
                    ))

            with open(str(output.fasta), 'w') as fh:
                fh.writelines(out_lines)

            with open(str(output.discarded_fasta), 'w') as fh:
                fh.writelines(discard_out_lines)

            with open(str(output.discarded_tsv), 'w', newline='') as fh:
                w = csv.writer(fh, delimiter='\t')
                w.writerow(DISCARD_TSV_COLUMNS)
                w.writerows(discard_rows)

            n_bins = len(bin_fastas(vbins_dir))
            counts = summarize(gate_results)
            with open(str(log[0]), 'w') as lf:
                lf.write(f'group={wildcards.group}\n')
                lf.write(f'vRhyme bins: {n_bins} ({len(binned)} binned sequences)\n')
                lf.write(f'composite gate: input={counts["total"]} '
                         f'kept_via_bin={counts["binned"]} '
                         f'kept_via_quality={counts["quality"]} '
                         f'kept_via_length>={params.min_length}={counts["length"]} '
                         f'dropped={counts["dropped"]}\n')
                lf.write(f'Discarded set: {len(discard_rows)} sequences -> '
                         f'{output.discarded_fasta} ({output.discarded_tsv})\n')
                if counts["total"] > 0 and counts["dropped"] == counts["total"]:
                    lf.write(f'WARNING: 0 of {counts["total"]} sequences kept by '
                             f'the composite gate — check upstream.\n')

# ── Group vOTU clustering + taxonomy view — both LR and SR ──────────────────
# Reopens the COASSEMBLY_VIRAL block (not the `not LONG_READS`-gated vRhyme
# block above) because skani_votu/skani_cluster/viral_votu_reps/viral_taxonomy
# run for BOTH tracks: LR groups have no group-level vRhyme (see the `if`
# just above), so their skani input below falls back to the pre-binning
# `coassembly_viral_trimmed` set — same LR behavior as before item (e)'s
# unification, preserved on purpose (see docs/ROADMAP_SIMPLIFICACAO.md (e)).
#
if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:

    # A cadeia skani LOCAL do grupo (coassembly_skani_votu /
    # coassembly_skani_cluster / coassembly_viral_votu_reps) foi APAGADA em
    # 2026-08-19. O catalogo global ja cobre estas sequencias: desde 18/08 o
    # `_catalog_sources()` (rules/votu_catalog.smk) pooleia exatamente o mesmo
    # FASTA que esta cadeia consumia, e a clusterizacao por amostra ja tinha
    # saido antes. Manter as duas significava duas definicoes de vOTU no mesmo
    # resultado -- uma global e uma por grupo, com IDs em espacos diferentes
    # (namespaced vs contig nu) e representantes possivelmente distintos.
    #
    # Quem lia a saida dela agora le o catalogo:
    #   coassembly_phist          -> votu_catalog_reps.mq_fasta (igual a
    #                                per-amostra: mesmo conjunto viral, MAGs
    #                                locais como candidatos a hospedeiro)
    #   relatorio (n_votus, vlen,
    #   curva de acumulacao)      -> vOTU_clusters.tsv do catalogo, filtrado
    #                                pela fonte do grupo
    #   final/viral/votu_*        -> final/votu_catalog/ (finalize_votu_catalog)

    # Taxonomia viral do grupo: passou a ser so uma VIEW da tabela global
    # (rule votu_taxonomy, rules/votu_catalog.smk) -- desde 2026-08-18 o
    # catalogo global ja inclui os grupos de coassembly como fonte (ver
    # _catalog_sources()), entao prodigal/mmseqs/genomad para o grupo ja
    # rodam la dentro; nao ha mais nada per-grupo para calcular aqui.
    # As gemeas coassembly_prodigal_viral / coassembly_mmseqs_taxonomy_viral /
    # coassembly_mmseqs_taxonomy_custom_viral foram removidas nesta mudanca
    # (redundantes com o catalogo global); coassembly_viral_taxonomy
    # continua existindo, mas agora so herda a VIEW (rule viral_taxonomy,
    # rules/taxonomy.smk), so com output/log/benchmark trocados para {group}.
    use rule viral_taxonomy as coassembly_viral_taxonomy with:
        output:
            tsv  = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/viral_taxonomy_merged.tsv",
            done = f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/taxonomy_done.txt",
        log:   f"{OUTDIR}/coassembly/{{group}}/logs/viral_taxonomy_view.log"
        benchmark: f"{OUTDIR}/coassembly/{{group}}/benchmarks/viral_taxonomy_view.tsv"

    # coassembly_pharokka / coassembly_phold removed 2026-08-18 (second half
    # of "(h)", docs/ROADMAP_SIMPLIFICACAO.md): pharokka/phold now run once
    # globally (votu_pharokka/votu_phold, rules/votu_catalog.smk) over the
    # whole vOTU catalog, which already includes every coassembly group as a
    # source (_catalog_sources() emits source_type 'group' too) -- a
    # per-group pharokka/phold pass would just re-annotate representatives
    # the global pass already covers.

    # coassembly_defensefinder_viral / coassembly_dbapis_viral removed
    # 2026-08-18 (second half of "(h)", docs/ROADMAP_SIMPLIFICACAO.md):
    # defensefinder_viral/dbapis_viral now run once globally
    # (votu_defensefinder_viral/votu_dbapis_viral, rules/votu_catalog.smk)
    # over the whole vOTU catalog, same move/rationale as
    # coassembly_pharokka/coassembly_phold above. Both per-group twins
    # already read the GLOBAL rules.votu_prodigal .faa (same fix applied
    # earlier the same day), so each per-group run was silently processing
    # the WHOLE catalog and writing it to a group-scoped path -- not just
    # wasted compute (N+G identical runs) but a meaning bug: coassembly's
    # finalize step (rule coassembly_organize_outputs, below) used to copy
    # that catalog-wide file into final/viral/defense_amr/ as if it were
    # the group's own defense systems.

# ── Group prok functional foundation: protein prediction (Plan 5) ──────────────
# `coassembly_prok_bin_proteins` e as seis gemeas de defesa/AMR do grupo
# (defensefinder, amrfinderplus, rgi_card, deeparg, abricate, argnorm,
# amr_consensus) foram APAGADAS em 2026-08-19: os MAGs do VAMB entram no
# catalogo global (rules/mag_catalog.smk) como qualquer outra fonte, as
# analises rodam uma vez nas representantes e `mag_views_group`
# (rules/defense_amr.smk) escreve os mesmos caminhos de antes.
if COASSEMBLY_ENABLED and COASSEMBLY_BINNING and not LONG_READS:

    # GUNC nos MAGs do co-binning VAMB do grupo. Herda `rule gunc`
    # (rules/prok_binning.smk); os bins do VAMB sao *.fna, entao bin_ext e
    # sobrescrito -- e a unica diferenca real de corpo entre as duas.
    use rule gunc as coassembly_gunc with:
        input:
            done = rules.vamb_cobinning.output.done,
        output:
            merged = f"{OUTDIR}/coassembly/{{group}}/bins/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/gunc.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/gunc.tsv"
        params:
            bins_dir = f"{OUTDIR}/coassembly/{{group}}/vamb/run/bins",
            bin_ext  = ".fna",
            outdir   = lambda wc, output: os.path.dirname(output.merged),
            db       = GUNC_DB,
            enabled  = GUNC_ENABLED,

# ── Group PHIST (viral vOTU host prediction vs group MAGs) ──────────────────────
# Needs BOTH the group's viral vOTU set (COASSEMBLY_VIRAL) and the group's
# prokaryotic MAGs (COASSEMBLY_BINNING) as candidate hosts, short reads only
# (co-binning/VAMB is short-read-only, see `if not LONG_READS:` at the top of
# this file).
if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL and COASSEMBLY_BINNING and not LONG_READS:

    # PHIST: hospedeiro dos vOTUs do grupo contra os MAGs do grupo. Herda
    # `rule phist` (rules/host_prediction.smk); bins do VAMB sao *.fna, logo
    # bin_ext e sobrescrito -- unica diferenca de corpo entre as duas.
    use rule phist as coassembly_phist with:
        input:
            # Mesmo conjunto viral do per-amostra (catalogo global), com os
            # MAGs DESTE grupo como candidatos a hospedeiro. Ate 2026-08-19
            # aqui entrava o votu_mq_reps do proprio grupo, o que deixava os
            # IDs de fago num espaco diferente (contig nu) do resto do
            # resultado -- e a documentacao afirmava o contrario do que o
            # codigo fazia.
            viral  = rules.votu_catalog_reps.output.mq_fasta,
            gtdbtk = rules.gtdbtk_group.output.done,
        output:
            done    = f"{OUTDIR}/coassembly/{{group}}/viral/phist/done.txt",
            results = f"{OUTDIR}/coassembly/{{group}}/viral/phist/phist_results.csv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/phist.log"
        benchmark:
            f"{OUTDIR}/coassembly/{{group}}/benchmarks/phist.tsv"
        params:
            bins_dir    = f"{OUTDIR}/coassembly/{{group}}/vamb/run/bins",
            bin_ext     = ".fna",
            vrhyme_dir  = f"{OUTDIR}/coassembly/{{group}}/bins/vrhyme",
            outdir      = lambda wc, output: os.path.dirname(output.done),
            scripts_dir = SCRIPTS_DIR,

# ══════════════════════════════════════════════════════════════════════
#  Organize final outputs per co-assembly group — mirrors the per-sample
#  `rule organize_outputs` (rules/finalize.smk). Builds a clean
#  coassembly/{group}/final/ tree with the group's main viral + prokaryotic
#  results. Inputs are given as explicit paths (not rules.<name>.output) so the
#  rule stays valid regardless of which coassembly_* rules are defined under
#  their own guards; the gating below mirrors `_t_coassembly` in the Snakefile,
#  so it never requests an output that was not scheduled to be produced.
# ══════════════════════════════════════════════════════════════════════
if COASSEMBLY_ENABLED and GROUPS:

    # Group prokaryotic co-binning layer runs on short reads only.
    _COAS_PROK = COASSEMBLY_BINNING and not LONG_READS

    def _coas_final_inputs(group):
        b = f"{OUTDIR}/coassembly/{group}"
        d = {}
        if _COAS_PROK:
            d["checkm2"] = f"{b}/checkm2/quality_report.tsv"
            d["gtdbtk"]  = f"{b}/gtdbtk/done.txt"
            d["bakta"]   = f"{b}/annotation/bakta/done.txt"
            if DEFENSE_AMR_ENABLED:
                d["pdef"]      = f"{b}/bins/defensefinder/done.txt"
                d["amrfinder"] = f"{b}/bins/amrfinderplus/done.txt"
                d["rgi"]       = f"{b}/bins/rgi/done.txt"
                d["deeparg"]   = f"{b}/bins/deeparg/done.txt"
            if ABRICATE_ENABLED:
                d["abricate"] = f"{b}/bins/abricate/done.txt"
            if ARGNORM_ENABLED:
                d["argnorm"] = f"{b}/bins/argnorm/done.txt"
        if COASSEMBLY_VIRAL:
            d["checkv"]   = f"{b}/viral/checkv/quality_summary.tsv"
            d["taxonomy"] = f"{b}/viral/taxonomy/taxonomy_done.txt"
            # vdef/dbapis removed 2026-08-18 (second half of "(h)"): viral
            # defense/anti-defense now lives only in the global vOTU
            # catalog (votu_defensefinder_viral/votu_dbapis_viral,
            # rules/votu_catalog.smk) -- there is no group-scoped file any
            # more to depend on here (see the removal note above
            # coassembly_organize_outputs).
            if not LONG_READS:
                d["vrhyme"] = f"{b}/bins/vrhyme/done.txt"
                # Item (e) unification: the group's "nonredundant" set in
                # final/ must be the POST-gate file (bins-first + composite
                # length/quality gate), the same one _catalog_sources()
                # (rules/votu_catalog.smk) feeds to the global catalog --
                # not the pre-binning coassembly_viral_trimmed set this
                # rule used to copy under that name.
                d["nonredundant"]   = f"{b}/viral/consensus/{group}_viral_nonredundant.fasta"
                d["discarded"]      = f"{b}/viral/consensus/{group}_viral_discarded.fasta"
                d["discarded_tsv"]  = f"{b}/viral/consensus/{group}_viral_discarded.tsv"
                if COASSEMBLY_BINNING:
                    d["phist"] = f"{b}/viral/phist/done.txt"
        return d

    rule coassembly_organize_outputs:
        """
        Assemble coassembly/{group}/final/ — the group analogue of the
        per-sample `rule organize_outputs`. Copies the group's main viral
        (vOTU reps, CheckV, taxonomy, host, defense) and prokaryotic (MAGs by
        domain, CheckM2, GTDB-Tk, AMR/defense, annotation) outputs, each with
        an existence guard. Group MAGs are VAMB `.fna` bins (not `.fa`).

        `viral_nonredundant.fasta` in final/ is short-read groups' POST-gate
        set (rule coassembly_viral_nonredundant) and their
        `viral_discarded.fasta`/`.tsv` sidecar (same schema as the per-sample
        `rule viral_nonredundant`, rules/viral_binning.smk); long-read groups
        have no group-level gate (no group vRhyme), so only the pre-binning
        `viral_nonredundant.fasta` is copied and there is no discard sidecar.
        """
        input:
            unpack(lambda wc: _coas_final_inputs(wc.group)),
        output:
            done = f"{OUTDIR}/coassembly/{{group}}/final/done.txt",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/organize_outputs.log"
        params:
            base       = f"{OUTDIR}/coassembly/{{group}}",
            bins_dir   = f"{OUTDIR}/coassembly/{{group}}/vamb/run/bins",
            vrhyme_dir = f"{OUTDIR}/coassembly/{{group}}/bins/vrhyme/vRhyme_best_bins_fasta",
        run:
            import os, glob, shutil

            b     = params.base
            final = f"{b}/final"
            for d in [
                f"{final}/viral/viral_bins",
                f"{final}/viral/taxonomy",
                f"{final}/viral/host_prediction",
                    f"{final}/bins/bacteria",
                f"{final}/bins/archaea",
                f"{final}/bins/unclassified",
                f"{final}/bins/taxonomy",
                f"{final}/bins/defense_amr",
                f"{final}/annotation",
            ]:
                os.makedirs(d, exist_ok=True)

            def cp(src, dst):
                if src and os.path.exists(str(src)) and os.path.getsize(str(src)) > 0:
                    shutil.copy(str(src), dst)

            with open(log[0], "w") as lf:

                # ── Viral core + vOTU reps ────────────────────────────
                if COASSEMBLY_VIRAL:
                    cp(f"{b}/viral/consensus/{wildcards.group}_viral_consensus.fasta",
                       f"{final}/viral/viral_consensus.fasta")
                    # "nonredundant" = the POST-gate set (item (e)
                    # unification): coassembly_viral_nonredundant for short
                    # reads (bins-first + composite length/quality gate,
                    # same source _catalog_sources() uses), or the
                    # pre-binning coassembly_viral_trimmed set for long
                    # reads (no group-level vRhyme there, so no gate to run
                    # -- see docs/ROADMAP_SIMPLIFICACAO.md item (e)).
                    if not LONG_READS:
                        cp(f"{b}/viral/consensus/{wildcards.group}_viral_nonredundant.fasta",
                           f"{final}/viral/viral_nonredundant.fasta")
                        cp(f"{b}/viral/consensus/{wildcards.group}_viral_discarded.fasta",
                           f"{final}/viral/viral_discarded.fasta")
                        cp(f"{b}/viral/consensus/{wildcards.group}_viral_discarded.tsv",
                           f"{final}/viral/viral_discarded.tsv")
                    else:
                        cp(f"{b}/viral/checkv/{wildcards.group}_viral_trimmed.fasta",
                           f"{final}/viral/viral_nonredundant.fasta")
                    cp(f"{b}/viral/checkv/quality_summary.tsv",
                       f"{final}/viral/checkv_quality.tsv")
                    # votu_all_reps/votu_annotation_reps sairam daqui em
                    # 2026-08-19 com a cadeia skani local: os representantes
                    # de vOTU sao do catalogo global e ja vao para
                    # final/votu_catalog/ (rule finalize_votu_catalog).
                    cp(f"{b}/viral/taxonomy/viral_taxonomy_merged.tsv",
                       f"{final}/viral/taxonomy/viral_taxonomy_merged.tsv")

                    # vRhyme vMAGs (short reads only)
                    if not LONG_READS:
                        vbins = glob.glob(f"{params.vrhyme_dir}/*.fasta") + \
                                glob.glob(f"{params.vrhyme_dir}/*.fna")
                        for vf in vbins:
                            shutil.copy(vf, f"{final}/viral/viral_bins/")
                        lf.write(f"vRhyme bins: {len(vbins)}\n")

                    # Host prediction
                    cp(f"{b}/viral/phist/phist_results.csv",
                       f"{final}/viral/host_prediction/phist_results.csv")

                    # Viral defense / anti-defense: no longer copied here.
                    # Removed 2026-08-18 (second half of "(h)") -- this used
                    # to copy coassembly/{group}/viral/defensefinder|dbapis,
                    # but since 2026-08-18 those rules read the GLOBAL
                    # rules.votu_prodigal .faa, so the file at that path
                    # held the WHOLE catalog's systems, not the group's.
                    # Copying it into final/ as if group-scoped was exactly
                    # the silent meaning bug (h) flagged; the only correct
                    # copy now is the single global one at
                    # {OUTDIR}/votu_catalog/{defensefinder,dbapis}/.

                # ── Prokaryotic MAGs — classify with GTDB-Tk ──────────
                archaea_bins, bacteria_bins = set(), set()
                if _COAS_PROK:
                    for tsv_path, is_arc in [
                        (f"{b}/gtdbtk/classify/gtdbtk.bac120.summary.tsv", False),
                        (f"{b}/gtdbtk/classify/gtdbtk.ar53.summary.tsv",   True),
                    ]:
                        if not os.path.exists(tsv_path):
                            continue
                        with open(tsv_path) as f:
                            hdr = None
                            for line in f:
                                parts = line.strip().split('\t')
                                if hdr is None:
                                    hdr = parts; continue
                                if not parts or len(parts) < 2:
                                    continue
                                (archaea_bins if is_arc else bacteria_bins).add(parts[0])

                    lf.write(f"GTDB-Tk — Bacteria: {len(bacteria_bins)}, "
                             f"Archaea: {len(archaea_bins)}\n")

                    copied = {'bacteria': 0, 'archaea': 0, 'unclassified': 0}
                    for bf in glob.glob(f"{params.bins_dir}/*.fna"):
                        name = os.path.splitext(os.path.basename(bf))[0]
                        if name in archaea_bins:
                            shutil.copy(bf, f"{final}/bins/archaea/"); copied['archaea'] += 1
                        elif name in bacteria_bins:
                            shutil.copy(bf, f"{final}/bins/bacteria/"); copied['bacteria'] += 1
                        else:
                            shutil.copy(bf, f"{final}/bins/unclassified/"); copied['unclassified'] += 1
                    lf.write(f"Bins copied: {copied}\n")

                    # Prok QC + taxonomy
                    cp(f"{b}/checkm2/quality_report.tsv",
                       f"{final}/bins/all_bins_checkm2.tsv")
                    cp(f"{b}/gtdbtk/classify/gtdbtk.bac120.summary.tsv",
                       f"{final}/bins/taxonomy/gtdbtk_bacteria.tsv")
                    cp(f"{b}/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
                       f"{final}/bins/taxonomy/gtdbtk_archaea.tsv")

                    # Prok defense / anti-defense + AMR
                    cp(f"{b}/bins/defensefinder/defensefinder_systems.tsv",
                       f"{final}/bins/defense_amr/defensefinder_systems.tsv")
                    cp(f"{b}/bins/defensefinder/antidefensefinder_systems.tsv",
                       f"{final}/bins/defense_amr/antidefensefinder_systems.tsv")
                    cp(f"{b}/bins/amrfinderplus/amrfinder_results.tsv",
                       f"{final}/bins/defense_amr/amrfinder_results.tsv")
                    cp(f"{b}/bins/rgi/rgi_results.txt",
                       f"{final}/bins/defense_amr/rgi_results.txt")
                    cp(f"{b}/bins/deeparg/deeparg_results.mapping.ARG",
                       f"{final}/bins/defense_amr/deeparg_results.tsv")
                    cp(f"{b}/bins/abricate/vfdb_results.tsv",
                       f"{final}/bins/defense_amr/vfdb_results.tsv")
                    cp(f"{b}/bins/abricate/plasmidfinder_results.tsv",
                       f"{final}/bins/defense_amr/plasmidfinder_results.tsv")
                    cp(f"{b}/bins/argnorm/amrfinderplus_normed.tsv",
                       f"{final}/bins/defense_amr/amrfinderplus_normed.tsv")
                    cp(f"{b}/bins/argnorm/deeparg_normed.tsv",
                       f"{final}/bins/defense_amr/deeparg_normed.tsv")

                    # Prok annotation summaries
                    cp(f"{b}/annotation/bakta/bakta_summary.tsv",
                       f"{final}/annotation/bakta_summary.tsv")

            with open(output.done, 'w') as f:
                f.write('ok\n')
