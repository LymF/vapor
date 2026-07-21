# co-assembly / co-binning rules (Plano 2)
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

wildcard_constraints:
    group  = "|".join(re.escape(g) for g in GROUPS) if GROUPS else "^$",
    sample = "|".join(re.escape(s) for s in _group_samples) if _group_samples else "^$",


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
