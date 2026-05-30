# ══════════════════════════════════════════════════════════════════════
# rules/viral_binning.smk — BLOCK 7: Viral Binning and QC
#
# checkv       — quality assessment on viral_consensus output
# vrhyme       — groups viral contigs into vMAGs (coverage + protein homology)
# checkv_vrhyme — quality assessment on vRhyme vMAGs
# ══════════════════════════════════════════════════════════════════════


rule checkv:
    """
    CheckV — quality assessment of consensus viral contigs.
    Classifies: Complete / High-quality / Medium-quality / Low-quality / Not-determined.
    NOTE: always removes output dir before running — CheckV skips gene calling
    if it finds existing files, causing KeyError when contig names changed.
    """
    input:
        viral = rules.viral_consensus.output.fasta,
    output:
        summary = f"{OUTDIR}/{{sample}}/viral/checkv/quality_summary.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/checkv.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/checkv.tsv"
    conda: "envs/env_viral.yaml"
    threads: THREADS
    shell:
        """
        rm -rf {OUTDIR}/{wildcards.sample}/viral/checkv
        checkv end_to_end \
            {input.viral} \
            {OUTDIR}/{wildcards.sample}/viral/checkv \
            -d {CHECKV_DB} \
            -t {threads} \
            > {log} 2>&1
        """


rule vrhyme:
    """
    vRhyme — group viral contigs into vMAGs using coverage + protein homology.
    NOTE: vRhyme creates its output dir itself — fails if it already exists.
          rm -rf before run; mkdir -p after ensures pipeline continues.
    NOTE: exits silently with no output if no contigs pass coverage threshold.
    """
    input:
        viral    = rules.viral_consensus.output.fasta,
        bam      = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bwa_done = (f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt"
                    if not LONG_READS else []),
    output:
        done = f"{OUTDIR}/{{sample}}/bins/vrhyme/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/vrhyme.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/vrhyme.tsv"
    conda: "envs/env_vrhyme.yaml"
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/bins/vrhyme",
    shell:
        """
        rm -rf {params.outdir}
        vRhyme \
            -i {input.viral} \
            -b {input.bam} \
            -o {params.outdir} \
            -t {threads} \
            -l {MIN_CONTIG} \
            > {log} 2>&1 || true
        mkdir -p {params.outdir}
        touch {output.done}
        """


rule checkv_vrhyme:
    """
    CheckV on vRhyme vMAGs.
    Concatenates only .fasta files (not .faa/.ffn) from vRhyme_best_bins_fasta/.
    Creates an empty summary if no bins were produced by vRhyme.
    """
    input:
        done = rules.vrhyme.output.done,
    output:
        summary = f"{OUTDIR}/{{sample}}/viral/checkv_vrhyme/quality_summary.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/checkv_vrhyme.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/checkv_vrhyme.tsv"
    conda: "envs/env_viral.yaml"
    threads: THREADS
    params:
        bin_dir  = f"{OUTDIR}/{{sample}}/bins/vrhyme/vRhyme_best_bins_fasta",
        out_dir  = f"{OUTDIR}/{{sample}}/viral/checkv_vrhyme",
        combined = f"{OUTDIR}/{{sample}}/viral/checkv_vrhyme/vrhyme_combined.fasta",
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


rule viral_nonredundant:
    """
    Non-redundant viral genome set for all downstream analyses.
    Strategy (bins-first):
      1. vRhyme bins  — each bin = one assembled viral genome (may span multiple contigs)
      2. Unbinned     — contigs from viral_consensus NOT in any vRhyme bin
    This replaces viral_consensus.fasta in taxonomy, host prediction and vConTACT3
    so that each viral genome is analysed exactly once.
    CheckV quality from rule checkv (run on viral_consensus.fasta) covers every
    contig in this set since viral_consensus contains all contig sequences.
    """
    input:
        viral = rules.viral_consensus.output.fasta,
        done  = rules.vrhyme.output.done,
    output:
        fasta = f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_viral_nonredundant.fasta",
    log:
        f"{OUTDIR}/{{sample}}/logs/viral_nonredundant.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/viral_nonredundant.tsv"
    params:
        vrhyme_dir = lambda wc: f"{OUTDIR}/{wc.sample}/bins/vrhyme",
    run:
        import glob, os
        vbins_dir = params.vrhyme_dir
        binned    = set()
        out_lines = []

        for vbin in sorted(glob.glob(os.path.join(vbins_dir, 'vRhyme_best_bins.*.fasta'))):
            with open(vbin) as f:
                for line in f:
                    if line.startswith('>'):
                        binned.add(line[1:].split()[0])
                    out_lines.append(line)

        skip = False
        with open(str(input.viral)) as f:
            for line in f:
                if line.startswith('>'):
                    name = line[1:].split()[0]
                    skip = name in binned
                if not skip:
                    out_lines.append(line)

        with open(str(output.fasta), 'w') as f:
            f.writelines(out_lines)

        n_bins  = len(glob.glob(os.path.join(vbins_dir, 'vRhyme_best_bins.*.fasta')))
        n_total = sum(1 for l in out_lines if l.startswith('>'))
        with open(str(log[0]), 'w') as lf:
            lf.write(f'vRhyme bins: {n_bins} ({len(binned)} binned contigs)\n')
            lf.write(f'Unbinned contigs: {n_total - len(binned)}\n')
            lf.write(f'Total non-redundant set: {n_total}\n')


rule make_votu_table:
    """
    Build the per-sample vOTU summary table.

    Crosses:
      - viral_nonredundant.fasta (vOTU representatives)
      - mmseqs cluster.tsv       (cluster size / n_members)
      - CheckV quality_summary   (completeness, quality tier, genome type)
      - VIBRANT output           (lifestyle: lytic/lysogenic; AMG count)
      - viral_taxonomy_merged    (family, genus, taxonomy source)
      - PHIST predictions        (host bin ID, soft-fail if absent)

    Output: viral/votu/{sample}_vOTU_table.tsv
    """
    input:
        viral_nr = rules.viral_nonredundant.output.fasta,
        cluster  = rules.mmseqs2.output.cluster,
        checkv   = rules.checkv.output.summary,
        vibrant  = rules.vibrant.output.done,
        taxonomy = f"{OUTDIR}/{{sample}}/viral/taxonomy/viral_taxonomy_merged.tsv",
        phist    = f"{OUTDIR}/{{sample}}/viral/phist/done.txt",
    output:
        tsv = f"{OUTDIR}/{{sample}}/viral/votu/{{sample}}_vOTU_table.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/make_votu_table.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/make_votu_table.tsv"
    params:
        vibrant_dir = f"{OUTDIR}/{{sample}}/viral/vibrant",
        phist_csv   = f"{OUTDIR}/{{sample}}/viral/phist/phist_results.csv",
    conda: "envs/env_viral.yaml"
    script:
        "../scripts/make_votu_table.py"
