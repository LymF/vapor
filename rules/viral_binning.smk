# ══════════════════════════════════════════════════════════════════════
# rules/viral_binning.smk — BLOCK 7: Viral Binning and QC
#
# cobra        — extend/join viral contigs using metaSPAdes assembly graph
# checkv       — quality assessment on (cobra_extended | viral_consensus)
# vrhyme       — groups viral contigs into vMAGs (coverage + protein homology)
# checkv_vrhyme — quality assessment on vRhyme vMAGs
# skani_votu   — ICTV vOTU clustering (95% ANI + 85% AF)
# ══════════════════════════════════════════════════════════════════════


def _viral_qc_input(wc):
    if COBRA_ENABLED:
        return f"{OUTDIR}/{wc.sample}/cobra/{wc.sample}_cobra_viral.fasta"
    return f"{OUTDIR}/{wc.sample}/viral/consensus/{wc.sample}_viral_consensus.fasta"


rule checkv:
    """
    CheckV — quality assessment of consensus viral contigs.
    Classifies: Complete / High-quality / Medium-quality / Low-quality / Not-determined.
    Input is COBRA-extended FASTA when cobra_enabled, else viral_consensus.
    NOTE: always removes output dir before running — CheckV skips gene calling
    if it finds existing files, causing KeyError when contig names changed.
    """
    input:
        viral = _viral_qc_input,
    output:
        summary = f"{OUTDIR}/{{sample}}/viral/checkv/quality_summary.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/checkv.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/checkv.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("checkv")
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
    Input is COBRA-extended FASTA when cobra_enabled, else viral_consensus.
    NOTE: vRhyme creates its output dir itself — fails if it already exists.
          rm -rf before run; mkdir -p after ensures pipeline continues.
    NOTE: exits silently with no output if no contigs pass coverage threshold.
    """
    input:
        viral    = _viral_qc_input,
        bam      = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bwa_done = (f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt"
                    if not LONG_READS else []),
    output:
        done = f"{OUTDIR}/{{sample}}/bins/vrhyme/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/vrhyme.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/vrhyme.tsv"
    conda: "../envs/env_vrhyme.yaml"
    container:  CONTAINERS.get("vrhyme")
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
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("checkv")
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
      2. Unbinned     — contigs from the viral QC input NOT in any vRhyme bin
    This replaces viral_consensus.fasta in taxonomy, host prediction and vConTACT3
    so that each viral genome is analysed exactly once.
    CheckV quality from rule checkv covers every contig in this set since the
    QC input contains all contig sequences.
    """
    input:
        viral = _viral_qc_input,
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


rule skani_votu:
    """
    skani triangle — pairwise ANI matrix for viral genomes.
    Clustering is done by the downstream skani_cluster rule (pure Python,
    runs in Snakemake's own interpreter — no container needed).
    """
    input:
        fasta = rules.viral_nonredundant.output.fasta,
    output:
        ani = f"{OUTDIR}/{{sample}}/viral/votu/skani_ani.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/skani_votu.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/skani_votu.tsv"
    conda: "../envs/env_derep.yaml"
    container: CONTAINERS.get("skani")
    threads: THREADS
    params:
        outdir  = f"{OUTDIR}/{{sample}}/viral/votu",
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


rule skani_cluster:
    """
    Greedy single-linkage vOTU clustering from the skani ANI matrix.
    Pure Python (stdlib only) — runs in Snakemake's interpreter, no container needed.
    ICTV / Roux 2019 definition: ANI >= VOTU_ANI AND max(af_q, af_r) >= VOTU_AF.
    """
    input:
        ani   = rules.skani_votu.output.ani,
        fasta = rules.viral_nonredundant.output.fasta,
    output:
        clusters = f"{OUTDIR}/{{sample}}/viral/votu/vOTU_clusters.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/skani_cluster.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/skani_cluster.tsv"
    params:
        ani_min = VOTU_ANI,
        af_min  = VOTU_AF,
        enabled = VOTU_CLUSTERING_ENABLED,
    run:
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
                        rep = comp[0]
                        for m in comp:
                            f.write(f"{rep}\t{m}\n")

                msg = (f"[skani_cluster] genomes={len(ids)} clusters={len(clusters)} "
                       f"ani>={params.ani_min} af>={params.af_min}\n")
                _lf.write(msg)
                print(msg, end="")


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
    # run: executes in the Snakemake Python process — bypasses container/conda
    # issues with script: calling `python` (not `python3`) inside containers.
    # The script uses only stdlib so no special environment is needed.
    run:
        import importlib.util, os as _os
        _path = _os.path.join(workflow.basedir, "scripts", "make_votu_table.py")
        _spec = importlib.util.spec_from_file_location("make_votu_table", _path)
        _mod  = importlib.util.module_from_spec(_spec)
        _spec.loader.exec_module(_mod)
        _mod.snakemake = snakemake
        _mod.main()
