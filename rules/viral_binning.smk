# ══════════════════════════════════════════════════════════════════════
# rules/viral_binning.smk — BLOCK 7: Viral Binning and QC
#
# cobra        — extend/join viral contigs using metaSPAdes assembly graph
# checkv       — quality assessment on (cobra_extended | viral_consensus)
# vrhyme       — groups viral contigs into vMAGs (coverage + protein homology)
# checkv_vrhyme — quality assessment on vRhyme vMAGs
# make_votu_table — per-sample vOTU membership table, sourced from the
#                    global catalog built in rules/votu_catalog.smk
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
    Outputs viruses.fna (complete viruses) and proviruses.fna (trimmed provirus
    regions — host DNA flanks removed) for use in viral_nonredundant.
    """
    input:
        viral = _viral_qc_input,
    output:
        summary   = f"{OUTDIR}/{{sample}}/viral/checkv/quality_summary.tsv",
        viruses   = f"{OUTDIR}/{{sample}}/viral/checkv/viruses.fna",
        proviruses = f"{OUTDIR}/{{sample}}/viral/checkv/proviruses.fna",
    log:
        f"{OUTDIR}/{{sample}}/logs/checkv.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/checkv.tsv"
    params:
        # derivado do output, nao de wildcards.sample: esta regra e herdada
        # por coassembly.smk via `use rule ... as ... with:`, e a copia
        # herdada usa o wildcard {group}.
        outdir = lambda wc, output: os.path.dirname(output.summary),
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("checkv")
    threads: THREADS
    shell:
        """
        rm -rf {params.outdir}
        checkv end_to_end \
            {input.viral} \
            {params.outdir} \
            -d {CHECKV_DB} \
            -t {threads} \
            > {log} 2>&1
        touch {output.viruses} {output.proviruses}
        """


rule vrhyme:
    """
    vRhyme — group viral contigs into vMAGs using coverage + protein homology.
    Input is COBRA-extended FASTA when cobra_enabled, else viral_consensus.
    NOTE: vRhyme creates its output dir itself — fails if it already exists.
          rm -rf before run; mkdir -p after ensures pipeline continues.
    NOTE: exits silently with no output if no contigs pass coverage threshold.
    NOTE: vRhyme refuses -l below 2000 ("minimum scaffold length cannot be set
          below 2000 for binning"), so MIN_CONTIG is clamped to that floor —
          same pattern as MetaBAT2's 1500 bp clamp (prok_binning.smk).
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
        VRHYME_MIN=$(( {MIN_CONTIG} > 2000 ? {MIN_CONTIG} : 2000 ))
        vRhyme \
            -i {input.viral} \
            -b {input.bam} \
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
        bin_dir  = lambda wc, input: os.path.join(os.path.dirname(str(input.done)), "vRhyme_best_bins_fasta"),
        # derivados do output (requisito da heranca por coassembly.smk).
        out_dir  = lambda wc, output: os.path.dirname(output.summary),
        combined = lambda wc, output: os.path.join(os.path.dirname(output.summary), "vrhyme_combined.fasta"),
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
    Strategy (bins-first, CheckV-trimmed):
      1. vRhyme bins  — each bin = one assembled viral genome (may span multiple contigs)
      2. Unbinned     — contigs from the viral QC input NOT in any vRhyme bin
    For every contig (binned or unbinned): if CheckV trimmed it (provirus), the
    trimmed sequence from viruses.fna / proviruses.fna replaces the original so
    that host DNA flanking proviruses is always removed. Sequences not processed
    by CheckV (rare, e.g. below CheckV internal length threshold) fall back to
    the original. One contig may yield >1 trimmed entry when CheckV found
    multiple provirus regions within it.
    """
    input:
        viral             = _viral_qc_input,
        done              = rules.vrhyme.output.done,
        checkv_viruses    = rules.checkv.output.viruses,
        checkv_proviruses = rules.checkv.output.proviruses,
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
        with open(str(input.viral)) as fh:
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

        # Process vRhyme bins first (bins-first strategy)
        vbins_dir = params.vrhyme_dir
        binned    = set()
        out_lines = []

        for vbin in sorted(glob.glob(os.path.join(vbins_dir, 'vRhyme_best_bins.*.fasta'))):
            with open(vbin) as fh:
                for line in fh:
                    if line.startswith('>'):
                        name = line[1:].split()[0]
                        if name not in binned:
                            binned.add(name)
                            emit(name, out_lines)

        # Unbinned: emit CheckV-trimmed version or original
        for name in orig_seqs:
            if name not in binned:
                emit(name, out_lines)

        with open(str(output.fasta), 'w') as fh:
            fh.writelines(out_lines)

        n_bins  = len(glob.glob(os.path.join(vbins_dir, 'vRhyme_best_bins.*.fasta')))
        n_total = sum(1 for l in out_lines if l.startswith('>'))
        n_trimmed = sum(1 for n in binned | set(orig_seqs) if n in trimmed)
        with open(str(log[0]), 'w') as lf:
            lf.write(f'vRhyme bins: {n_bins} ({len(binned)} binned contigs)\n')
            lf.write(f'CheckV-trimmed sequences: {n_trimmed}\n')
            lf.write(f'Total non-redundant set: {n_total}\n')


rule make_votu_table:
    """
    Build the per-sample vOTU membership table — one row per cluster member
    that belongs to THIS sample.

    vOTU identity (votu_id, representative) is global and stays namespaced
    ("{source_id}|{contig_id}") since a vOTU's representative may belong to
    a different sample than the one this table is built for. Members are
    filtered to this sample's own contigs and their IDs are stripped back
    to bare form so they match this sample's CheckV/taxonomy/PHIST
    tables; each member's annotations come from its own bare ID, not the
    representative's. A vOTU with no member in this sample produces no row.

    Crosses:
      - vOTU_clusters.tsv        (global catalog votu_id/representative/member triples)
      - catalog_all_reps.fasta   (representative lengths)
      - CheckV quality_summary   (completeness, quality tier, genome type)
      - viral_taxonomy_merged    (family, genus, taxonomy source)
      - PHIST predictions        (host bin ID, soft-fail if absent)

    Output: viral/votu/{sample}_vOTU_table.tsv
    """
    input:
        votu_clusters = rules.votu_catalog_cluster.output.clusters,
        votu_reps     = rules.votu_catalog_reps.output.all_fasta,
        checkv        = rules.checkv.output.summary,
        taxonomy      = f"{OUTDIR}/{{sample}}/viral/taxonomy/viral_taxonomy_merged.tsv",
        # PHIST (host prediction) requires the prok bins/GTDB-Tk chain — only
        # pull it in when the integration between viral+prok tracks is on,
        # otherwise this single hard input drags the whole prok pipeline
        # into a viral-only run. Script already soft-fails when absent.
        **({
            "phist": f"{OUTDIR}/{{sample}}/viral/phist/done.txt",
        } if INTEGRATION_ENABLED else {}),
    output:
        tsv = f"{OUTDIR}/{{sample}}/viral/votu/{{sample}}_vOTU_table.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/make_votu_table.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/make_votu_table.tsv"
    params:
        phist_csv   = f"{OUTDIR}/{{sample}}/viral/phist/phist_results.csv",
    # run: executes in the Snakemake Python process — bypasses container/conda
    # issues with script: calling `python` (not `python3`) inside containers.
    # Note: run: blocks expose input/output/params/wildcards directly (no
    # `snakemake` wrapper), so we build a SimpleNamespace to satisfy the script.
    run:
        import importlib.util, os as _os, types, traceback
        _path = _os.path.join(workflow.basedir, "scripts", "make_votu_table.py")
        _spec = importlib.util.spec_from_file_location("make_votu_table", _path)
        _mod  = importlib.util.module_from_spec(_spec)
        _spec.loader.exec_module(_mod)
        _mod.snakemake = types.SimpleNamespace(
            input=input, output=output, params=params,
            wildcards=wildcards, log=log, threads=threads,
        )
        try:
            _mod.main()
        except Exception:
            _os.makedirs(_os.path.dirname(str(log[0])), exist_ok=True)
            with open(str(log[0]), "a") as _lf:
                traceback.print_exc(file=_lf)
            raise
