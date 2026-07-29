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
        touch {output.viruses} {output.proviruses}
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
    The cluster representative is the member with the highest CheckV
    completeness (ties broken by FASTA order) — not simply the first contig
    encountered, since that's an arbitrary assembly-order artifact and the
    representative's sequence/length is what downstream genome maps and
    vOTU abundance use.
    """
    input:
        ani    = rules.skani_votu.output.ani,
        fasta  = rules.viral_nonredundant.output.fasta,
        checkv = rules.checkv.output.summary,
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


rule viral_votu_reps:
    """
    Extract vOTU representative sequences from viral_nonredundant.fasta and
    produce quality-filtered subsets used by all downstream analyses.

    Outputs:
      all_fasta      — all representatives (one per vOTU cluster), for reporting
      mq_fasta       — MQ+ (Complete/HQ/MQ or completeness>=50%) representatives,
                       for prodigal_viral, PHIST, Pharokka, genome maps
      hq_10kb_fasta  — HQ+/Complete AND >= 10 kb representatives, for vConTACT3

    Why only representatives?
      skani_cluster picks the highest-completeness member per cluster as the
      canonical genome. Running taxonomy/PHIST/etc. on all members would inflate
      compute by ~20–40% and produce duplicate annotations for near-identical
      sequences, with no scientific gain. The make_votu_table rule propagates
      representative annotations to all cluster members.
    """
    input:
        fasta    = rules.viral_nonredundant.output.fasta,
        clusters = rules.skani_cluster.output.clusters,
        checkv   = rules.checkv.output.summary,
    output:
        all_fasta     = f"{OUTDIR}/{{sample}}/viral/votu/votu_all_reps.fasta",
        mq_fasta      = f"{OUTDIR}/{{sample}}/viral/votu/votu_mq_reps.fasta",
        hq_10kb_fasta = f"{OUTDIR}/{{sample}}/viral/votu/votu_hq_10kb_reps.fasta",
    log:
        f"{OUTDIR}/{{sample}}/logs/viral_votu_reps.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/viral_votu_reps.tsv"
    params:
        hq_min_len = 10000,
    run:
        import csv, os

        # Collect unique representatives from cluster TSV
        reps = set()
        with open(str(input.clusters)) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                reps.add(row['representative'])

        # CheckV quality sets for filtering.
        #   keep_mq — downstream annotation subset (taxonomy/PHIST/Pharokka/maps).
        #             Tier gate controlled by config `viral_min_quality`
        #             (VIRAL_KEEP_TIERS). At the default "medium" threshold,
        #             completeness >= 50% is also admitted (parity with CheckV's
        #             own MQ definition).
        #   hq_plus — Complete / HQ only, for vConTACT3. Config-independent.
        keep_mq = set()
        hq_plus = set()
        _comp_fallback = VIRAL_MIN_QUALITY_RANK == 2
        with open(str(input.checkv)) as fh:
            for row in csv.DictReader(fh, delimiter='\t'):
                cid = row.get('contig_id', '').strip()
                q   = row.get('checkv_quality', '').strip()
                try:
                    comp = float(row.get('completeness', '0') or 0)
                except (ValueError, TypeError):
                    comp = 0.0
                if q in VIRAL_KEEP_TIERS or (_comp_fallback and comp >= 50):
                    keep_mq.add(cid)
                if q in ('Complete', 'High-quality'):
                    hq_plus.add(cid)

        # Parse viral_nonredundant.fasta into dict
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
        mq_reps    = reps & keep_mq
        n_mq       = write_filtered(mq_reps, str(output.mq_fasta))
        hq_reps    = reps & hq_plus
        n_hq_10kb  = write_filtered(hq_reps, str(output.hq_10kb_fasta),
                                    length_min=params.hq_min_len)

        with open(str(log[0]), 'w') as lf:
            lf.write(f'[viral_votu_reps] Total vOTU reps: {n_all}\n')
            lf.write(f'[viral_votu_reps] min quality gate: {VIRAL_MIN_QUALITY} '
                     f'(tiers kept: {sorted(VIRAL_KEEP_TIERS)})\n')
            lf.write(f'[viral_votu_reps] annotation subset (taxonomy/PHIST/Pharokka): {n_mq}\n')
            lf.write(f'[viral_votu_reps] HQ+/>=10kb reps (vConTACT3): {n_hq_10kb}\n')


rule make_votu_table:
    """
    Build the per-sample vOTU membership table — one row per cluster member.

    Representative-level annotations (CheckV, taxonomy, lifestyle, host) are
    propagated to all members; cluster_size and is_rep columns let the report
    reconstruct full cluster structure without a separate join.

    Crosses:
      - vOTU_clusters.tsv        (skani representative/member pairs)
      - votu_all_reps.fasta      (representative lengths)
      - CheckV quality_summary   (completeness, quality tier, genome type)
      - VIBRANT output           (lifestyle: lytic/lysogenic; AMG count)
      - viral_taxonomy_merged    (family, genus, taxonomy source)
      - PHIST predictions        (host bin ID, soft-fail if absent)

    Output: viral/votu/{sample}_vOTU_table.tsv
    """
    input:
        votu_clusters = rules.skani_cluster.output.clusters,
        votu_reps     = rules.viral_votu_reps.output.all_fasta,
        checkv        = rules.checkv.output.summary,
        vibrant       = rules.vibrant.output.done,
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
        vibrant_dir = f"{OUTDIR}/{{sample}}/viral/vibrant",
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
