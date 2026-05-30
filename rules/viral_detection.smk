# ══════════════════════════════════════════════════════════════════════
# rules/viral_detection.smk — BLOCK 5: Viral Detection
#
# Tools (all receive rules.mmseqs2.output.rep as input):
#   virsorter2 — HMM/ML, dsDNA + ssDNA + NCLDV + RNA + lavidaviridae
#   genomad    — marker genes + NN (NN disabled: kernel execstack issue)
#   vibrant    — HMM metabolic/structural, highest precision tool
#
# viral_consensus — integrates all three tools via configurable strategy:
#   "count"  : >= MIN_VIRAL_TOOLS tools agree (default: 2 of 3)
#   "score"  : any tool score >= threshold (SCORE_*_MIN)
#   "hybrid" : count OR single high-confidence tool
# ══════════════════════════════════════════════════════════════════════


rule virsorter2:
    """
    VirSorter2 viral detection.
    --include-groups: all major viral groups (default is only dsDNA + ssDNA).
    --hallmark-required-on-short: reduces FP on short contigs.
    """
    input:
        contigs = rules.mmseqs2.output.rep,
    output:
        viral = f"{OUTDIR}/{{sample}}/viral/virsorter2/final-viral-combined.fa",
    log:
        f"{OUTDIR}/{{sample}}/logs/virsorter2.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/virsorter2.tsv"
    conda: "envs/env_viral.yaml"
    container:  CONTAINERS.get("virsorter")
    threads: THREADS
    shell:
        """
        virsorter run \
            -i {input.contigs} \
            -w {OUTDIR}/{wildcards.sample}/viral/virsorter2 \
            --db-dir {VS2_DB} \
            --include-groups dsDNAphage,ssDNA,NCLDV,RNA,lavidaviridae \
            --min-score 0.5 \
            --hallmark-required-on-short \
            -j {threads} all \
            > {log} 2>&1
        """


rule genomad:
    """
    GeNomad viral (and plasmid) detection.
    Output: *_summary/*_virus_summary.tsv — original contig names, no renaming.
    """
    input:
        contigs = rules.mmseqs2.output.rep,
    output:
        done = f"{OUTDIR}/{{sample}}/viral/genomad/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/genomad.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/genomad.tsv"
    conda: "envs/env_genomad.yaml"
    container:  CONTAINERS.get("genomad")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/viral/genomad",
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


rule vibrant:
    """
    VIBRANT: 5th viral detection tool. Strong on integrated proviruses
    and divergent dsDNA/ssDNA phages missed by marker-based tools.
    VIBRANT outputs to CWD — cd to outdir first, using absolute paths.
    """
    input:
        contigs = rules.mmseqs2.output.rep,
    output:
        done = f"{OUTDIR}/{{sample}}/viral/vibrant/done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/vibrant.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/vibrant.tsv"
    conda: "envs/phage_vibrant.yaml"
    container:  CONTAINERS.get("vibrant")
    threads: THREADS
    params:
        outdir   = f"{OUTDIR}/{{sample}}/viral/vibrant",
        minlen   = MIN_CONTIG,
        hmm_kegg = f"{_VIBRANT_BASE}/databases/KEGG_profiles_prokaryotes.HMM",
        hmm_pfam = f"{_VIBRANT_BASE}/databases/Pfam-A_v32.HMM",
        hmm_vog  = f"{_VIBRANT_BASE}/databases/VOGDB94_phage.HMM",
        hmm_plas = f"{_VIBRANT_BASE}/databases/Pfam-A_plasmid_v32.HMM",
        hmm_phag = f"{_VIBRANT_BASE}/databases/Pfam-A_phage_v32.HMM",
        f_cat    = f"{_VIBRANT_BASE}/files/VIBRANT_categories.tsv",
        f_names  = f"{_VIBRANT_BASE}/files/VIBRANT_names.tsv",
        f_kegg   = f"{_VIBRANT_BASE}/files/VIBRANT_KEGG_pathways_summary.tsv",
        f_model  = f"{_VIBRANT_BASE}/files/VIBRANT_machine_model.sav",
        f_amgs   = f"{_VIBRANT_BASE}/files/VIBRANT_AMGs.tsv",
    shell:
        """
        mkdir -p {params.outdir}
        cp {input.contigs} {params.outdir}/input.fasta
        ABS_LOG=$(realpath {log})
        ABS_DONE=$(realpath {output.done})
        cd {params.outdir}
        VIBRANT_run.py \
            -i input.fasta \
            -f nucl \
            -t {threads} \
            -l {params.minlen} \
            -no_plot \
            -k {params.hmm_kegg} \
            -p {params.hmm_pfam} \
            -v {params.hmm_vog} \
            -e {params.hmm_plas} \
            -a {params.hmm_phag} \
            -c {params.f_cat} \
            -n {params.f_names} \
            -s {params.f_kegg} \
            -m {params.f_model} \
            -g {params.f_amgs} \
            > "$ABS_LOG" 2>&1 || true
        touch "$ABS_DONE"
        """


rule viral_consensus:
    """
    Integrate results from all 3 viral detection tools.

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
        contigs      = rules.mmseqs2.output.rep,
        vs2_done     = rules.virsorter2.output.viral,
        genomad_done = rules.genomad.output.done,
        vibrant_done = rules.vibrant.output.done,
    output:
        fasta   = f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_viral_consensus.fasta",
        support = f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_tool_support.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/viral_consensus.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/viral_consensus.tsv"
    params:
        vs2_dir     = f"{OUTDIR}/{{sample}}/viral/virsorter2",
        genomad_dir = f"{OUTDIR}/{{sample}}/viral/genomad",
        vibrant_dir = f"{OUTDIR}/{{sample}}/viral/vibrant",
        outdir      = f"{OUTDIR}/{{sample}}/viral/consensus",
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

        print(f"[viral_consensus] {wildcards.sample}: {kept} consensus viral contigs")
