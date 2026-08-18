# ══════════════════════════════════════════════════════════════════════
# rules/viral_detection.smk — BLOCK 5: Viral Detection
#
# Tools (all receive the sample assembly -- _sample_contigs -- as input):
#   virsorter2 — HMM/ML, dsDNA + ssDNA + NCLDV + RNA + lavidaviridae
#   genomad    — marker genes + NN (NN disabled: kernel execstack issue)
#
# viral_consensus — integrates both tools via configurable strategy:
#   "count"  : >= MIN_VIRAL_TOOLS tools agree (default: 2 of 2)
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
        contigs = _sample_contigs,
    output:
        viral = f"{OUTDIR}/{{sample}}/viral/virsorter2/final-viral-combined.fa",
    log:
        f"{OUTDIR}/{{sample}}/logs/virsorter2.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/virsorter2.tsv"
    params:
        # Derived from output, never from wildcards.sample: this rule is
        # inherited by coassembly.smk via `use rule ... as ... with:`, and the
        # inherited copy has a {group} wildcard instead of {sample}. Any
        # wildcards.<name> in the body would break there.
        wdir = lambda wc, output: os.path.dirname(output.viral),
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("virsorter")
    threads: THREADS
    shell:
        """
        virsorter run \
            -i {input.contigs} \
            -w {params.wdir} \
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
        contigs = _sample_contigs,
    output:
        done = f"{OUTDIR}/{{sample}}/viral/genomad/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/genomad.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/genomad.tsv"
    conda: "../envs/env_genomad.yaml"
    container:  CONTAINERS.get("genomad")
    threads: THREADS
    params:
        # derivado do output, nao de wildcards.sample: esta regra e herdada
        # por coassembly.smk via `use rule ... as ... with:` e a copia herdada
        # usa o wildcard {group}.
        outdir = lambda wc, output: os.path.dirname(output.done),
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



rule viral_consensus:
    """
    Integrate results from all viral detection tools.

    CONTIG NAME NORMALIZATION:
    - VirSorter2 : appends ||full or ||lt0.5 — stripped with split("||")[0]
    - GeNomad    : uses original names directly; provirus "contig|prov_X_Y" → "contig"

    CONSENSUS STRATEGIES (VIRAL_CONSENSUS_MODE):
    - "count"  : keep contigs called by >= MIN_VIRAL_TOOLS tools (default: 2 of 2)
    - "score"  : keep contigs with any tool score >= threshold
    - "hybrid" : count OR single high-confidence tool (score mode)

    OUTPUT:
      *_viral_consensus.fasta — confirmed viral contigs
      *_tool_support.tsv      — all contigs × tools matrix
    """
    input:
        contigs      = _sample_contigs,
        vs2_done     = rules.virsorter2.output.viral,
        genomad_done = rules.genomad.output.done,
    output:
        fasta   = f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_viral_consensus.fasta",
        support = f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_tool_support.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/viral_consensus.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/viral_consensus.tsv"
    params:
        # todos derivados do output (ver nota em `rule genomad`): a copia
        # herdada em coassembly.smk usa {group} no lugar de {sample}.
        vs2_dir     = lambda wc, output: os.path.join(os.path.dirname(os.path.dirname(output.fasta)), "virsorter2"),
        genomad_dir = lambda wc, output: os.path.join(os.path.dirname(os.path.dirname(output.fasta)), "genomad"),
        outdir      = lambda wc, output: os.path.dirname(output.fasta),
        unit_label  = lambda wc, output: os.path.basename(output.fasta).replace("_viral_consensus.fasta", ""),
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

        # ── Count tool support per contig ─────────────────────────────
        tool_hits = defaultdict(list)
        for n in vs2_names:     tool_hits[n].append("VirSorter2")
        for n in genomad_names: tool_hits[n].append("GeNomad")

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
            lf.write(f"Union total: {len(tool_hits)}\n")
            lf.write(f"Consensus mode={VIRAL_CONSENSUS_MODE} (min_tools={MIN_VIRAL_TOOLS}): {len(consensus)}\n")
            lf.write(f"FASTA output: {kept} → {output.fasta}\n")

        print(f"[viral_consensus] {params.unit_label}: {kept} consensus viral contigs")
