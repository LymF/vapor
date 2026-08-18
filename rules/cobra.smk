# ══════════════════════════════════════════════════════════════════════
# rules/cobra.smk — BLOCK 5.5: COBRA Viral Contig Extension
#
# cobra_coverage   — convert jgi depth → 2-col TSV (run: block, no container)
# cobra_megahit    — COBRA with MEGAHIT k-mer params (maxK=141)
# cobra_spades     — COBRA with metaSPAdes k-mer params (maxK=127)
# cobra_merge      — merge both runs; longest extension wins (run: block)
#
# Only active when cobra_enabled: true in config.yaml.
# When disabled, viral_binning rules consume viral_consensus.fasta directly
# via the _viral_qc_input() hook (defined in viral_binning.smk).
#
# COMPATIBILITY:
#   COBRA extends contigs using k-mer overlap AND paired-end read linkage.
#   Only useful for paired-end Illumina (PE) short reads.
#   Set cobra_enabled: false for long-read or single-end samples.
# ══════════════════════════════════════════════════════════════════════


rule cobra_coverage:
    """
    Convert jgi_summarize_bam_contig_depths output → 2-column TSV (name, depth)
    as required by COBRA's -c flag.
    Pure Python run: block — no container needed.
    """
    input:
        depth = f"{OUTDIR}/{{sample}}/mapping/{{sample}}_depth.txt",
    output:
        cov = f"{OUTDIR}/{{sample}}/cobra/coverage.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/cobra_coverage.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/cobra_coverage.tsv"
    run:
        import os
        os.makedirs(os.path.dirname(output.cov), exist_ok=True)
        n = 0
        with open(input.depth) as fin, \
             open(output.cov, "w") as fout, \
             open(log[0], "w") as lf:
            fin.readline()  # skip jgi header row
            for line in fin:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 4:
                    fout.write(f"{parts[0]}\t{parts[3]}\n")
                    n += 1
            lf.write(f"Wrote {n} contig coverage entries\n")


rule cobra_megahit:
    """
    COBRA with MEGAHIT k-mer params (assembler=megahit, default maxK=141).
    Extends viral contigs using k-mer overlap from the MEGAHIT assembly graph.
    NOTE: uses || true — COBRA exits non-zero when no contigs are extendable,
          which is normal for highly diverse communities.
    """
    input:
        all_contigs = f"{OUTDIR}/{{sample}}/mmseqs/{{sample}}_rep_seq.fasta",
        query       = f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_viral_consensus.fasta",
        coverage    = f"{OUTDIR}/{{sample}}/cobra/coverage.tsv",
        bam         = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bwa_done    = (f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt"
                       if not LONG_READS else []),
    output:
        done = f"{OUTDIR}/{{sample}}/cobra/megahit/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/cobra_megahit.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/cobra_megahit.tsv"
    conda: "../envs/env_cobra.yaml"
    container: CONTAINERS.get("cobra_meta")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/cobra/megahit",
        mink   = COBRA_MEGAHIT_MINK,
        maxk   = COBRA_MEGAHIT_MAXK,
    shell:
        """
        rm -rf {params.outdir}
        cobra-meta \
            -f {input.all_contigs} \
            -q {input.query} \
            -m {input.bam} \
            -c {input.coverage} \
            -a megahit \
            -mink {params.mink} \
            -maxk {params.maxk} \
            -t {threads} \
            -o {params.outdir} \
            > {log} 2>&1 && RC=0 || RC=$?
        if [ "$RC" -ne 0 ]; then
            echo "failed: cobra-meta exit $RC" > {output.done}
        else
            echo "ok" > {output.done}
        fi
        """


rule cobra_spades:
    """
    COBRA with metaSPAdes k-mer params (assembler=metaspades, default maxK=127).
    Extends viral contigs using k-mer overlap from the SPAdes assembly graph.
    NOTE: uses || true for the same reason as cobra_megahit.
    """
    input:
        all_contigs = f"{OUTDIR}/{{sample}}/mmseqs/{{sample}}_rep_seq.fasta",
        query       = f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_viral_consensus.fasta",
        coverage    = f"{OUTDIR}/{{sample}}/cobra/coverage.tsv",
        bam         = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bwa_done    = (f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt"
                       if not LONG_READS else []),
    output:
        done = f"{OUTDIR}/{{sample}}/cobra/spades/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/cobra_spades.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/cobra_spades.tsv"
    conda: "../envs/env_cobra.yaml"
    container: CONTAINERS.get("cobra_meta")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/cobra/spades",
        mink   = COBRA_SPADES_MINK,
        maxk   = COBRA_SPADES_MAXK,
    shell:
        """
        rm -rf {params.outdir}
        cobra-meta \
            -f {input.all_contigs} \
            -q {input.query} \
            -m {input.bam} \
            -c {input.coverage} \
            -a metaspades \
            -mink {params.mink} \
            -maxk {params.maxk} \
            -t {threads} \
            -o {params.outdir} \
            > {log} 2>&1 && RC=0 || RC=$?
        if [ "$RC" -ne 0 ]; then
            echo "failed: cobra-meta exit $RC" > {output.done}
        else
            echo "ok" > {output.done}
        fi
        """


rule cobra_merge:
    """
    Merge COBRA-extended sequences from MEGAHIT and SPAdes runs.
    Strategy per contig:
      1. Collect all non-orphan sequences from both run output dirs
      2. Map back to original query names via the joining summary
      3. Keep the longest extension across both runs
      4. Fall back to original sequence if not extended
    Output: {sample}_cobra_viral.fasta — used by viral_binning rules when enabled.
    Pure Python run: block — no container needed.
    """
    input:
        viral   = f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_viral_consensus.fasta",
        done_mh = f"{OUTDIR}/{{sample}}/cobra/megahit/done.txt",
        done_sp = f"{OUTDIR}/{{sample}}/cobra/spades/done.txt",
    output:
        fasta = f"{OUTDIR}/{{sample}}/cobra/{{sample}}_cobra_viral.fasta",
    log:
        f"{OUTDIR}/{{sample}}/logs/cobra_merge.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/cobra_merge.tsv"
    params:
        mh_dir = f"{OUTDIR}/{{sample}}/cobra/megahit",
        sp_dir = f"{OUTDIR}/{{sample}}/cobra/spades",
    run:
        import os, glob

        def read_fasta(path):
            if not os.path.exists(path):
                return {}
            seqs, name, buf = {}, None, []
            with open(path) as f:
                for line in f:
                    line = line.rstrip()
                    if line.startswith(">"):
                        if name:
                            seqs[name] = "".join(buf)
                        name, buf = line[1:].split()[0], []
                    else:
                        buf.append(line)
            if name:
                seqs[name] = "".join(buf)
            return seqs

        def parse_summary(path):
            """Return {query_name: category} from a COBRA joining summary TSV."""
            result = {}
            if not os.path.exists(path):
                return result
            with open(path) as f:
                for line in f:
                    line = line.strip()
                    if not line or line.startswith("#") or line.startswith("QueryName"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 3:
                        result[parts[0]] = parts[2]
            return result

        orig = read_fasta(str(input.viral))
        # best[orig_name] = (header_in_output, sequence)
        best = {}

        for run_dir in [params.mh_dir, params.sp_dir]:
            # --- find joining summary ---
            summaries = glob.glob(os.path.join(run_dir, "*joining_summary*"))
            categories = parse_summary(summaries[0]) if summaries else {}

            # --- collect non-orphan sequences from this run ---
            # Primary: look for contigs.new.fa (COBRA v1.3+)
            # Fallback: glob all *.fasta excluding orphan_end
            ext_seqs = {}
            new_fa = os.path.join(run_dir, "contigs.new.fa")
            if os.path.exists(new_fa):
                ext_seqs = read_fasta(new_fa)
                # contigs.new.fa contains ALL queries including orphan_end;
                # filter out those that are categorised as orphan_end
                ext_seqs = {
                    k: v for k, v in ext_seqs.items()
                    if categories.get(k, "extended") not in ("orphan_end",)
                       and len(v) > len(orig.get(k, ""))
                }
            else:
                for fpath in sorted(glob.glob(os.path.join(run_dir, "*.fasta"))):
                    if "orphan" in os.path.basename(fpath):
                        continue
                    ext_seqs.update(read_fasta(fpath))

            # --- map original name → best extended sequence ---
            for ext_name, ext_seq in ext_seqs.items():
                # The original query name is the first token before any COBRA suffix
                orig_name = ext_name
                for sep in ("_join_", "_extended_", "_circular", "_self_circular"):
                    if sep in ext_name:
                        orig_name = ext_name.split(sep)[0]
                        break
                if orig_name not in orig:
                    # Try direct match (COBRA may keep original name unchanged)
                    orig_name = ext_name
                if orig_name in orig:
                    if orig_name not in best or len(ext_seq) > len(best[orig_name][1]):
                        best[orig_name] = (ext_name, ext_seq)

        n_ext = 0
        out_lines = []
        for orig_name, orig_seq in orig.items():
            if orig_name in best:
                new_hdr, new_seq = best[orig_name]
                out_lines.append(f">{new_hdr}\n{new_seq}\n")
                n_ext += 1
            else:
                out_lines.append(f">{orig_name}\n{orig_seq}\n")

        with open(str(output.fasta), "w") as f:
            f.writelines(out_lines)

        with open(str(log[0]), "w") as lf:
            lf.write(f"Original sequences  : {len(orig)}\n")
            lf.write(f"COBRA-extended      : {n_ext}\n")
            lf.write(f"Output sequences    : {len(orig)}\n")
