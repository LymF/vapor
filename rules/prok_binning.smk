# ══════════════════════════════════════════════════════════════════════
# rules/prok_binning.smk — BLOCK 8: Prokaryote Binning
#
# Pre-binning:
#   filter_viral_for_prok — remove free-living viral contigs (keep provirus)
#
# Binners:
#   metabat2          — tetranucleotide + coverage
#   vamb              — Variational Autoencoder (v5 interface)
#   semibin2          — semi-supervised with environment priors
#   comebin           — transformer + contrastive learning (optional)
#
# Integration:
#   prepare_scaffold2bin — convert binner outputs to scaffold2bin.tsv
#   binette              — best-bin selection (successor to DAS Tool)
#
# QC + refinement + taxonomy:
#   checkm2           — completeness / contamination (MIMAG standards)
#   gunc              — chimera detection (optional)
#   galah_derep       — CheckM2-aware MAG dereplication by ANI (optional)
#   gtdbtk            — GTDB taxonomy for MAGs
# ══════════════════════════════════════════════════════════════════════


def _prok_input_contigs(wc):
    """FASTA used as input to all prok binners.
    Returns filtered (viral-free, provirus-aware) FASTA when PROK_FILTER_VIRAL
    is enabled; otherwise falls back to the raw MMseqs2 representative set."""
    if PROK_FILTER_VIRAL:
        return f"{OUTDIR}/{wc.sample}/prok_input/{wc.sample}_rep_seq_nonviral.fasta"
    return f"{OUTDIR}/{wc.sample}/mmseqs/{wc.sample}_rep_seq.fasta"


def _gtdbtk_bins_dir(wc):
    """Bin set passed to GTDB-Tk: dereplicated when galah is enabled."""
    if MAG_DEREP_ENABLED:
        return f"{OUTDIR}/{wc.sample}/bins/derep/derep_bins"
    return f"{OUTDIR}/{wc.sample}/bins/binette/final_bins"


rule filter_viral_for_prok:
    """
    Remove free-living viral contigs from the prokaryotic binning input.

    Strategy:
      1. viral_consensus.fasta              → set of viral contigs
      2. CheckV `provirus=Yes` + GeNomad `|provirus_` suffix → set of provirus contigs
      3. remove = viral_consensus MINUS provirus     (when keep_provirus=True)
                = viral_consensus                    (when keep_provirus=False)
      4. {sample}_rep_seq.fasta MINUS remove → {sample}_rep_seq_nonviral.fasta

    Provirus-bearing contigs stay in the prok input because the provirus
    region is integrated within a bacterial host contig — removing it would
    drop the host genome with the prophage. Free-living viruses (no host
    chromosome context) are removed to clean up MAGs.

    NOTE: mapping/calc_depth STILL uses the unfiltered rep_seq.fasta for
    statistically correct coverage estimation; only contig→bin assignment
    uses the filtered FASTA.
    """
    input:
        contigs   = rules.mmseqs2.output.rep,
        viral     = rules.viral_consensus.output.fasta,
        checkv    = rules.checkv.output.summary,
        genomad   = rules.genomad.output.done,
    output:
        filtered  = f"{OUTDIR}/{{sample}}/prok_input/{{sample}}_rep_seq_nonviral.fasta",
        stats     = f"{OUTDIR}/{{sample}}/prok_input/filter_stats.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/filter_viral_for_prok.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/filter_viral_for_prok.tsv"
    params:
        keep_provirus = PROK_FILTER_KEEP_PROVIRUS,
        genomad_dir   = lambda wc: f"{OUTDIR}/{wc.sample}/viral/genomad",
    run:
        import os, glob, csv

        viral_set = set()
        with open(input.viral) as f:
            for line in f:
                if line.startswith(">"):
                    viral_set.add(line[1:].strip().split()[0])

        provirus_set = set()
        if params.keep_provirus:
            # CheckV column: provirus=Yes
            try:
                with open(input.checkv) as f:
                    rdr = csv.DictReader(f, delimiter="\t")
                    for row in rdr:
                        if (row.get("provirus", "") or "").strip().lower() == "yes":
                            cid = (row.get("contig_id", "") or "").strip()
                            if cid:
                                provirus_set.add(cid)
            except Exception:
                pass

            # GeNomad: virus_summary entries with "|provirus_" suffix in seq_name
            for gf in glob.glob(os.path.join(str(params.genomad_dir),
                                             "**", "*_virus_summary.tsv"),
                                recursive=True):
                try:
                    with open(gf) as f:
                        rdr = csv.DictReader(f, delimiter="\t")
                        for row in rdr:
                            seq = (row.get("seq_name", "") or "").strip()
                            if "|provirus_" in seq:
                                provirus_set.add(seq.split("|")[0])
                except Exception:
                    pass

        # Contigs we will remove from the prok input
        remove_set = viral_set - provirus_set

        os.makedirs(os.path.dirname(output.filtered), exist_ok=True)
        kept = removed = total = 0
        with open(input.contigs) as fin, open(output.filtered, "w") as fout:
            write = False
            header = None
            for line in fin:
                if line.startswith(">"):
                    if header is not None:
                        # finalize previous record (already written if write=True)
                        pass
                    header = line[1:].strip().split()[0]
                    total += 1
                    if header in remove_set:
                        write = False
                        removed += 1
                    else:
                        write = True
                        kept += 1
                        fout.write(line)
                else:
                    if write:
                        fout.write(line)

        with open(output.stats, "w") as s:
            s.write("metric\tcount\n")
            s.write(f"total_contigs\t{total}\n")
            s.write(f"viral_total\t{len(viral_set)}\n")
            s.write(f"provirus_kept\t{len(provirus_set)}\n")
            s.write(f"removed\t{removed}\n")
            s.write(f"kept_for_prok\t{kept}\n")

        with open(log[0], "w") as lf:
            lf.write(f"Total contigs        : {total}\n")
            lf.write(f"Viral consensus      : {len(viral_set)}\n")
            lf.write(f"Provirus (kept)      : {len(provirus_set)}\n")
            lf.write(f"Removed (viral-only) : {removed}\n")
            lf.write(f"Kept for prok binning: {kept}\n")


rule metabat2:
    """
    MetaBAT2 — tetranucleotide + coverage binning.
    --minContig enforced minimum is 1500 (tool limit) — we use max(MIN_CONTIG, 1500).
    Input FASTA is viral-filtered when PROK_FILTER_VIRAL is enabled.
    """
    input:
        contigs = _prok_input_contigs,
        depth   = rules.calc_depth.output.depth,
    output:
        done = f"{OUTDIR}/{{sample}}/bins/metabat2/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/metabat2.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/metabat2.tsv"
    conda: "../envs/env_binning.yaml"
    container:  CONTAINERS.get("metabat2")
    threads: THREADS
    params:
        outdir     = f"{OUTDIR}/{{sample}}/bins/metabat2",
        low_depth  = LOW_DEPTH_MODE,
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.low_depth}" = "True" ]; then
            echo "[MetaBAT2] Skipped -- low_depth_mode enabled" | tee {log}
            touch {output.done}; exit 0
        fi
        METABAT_MIN=$(( {MIN_CONTIG} > 1500 ? {MIN_CONTIG} : 1500 ))
        metabat2 \
            -i {input.contigs} \
            -a {input.depth} \
            -o {params.outdir}/bin \
            -t {threads} \
            --minContig $METABAT_MIN \
            --unbinned \
            > {log} 2>&1
        touch {output.done}
        """


rule vamb:
    """
    VAMB v5 — Variational Autoencoder binning.
    IMPORTANT v5 changes from v4:
      - Subcommand: vamb bin default
      - No --bamfiles; requires --abundance_tsv
      - TSV header must be: contigname <TAB> sample_name
    NOTE: VAMB creates its output dir itself — rm -rf before run.
    Input FASTA is viral-filtered when PROK_FILTER_VIRAL is enabled.
    """
    input:
        contigs = _prok_input_contigs,
        depth   = rules.calc_depth.output.depth,
    output:
        done = f"{OUTDIR}/{{sample}}/bins/vamb/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/vamb.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/vamb.tsv"
    conda: "../envs/env_binning.yaml"
    container:  CONTAINERS.get("vamb")
    threads: THREADS
    params:
        outdir     = f"{OUTDIR}/{{sample}}/bins/vamb",
        low_depth  = LOW_DEPTH_MODE,
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.low_depth}" = "True" ]; then
            echo "[VAMB] Skipped -- low_depth_mode enabled" | tee {log}
            touch {output.done}; exit 0
        fi
        rm -rf {params.outdir}
        mkdir -p $(dirname {params.outdir}_abundance.tsv)
        # Build abundance TSV filtered to only contigs present in the (viral-filtered) FASTA.
        # The depth file covers all contigs (viral + non-viral); VAMB requires an exact match
        # with the FASTA sequence set — any surplus rows cause a ValueError.
        grep "^>" {input.contigs} | sed 's/^>//' | awk '{{print $1}}' \
            > {params.outdir}_contig_names.txt
        echo -e "contigname\\t{wildcards.sample}" > {params.outdir}_abundance.tsv
        awk 'NR==FNR{{keep[$1]=1; next}} FNR>1 && ($1 in keep){{print $1"\\t"$3}}' \
            {params.outdir}_contig_names.txt {input.depth} \
            >> {params.outdir}_abundance.tsv
        rm -f {params.outdir}_contig_names.txt
        CUDA_FLAG=""
        if [ "{USE_GPU}" = "True" ]; then CUDA_FLAG="--cuda"; fi
        vamb bin default \
            --outdir {params.outdir} \
            --fasta {input.contigs} \
            --abundance_tsv {params.outdir}_abundance.tsv \
            --minfasta 200000 \
            -p {threads} \
            $CUDA_FLAG \
            > {log} 2>&1
        touch {output.done}
        """


rule semibin2:
    """
    SemiBin2 — semi-supervised binning with environment-specific priors.
    Set SEMIBIN_ENV to match your sample type:
    soil / ocean / gut / wastewater / global (use global if unsure).
    Input FASTA is viral-filtered when PROK_FILTER_VIRAL is enabled.
    """
    input:
        contigs  = _prok_input_contigs,
        bam      = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bwa_done = (f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt"
                    if not LONG_READS else []),
    output:
        done = f"{OUTDIR}/{{sample}}/bins/semibin2/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/semibin2.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/semibin2.tsv"
    conda: "../envs/env_binning.yaml"
    container:  CONTAINERS.get("semibin")
    threads: THREADS
    params:
        outdir     = f"{OUTDIR}/{{sample}}/bins/semibin2",
        low_depth  = LOW_DEPTH_MODE,
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.low_depth}" = "True" ]; then
            echo "[SemiBin2] Skipped -- low_depth_mode enabled" | tee {log}
            touch {output.done}; exit 0
        fi
        ENGINE="cpu"
        if [ "{USE_GPU}" = "True" ]; then ENGINE="gpu"; fi
        SemiBin2 single_easy_bin \
            -i {input.contigs} \
            -b {input.bam} \
            -o {params.outdir} \
            --environment {SEMIBIN_ENV} \
            --engine $ENGINE \
            -t {threads} \
            > {log} 2>&1 || echo "[SemiBin2] WARNING: binning failed (e.g. too few long contigs) — skipping" | tee -a {log}
        touch {output.done}
        """


rule comebin:
    """
    COMEBin — transformer embedding + contrastive learning binner.
    Rank 1 in 2025 metagenome binning benchmark (Nature Communications).
    Added as 4th binner alongside MetaBAT2 + VAMB + SemiBin2.
    Binette uses all 4 scaffold2bin files for refined bin selection.
    Set comebin_enabled: false in config.yaml to skip.
    Input FASTA is viral-filtered when PROK_FILTER_VIRAL is enabled.
    """
    input:
        contigs  = _prok_input_contigs,
        bam      = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam",
        bai      = f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam.bai",
        bwa_done = (f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt"
                    if not LONG_READS else []),
    output:
        done = f"{OUTDIR}/{{sample}}/bins/comebin/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/comebin.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/comebin.tsv"
    conda: "../envs/env_comebin.yaml"
    container:  CONTAINERS.get("comebin")
    threads: THREADS
    params:
        outdir  = f"{OUTDIR}/{{sample}}/bins/comebin",
        bam_dir = f"{OUTDIR}/{{sample}}/mapping",
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{COMEBIN_ENABLED}" != "True" ]; then
            echo "[COMEBin] Disabled via config (comebin_enabled: false)" | tee {log}
            touch {output.done}; exit 0
        fi
        if [ "{LOW_DEPTH_MODE}" = "True" ]; then
            echo "[COMEBin] Skipped -- low_depth_mode enabled" | tee {log}
            touch {output.done}; exit 0
        fi
        # COMEBin internally runs 'samtools depth -aa' on the BAM.
        # -aa outputs ALL reference sequences listed in the BAM @SQ header,
        # even those with zero coverage. If the BAM was mapped to all contigs
        # (viral + non-viral) the depth file will contain viral contig names
        # that are absent from aug_seq_info_dict → KeyError in gen_cov.py.
        # Fix: build a BAM that has (1) only reads on non-viral contigs and
        # (2) @SQ header lines stripped of viral references.
        FILT_DIR="{params.outdir}/bam_nonviral"
        mkdir -p "$FILT_DIR"
        # Non-viral contig name list + BED for read filtering
        python3 -c "
with open('{input.contigs}') as fh:
    for line in fh:
        if line.startswith('>'):
            print(line[1:].split()[0])
" > "$FILT_DIR/names.txt"
        awk '{{print $1"\t0\t999999999"}}' "$FILT_DIR/names.txt" \
            > "$FILT_DIR/nonviral.bed"
        # Build filtered header (SAM text, non-viral @SQ only)
        python3 -c "
import subprocess
names = set(open('$FILT_DIR/names.txt').read().split())
hdr = subprocess.check_output(['samtools','view','-H','{input.bam}']).decode()
for l in hdr.splitlines():
    if l.startswith('@SQ') and not any(f.startswith('SN:') and f[3:] in names for f in l.split('\t')):
        continue
    print(l)
" > "$FILT_DIR/new_header.sam" 2>> {log}
        # Stream: new header (SAM text) + reads on non-viral contigs (SAM text)
        # → samtools view -bS maps contig NAMES → correct new indices → valid BAM
        # (samtools reheader only replaces text but not binary indices → corrupt BAM)
        {{ cat "$FILT_DIR/new_header.sam"; \
           samtools view -L "$FILT_DIR/nonviral.bed" {input.bam}; }} \
        | samtools view -bS -o "$FILT_DIR/nonviral.bam" 2>> {log}
        samtools index "$FILT_DIR/nonviral.bam" 2>> {log}
        if [ "{USE_GPU}" = "True" ]; then
            export CUDA_VISIBLE_DEVICES=0
        else
            export CUDA_VISIBLE_DEVICES=""
        fi
        run_comebin.sh \
            -a {input.contigs} \
            -o {params.outdir} \
            -p "$FILT_DIR" \
            -t {threads} \
            >> {log} 2>&1 || true
        touch {output.done}
        N=$(find {params.outdir}/comebin_res/comebin_res_bins -name "*.fa" 2>/dev/null | wc -l) || true
        echo "[COMEBin] $N bins produced" | tee -a {log}
        """


rule prepare_scaffold2bin:
    """
    Convert each binner output to scaffold2bin.tsv (contig <TAB> bin).
    Each binner has a different format:
      MetaBAT2 : bin.X.fa files
      VAMB v5  : clusters.tsv — REVERSED columns (bin_id TAB contig_id)
      SemiBin2 : contig_bins.tsv — contig TAB bin_number
      COMEBin  : bins/*.fa files (if comebin_enabled)
    """
    input:
        mb2     = rules.metabat2.output.done,
        vamb    = rules.vamb.output.done,
        sb2     = rules.semibin2.output.done,
        comebin = rules.comebin.output.done,
    output:
        done = f"{OUTDIR}/{{sample}}/bins/scaffold2bin/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/scaffold2bin.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/prepare_scaffold2bin.tsv"
    params:
        s = f"{OUTDIR}/{{sample}}",
    run:
        import os, glob

        outdir = f"{params.s}/bins/scaffold2bin"
        os.makedirs(outdir, exist_ok=True)

        def write_s2b(bin_files, outfile, strip_ext):
            with open(outfile, "w") as out:
                for bf in bin_files:
                    bin_name = os.path.basename(bf).replace(strip_ext, "")
                    with open(bf) as f:
                        for line in f:
                            if line.startswith(">"):
                                contig = line[1:].strip().split()[0]
                                out.write(f"{contig}\t{bin_name}\n")

        mb_bins = glob.glob(f"{params.s}/bins/metabat2/bin.*.fa")
        if mb_bins:
            write_s2b(mb_bins, f"{outdir}/metabat2_s2b.tsv", ".fa")

        vamb_clusters = f"{params.s}/bins/vamb/clusters.tsv"
        if os.path.exists(vamb_clusters):
            with open(vamb_clusters) as fin, open(f"{outdir}/vamb_s2b.tsv", "w") as fout:
                for line in fin:
                    parts = line.strip().split("\t")
                    if len(parts) == 2 and not parts[0].startswith("clusterid"):
                        fout.write(f"{parts[1]}\t{parts[0]}\n")

        # SemiBin2: use contig_bins.tsv directly (bins are SemiBin_N.fa.gz)
        sb_tsv = f"{params.s}/bins/semibin2/contig_bins.tsv"
        if os.path.exists(sb_tsv):
            with open(sb_tsv) as fin, open(f"{outdir}/semibin2_s2b.tsv", "w") as fout:
                header = True
                for line in fin:
                    if header:
                        header = False
                        continue
                    parts = line.strip().split("\t")
                    if len(parts) == 2:
                        contig, bin_num = parts
                        fout.write(f"{contig}\tSemiBin_{bin_num}\n")

        # COMEBin: comebin_res/comebin_res_bins/*.fa (run_comebin.sh output layout)
        cb_bins = glob.glob(f"{params.s}/bins/comebin/comebin_res/comebin_res_bins/*.fa")
        if cb_bins:
            write_s2b(cb_bins, f"{outdir}/comebin_s2b.tsv", ".fa")

        with open(output.done, "w") as f:
            f.write("ok\n")


rule binette:
    """
    Binette — best-bin selection across multiple binners (2023).
    Successor to DAS Tool: faster, better scoring, same scaffold2bin input format.
    Outputs: final_bins/ (FASTA) + binette_results.tsv (quality table).
    Filters scaffold2bin to only contigs present in the assembly (safety check).
    Uses the viral-filtered FASTA when PROK_FILTER_VIRAL is enabled to keep
    bin assignment consistent with the binners.
    """
    input:
        contigs = _prok_input_contigs,
        s2b     = rules.prepare_scaffold2bin.output.done,
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/binette/done.txt",
        summary = f"{OUTDIR}/{{sample}}/bins/binette/binette_results.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/binette.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/binette.tsv"
    conda: "../envs/env_binette.yaml"
    container:  CONTAINERS.get("binette")
    threads: THREADS
    params:
        s2b_dir = f"{OUTDIR}/{{sample}}/bins/scaffold2bin",
        outdir  = f"{OUTDIR}/{{sample}}/bins/binette",
    shell:
        """
        mkdir -p {params.outdir}

        VALID_CONTIGS=$(mktemp)
        grep "^>" {input.contigs} | sed 's/^>//' | awk '{{print $1}}' > "$VALID_CONTIGS"
        echo "[binette] $(wc -l < $VALID_CONTIGS) assembly contigs" | tee {log}

        FILTERED_DIR={params.s2b_dir}/filtered
        mkdir -p "$FILTERED_DIR"
        rm -f "$FILTERED_DIR"/*.tsv
        for S2B in {params.s2b_dir}/*_s2b.tsv; do
            [ -f "$S2B" ] || continue
            BNAME=$(basename "$S2B")
            awk 'NR==FNR{{valid[$1]=1;next}} ($1 in valid){{print}}' \
                "$VALID_CONTIGS" "$S2B" > "$FILTERED_DIR/$BNAME"
            [ -s "$FILTERED_DIR/$BNAME" ] || rm -f "$FILTERED_DIR/$BNAME"
        done
        rm -f "$VALID_CONTIGS"

        S2B_FILES=$(ls "$FILTERED_DIR"/*_s2b.tsv 2>/dev/null | tr '\n' ' ')
        if [ -z "$(echo $S2B_FILES | tr -d ' ')" ]; then
            echo "[WARNING] No valid s2b files — skipping Binette" | tee -a {log}
            mkdir -p {params.outdir}/final_bins
            printf "Name\tCompleteness\tContamination\tScore\n" > {output.summary}
            touch {output.done}; exit 0
        fi

        binette -b $S2B_FILES -c {input.contigs} -o {params.outdir} -t {threads} --checkm2_db {CHECKM2_DB} >> {log} 2>&1

        # Locate summary (Binette may name it differently across versions)
        for candidate in {params.outdir}/binette_results.tsv \
                          {params.outdir}/final_bins_quality.tsv \
                          {params.outdir}/results.tsv; do
            [ -f "$candidate" ] && cp "$candidate" {output.summary} && break
        done
        [ -f {output.summary} ] || \
            printf "Name\tCompleteness\tContamination\tScore\n" > {output.summary}
        touch {output.done}
        """


rule prok_bin_proteins:
    """
    Per-genome protein prediction (Prodigal) shared by every defense/AMR
    annotation tool in rules/defense_amr.smk plus diamond_custom_prok in
    rules/taxonomy.smk -- runs once regardless of how many downstream tools
    consume it (same reuse pattern as rules.prodigal_viral.output.faa for
    viral taxonomy). Lives here (not in defense_amr.smk, where it
    conceptually originated) because both defense_amr.smk and taxonomy.smk
    reference rules.prok_bin_proteins, and Snakemake requires a rule to be
    defined (via include: order in the Snakefile) before anything reuses it
    through `rules.<name>`.

    Normal path: one Prodigal call per Binette final bin (single mode).
    Low-depth fallback: one Prodigal call (meta mode) on the whole
    viral-filtered contig set used as binning input, treated as a single
    pseudo-genome -- triggered either when binning produced zero bins, or
    explicitly by low_depth_mode (config.yaml; binners are already skipped
    upstream in that case, but the flag is also checked here directly).
    Gated by defense_amr_contig_fallback; every downstream rule just
    iterates whatever the manifest contains, with no separate low-depth
    codepath of its own.
    """
    input:
        contigs = _prok_input_contigs,
        done    = rules.binette.output.done,
    output:
        manifest = f"{OUTDIR}/{{sample}}/bins/proteins/manifest.txt",
        done     = f"{OUTDIR}/{{sample}}/bins/proteins/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/prok_bin_proteins.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/prok_bin_proteins.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("prodigal")
    threads: 1
    params:
        bins_dir   = lambda wc: f"{OUTDIR}/{wc.sample}/bins/binette/final_bins",
        outdir     = f"{OUTDIR}/{{sample}}/bins/proteins",
        fallback   = DEFENSE_AMR_CONTIG_FALLBACK,
        low_depth  = LOW_DEPTH_MODE,
        enabled    = DEFENSE_AMR_ENABLED,
    run:
        import glob, os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        manifest_rows = []

        with open(str(log[0]), "w") as lf:
            if not params.enabled:
                lf.write("[prok_bin_proteins] defense_amr_enabled=False -- skipping\n")
            else:
                # low_depth_mode forces the single-pseudo-genome path regardless
                # of whether bins exist -- the binners themselves are already
                # skipped (prok_binning.smk) when this flag is set, but checking
                # it explicitly here documents the intent rather than relying on
                # bins_dir happening to be empty.
                bins = ([] if params.low_depth else
                        sorted(glob.glob(os.path.join(params.bins_dir, "*.fa"))))
                if bins:
                    lf.write(f"[prok_bin_proteins] {len(bins)} bins -- per-genome protein prediction\n")
                    for bin_fa in bins:
                        name = os.path.splitext(os.path.basename(bin_fa))[0]
                        faa  = os.path.join(params.outdir, f"{name}.faa")
                        gff  = os.path.join(params.outdir, f"{name}.gff")
                        shell(
                            "prodigal -i {bin_fa} -a {faa} -f gff -o {gff} -p single -q "
                            ">> {log} 2>&1 || true"
                        )
                        if os.path.exists(faa) and os.path.getsize(faa) > 0:
                            manifest_rows.append((name, "bins", bin_fa, faa, gff))
                elif ((params.fallback or params.low_depth) and os.path.exists(str(input.contigs))
                      and os.path.getsize(str(input.contigs)) > 0):
                    lf.write("[prok_bin_proteins] No bins (low depth) -- "
                             "fallback: contigs as pseudo-genome\n")
                    name = "contigs_pseudogenome"
                    faa  = os.path.join(params.outdir, f"{name}.faa")
                    gff  = os.path.join(params.outdir, f"{name}.gff")
                    contigs = str(input.contigs)
                    shell(
                        "prodigal -i {contigs} -a {faa} -f gff -o {gff} -p meta -q "
                        ">> {log} 2>&1 || true"
                    )
                    if os.path.exists(faa) and os.path.getsize(faa) > 0:
                        manifest_rows.append((name, "contigs", contigs, faa, gff))
                else:
                    lf.write("[prok_bin_proteins] No bins and fallback disabled/no contigs "
                             "-- skipping\n")

            with open(str(output.manifest), "w") as mf:
                for name, mode, fna, faa, gff in manifest_rows:
                    mf.write(f"{name}\t{mode}\t{fna}\t{faa}\t{gff}\n")

            lf.write(f"[prok_bin_proteins] {len(manifest_rows)} genome unit(s) in manifest\n")

        Path(str(output.done)).touch()


def _read_manifest(path):
    """Yield (name, mode, fna, faa, gff) tuples from a prok_bin_proteins manifest."""
    import os
    if not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 5:
                continue
            yield parts[0], parts[1], parts[2], parts[3], parts[4]


def _concat_proteins(manifest_path, out_faa):
    """Concatenate every genome unit's proteins into one FASTA, prefixing
    each header with its genome/bin name so gene-level hits (AMRFinderPlus,
    RGI, DeepARG) can be attributed back to a bin or the contig fallback."""
    import os
    wrote_any = False
    with open(out_faa, "w") as out_f:
        for name, mode, fna, faa, gff in _read_manifest(manifest_path):
            if not os.path.exists(faa) or os.path.getsize(faa) == 0:
                continue
            with open(faa) as f:
                for line in f:
                    if line.startswith(">"):
                        out_f.write(f">{name}__{line[1:]}")
                    else:
                        out_f.write(line)
            wrote_any = True
    return wrote_any


rule checkm2:
    """
    CheckM2 — completeness and contamination of Binette final bins.
    Uses universal ML model — no lineage-specific marker database needed.
    MIMAG thresholds:
      High-quality   : completeness >= 90%, contamination <= 5%
      Medium-quality : completeness >= 50%, contamination <= 10%
    """
    input:
        done = rules.binette.output.done,
    output:
        report = f"{OUTDIR}/{{sample}}/bins/checkm2/quality_report.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/checkm2.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/checkm2.tsv"
    conda: "../envs/env_checkm2.yaml"
    container:  CONTAINERS.get("checkm2")
    threads: THREADS
    params:
        bins_dir = lambda wc: f"{OUTDIR}/{wc.sample}/bins/binette/final_bins",
        outdir   = f"{OUTDIR}/{{sample}}/bins/checkm2",
    shell:
        """
        rm -rf {params.outdir}
        mkdir -p {params.outdir}

        N_BINS=$(find {params.bins_dir} -maxdepth 1 -name "*.fa" 2>/dev/null | wc -l)
        if [ "$N_BINS" -eq 0 ]; then
            echo "[checkm2] No bins found — skipping" | tee {log}
            printf "Name\tCompleteness\tContamination\tGenome_Size\n" > {output.report}
            exit 0
        fi

        checkm2 predict \
            --threads {threads} \
            --input {params.bins_dir} \
            --output-directory {params.outdir} \
            -x fa \
            --database_path {CHECKM2_DB} \
            > {log} 2>&1 || echo "[checkm2] WARNING: predict failed" | tee -a {log}

        # Ensure output exists
        [ -f {output.report} ] || \
            printf "Name\tCompleteness\tContamination\tGenome_Size\n" > {output.report}
        """


rule gunc:
    """
    GUNC — detect chimeric MAGs by checking taxonomic consistency across
    Diamond-annotated genes (Orakov et al. 2021).
    Runs on Binette final_bins in parallel with CheckM2. The report flags
    chimeras for manual review; bins are NOT auto-excluded here.

    Database (~13 GB) must be downloaded once — see INSTALL.md:
        gunc download_db -db progenomes <dir>
    """
    input:
        done = rules.binette.output.done,
    output:
        merged = f"{OUTDIR}/{{sample}}/bins/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/gunc.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/gunc.tsv"
    conda: "../envs/env_gunc.yaml"
    container: CONTAINERS.get("gunc")
    threads: THREADS
    params:
        bins_dir = lambda wc: f"{OUTDIR}/{wc.sample}/bins/binette/final_bins",
        outdir   = f"{OUTDIR}/{{sample}}/bins/gunc",
        db       = GUNC_DB,
        enabled  = GUNC_ENABLED,
    shell:
        """
        mkdir -p {params.outdir}
        if [ "{params.enabled}" != "True" ]; then
            echo "[gunc] Disabled via config (gunc_enabled: false)" | tee {log}
            printf "genome\tpass.GUNC\tn_genes_called\n" > {output.merged}
            exit 0
        fi
        if [ -z "{params.db}" ] || [ ! -e "{params.db}" ]; then
            echo "[gunc] WARNING: gunc_db not set or missing — skipping" | tee {log}
            printf "genome\tpass.GUNC\tn_genes_called\n" > {output.merged}
            exit 0
        fi
        N_BINS=$(find {params.bins_dir} -maxdepth 1 -name "*.fa" 2>/dev/null | wc -l)
        if [ "$N_BINS" -eq 0 ]; then
            echo "[gunc] No bins to evaluate" | tee {log}
            printf "genome\tpass.GUNC\tn_genes_called\n" > {output.merged}
            exit 0
        fi
        gunc run \
            --input_dir {params.bins_dir} \
            --out_dir   {params.outdir} \
            --db_file   {params.db} \
            --threads   {threads} \
            --file_suffix .fa \
            > {log} 2>&1 || echo "[gunc] WARNING: gunc run failed" | tee -a {log}
        # Locate produced TSV (filename includes DB version)
        for cand in {params.outdir}/GUNC.progenomes_2.1.maxCSS_level.tsv \
                    {params.outdir}/GUNC.*.maxCSS_level.tsv; do
            [ -f "$cand" ] && cp "$cand" {output.merged} && break
        done
        [ -f {output.merged} ] || \
            printf "genome\tpass.GUNC\tn_genes_called\n" > {output.merged}
        """


rule galah_derep:
    """
    galah — CheckM2-aware MAG dereplication by ANI (Woodcroft).
    Clusters Binette final_bins at MAG_DEREP_ANI; keeps the bin with the
    highest CheckM2 quality score per cluster (completeness − 5×contamination).
    Output dir feeds GTDB-Tk so taxonomy runs only on representatives.
    """
    input:
        bins_done   = rules.binette.output.done,
        checkm2_tsv = rules.checkm2.output.report,
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/derep/done.txt",
        cluster = f"{OUTDIR}/{{sample}}/bins/derep/galah_clusters.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/galah_derep.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/galah_derep.tsv"
    conda: "../envs/env_derep.yaml"
    container: CONTAINERS.get("galah")
    threads: THREADS
    params:
        bins_dir = lambda wc: f"{OUTDIR}/{wc.sample}/bins/binette/final_bins",
        outdir   = f"{OUTDIR}/{{sample}}/bins/derep",
        repdir   = f"{OUTDIR}/{{sample}}/bins/derep/derep_bins",
        ani      = MAG_DEREP_ANI,
        enabled  = MAG_DEREP_ENABLED,
    shell:
        """
        mkdir -p {params.outdir} {params.repdir}
        if [ "{params.enabled}" != "True" ]; then
            echo "[galah] Disabled via config — symlinking original bins" | tee {log}
            # Mirror binette/final_bins into derep_bins for a uniform GTDB-Tk input
            shopt -s nullglob
            for fa in {params.bins_dir}/*.fa; do
                ln -sf "$fa" {params.repdir}/
            done
            printf "representative\tmember\n" > {output.cluster}
            touch {output.done}; exit 0
        fi
        shopt -s nullglob
        BINS=({params.bins_dir}/*.fa)
        if [ ${{#BINS[@]}} -eq 0 ]; then
            echo "[galah] No bins to dereplicate" | tee {log}
            printf "representative\tmember\n" > {output.cluster}
            touch {output.done}; exit 0
        fi
        galah cluster \
            --genome-fasta-files "${{BINS[@]}}" \
            --ani {params.ani} \
            --checkm2-quality-report {input.checkm2_tsv} \
            --output-cluster-definition {output.cluster} \
            --output-representative-fasta-directory {params.repdir} \
            --threads {threads} \
            > {log} 2>&1 || echo "[galah] WARNING: cluster failed — falling back to original bins" | tee -a {log}
        # Fallback: if galah produced nothing usable, mirror original bins
        if [ -z "$(ls {params.repdir}/*.fa 2>/dev/null)" ]; then
            for fa in {params.bins_dir}/*.fa; do
                ln -sf "$fa" {params.repdir}/
            done
            [ -s {output.cluster} ] || printf "representative\tmember\n" > {output.cluster}
        fi
        touch {output.done}
        """


rule gtdbtk:
    """
    GTDB-Tk classify_wf: assign GTDB taxonomy to MAGs.
    Replaces NCBI taxonomy for metagenome-assembled genomes.
    Outputs: bac120.summary.tsv + ar53.summary.tsv
    Creates empty outputs if no bins are available (safe fallback).

    Input bins:
      - MAG_DEREP_ENABLED=True  → bins/derep/derep_bins/   (galah representatives)
      - MAG_DEREP_ENABLED=False → bins/binette/final_bins/ (all Binette bins)
    """
    input:
        checkm2_done = rules.checkm2.output.report,
        derep_done   = (rules.galah_derep.output.done if MAG_DEREP_ENABLED else []),
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/gtdbtk/done.txt",
        bac_tsv = f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.bac120.summary.tsv",
        ar_tsv  = f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/gtdbtk.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/gtdbtk.tsv"
    conda: "../envs/env_gtdbtk.yaml"
    container:  CONTAINERS.get("gtdbtk")
    threads: THREADS
    params:
        bins_dir = _gtdbtk_bins_dir,
        outdir   = f"{OUTDIR}/{{sample}}/bins/gtdbtk",
    shell:
        """
        mkdir -p {params.outdir}

        # If no bins, create empty outputs
        N_BINS=$(find {params.bins_dir} -maxdepth 1 -name "*.fa" 2>/dev/null | wc -l)
        if [ "$N_BINS" -eq 0 ]; then
            echo "[gtdbtk] No bins found — skipping" | tee {log}
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
            --extension  fa \
            >> {log} 2>&1 || echo "[gtdbtk] WARNING: classify_wf failed — creating empty outputs" | tee -a {log}

        # Always ensure output files exist — gtdbtk may fail if bins are low quality
        mkdir -p {params.outdir}/classify
        [ -f {output.bac_tsv} ] || printf "user_genome\tclassification\n" > {output.bac_tsv}
        [ -f {output.ar_tsv}  ] || printf "user_genome\tclassification\n" > {output.ar_tsv}
        touch {output.done}
        """
