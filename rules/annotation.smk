# ══════════════════════════════════════════════════════════════════════
# rules/annotation.smk — BLOCK 6: Genome Annotation
#
# pharokka    — bacteriophage annotation (PHROGS) on HQ CheckV phages
# phold       — structure-based annotation of viral hypothetical proteins
# bakta       — prokaryotic MAG annotation (HQ/MQ bins from CheckM2)
# eggnog_prok — COG/KEGG/CAZy/GO functional annotation of MAG proteins
# extract_kegg_kos — KO extraction per MAG from EggNOG output (KOALA-format TSV)
#
# All rules soft-fail (touch output) when their database is not configured.
# Config keys: pharokka_db, phold_db, pharokka_min_completeness,
#              bakta_db, bakta_min_completeness,
#              bakta_max_contamination, eggnog_db
# ══════════════════════════════════════════════════════════════════════


rule pharokka:
    """
    Bacteriophage genome annotation with Pharokka (PHROGS database).
    Runs in --meta mode on ALL viral sequences >= PHAROKKA_MIN_COMPLETENESS.

    Selection criteria:
      - completeness >= PHAROKKA_MIN_COMPLETENESS  (quality tier NOT required —
        novel phages often receive "Not-determined" from CheckV despite high
        completeness because no close reference cluster exists)
      - No cap on genome count: every qualifying sequence is annotated, since
        the resulting PHROGS hallmark-gene hits are what genome_map_universal.py
        uses to decide phage vs. non-phage virus (see genome_map_phage/_virus
        below) — taxonomy alone is unreliable for novel/divergent sequences.

    Output GBK/TSV are later used by phold and genome_map_universal.py.
    Skipped if PHAROKKA_DB is not configured (empty string).
    """
    input:
        viral_nr = rules.viral_votu_reps.output.mq_fasta,
        checkv   = rules.checkv.output.summary,
    output:
        done = f"{OUTDIR}/{{sample}}/annotation/pharokka/done.txt",
        gbk  = f"{OUTDIR}/{{sample}}/annotation/pharokka/pharokka.gbk",
        tsv  = f"{OUTDIR}/{{sample}}/annotation/pharokka/pharokka_cds_final_merged_output.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/pharokka.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/pharokka.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("pharokka")
    threads: THREADS
    params:
        outdir   = f"{OUTDIR}/{{sample}}/annotation/pharokka",
        db       = PHAROKKA_DB,
        min_comp = PHAROKKA_MIN_COMPLETENESS,
        # NOT inside outdir: pharokka.py --force deletes/recreates its own
        # -o directory on startup, which would delete this -i input too if
        # it lived underneath it.
        hq_fa    = f"{OUTDIR}/{{sample}}/annotation/pharokka_hq_phages.fasta",
    run:
        import csv, os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        log_path = str(log[0])

        def touch_empty(path):
            Path(path).touch()

        if not params.db or not os.path.isdir(str(params.db)):
            with open(log_path, "w") as lf:
                lf.write("[pharokka] PHAROKKA_DB not configured — skipping\n")
            touch_empty(output.done)
            touch_empty(output.gbk)
            touch_empty(output.tsv)
            return

        # Filter CheckV by completeness only.
        # Novel phages (Caudovirales, etc.) often receive "Not-determined" quality
        # from CheckV even when completeness is high, because CheckV cannot assign
        # a quality tier without a close reference cluster.  Filtering by quality
        # tier would silently drop these bona-fide phages from Pharokka annotation.
        hq_ids = []
        with open(str(input.checkv)) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                try:
                    comp = float(row.get("completeness", "0") or 0)
                except ValueError:
                    # CheckV writes "NA" for completeness on some
                    # "Not-determined" quality contigs -- treat as 0.
                    comp = 0.0
                if comp >= float(params.min_comp):
                    hq_ids.append((row["contig_id"], comp))
        hq_ids.sort(key=lambda x: x[1], reverse=True)
        hq_set = {cid for cid, _ in hq_ids}

        # Extract sequences
        with open(str(params.hq_fa), "w") as out_fa, \
             open(str(input.viral_nr)) as in_fa, \
             open(log_path, "w") as lf:
            lf.write(f"[pharokka] Phages selected (completeness >= {params.min_comp}%): {len(hq_set)}\n")
            write = False
            for line in in_fa:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    write = name in hq_set
                if write:
                    out_fa.write(line)

        if not hq_set or os.path.getsize(str(params.hq_fa)) == 0:
            with open(log_path, "a") as lf:
                lf.write("[pharokka] No HQ phages found — skipping\n")
            touch_empty(output.done)
            touch_empty(output.gbk)
            touch_empty(output.tsv)
            return

        # --meta/--meta_hmm are multi-FASTA only; pharokka.py refuses to run
        # with --meta on a single-contig input ("ERROR: -m meta mode
        # specified when the input file only contains 1 contig").
        with open(str(params.hq_fa)) as f:
            n_seqs = sum(1 for line in f if line.startswith(">"))
        meta_flags = "--meta --meta_hmm" if n_seqs > 1 else ""
        with open(log_path, "a") as lf:
            lf.write(f"[pharokka] {n_seqs} sequence(s) in input — "
                      f"{'meta' if n_seqs > 1 else 'single-genome'} mode\n")

        shell(
            "pharokka.py"
            " -i {params.hq_fa}"
            " -o {params.outdir}"
            " -d {params.db}"
            " -t {threads}"
            " {meta_flags}"
            " --dnaapler"
            " --force"
            " >> {log} 2>&1"
        )

        # Standardize output filenames — Pharokka may name them differently in meta mode
        for candidate in [
            f"{params.outdir}/pharokka.gbk",
            f"{params.outdir}/pharokka_meta.gbk",
        ]:
            if os.path.exists(candidate):
                if candidate != str(output.gbk):
                    os.rename(candidate, str(output.gbk))
                break
        else:
            touch_empty(output.gbk)

        for candidate in [
            f"{params.outdir}/pharokka_cds_final_merged_output.tsv",
            f"{params.outdir}/pharokka_annotations.tsv",
        ]:
            if os.path.exists(candidate):
                if candidate != str(output.tsv):
                    os.rename(candidate, str(output.tsv))
                break
        else:
            touch_empty(output.tsv)

        shell("touch {output.done}")


rule phold:
    """
    Phold — structure-based annotation of viral hypothetical proteins.
    Takes Pharokka GBK as input; updates annotations via ProstT5 protein
    structure similarity. Annotates ~40–60% of hypothetical ORFs missed
    by PHROGS sequence similarity.
    Produces updated GBK used by genome_map_universal.py.
    Skipped if Pharokka found no HQ phages (gbk is empty).
    """
    input:
        pharokka_done = rules.pharokka.output.done,
        pharokka_gbk  = rules.pharokka.output.gbk,
    output:
        done = f"{OUTDIR}/{{sample}}/annotation/phold/done.txt",
        gbk  = f"{OUTDIR}/{{sample}}/annotation/phold/phold.gbk",
    log:
        f"{OUTDIR}/{{sample}}/logs/phold.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/phold.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("phold")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/{{sample}}/annotation/phold",
        db     = PHOLD_DB,
    shell:
        """
        mkdir -p {params.outdir}
        if [ ! -s {input.pharokka_gbk} ]; then
            echo "[phold] Pharokka GBK empty — skipping" | tee {log}
            touch {output.gbk} {output.done}; exit 0
        fi
        DB_FLAG=""
        [ -n "{params.db}" ] && [ -d "{params.db}" ] && DB_FLAG="-d {params.db}"
        phold run \
            -i {input.pharokka_gbk} \
            -o {params.outdir} \
            -t {threads} \
            --cpu \
            --hyps \
            --force \
            $DB_FLAG \
            > {log} 2>&1 || true
        # phold outputs phold.gbk in the output dir
        [ -f {params.outdir}/phold.gbk ] && \
            cp {params.outdir}/phold.gbk {output.gbk} || \
            touch {output.gbk}
        touch {output.done}
        """


rule bakta:
    """
    Bakta — prokaryotic MAG annotation (replaces Prokka).
    Annotates HQ/MQ bins from CheckM2/Binette that pass quality thresholds.
    Loops over all qualifying MAGs within a single rule execution.
    Generates GBK (for genome maps) and TSV (for EggNOG/functional analysis).

    Thresholds (config.yaml):
      bakta_min_completeness:  70.0  (%)
      bakta_max_contamination: 10.0  (%)
    Skipped if BAKTA_DB is not configured (empty string).
    """
    input:
        checkm2 = rules.checkm2.output.report,
        binette = rules.binette.output.done,
    output:
        done    = f"{OUTDIR}/{{sample}}/annotation/bakta/done.txt",
        summary = f"{OUTDIR}/{{sample}}/annotation/bakta/bakta_summary.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/bakta.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/bakta.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("bakta")
    threads: THREADS
    params:
        outdir    = f"{OUTDIR}/{{sample}}/annotation/bakta",
        bins_dir  = f"{OUTDIR}/{{sample}}/bins/binette/final_bins",
        min_comp  = BAKTA_MIN_COMPLETENESS,
        max_cont  = BAKTA_MAX_CONTAMINATION,
    shell:
        """
        mkdir -p {params.outdir}
        export BAKTA_DB="{BAKTA_DB}"

        if [ -z "{BAKTA_DB}" ] || [ ! -d "{BAKTA_DB}" ]; then
            echo "[bakta] BAKTA_DB not configured — skipping" | tee {log}
            echo -e "bin\tstatus" > {output.summary}
            touch {output.done}; exit 0
        fi

        # Identify qualifying MAGs from CheckM2 report
        python3 - <<'PYEOF'
import csv, sys
qualifying = []
with open("{input.checkm2}") as f:
    for row in csv.DictReader(f, delimiter="\\t"):
        name = row.get("Name", "").strip()
        comp = float(row.get("Completeness", 0) or 0)
        cont = float(row.get("Contamination", 100) or 100)
        if comp >= {params.min_comp} and cont <= {params.max_cont}:
            qualifying.append(name)
with open("{params.outdir}/qualifying_bins.txt", "w") as f:
    f.write("\\n".join(qualifying) + "\\n")
print(f"[bakta] Qualifying MAGs: {{len(qualifying)}}", file=sys.stderr)
PYEOF

        N_QUAL=$(wc -l < {params.outdir}/qualifying_bins.txt 2>/dev/null || echo 0)
        echo "[bakta] $N_QUAL qualifying MAGs (>={params.min_comp}% comp, <={params.max_cont}% cont)" \
            | tee -a {log}

        if [ "$N_QUAL" -eq 0 ]; then
            echo "[bakta] No qualifying MAGs — skipping" | tee -a {log}
            echo -e "bin\tstatus" > {output.summary}
            touch {output.done}; exit 0
        fi

        echo -e "bin\tstatus" > {output.summary}

        while IFS= read -r BIN_NAME || [ -n "$BIN_NAME" ]; do
            [ -z "$BIN_NAME" ] && continue
            BIN_FA="{params.bins_dir}/$BIN_NAME.fa"
            [ -f "$BIN_FA" ] || {{ printf "%s\tmissing\n" "$BIN_NAME" >> {output.summary}; continue; }}

            BIN_OUT="{params.outdir}/$BIN_NAME"
            mkdir -p "$BIN_OUT"
            bakta \
                --db {BAKTA_DB} \
                --output "$BIN_OUT" \
                --prefix "$BIN_NAME" \
                --threads {threads} \
                --meta \
                --force \
                --skip-plot \
                "$BIN_FA" \
                >> {log} 2>&1 && \
                echo "$BIN_NAME\tok" >> {output.summary} || \
                echo "$BIN_NAME\tfailed" >> {output.summary}
        done < {params.outdir}/qualifying_bins.txt

        touch {output.done}
        echo "[bakta] Done — $(grep -c 'ok' {output.summary}) MAGs annotated" | tee -a {log}
        """


rule eggnog_prok:
    """
    EggNOG-mapper v2 — COG/KEGG/CAZy/GO functional annotation of MAG proteins.
    Concatenates FAA files from all Bakta-annotated MAGs, then runs emapper.py.
    Standard for MAG functional characterisation; required for publication in
    high-impact journals (ATLAS, Aviary, SqueezeMeta all include this step).
    Skipped if EGGNOG_DB is not configured or no Bakta FAA files exist.
    """
    input:
        bakta_done = rules.bakta.output.done,
    output:
        done      = f"{OUTDIR}/{{sample}}/annotation/eggnog/done.txt",
        annot_tsv = f"{OUTDIR}/{{sample}}/annotation/eggnog/eggnog_annotations.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/eggnog_prok.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/eggnog_prok.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("eggnog_mapper")
    threads: THREADS
    params:
        outdir      = f"{OUTDIR}/{{sample}}/annotation/eggnog",
        bakta_dir   = f"{OUTDIR}/{{sample}}/annotation/bakta",
        all_faa     = f"{OUTDIR}/{{sample}}/annotation/eggnog/all_mags.faa",
    shell:
        """
        mkdir -p {params.outdir}

        if [ -z "{EGGNOG_DB}" ] || [ ! -d "{EGGNOG_DB}" ]; then
            echo "[eggnog] EGGNOG_DB not configured — skipping" | tee {log}
            touch {output.annot_tsv} {output.done}; exit 0
        fi

        # Concatenate Bakta FAA files (one per MAG subdirectory)
        rm -f {params.all_faa}
        for FAA in {params.bakta_dir}/*/*.faa; do
            [ -f "$FAA" ] && [ -s "$FAA" ] && cat "$FAA" >> {params.all_faa}
        done

        if [ ! -s {params.all_faa} ]; then
            echo "[eggnog] No Bakta FAA files found — skipping" | tee -a {log}
            touch {output.annot_tsv} {output.done}; exit 0
        fi

        N_PROT=$(grep -c "^>" {params.all_faa} || echo 0)
        echo "[eggnog] Running on $N_PROT proteins" | tee -a {log}

        emapper.py \
            -m diamond \
            --itype proteins \
            -i {params.all_faa} \
            -o eggnog_annotations \
            --output_dir {params.outdir} \
            --cpu {threads} \
            --data_dir {EGGNOG_DB} \
            --override \
            >> {log} 2>&1

        [ -f {params.outdir}/eggnog_annotations.emapper.annotations ] && \
            cp {params.outdir}/eggnog_annotations.emapper.annotations {output.annot_tsv} || \
            touch {output.annot_tsv}

        touch {output.done}
        echo "[eggnog] Done" | tee -a {log}
        """


rule extract_kegg_kos:
    """
    Extrai KO numbers por MAG a partir do output do EggNOG-mapper.
    Produz ko_per_mag.tsv (gene_id TAB KO) pronto para KEGG-Decoder,
    IPATH3, ou qualquer análise downstream de vias metabólicas.
    Skipped if eggnog annotations are empty.
    """
    input:
        eggnog_done = rules.eggnog_prok.output.done,
        annot_tsv   = rules.eggnog_prok.output.annot_tsv,
    output:
        done     = f"{OUTDIR}/{{sample}}/annotation/kegg_decoder/done.txt",
        ko_table = f"{OUTDIR}/{{sample}}/annotation/kegg_decoder/ko_per_mag.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/extract_kegg_kos.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/extract_kegg_kos.tsv"
    conda: "../envs/env_annotation.yaml"
    threads: 1
    params:
        outdir = f"{OUTDIR}/{{sample}}/annotation/kegg_decoder",
    run:
        import re, os, sys
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        log_path = str(log[0])

        if not os.path.exists(str(input.annot_tsv)) or \
           os.path.getsize(str(input.annot_tsv)) == 0:
            with open(log_path, "w") as lf:
                lf.write("[extract_kegg_kos] EggNOG annotations empty — skipping\n")
            Path(str(output.ko_table)).touch()
            Path(str(output.done)).touch()
            return

        records = []
        with open(str(input.annot_tsv)) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 12:
                    continue
                query   = cols[0]
                kegg_ko = cols[11]
                if kegg_ko == "-" or not kegg_ko.strip():
                    continue
                mag = re.sub(r"_CDS_\d+$", "", query)
                mag = re.sub(r"_\d+$", "", mag)
                for ko_entry in kegg_ko.split(","):
                    ko = ko_entry.strip().replace("ko:", "")
                    if re.match(r"K\d{5}", ko):
                        records.append((mag, ko))

        with open(str(output.ko_table), "w") as f:
            f.write("mag\tko\n")
            for mag, ko in records:
                f.write(f"{mag}\t{ko}\n")

        n_mags = len({r[0] for r in records})
        with open(log_path, "w") as lf:
            lf.write(f"[extract_kegg_kos] {len(records)} KO entries, {n_mags} MAGs\n")
            lf.write(f"[extract_kegg_kos] Output: {output.ko_table}\n")

        Path(str(output.done)).touch()


rule genome_map_phage:
    """
    Circular genome maps for confirmed phages (PHROGS color scheme).
    Reads Phold GBK (multi-record, every sequence pharokka annotated) and the
    Pharokka CDS TSV; a genome is only drawn as a phage if it has >=1 PHROGS
    structural/packaging hallmark gene (head-and-packaging, connector, tail,
    lysis) — this is what actually distinguishes phages from other high-
    completeness viral sequences that happened to pass through pharokka.
    The full confirmed-phage ID list (not just the displayed top-N) is written
    to confirmed_phage_ids.txt, consumed by genome_map_virus below to exclude
    real phages regardless of completeness ranking.
    Skipped if Phold GBK is empty (no qualifying sequences in sample).
    """
    input:
        phold_done   = rules.phold.output.done,
        phold_gbk    = rules.phold.output.gbk,
        pharokka_tsv = rules.pharokka.output.tsv,
        checkv       = rules.checkv.output.summary,
        viral_nr     = rules.viral_votu_reps.output.mq_fasta,
    output:
        done = f"{OUTDIR}/{{sample}}/annotation/genome_maps/phage_maps_done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/genome_map_phage.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/genome_map_phage.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("genome_map")
    threads: 1
    params:
        outdir      = f"{OUTDIR}/{{sample}}/annotation/genome_maps/phage",
        min_comp    = GENOME_MAP_MIN_COMP_VIRAL,
        top_n       = GENOME_MAP_TOP_N,
        scripts_dir = SCRIPTS_DIR,
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.scripts_dir}/genome_map_universal.py \
            --mode phage \
            --gbk {input.phold_gbk} \
            --pharokka-tsv {input.pharokka_tsv} \
            --fasta {input.viral_nr} \
            --checkv {input.checkv} \
            --outdir {params.outdir} \
            --min-completeness {params.min_comp} \
            --top-n {params.top_n} \
            > {log} 2>&1 || true
        touch {output.done}
        """


rule genome_map_virus:
    """
    Circular genome maps for non-phage viral sequences.
    Uses GeNomad gene predictions for gene positions and functional categories.
    Excludes confirmed phages — the FULL set written by genome_map_phage to
    confirmed_phage_ids.txt (every sequence with a PHROGS hallmark gene, not
    just the top-N that got drawn), so a real phage ranked below the display
    cutoff cannot leak into the virus map. Anything else >= min completeness
    qualifies, including sequences with no taxonomy assignment — pharokka
    already had the chance to find phage evidence and found none.
    Selects top-N viral genomes by CheckV completeness.
    """
    input:
        checkv     = rules.checkv.output.summary,
        viral_nr   = rules.viral_votu_reps.output.mq_fasta,
        phage_done = rules.genome_map_phage.output.done,
        genomad    = f"{OUTDIR}/{{sample}}/viral/genomad/done.txt",
    output:
        done = f"{OUTDIR}/{{sample}}/annotation/genome_maps/virus_maps_done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/genome_map_virus.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/genome_map_virus.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("genome_map")
    threads: 1
    params:
        outdir      = f"{OUTDIR}/{{sample}}/annotation/genome_maps/virus",
        phage_dir   = f"{OUTDIR}/{{sample}}/annotation/genome_maps/phage",
        genomad_dir = f"{OUTDIR}/{{sample}}/viral/genomad",
        tax_dir     = f"{OUTDIR}/{{sample}}/viral/taxonomy",
        min_comp    = GENOME_MAP_MIN_COMP_VIRAL,
        top_n       = GENOME_MAP_TOP_N,
        scripts_dir = SCRIPTS_DIR,
    run:
        import os, shutil
        from pathlib import Path

        # Wipe stale output: with a relative script path this rule used to fail
        # silently whenever Snakemake's cwd wasn't the pipeline directory (see
        # genome_map_virus.log), leaving old SVGs from a prior run/classification
        # logic (e.g. before the PHROGS-hallmark phage/virus split) sitting here
        # forever -- a contig later confirmed as a phage would still show its
        # outdated "Virus" map. Always start from a clean directory.
        shutil.rmtree(params.outdir, ignore_errors=True)
        os.makedirs(params.outdir, exist_ok=True)
        log_path = str(log[0])

        # Exclude-IDs file: the FULL confirmed-phage list from genome_map_phage
        # (confirmed_phage_ids.txt), not just the genomes it drew.
        exclude_file = os.path.join(params.outdir, "phage_ids_to_exclude.txt")
        confirmed_file = os.path.join(params.phage_dir, "confirmed_phage_ids.txt")
        confirmed_ids = []
        if os.path.exists(confirmed_file):
            with open(confirmed_file) as cf:
                confirmed_ids = [line.strip() for line in cf if line.strip()]

        with open(exclude_file, "w") as ef:
            ef.write("\n".join(confirmed_ids) + ("\n" if confirmed_ids else ""))

        # Locate GeNomad gene summary TSV (genes in GFF-like format)
        genomad_genes = ""
        for candidate in Path(params.genomad_dir).rglob("*_genes.tsv"):
            genomad_genes = str(candidate)
            break

        # Viral taxonomy TSV (optional — enriches genome map titles with class/family)
        tax_tsv = os.path.join(params.tax_dir, "viral_taxonomy_merged.tsv")

        shell(
            "python3 {params.scripts_dir}/genome_map_universal.py"
            " --mode virus"
            " --fasta {input.viral_nr}"
            + (f" --genomad-genes {genomad_genes}" if genomad_genes else "") +
            " --checkv {input.checkv}"
            " --exclude-ids " + exclude_file +
            " --outdir {params.outdir}"
            " --min-completeness {params.min_comp}"
            " --top-n {params.top_n}"
            + (f" --viral-taxonomy {tax_tsv}" if os.path.exists(tax_tsv) else "") +
            " >> {log} 2>&1 || true"
        )
        Path(str(output.done)).touch()


rule genome_map_prok:
    """
    Circular genome maps for prokaryotic MAGs (COG color scheme).
    Reads Bakta GBK files for each qualifying MAG (≥90% completeness, ≤5%
    contamination, ≤5 contigs). Multi-contig MAGs are drawn as multi-sector
    Circos plots. Genes colored by keyword-inferred COG functional category.
    Skipped if no Bakta outputs exist.
    """
    input:
        bakta_done = rules.bakta.output.done,
        checkm2    = rules.checkm2.output.report,
    output:
        done = f"{OUTDIR}/{{sample}}/annotation/genome_maps/prok_maps_done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/genome_map_prok.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/genome_map_prok.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("genome_map")
    threads: 1
    params:
        outdir      = f"{OUTDIR}/{{sample}}/annotation/genome_maps/prok",
        bakta_dir   = f"{OUTDIR}/{{sample}}/annotation/bakta",
        min_comp    = GENOME_MAP_MIN_COMP_PROK,
        max_cont    = GENOME_MAP_MAX_CONT_PROK,
        max_contigs = GENOME_MAP_MAX_CONTIGS_PROK,
        top_n       = GENOME_MAP_TOP_N,
        scripts_dir = SCRIPTS_DIR,
    shell:
        """
        mkdir -p {params.outdir}
        python3 {params.scripts_dir}/genome_map_universal.py \
            --mode prok \
            --bakta-dir {params.bakta_dir} \
            --checkm2 {input.checkm2} \
            --outdir {params.outdir} \
            --min-completeness {params.min_comp} \
            --max-contamination {params.max_cont} \
            --max-contigs {params.max_contigs} \
            --top-n {params.top_n} \
            > {log} 2>&1 || true
        touch {output.done}
        """
