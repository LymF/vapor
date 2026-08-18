# ══════════════════════════════════════════════════════════════════════
# rules/annotation.smk — BLOCK 6: Genome Annotation
#
# bakta       — prokaryotic MAG annotation (HQ/MQ bins from CheckM2)
# eggnog_prok — COG/KEGG/CAZy/GO functional annotation of MAG proteins
# extract_kegg_kos — KO extraction per MAG from EggNOG output (KOALA-format TSV)
# genome_map_prok — circular genome maps for prokaryotic MAGs (per sample:
#                    prokaryote bins are not deduplicated into a global
#                    catalog the way vOTUs are, so there is no global
#                    representative to compute this against)
#
# pharokka/phold/genome_map_phage/genome_map_virus moved to
# rules/votu_catalog.smk (votu_pharokka/votu_phold/votu_genome_map_phage/
# votu_genome_map_virus) on 2026-08-18 -- they depend only on the vOTU
# representative SEQUENCE, never on which sample it came from, so they now
# run once over the global catalog instead of once per sample. See
# docs/ROADMAP_SIMPLIFICACAO.md "(h) Princípio: computar no representante,
# herdar no membro".
#
# All rules soft-fail (touch output) when their database is not configured.
# Config keys: bakta_db, bakta_min_completeness,
#              bakta_max_contamination, eggnog_db
# ══════════════════════════════════════════════════════════════════════


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
        # derivado do output, nao de {{sample}} (requisito da heranca).
        outdir    = lambda wc, output: os.path.dirname(output.done),
        # bins_dir/bin_ext sao os dois pontos onde a trilha de grupo difere:
        # per-sample usa bins do Binette (*.fa), co-assembly usa VAMB (*.fna).
        # Parametrizados para coassembly.smk poder herdar esta regra via
        # `use rule ... as ... with:` em vez de manter uma segunda copia.
        bins_dir  = f"{OUTDIR}/{{sample}}/bins/binette/final_bins",
        bin_ext   = ".fa",
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
            BIN_FA="{params.bins_dir}/$BIN_NAME{params.bin_ext}"
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
        # derivados do output, nao de {{sample}}: regra herdada por
        # coassembly.smk via `use rule ... as ... with:` (wildcard {group}).
        outdir      = lambda wc, output: os.path.dirname(output.done),
        bakta_dir   = lambda wc, output: os.path.join(os.path.dirname(os.path.dirname(output.done)), "bakta"),
        all_faa     = lambda wc, output: os.path.join(os.path.dirname(output.done), "all_mags.faa"),
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
        # derivado do output, nao de {{sample}}: regra herdada por
        # coassembly.smk via `use rule ... as ... with:` (wildcard {group}).
        outdir = lambda wc, output: os.path.dirname(output.done),
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
            > {log} 2>&1 && RC=0 || RC=$?
        if [ "$RC" -ne 0 ]; then
            echo "failed: genome_map_universal.py exit $RC" > {output.done}
        else
            echo "ok" > {output.done}
        fi
        """
