# ══════════════════════════════════════════════════════════════════════
# rules/annotation.smk — BLOCK 6: Genome Annotation
#
# bakta       — prokaryotic MAG annotation (HQ/MQ bins from CheckM2)
# eggnog_prok — COG/KEGG/CAZy/GO functional annotation of MAG proteins
# extract_kegg_kos — KO extraction per MAG from EggNOG output (KOALA-format TSV)
#
# pharokka/phold moved to rules/votu_catalog.smk (votu_pharokka/votu_phold) on
# 2026-08-18 -- dependem so da SEQUENCIA do representante de vOTU, nunca de qual
# amostra ela veio, entao rodam uma vez sobre o catalogo global. Ver
# docs/ROADMAP_SIMPLIFICACAO.md "(h) Princípio: computar no representante,
# herdar no membro".
#
# Os genome maps (genome_map_prok/phage/virus e scripts/genome_map*.py) foram
# REMOVIDOS em 2026-08-19 a pedido do usuario. Levaram junto duas regras que so
# existiam para alimenta-los: votu_catalog_quality_summary e
# votu_catalog_genomad_genes.
#
# All rules soft-fail (touch output) when their database is not configured.
# Config keys: bakta_db, bakta_min_completeness,
#              bakta_max_contamination, eggnog_db
# ══════════════════════════════════════════════════════════════════════


rule mag_bakta:
    """
    Bakta — prokaryotic MAG annotation (replaces Prokka).

    GLOBAL desde 2026-08-19: anota as REPRESENTANTES do catálogo
    (rules/mag_catalog.smk) que passam no corte de qualidade, uma vez cada,
    e as vistas por amostra/grupo distribuem o sumário. O corte usa o
    `checkm2_quality_report.tsv` re-chaveado do catálogo, cuja coluna `Name`
    já são os IDs namespaced — os mesmos nomes dos FASTAs em
    `representatives/`.

    Thresholds (config.yaml):
      bakta_min_completeness:  70.0  (%)
      bakta_max_contamination: 10.0  (%)
    Skipped if BAKTA_DB is not configured (empty string).
    """
    input:
        checkm2 = rules.mag_catalog_quality.output.tsv,
        derep   = rules.mag_catalog_derep.output.done,
    output:
        done    = f"{MAG_CATALOG_DIR}/bakta/done.txt",
        summary = f"{MAG_CATALOG_DIR}/bakta/bakta_summary.tsv",
    log:
        f"{OUTDIR}/logs/mag_bakta.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_bakta.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("bakta")
    threads: THREADS
    params:
        outdir    = lambda wc, output: os.path.dirname(output.done),
        # O pool normaliza toda extensao para .fa, entao aqui nao ha mais a
        # bifurcacao Binette (*.fa) / VAMB (*.fna) que existia quando esta
        # regra era herdada pela trilha de grupo.
        bins_dir  = f"{MAG_CATALOG_DIR}/representatives",
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

        printf "bin\tstatus\n" > {output.summary}

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
                printf "%s\tok\n" "$BIN_NAME" >> {output.summary} || \
                printf "%s\tfailed\n" "$BIN_NAME" >> {output.summary}
        done < {params.outdir}/qualifying_bins.txt

        touch {output.done}
        echo "[bakta] Done — $(grep -c 'ok' {output.summary}) MAGs annotated" | tee -a {log}
        """


rule mag_eggnog_prok:
    """
    EggNOG-mapper v2 — COG/KEGG/CAZy/GO functional annotation of MAG proteins.
    Concatenates FAA files from all Bakta-annotated MAGs, then runs emapper.py.
    Standard for MAG functional characterisation; required for publication in
    high-impact journals (ATLAS, Aviary, SqueezeMeta all include this step).
    Skipped if EGGNOG_DB is not configured or no Bakta FAA files exist.
    """
    input:
        bakta_done = rules.mag_bakta.output.done,
    output:
        done      = f"{MAG_CATALOG_DIR}/eggnog/done.txt",
        annot_tsv = f"{MAG_CATALOG_DIR}/eggnog/eggnog_annotations.tsv",
    log:
        f"{OUTDIR}/logs/mag_eggnog_prok.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_eggnog_prok.tsv"
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

        # Concatena os FAA do Bakta (um subdiretorio por MAG) PREFIXANDO
        # cada proteina com o nome do genoma. Sem isso o ID e o locus tag do
        # Bakta ("LLOGBO_00001"), que nao carrega vinculo nenhum com o MAG:
        # a saida do eggNOG ficaria sem atribuicao de genoma e nem a vista
        # por amostra nem o ko_per_mag.tsv teriam como saber de quem e cada
        # linha. Mesma convencao do _concat_proteins do lado AMR.
        rm -f {params.all_faa}
        for FAA in {params.bakta_dir}/*/*.faa; do
            [ -f "$FAA" ] && [ -s "$FAA" ] || continue
            GENOME=$(basename $(dirname "$FAA"))
            awk -v g="$GENOME" '/^>/ {{ sub(/^>/, ">" g "__"); print; next }} {{ print }}' \
                "$FAA" >> {params.all_faa}
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

        # Grava um TSV de verdade: o .emapper.annotations vem com 4 linhas
        # "##" antes do cabecalho e mais algumas no fim. Qualquer leitor de
        # TSV que tome a primeira linha como cabecalho (o load_tsv do
        # relatorio, o csv.DictReader das vistas) leria "## <data>" como
        # cabecalho e devolveria None para TODA coluna -- foi o que fez a
        # aba de COG do relatorio contar tudo como "Function unknown".
        # O arquivo original fica ao lado, intacto.
        [ -f {params.outdir}/eggnog_annotations.emapper.annotations ] && \
            sed '/^##/d' {params.outdir}/eggnog_annotations.emapper.annotations \
                > {output.annot_tsv} || \
            touch {output.annot_tsv}

        touch {output.done}
        echo "[eggnog] Done" | tee -a {log}
        """


rule mag_extract_kegg_kos:
    """
    Extrai KO numbers e famílias CAZy por MAG a partir do output do EggNOG-mapper.

    Produz ko_per_mag.tsv (mag TAB KO), consumido pelo `mag_kegg_completeness`
    logo abaixo, e cazy_per_mag.tsv (mag TAB família) — a coluna `CAZy` do
    emapper existia desde sempre no arquivo e nenhum consumidor a lia até
    2026-08-19. Sai daqui, e não de uma regra nova, porque é o MESMO parse do
    MESMO arquivo: +0 jobs no DAG e um banco a menos para instalar.

    A coluna `mag` é o nome do genoma, recuperado do prefixo `{genome}__`
    que o `mag_eggnog_prok` põe em cada proteína. Antes de 2026-08-19 ela era
    derivada por regex do ID (`LLOGBO_00001` -> `LLOGBO`), ou seja, o
    prefixo de locus tag que o Bakta sorteia — um proxy ilegível do MAG,
    apesar do nome do arquivo. Skipped if eggnog annotations are empty.
    """
    input:
        eggnog_done = rules.mag_eggnog_prok.output.done,
        annot_tsv   = rules.mag_eggnog_prok.output.annot_tsv,
        bakta       = rules.mag_bakta.output.summary,
    output:
        done       = f"{MAG_CATALOG_DIR}/kegg/done.txt",
        ko_table   = f"{MAG_CATALOG_DIR}/kegg/ko_per_mag.tsv",
        cazy_table = f"{MAG_CATALOG_DIR}/kegg/cazy_per_mag.tsv",
    log:
        f"{OUTDIR}/logs/mag_extract_kegg_kos.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_extract_kegg_kos.tsv"
    conda: "../envs/env_annotation.yaml"
    threads: 1
    params:
        # derivado do output, nao de {{sample}}: regra herdada por
        # coassembly.smk via `use rule ... as ... with:` (wildcard {group}).
        outdir = lambda wc, output: os.path.dirname(output.done),
    run:
        import csv as _csv
        import re, os, sys
        from pathlib import Path

        sys.path.insert(0, SCRIPTS_DIR)
        from mag_catalog import resolve_prefixed_id
        from annotation_tables import CAZY_COL, KEGG_KO_COL, parse_cazy_field

        os.makedirs(params.outdir, exist_ok=True)
        log_path = str(log[0])

        if not os.path.exists(str(input.annot_tsv)) or \
           os.path.getsize(str(input.annot_tsv)) == 0:
            with open(log_path, "w") as lf:
                lf.write("[extract_kegg_kos] EggNOG annotations empty — skipping\n")
            Path(str(output.ko_table)).touch()
            Path(str(output.cazy_table)).touch()
            Path(str(output.done)).touch()
            return

        # Genomas conhecidos = as representantes que o Bakta anotou. O ID do
        # catalogo ja contem "__" ({source}__{bin}), entao cortar a proteina
        # no primeiro separador devolveria a AMOSTRA em vez do MAG.
        genomes = set()
        with open(str(input.bakta), newline="") as bf:
            for row in _csv.DictReader(bf, delimiter="\t"):
                name = (row.get("bin") or "").strip()
                if name:
                    genomes.add(name)

        records, cazy_records, unresolved = [], [], 0
        with open(str(input.annot_tsv)) as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                cols = line.rstrip("\n").split("\t")
                if len(cols) < KEGG_KO_COL:
                    continue
                query   = cols[0]
                kegg_ko = cols[KEGG_KO_COL - 1]
                # CAZy é a última coluna que interessa aqui e nem toda linha
                # chega até ela; sem o guarda, uma linha truncada mataria a
                # extração de KO junto.
                cazy    = cols[CAZY_COL - 1] if len(cols) >= CAZY_COL else ""

                has_ko   = bool(kegg_ko.strip()) and kegg_ko != "-"
                families = parse_cazy_field(cazy)
                if not has_ko and not families:
                    continue

                mag, _rest = resolve_prefixed_id(query, genomes)
                if mag is None:
                    unresolved += 1
                    continue

                if has_ko:
                    for ko_entry in kegg_ko.split(","):
                        ko = ko_entry.strip().replace("ko:", "")
                        if re.match(r"K\d{5}", ko):
                            records.append((mag, ko))
                for fam in families:
                    cazy_records.append((mag, fam))

        with open(str(output.ko_table), "w") as f:
            f.write("mag\tko\n")
            for mag, ko in records:
                f.write(f"{mag}\t{ko}\n")

        with open(str(output.cazy_table), "w") as f:
            f.write("mag\tcazy_family\n")
            for mag, fam in cazy_records:
                f.write(f"{mag}\t{fam}\n")

        n_mags = len({r[0] for r in records})
        n_cazy_mags = len({r[0] for r in cazy_records})
        with open(log_path, "w") as lf:
            lf.write(f"[extract_kegg_kos] {len(records)} KO entries, {n_mags} MAGs\n")
            lf.write(f"[extract_kegg_kos] {len(cazy_records)} CAZy entries, "
                     f"{n_cazy_mags} MAGs\n")
            if unresolved:
                lf.write(f"[extract_kegg_kos] AVISO: {unresolved} proteinas sem "
                         "prefixo de genoma conhecido -- descartadas. Se for a "
                         "maioria, o mag_eggnog_prok rodou sem prefixar.\n")
            lf.write(f"[extract_kegg_kos] Output: {output.ko_table}\n")

        Path(str(output.done)).touch()


rule mag_kegg_completeness:
    """
    Completude de módulo KEGG por MAG (kegg-pathways-completeness, EBI/MGnify).

    O `ko_per_mag.tsv` diz desde sempre que é "pronto para KEGG-Decoder", mas
    ninguém rodava o passo — a pipeline parava em "este MAG tem estes KOs" sem
    chegar em "este MAG fecha esta via". Esta regra fecha a lacuna.

    Escolhida em vez do KEGG-Decoder por três razões medidas em 2026-08-19:
      1. o KEGG-Decoder NÃO existe no bioconda (`conda search` não acha);
      2. esta traz `modules_table.tsv` + `graphs.pkl` embutidos (573 módulos,
         jan/2026) — nenhum banco a baixar, nenhuma licença KEGG;
      3. o genoma é uma COLUNA do input (`nome<TAB>KO<TAB>KO...`), não um
         prefixo cortado no primeiro underscore. Testado com a v1.4.4 real:
         `S1__binette_bin1` volta intacto na coluna `contig`. O KEGG-Decoder
         corta no underscore e teria atribuído os módulos à AMOSTRA (`S1`) —
         a mesma família de bug que docs/ROADMAP_SIMPLIFICACAO.md persegue.

    Duas saídas: uma por MAG (que a vista distribui) e uma agregada do
    catálogo inteiro. Ambas trazem `missing_ko`, não só o percentual — para
    interpretar uma via incompleta é preciso saber o que falta nela.
    """
    input:
        kegg_done = rules.mag_extract_kegg_kos.output.done,
        ko_table  = rules.mag_extract_kegg_kos.output.ko_table,
    output:
        done     = f"{MAG_CATALOG_DIR}/kegg_modules/done.txt",
        per_mag  = f"{MAG_CATALOG_DIR}/kegg_modules/module_completeness.tsv",
        catalog  = f"{MAG_CATALOG_DIR}/kegg_modules/module_completeness_catalog.tsv",
    log:
        f"{OUTDIR}/logs/mag_kegg_completeness.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_kegg_completeness.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("kegg_pathways_completeness")
    threads: 1
    params:
        outdir = lambda wc, output: os.path.dirname(output.done),
        wide   = lambda wc, output: os.path.join(os.path.dirname(output.done), "ko_wide.tsv"),
        script = os.path.join(SCRIPTS_DIR, "ko_to_wide.py"),
    shell:
        """
        mkdir -p {params.outdir}

        if [ ! -s {input.ko_table} ]; then
            echo "[kegg_completeness] ko_per_mag.tsv vazio -- skipping" | tee {log}
            : > {output.per_mag}; : > {output.catalog}
            printf "skipped: no KOs\n" > {output.done}; exit 0
        fi

        python {params.script} {input.ko_table} {params.wide} >> {log} 2>&1

        if [ ! -s {params.wide} ]; then
            echo "[kegg_completeness] nenhum MAG com KO -- skipping" | tee -a {log}
            : > {output.per_mag}; : > {output.catalog}
            printf "skipped: no MAG with KOs\n" > {output.done}; exit 0
        fi

        give_completeness -i {params.wide} -o {params.outdir} -m >> {log} 2>&1

        # A coluna do genoma sai como `contig` (a ferramenta é agnóstica ao
        # que a linha representa); aqui ela é um MAG, e a vista casa por
        # `mag` como em todas as outras tabelas do catálogo.
        sed '1s/^contig\t/mag\t/' {params.outdir}/summary.kegg_contigs.tsv > {output.per_mag}
        cp {params.outdir}/summary.kegg_pathways.tsv {output.catalog}

        N_MAG=$(tail -n +2 {output.per_mag} | cut -f1 | sort -u | wc -l)
        N_MOD=$(tail -n +2 {output.per_mag} | cut -f2 | sort -u | wc -l)
        echo "[kegg_completeness] $N_MAG MAGs, $N_MOD modulos" | tee -a {log}
        printf "ok\n" > {output.done}
        """
