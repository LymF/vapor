# ══════════════════════════════════════════════════════════════════════
# rules/defense_amr.smk — BLOCK 10.5: Defense Systems + AMR (MAGs)
#
# TODAS as regras deste arquivo são GLOBAIS desde 2026-08-19: rodam uma vez
# sobre os MAGs REPRESENTANTES do catálogo (`mag_catalog_proteins`), não uma
# vez por amostra sobre todos os bins. É o item "(h)" do
# docs/ROADMAP_SIMPLIFICACAO.md aplicado ao lado procariótico: um MAG da
# mesma espécie recuperado em três amostras é a mesma entidade biológica,
# então anotá-lo três vezes é desperdício — e pior, nada garantia que as três
# anotações coincidissem. As vistas por amostra/grupo (no fim deste arquivo)
# distribuem o resultado do representante para cada bin membro, escrevendo
# nos MESMOS caminhos de antes, para que finalize/relatório não mudem.
#
# `_read_manifest`/`_concat_proteins` continuam em rules/prok_binning.smk.
#
# defensefinder — MacSyFinder anti-phage defense systems + built-in
#                 AntiDefenseFinder (-a flag); per-genome (architecture-aware)
#                 PADLOC was evaluated as a complementary 2nd detector but
#                 dropped: its biocontainers image ships a BusyBox `rm`
#                 lacking `-d`, which the bundled --db-update script needs,
#                 so the database can never be (re)built inside the container
#                 (confirmed broken on litrp4 2026-06-19) -- conda works fine,
#                 but every other tool in this pipeline must work via Docker.
# amrfinderplus — curated AMR genes/point mutations (NCBI); gene-level,
#                 single batch call on the concatenated protein set
# rgi_card      — curated AMR genes (CARD); gene-level, single batch call
# deeparg       — deep-learning AMR genes; exploratory/sensitivity-oriented
#                 complement to amrfinderplus+rgi's curated calls
#                 (Serrana et al. 2026) -- reported separately, never merged
# abricate      — BLASTN screen for VFDB (virulence factors) + PlasmidFinder
#                 (plasmid replicons) only -- NOT AMR; its bundled AMR
#                 databases are flat BLASTN screens with no point-mutation
#                 detection or SNP/variant models, a downgrade vs the 3
#                 tools above if used for that purpose
# argnorm_normalize — maps AMRFinderPlus/DeepARG gene calls onto the
#                 Antibiotic Resistance Ontology (ARO), so they can be
#                 compared with each other and with RGI (which already
#                 speaks ARO natively via CARD -- not routed through argNorm;
#                 a hAMRonization bridge was tried and dropped, argNorm
#                 1.1.0 has no working RGI support despite the docs implying
#                 otherwise, see argnorm_normalize docstring). Normalized
#                 tables stay separate per tool -- same "never merged" rule
#                 as the raw outputs.
#
# All tools soft-fail (header-only TSV + done.txt) when disabled
# (defense_amr_enabled / abricate_enabled / argnorm_enabled: false) or no
# genome units / upstream hits exist.
# Config keys: defense_amr_enabled, card_db,
# abricate_enabled, argnorm_enabled
# ══════════════════════════════════════════════════════════════════════


rule mag_defensefinder:
    """
    DefenseFinder -- systematic detection of anti-phage defense systems
    (MacSyFinder + HMM models), with built-in AntiDefenseFinder
    (--antidefensefinder flag) for anti-defense proteins in the same pass.
    Runs once per genome unit (manifesto do mag_catalog_proteins, ou seja,
    uma vez por MAG REPRESENTANTE do catálogo global): MacSyFinder
    needs gene order within a replicon, so genomes cannot be concatenated
    (unlike the gene-level AMR tools below).
    """
    input:
        manifest = rules.mag_catalog_proteins.output.manifest,
        done     = rules.mag_catalog_proteins.output.done,
    output:
        done        = f"{MAG_CATALOG_DIR}/defensefinder/done.txt",
        systems     = f"{MAG_CATALOG_DIR}/defensefinder/defensefinder_systems.tsv",
        antisystems = f"{MAG_CATALOG_DIR}/defensefinder/antidefensefinder_systems.tsv",
        # Sumario por genoma (bin TAB status, mesma convencao de
        # bakta_summary.tsv em rules/annotation.smk): o laco abaixo e
        # per-genoma e so loga um WARNING quando falha, sem derrubar a
        # regra -- entao um genoma ausente de defensefinder_systems.tsv e
        # ambiguo por construcao (rodou-e-nao-achou-nada vs quebrou). Este
        # arquivo e a unica fonte que distingue os dois casos; ver
        # scripts/pangenome_matrix.py.
        summary     = f"{MAG_CATALOG_DIR}/defensefinder/defensefinder_summary.tsv",
    log:
        f"{OUTDIR}/logs/mag_defensefinder.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_defensefinder.tsv"
    conda: "../envs/env_defense.yaml"
    container:  CONTAINERS.get("defense_finder")
    threads: THREADS
    params:
        # derivado do output: a regra e global (uma execucao para todo o
        # catalogo), entao nao ha wildcard nenhum a citar.
        outdir     = lambda wc, output: os.path.dirname(output.done),
        models_dir = DEFENSE_FINDER_MODELS_DB,
        enabled    = DEFENSE_AMR_ENABLED,
    run:
        import csv, glob, os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)

        def write_empty(msg, status):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.systems)).write_text("genome\n")
            Path(str(output.antisystems)).write_text("genome\n")
            Path(str(output.summary)).write_text("genome\tstatus\n")
            write_status(str(output.done), status)

        if (not params.enabled or not os.path.exists(str(input.manifest))
                or os.path.getsize(str(input.manifest)) == 0):
            write_empty("[defensefinder] Disabled or no genome units -- skipping",
                        "skipped: disabled or no genome units")
            return

        # Without --models-dir, defense-finder caches models under
        # $HOME/.macsyfinder/models -- but $HOME isn't stable across how
        # this rule gets invoked (conda vs apptainer vs whatever sets the
        # working directory), so different runs ended up with different,
        # sometimes only-partially-populated caches. Confirmed on litrp4
        # 2026-06-20: a batch run left every sample with 0 systems despite
        # macsyfinder completing its HMMER passes (DefenseFinder/RM tmp
        # dirs were fully populated) -- the final consolidated
        # *_defense_finder_systems.tsv was never written, while the exact
        # same command against an explicit models dir, run standalone,
        # worked and produced it immediately. Pinning --models-dir removes
        # the ambiguity: one shared, always-consistent location, same as
        # card_db/deeparg_db below.
        models_dir = params.models_dir or os.path.join(params.outdir, "models")
        os.makedirs(models_dir, exist_ok=True)

        # defense-finder update always pings GitHub for the latest model
        # release, even when models are already cached -- across a batch of
        # many samples that means one GitHub API call per sample, which
        # exhausts the unauthenticated 60-req/hour limit fast (confirmed on
        # litrp4 2026-06-20: every sample after the first few failed with
        # "maximum number of request per hour"). Only call update when this
        # models dir is empty -- once any one sample populates it, every
        # later sample (in this or any future batch) skips the network call.
        if os.listdir(models_dir):
            with open(str(log[0]), "a") as lf:
                lf.write(f"[defensefinder] Models already cached in {models_dir} -- skipping update\n")
        else:
            shell("defense-finder update --models-dir {models_dir} >> {log} 2>&1 || "
                  "echo '[defensefinder] WARNING: model update failed -- "
                  "GitHub rate limit? set GITHUB_TOKEN or retry later' >> {log}")
        # NOTE: env_defense.yaml/containers.yaml pin defense-finder=3.0.0, not
        # 2.0.0/2.0.1 -- those fail against current CasFinder releases with
        # "has not the right version" (macsypy.error.MacsypyError), see
        # https://github.com/mdmparis/defense-finder/issues/95. 3.0.0 resolves
        # a compatible CasFinder (3.1.0) automatically. If this still errors
        # on a fresh CasFinder release, the per-genome loop below already
        # degrades gracefully (warns + 0 rows, doesn't fail the rule).

        genome_status = {}
        for name, mode, fna, faa, gff in _read_manifest(str(input.manifest)):
            if not os.path.exists(faa) or os.path.getsize(faa) == 0:
                genome_status[name] = "failed"
                continue
            genome_out = os.path.join(params.outdir, name)
            os.makedirs(genome_out, exist_ok=True)
            shell(
                "defense-finder run -o {genome_out} --models-dir {models_dir} "
                "--antidefensefinder {faa} "
                ">> {log} 2>&1 || echo '[defensefinder] WARNING: failed on {name}' >> {log}"
            )
            # O `|| echo WARNING` acima nao derruba a regra de proposito
            # (um genoma ruim nao pode travar o catalogo inteiro) -- por
            # isso o status real e lido do arquivo que o defense-finder
            # deveria ter escrito, nao do exit code do shell().
            out_tsv = os.path.join(genome_out, f"{name}_defense_finder_systems.tsv")
            genome_status[name] = "ok" if os.path.exists(out_tsv) else "failed"

        # defense-finder 3.0.0 (--antidefensefinder) writes ONE file per
        # genome -- {name}_defense_finder_systems.tsv -- with defense AND
        # anti-defense systems both in it, distinguished by the "activity"
        # column ("Defense" vs "Anti-defense"), not by a separate filename
        # (confirmed on litrp4 2026-06-20: no *anti*systems.tsv file is ever
        # produced). Split rows by that column instead of by filename.
        def_rows, anti_rows, header = [], [], None
        for tsv in sorted(glob.glob(os.path.join(params.outdir, "*", "*_defense_finder_systems.tsv"))):
            genome = os.path.basename(os.path.dirname(tsv))
            with open(tsv) as f:
                r = csv.reader(f, delimiter="\t")
                h = next(r, None)
                if h is None:
                    continue
                if header is None:
                    header = h
                activity_idx = h.index("activity") if "activity" in h else None
                for row in r:
                    is_anti = (activity_idx is not None and len(row) > activity_idx
                               and "anti" in row[activity_idx].lower())
                    (anti_rows if is_anti else def_rows).append([genome] + row)

        def write(path, rows):
            with open(path, "w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["genome"] + (header or []))
                w.writerows(rows)

        write(str(output.systems), def_rows)
        write(str(output.antisystems), anti_rows)

        with open(str(output.summary), "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["genome", "status"])
            for name, mode, fna, faa, gff in _read_manifest(str(input.manifest)):
                w.writerow([name, genome_status.get(name, "failed")])

        with open(str(log[0]), "a") as lf:
            lf.write(f"[defensefinder] Done -- {len(def_rows)} defense, "
                      f"{len(anti_rows)} antidefense system rows\n")
        write_status(str(output.done), "ok")

# rule defensefinder_viral / rule dbapis_viral removed 2026-08-18 (second
# half of "(h)", docs/ROADMAP_SIMPLIFICACAO.md): both moved to the global
# vOTU catalog as rule votu_defensefinder_viral / rule votu_dbapis_viral
# (rules/votu_catalog.smk). Both already read the GLOBAL rules.votu_prodigal
# output, so the per-sample/per-group fan-out here was pure waste (N+G
# runs on byte-identical input), and worse: each per-group run silently
# processed the WHOLE catalog, so coassembly finalize copied catalog-wide
# systems into a group-scoped path. See rules/votu_catalog.smk for the
# moved rules.

rule mag_amrfinderplus:
    """
    AMRFinderPlus -- curated, alignment-based AMR gene + point-mutation
    detection (NCBI Reference Gene/Point Mutation databases).
    Gene-level: runs once on the concatenated protein set (all genome
    units from mag_catalog_proteins), unlike DefenseFinder.
    Reuses env_annotation.yaml, which already pins ncbi-amrfinderplus.
    """
    input:
        manifest = rules.mag_catalog_proteins.output.manifest,
        done     = rules.mag_catalog_proteins.output.done,
    output:
        done    = f"{MAG_CATALOG_DIR}/amrfinderplus/done.txt",
        results = f"{MAG_CATALOG_DIR}/amrfinderplus/amrfinder_results.tsv",
    log:
        f"{OUTDIR}/logs/mag_amrfinderplus.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_amrfinderplus.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("ncbi-amrfinderplus")
    threads: THREADS
    params:
        # derivado do output: a regra e global (uma execucao para todo o
        # catalogo), entao nao ha wildcard nenhum a citar.
        outdir  = lambda wc, output: os.path.dirname(output.done),
        enabled = DEFENSE_AMR_ENABLED,
    run:
        import os
        import subprocess as sp
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        all_faa = os.path.join(params.outdir, "all_genomes.faa")

        def write_empty(msg, status):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.results)).write_text("Protein identifier\tGene symbol\n")
            write_status(str(output.done), status)

        if (not params.enabled or not os.path.exists(str(input.manifest))
                or os.path.getsize(str(input.manifest)) == 0):
            write_empty("[amrfinderplus] Disabled or no genome units -- skipping",
                        "skipped: disabled or no genome units")
            return

        if not _concat_proteins(str(input.manifest), all_faa):
            write_empty("[amrfinderplus] No proteins -- skipping", "skipped: no proteins")
            return

        shell("amrfinder -u >> {log} 2>&1 || "
              "echo '[amrfinderplus] WARNING: database update failed (may already be current)' >> {log}")
        amr_err = None
        try:
            shell(
                "amrfinder -p {all_faa} --plus -o {output.results} --threads {threads} "
                ">> {log} 2>&1"
            )
        except sp.CalledProcessError as exc:
            amr_err = exc.returncode

        if not os.path.exists(str(output.results)) or os.path.getsize(str(output.results)) == 0:
            Path(str(output.results)).write_text("Protein identifier\tGene symbol\n")

        if amr_err is not None:
            with open(str(log[0]), "a") as lf:
                lf.write(f"[amrfinderplus] WARNING: amrfinder failed (exit {amr_err})\n")
            write_status(str(output.done), f"failed: amrfinder exit {amr_err}")
        else:
            write_status(str(output.done), "ok")


rule mag_rgi_card:
    """
    RGI (CARD) -- curated AMR detection via the Comprehensive Antibiotic
    Resistance Database (homology + SNP models). Same concatenated batch
    input as AMRFinderPlus; reported alongside it as "curated" AMR.
    CARD database path is configurable (card_db) -- auto-downloaded and
    loaded into that directory on first use (public CARD download URL),
    so a fresh clone is reproducible without a manual prep step.
    """
    input:
        manifest = rules.mag_catalog_proteins.output.manifest,
        done     = rules.mag_catalog_proteins.output.done,
    output:
        done    = f"{MAG_CATALOG_DIR}/rgi/done.txt",
        results = f"{MAG_CATALOG_DIR}/rgi/rgi_results.txt",
    log:
        f"{OUTDIR}/logs/mag_rgi.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_rgi.tsv"
    conda: "../envs/env_rgi.yaml"
    container:  CONTAINERS.get("rgi")
    threads: THREADS
    params:
        # derivado do output: a regra e global (uma execucao para todo o
        # catalogo), entao nao ha wildcard nenhum a citar.
        outdir  = lambda wc, output: os.path.dirname(output.done),
        card_db = CARD_DB,
        enabled = DEFENSE_AMR_ENABLED,
    run:
        import os
        import subprocess as sp
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        all_faa = os.path.join(params.outdir, "all_genomes.faa")

        def write_empty(msg, status):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.results)).write_text("ORF_ID\tBest_Hit_ARO\n")
            write_status(str(output.done), status)

        if (not params.enabled or not os.path.exists(str(input.manifest))
                or os.path.getsize(str(input.manifest)) == 0):
            write_empty("[rgi] Disabled or no genome units -- skipping",
                        "skipped: disabled or no genome units")
            return

        if not _concat_proteins(str(input.manifest), all_faa):
            write_empty("[rgi] No proteins -- skipping", "skipped: no proteins")
            return

        card_dir = params.card_db or os.path.join(params.outdir, "card_db")
        os.makedirs(card_dir, exist_ok=True)
        card_json = os.path.join(card_dir, "card.json")
        if not os.path.exists(card_json) or os.path.getsize(card_json) == 0:
            shell(
                "curl -sL https://card.mcmaster.ca/latest/data -o {card_dir}/card_data.tar.bz2 "
                ">> {log} 2>&1 && tar -xjf {card_dir}/card_data.tar.bz2 -C {card_dir} "
                # archive members are stored as "./card.json" etc -- an exact
                # "card.json" filter doesn't match, so extract everything.
                ">> {log} 2>&1 || true"
            )

        if not os.path.exists(card_json) or os.path.getsize(card_json) == 0:
            write_empty("[rgi] WARNING: CARD database unavailable -- skipping",
                        "failed: CARD database unavailable")
            return

        # RGI rejects the whole input with "invalid protein fasta due to: '*'"
        # when Prodigal's stop codons are present, so it silently ran on a
        # subset before this. Strip them into a separate file -- all_genomes.faa
        # is shared with AMRFinderPlus and DeepARG and must stay untouched.
        rgi_faa = os.path.join(params.outdir, "all_genomes_nostop.faa")
        with open(all_faa) as fin, open(rgi_faa, "w") as fout:
            for line in fin:
                fout.write(line if line.startswith(">") else line.replace("*", ""))

        # --local makes RGI create/read ./localDB relative to the CURRENT
        # working directory, not the --card_json path -- 'load' and 'main'
        # must run from the exact same directory or 'main' can't find the
        # database 'load' just built.
        shell("cd {params.outdir} && rgi load --card_json {card_json} --local >> {log} 2>&1")
        rgi_err = None
        try:
            shell(
                "cd {params.outdir} && rgi main -i {rgi_faa} -t protein --include_loose "
                "--local -o rgi_results -n {threads} >> {log} 2>&1"
            )
        except sp.CalledProcessError as exc:
            rgi_err = exc.returncode

        # 'rgi main -o rgi_results' inside params.outdir already writes
        # directly to output.results (rgi_results.txt) -- same path, no
        # copy needed (shutil.copy onto an identical path raises
        # SameFileError).
        produced = os.path.join(params.outdir, "rgi_results.txt")
        if not os.path.exists(produced) or os.path.getsize(produced) == 0:
            Path(str(output.results)).write_text("ORF_ID\tBest_Hit_ARO\n")

        if rgi_err is not None:
            with open(str(log[0]), "a") as lf:
                lf.write(f"[rgi] WARNING: rgi main failed (exit {rgi_err})\n")
            write_status(str(output.done), f"failed: rgi main exit {rgi_err}")
        else:
            write_status(str(output.done), "ok")


rule mag_deeparg:
    """
    DeepARG -- deep-learning AMR gene prediction (CNN trained on
    CARD/ARDB/UniProt-derived sequences). Exploratory/sensitivity-oriented
    complement to AMRFinderPlus+RGI's curated calls (Serrana et al. 2026,
    iMetaOmics): higher recall on divergent/novel environmental ARGs that
    curated homology search misses, but with no accuracy benchmark --
    kept and reported separately, never merged into the curated AMR count.
    Same concatenated batch input as AMRFinderPlus/RGI.

    NOTE: bioconda's deeparg=1.0.4 is the classic Python2/Theano codebase,
    not the newer PyTorch/HuggingFace rewrite some docs describe -- there is
    no auto-download and no --threads flag. `deeparg predict` requires an
    explicit --data-path (fetched once via `deeparg download_data`).
    """
    input:
        manifest = rules.mag_catalog_proteins.output.manifest,
        done     = rules.mag_catalog_proteins.output.done,
    output:
        done    = f"{MAG_CATALOG_DIR}/deeparg/done.txt",
        results = f"{MAG_CATALOG_DIR}/deeparg/deeparg_results.mapping.ARG",
    log:
        f"{OUTDIR}/logs/mag_deeparg.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_deeparg.tsv"
    conda: "../envs/env_deeparg.yaml"
    container:  CONTAINERS.get("deeparg")
    threads: THREADS
    params:
        # derivado do output: a regra e global (uma execucao para todo o
        # catalogo), entao nao ha wildcard nenhum a citar.
        outdir   = lambda wc, output: os.path.dirname(output.done),
        data_dir = DEEPARG_DB,
        enabled  = DEFENSE_AMR_ENABLED,
    run:
        import os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        all_faa = os.path.join(params.outdir, "all_genomes.faa")

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.results)).write_text(
                "#ARG\tquery-start\tquery-end\tread_id\n")
            Path(str(output.done)).touch()

        if (not params.enabled or not os.path.exists(str(input.manifest))
                or os.path.getsize(str(input.manifest)) == 0):
            write_empty("[deeparg] Disabled or no genome units -- skipping")
            return

        if not _concat_proteins(str(input.manifest), all_faa):
            write_empty("[deeparg] No proteins -- skipping")
            return

        data_dir = params.data_dir or os.path.join(params.outdir, "deeparg_data")
        os.makedirs(data_dir, exist_ok=True)
        if not os.listdir(data_dir):
            shell("deeparg download_data -o {data_dir} >> {log} 2>&1 || "
                  "echo '[deeparg] WARNING: download_data failed' >> {log}")

        if not os.listdir(data_dir):
            write_empty("[deeparg] WARNING: deeparg data unavailable -- skipping")
            return

        out_prefix = os.path.join(params.outdir, "deeparg_results")
        shell(
            "deeparg predict --model LS --type prot -i {all_faa} -o {out_prefix} "
            "-d {data_dir} >> {log} 2>&1 || "
            "echo '[deeparg] WARNING: deeparg predict failed' >> {log}"
        )

        if not os.path.exists(str(output.results)) or os.path.getsize(str(output.results)) == 0:
            Path(str(output.results)).write_text(
                "#ARG\tquery-start\tquery-end\tread_id\n")
        Path(str(output.done)).touch()


def _has_data_rows(path):
    """True if a TSV has more than just a header line -- distinguishes a
    real result set from the header-only stub files the soft-fail paths
    above write when a tool was disabled or had nothing to call."""
    import os
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return False
    with open(path) as f:
        next(f, None)
        return next(f, None) is not None


rule mag_abricate:
    """
    ABRicate -- BLASTN mass screening of contigs/bins for gene presence.
    Used here only for the two databases AMRFinderPlus/RGI/DeepARG do not
    cover -- virulence factors (VFDB) and plasmid replicons (PlasmidFinder)
    -- NOT for AMR gene calling: ABRicate's bundled AMR databases are
    static flat-file BLASTN screens with no point-mutation detection
    (AMRFinderPlus) or SNP/variant models (RGI/CARD), so they would be a
    strict downgrade if used to replace either tool.
    Runs per genome unit (manifest from mag_catalog_proteins) on the
    nucleotide sequence -- ABRicate works on DNA via BLASTN, unlike the
    protein-level AMR/defense tools above.
    """
    input:
        manifest = rules.mag_catalog_proteins.output.manifest,
        done     = rules.mag_catalog_proteins.output.done,
    output:
        done          = f"{MAG_CATALOG_DIR}/abricate/done.txt",
        vfdb          = f"{MAG_CATALOG_DIR}/abricate/vfdb_results.tsv",
        plasmidfinder = f"{MAG_CATALOG_DIR}/abricate/plasmidfinder_results.tsv",
    log:
        f"{OUTDIR}/logs/mag_abricate.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_abricate.tsv"
    conda: "../envs/env_abricate.yaml"
    container:  CONTAINERS.get("abricate")
    threads: THREADS
    params:
        # derivado do output: a regra e global (uma execucao para todo o
        # catalogo), entao nao ha wildcard nenhum a citar.
        outdir  = lambda wc, output: os.path.dirname(output.done),
        dbs     = ["vfdb", "plasmidfinder"],
        enabled = ABRICATE_ENABLED,
    run:
        import csv, os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        out_paths = {"vfdb": str(output.vfdb), "plasmidfinder": str(output.plasmidfinder)}

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            for db in params.dbs:
                Path(out_paths[db]).write_text("genome\n")
            Path(str(output.done)).touch()

        if (not params.enabled or not os.path.exists(str(input.manifest))
                or os.path.getsize(str(input.manifest)) == 0):
            write_empty("[abricate] Disabled or no genome units -- skipping")
            return

        shell("abricate --setupdb >> {log} 2>&1 || "
              "echo '[abricate] WARNING: setupdb failed (may already be set up)' >> {log}")

        for db in params.dbs:
            rows, header = [], None
            for name, mode, fna, faa, gff in _read_manifest(str(input.manifest)):
                if not os.path.exists(fna) or os.path.getsize(fna) == 0:
                    continue
                genome_tsv = os.path.join(params.outdir, f"{name}.{db}.tsv")
                shell(
                    "abricate --db {db} --quiet {fna} > {genome_tsv} "
                    "2>> {log} || echo '[abricate] WARNING: failed on {name} ({db})' >> {log}"
                )
                if os.path.exists(genome_tsv) and os.path.getsize(genome_tsv) > 0:
                    with open(genome_tsv) as f:
                        r = csv.reader(f, delimiter="\t")
                        h = next(r, None)
                        if h is None:
                            continue
                        if header is None:
                            header = h
                        for row in r:
                            rows.append([name] + row)

            with open(out_paths[db], "w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["genome"] + (header or []))
                w.writerows(rows)

        with open(str(log[0]), "a") as lf:
            lf.write(f"[abricate] Done -- {len(params.dbs)} database(s) screened\n")
        Path(str(output.done)).touch()


rule mag_argnorm_normalize:
    """
    argNorm -- maps AMRFinderPlus/DeepARG gene calls onto the Antibiotic
    Resistance Ontology (ARO), so the same gene reported under different
    names by the two tools can be compared with each other and with RGI.
    RGI is deliberately NOT routed through argNorm: it is built directly
    on CARD, so its native output (Best_Hit_ARO) already speaks ARO --
    there is nothing to normalize. A hAMRonization+argnorm bridge for RGI
    was tried and dropped (confirmed in argNorm 1.1.0 source,
    normalizers.py's input_id_lookup has no 'rgi' key at all -- the
    documented hAMRonization pathway for RGI does not actually work,
    confirmed on litrp4 2026-06-20 by a hamronized RGI file with a
    correct analysis_software_name="rgi" still getting rejected as
    "abricate is not a supported ARG annotation tool").
    The 2 normalized tables stay separate -- same "never merged" rule as
    the raw AMR outputs above, this only adds a common ARO/drug-class
    vocabulary to each.
    """
    input:
        amrfinder      = rules.mag_amrfinderplus.output.results,
        amrfinder_done = rules.mag_amrfinderplus.output.done,
        deeparg        = rules.mag_deeparg.output.results,
        deeparg_done   = rules.mag_deeparg.output.done,
    output:
        done             = f"{MAG_CATALOG_DIR}/argnorm/done.txt",
        amrfinder_normed = f"{MAG_CATALOG_DIR}/argnorm/amrfinderplus_normed.tsv",
        deeparg_normed   = f"{MAG_CATALOG_DIR}/argnorm/deeparg_normed.tsv",
    log:
        f"{OUTDIR}/logs/mag_argnorm.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_argnorm.tsv"
    conda: "../envs/env_argnorm.yaml"
    container:  CONTAINERS.get("argnorm")
    threads: 1
    params:
        enabled = ARGNORM_ENABLED,
    run:
        import os
        from pathlib import Path

        os.makedirs(os.path.dirname(str(output.done)), exist_ok=True)

        def stub(path):
            Path(path).write_text("ARO\n")

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            stub(output.amrfinder_normed)
            stub(output.deeparg_normed)
            Path(str(output.done)).touch()

        if not params.enabled:
            write_empty("[argnorm] argnorm_enabled=False -- skipping")
            return

        if _has_data_rows(str(input.amrfinder)):
            shell(
                "argnorm amrfinderplus -i {input.amrfinder} -o {output.amrfinder_normed} "
                ">> {log} 2>&1 || echo '[argnorm] WARNING: amrfinderplus normalization failed' >> {log}"
            )
        if not os.path.exists(str(output.amrfinder_normed)) or os.path.getsize(str(output.amrfinder_normed)) == 0:
            stub(output.amrfinder_normed)

        if _has_data_rows(str(input.deeparg)):
            shell(
                "argnorm deeparg -i {input.deeparg} -o {output.deeparg_normed} "
                ">> {log} 2>&1 || echo '[argnorm] WARNING: deeparg normalization failed' >> {log}"
            )
        if not os.path.exists(str(output.deeparg_normed)) or os.path.getsize(str(output.deeparg_normed)) == 0:
            stub(output.deeparg_normed)

        with open(str(log[0]), "a") as lf:
            lf.write("[argnorm] Done\n")
        Path(str(output.done)).touch()


rule mag_amr_consensus:
    """
    AMR consensus -- merge AMRFinderPlus + RGI/CARD + DeepARG hits by CDS
    locus using ARO as the common vocabulary.  Consensus score = number of
    tools that detected the locus divided by 3.  The curated tools
    (RGI > AMRFinderPlus) supply the canonical annotation fields; DeepARG
    contributes to the score but is not used as the preferred annotation
    source.  argNorm-normalized files are used for AMRFinderPlus/DeepARG so
    all three tools share a common ARO namespace.
    """
    input:
        argnorm_done     = rules.mag_argnorm_normalize.output.done,
        rgi_done         = rules.mag_rgi_card.output.done,
        amrfinder_normed = rules.mag_argnorm_normalize.output.amrfinder_normed,
        deeparg_normed   = rules.mag_argnorm_normalize.output.deeparg_normed,
        rgi_results      = rules.mag_rgi_card.output.results,
    output:
        done      = f"{MAG_CATALOG_DIR}/amr_consensus/done.txt",
        consensus = f"{MAG_CATALOG_DIR}/amr_consensus/amr_consensus.tsv",
    log:
        f"{OUTDIR}/logs/mag_amr_consensus.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_amr_consensus.tsv"
    conda: "../envs/env_argnorm.yaml"
    threads: 1
    params:
        enabled = AMR_CONSENSUS_ENABLED,
        script  = os.path.join(workflow.basedir, "scripts", "consolidate_amr.py"),
    run:
        import os
        from pathlib import Path

        os.makedirs(os.path.dirname(str(output.done)), exist_ok=True)

        def write_empty(msg, status):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            cols = "\t".join([
                "locus", "aro_accession", "gene_name", "drug_class",
                "resistance_mechanism", "n_tools", "consensus_score", "tools_detected",
            ])
            Path(str(output.consensus)).write_text(cols + "\n")
            write_status(str(output.done), status)

        if not params.enabled:
            write_empty("[amr_consensus] use_amr_consensus=False -- skipping",
                        "skipped: use_amr_consensus=False")
            return

        shell(
            "python {params.script} "
            "--amrfinder-normed {input.amrfinder_normed} "
            "--rgi-results {input.rgi_results} "
            "--deeparg-normed {input.deeparg_normed} "
            "-o {output.consensus} >> {log} 2>&1"
        )

        if not os.path.exists(str(output.consensus)) or os.path.getsize(str(output.consensus)) == 0:
            write_empty("[amr_consensus] WARNING: script produced no output",
                        "failed: script produced no output")
            return

        with open(str(log[0]), "a") as lf:
            lf.write("[amr_consensus] Done\n")
        write_status(str(output.done), "ok")


# ══════════════════════════════════════════════════════════════════════
# VISTAS por fonte
#
# Ficam AQUI, e não em rules/mag_catalog.smk onde os helpers vivem, apenas
# por ordem de include: elas referenciam `rules.mag_*`, que só existem
# depois deste arquivo e de rules/taxonomy.smk.
#
# São a ÚNICA coisa que escreve em `{sample}/bins/...` e
# `coassembly/{group}/bins/...` desde 2026-08-19 — as regras por amostra
# que escreviam ali foram APAGADAS, não gateadas. É por isso que nenhum
# `ruleorder` é necessário: não há dois produtores do mesmo caminho.
# ══════════════════════════════════════════════════════════════════════


rule mag_views_sample:
    """Distribui para os bins da amostra o que foi computado nas representantes."""
    input:
        membership = rules.mag_catalog_membership.output.tsv,
        **{k: v for k, v in _MAG_VIEW_GLOBAL().items()}
    output:
        **{k: v for k, v in _mag_view_io(f"{OUTDIR}/{{sample}}").items()},
        mmseqs      = f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/taxonomy.tsv",
        mmseqs_done = f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/mag_views.log"
    run:
        _mag_write_views(wildcards.sample, str(input.membership),
                         dict(input.items()), dict(output.items()), str(log[0]))


if COASSEMBLY_ENABLED and COASSEMBLY_BINNING and not LONG_READS:

    rule mag_views_group:
        """Gêmea de `mag_views_sample` para os MAGs do co-binning VAMB.

        Sem a vista de MMseqs2: a trilha de grupo nunca teve
        `bins/mmseqs_taxonomy_prok` (não havia gêmea `coassembly_` dela), e
        a vista não é lugar de inventar saída nova.
        """
        input:
            membership = rules.mag_catalog_membership.output.tsv,
            **{k: v for k, v in _MAG_VIEW_GLOBAL().items() if not k.startswith("mmseqs")}
        output:
            **_mag_view_io(f"{OUTDIR}/coassembly/{{group}}")
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/mag_views.log"
        run:
            _mag_write_views(wildcards.group, str(input.membership),
                             dict(input.items()), dict(output.items()), str(log[0]))
