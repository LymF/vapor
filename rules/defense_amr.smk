# ══════════════════════════════════════════════════════════════════════
# rules/defense_amr.smk — BLOCK 10.5: Defense Systems + AMR (prok bins)
#
# prok_bin_proteins (shared per-genome Prodigal proteins, feeds all 4 tools
# below + mmseqs_taxonomy_prok in taxonomy.smk) now lives in
# rules/prok_binning.smk -- it must be defined before BOTH taxonomy.smk and
# this file in the Snakefile's include: order, since both reference
# rules.prok_bin_proteins. _read_manifest/_concat_proteins moved with it.
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
# Config keys: defense_amr_enabled, low_depth_mode, card_db,
# abricate_enabled, argnorm_enabled
# ══════════════════════════════════════════════════════════════════════


rule defensefinder:
    """
    DefenseFinder -- systematic detection of anti-phage defense systems
    (MacSyFinder + HMM models), with built-in AntiDefenseFinder
    (--antidefensefinder flag) for anti-defense proteins in the same pass.
    Runs once per genome unit (manifest from prok_bin_proteins): MacSyFinder
    needs gene order within a replicon, so genomes cannot be concatenated
    (unlike the gene-level AMR tools below).
    """
    input:
        manifest = rules.prok_bin_proteins.output.manifest,
        done     = rules.prok_bin_proteins.output.done,
    output:
        done        = f"{OUTDIR}/{{sample}}/bins/defensefinder/done.txt",
        systems     = f"{OUTDIR}/{{sample}}/bins/defensefinder/defensefinder_systems.tsv",
        antisystems = f"{OUTDIR}/{{sample}}/bins/defensefinder/antidefensefinder_systems.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/defensefinder.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/defensefinder.tsv"
    conda: "../envs/env_defense.yaml"
    container:  CONTAINERS.get("defense_finder")
    threads: THREADS
    params:
        outdir     = f"{OUTDIR}/{{sample}}/bins/defensefinder",
        models_dir = DEFENSE_FINDER_MODELS_DB,
        enabled    = DEFENSE_AMR_ENABLED,
    run:
        import csv, glob, os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.systems)).write_text("genome\n")
            Path(str(output.antisystems)).write_text("genome\n")
            Path(str(output.done)).touch()

        if (not params.enabled or not os.path.exists(str(input.manifest))
                or os.path.getsize(str(input.manifest)) == 0):
            write_empty("[defensefinder] Disabled or no genome units -- skipping")
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

        for name, mode, fna, faa, gff in _read_manifest(str(input.manifest)):
            if not os.path.exists(faa) or os.path.getsize(faa) == 0:
                continue
            genome_out = os.path.join(params.outdir, name)
            os.makedirs(genome_out, exist_ok=True)
            shell(
                "defense-finder run -o {genome_out} --models-dir {models_dir} "
                "--antidefensefinder {faa} "
                ">> {log} 2>&1 || echo '[defensefinder] WARNING: failed on {name}' >> {log}"
            )

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

        with open(str(log[0]), "a") as lf:
            lf.write(f"[defensefinder] Done -- {len(def_rows)} defense, "
                      f"{len(anti_rows)} antidefense system rows\n")
        Path(str(output.done)).touch()


rule defensefinder_viral:
    """
    Anti-defense systems on VIRAL proteins (Han et al. 2026, Nat Commun
    cold seep defensome paper, Fig 6b/6c) -- DefenseFinder's
    --antidefensefinder side of the SAME tool/models/container as
    `defensefinder` above, just pointed at viral ORFs instead of bin
    proteins. Today DefenseFinder/AntiDefenseFinder only runs on the host
    (bin) side, so every Host<->Virus defense/anti-defense cross-link is
    one-sided (we know which bins defend, never which phages counter-defend).
    One call across the whole sample's viral protein set (rules.prodigal_viral
    ran in meta mode on the multi-genome viral FASTA already, same way the
    low-depth "contigs_pseudogenome" fallback above is tolerated without
    per-genome splitting).
    """
    input:
        faa  = rules.prodigal_viral.output.faa,
        done = rules.prodigal_viral.output.done,
    output:
        done        = f"{OUTDIR}/{{sample}}/viral/defensefinder/done.txt",
        systems     = f"{OUTDIR}/{{sample}}/viral/defensefinder/viral_defense_systems.tsv",
        antisystems = f"{OUTDIR}/{{sample}}/viral/defensefinder/viral_antidefense_systems.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/defensefinder_viral.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/defensefinder_viral.tsv"
    conda: "../envs/env_defense.yaml"
    container:  CONTAINERS.get("defense_finder")
    threads: THREADS
    params:
        outdir     = f"{OUTDIR}/{{sample}}/viral/defensefinder",
        models_dir = DEFENSE_FINDER_MODELS_DB,
        enabled    = DEFENSE_AMR_VIRAL_ENABLED,
    run:
        import csv, glob, os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.systems)).write_text("genome\n")
            Path(str(output.antisystems)).write_text("genome\n")
            Path(str(output.done)).touch()

        if (not params.enabled or not os.path.exists(str(input.faa))
                or os.path.getsize(str(input.faa)) == 0):
            write_empty("[defensefinder_viral] Disabled or no viral proteins -- skipping")
            return

        models_dir = params.models_dir or os.path.join(params.outdir, "models")
        os.makedirs(models_dir, exist_ok=True)
        if os.listdir(models_dir):
            with open(str(log[0]), "a") as lf:
                lf.write(f"[defensefinder_viral] Models already cached in {models_dir}\n")
        else:
            shell("defense-finder update --models-dir {models_dir} >> {log} 2>&1 || "
                  "echo '[defensefinder_viral] WARNING: model update failed' >> {log}")

        run_out = os.path.join(params.outdir, "run")
        os.makedirs(run_out, exist_ok=True)
        shell(
            "defense-finder run -o {run_out} --models-dir {models_dir} "
            "--antidefensefinder {input.faa} "
            ">> {log} 2>&1 || echo '[defensefinder_viral] WARNING: run failed' >> {log}"
        )

        # Same split-by-"activity" logic as the host-side rule above (3.0.0
        # writes one consolidated *_defense_finder_systems.tsv, Defense and
        # Anti-defense rows distinguished by the "activity" column).
        def_rows, anti_rows, header = [], [], None
        for tsv in sorted(glob.glob(os.path.join(run_out, "*_defense_finder_systems.tsv"))):
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
                    # "genome" column here is the sample (single viral protein
                    # set), not a per-bin name -- the protein_in_syst column
                    # (already part of `header`) carries the actual viral
                    # contig/ORF IDs needed to attribute hits downstream.
                    (anti_rows if is_anti else def_rows).append([wildcards.sample] + row)

        def write(path, rows):
            with open(path, "w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["genome"] + (header or []))
                w.writerows(rows)

        write(str(output.systems), def_rows)
        write(str(output.antisystems), anti_rows)

        with open(str(log[0]), "a") as lf:
            lf.write(f"[defensefinder_viral] Done -- {len(def_rows)} defense, "
                      f"{len(anti_rows)} antidefense system rows\n")
        Path(str(output.done)).touch()


rule dbapis_viral:
    """
    Anti-defense systems on VIRAL proteins via dbAPIS (Yan et al. 2023, NAR
    -- pro.unl.edu/dbAPIS, moved from the old bcb.unl.edu host; bcb.unl.edu
    now 302-redirects to the bare /dbAPIS homepage instead of the requested
    file, so the old direct /downloads/<file> URLs would have silently
    saved an HTML page as if it were the .pep/.tsv -- confirmed live by
    fetching both old and new URLs), DIAMOND blastp only (no HMMER pass --
    keeps this lightweight, matches the DIAMOND command from the dbAPIS
    README).
    Sequence-similarity-based (no genetic-architecture rule), unlike
    DefenseFinder/MacSyFinder -- complementary detector for single scattered
    anti-defense genes in small phage genomes, kept separate and never
    merged with defensefinder_viral (same "never merge tiers" rule as AMR
    curated/exploratory, see module docstring).
    DB is tiny (~4.4k curated proteins, a few MB) -- downloaded once into a
    shared cache dir, same auto-populate pattern as card_db/deeparg_db.

    Also downloads seed_and_familyrep_all_infor.tsv -- one row per APIS
    family (APIS001, APIS002, ...) with a short characterized gene name
    ("APIS genes", e.g. "Apyc1") and a readable inhibited-defense-system
    label ("Defense systems", e.g. "pyrimidine cyclase system for
    antiphage resistance (Pycsar)") -- confirmed against a real download
    2026-06-23. Used by load_dbapis_viral (scripts/report/data_loaders.py)
    to translate the bare family/gene ID a dbAPIS hit reports into that
    readable gene name + defense-system label.
    """
    input:
        faa  = rules.prodigal_viral.output.faa,
        done = rules.prodigal_viral.output.done,
    output:
        done = f"{OUTDIR}/{{sample}}/viral/dbapis/done.txt",
        hits = f"{OUTDIR}/{{sample}}/viral/dbapis/dbapis_hits.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/dbapis_viral.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/dbapis_viral.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("diamond")
    threads: THREADS
    params:
        outdir   = f"{OUTDIR}/{{sample}}/viral/dbapis",
        apis_dir = APIS_DB or f"{OUTDIR}/dbapis_db",
        enabled  = DEFENSE_AMR_VIRAL_ENABLED,
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        if [ "{params.enabled}" != "True" ] || [ ! -s {input.faa} ]; then
            echo "[dbapis_viral] Disabled or no viral proteins -- skipping" | tee {log}
            printf "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\n" > {output.hits}
            touch {output.done}
            exit 0
        fi

        APIS_DIR="{params.apis_dir}"
        mkdir -p "$APIS_DIR"

        if [ ! -s "$APIS_DIR/APIS_db.dmnd" ]; then
            echo "[dbapis_viral] Building dbAPIS Diamond DB in $APIS_DIR" | tee -a {log}
            wget -q -O "$APIS_DIR/anti_defense.pep" \
                "https://pro.unl.edu/dbAPIS/download_file.php?file=anti_defense.pep" \
                >> {log} 2>&1 || echo "[dbapis_viral] WARNING: wget anti_defense.pep failed" >> {log}
            # Family -> (gene name, inhibited defense-type) mapping, one row
            # per APIS family -- confirmed against a real download 2026-06-23
            # (columns: "APIS families", "APIS genes", "Defense systems", ...).
            wget -q -O "$APIS_DIR/seed_and_familyrep_all_infor.tsv" \
                "https://pro.unl.edu/dbAPIS/download_file.php?file=seed_and_familyrep_all_infor.tsv" \
                >> {log} 2>&1 || echo "[dbapis_viral] WARNING: wget seed_and_familyrep_all_infor.tsv failed (mapping disabled, family IDs still reported)" >> {log}
            if [ -s "$APIS_DIR/anti_defense.pep" ]; then
                diamond makedb --in "$APIS_DIR/anti_defense.pep" -d "$APIS_DIR/APIS_db" >> {log} 2>&1
            fi
        fi

        if [ -s "$APIS_DIR/APIS_db.dmnd" ]; then
            diamond blastp --db "$APIS_DIR/APIS_db" -q {input.faa} \
                -f 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen \
                --max-target-seqs 25 --evalue 1e-10 --id 30 \
                -o {output.hits}.raw --threads {threads} >> {log} 2>&1
            printf "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\n" > {output.hits}
            cat {output.hits}.raw >> {output.hits}
            rm -f {output.hits}.raw
        else
            echo "[dbapis_viral] No dbAPIS DB available -- writing empty hits" | tee -a {log}
            printf "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\n" > {output.hits}
        fi
        touch {output.done}
        echo "[dbapis_viral] Done -- $(( $(wc -l < {output.hits}) - 1 )) hits" | tee -a {log}
        """


rule amrfinderplus:
    """
    AMRFinderPlus -- curated, alignment-based AMR gene + point-mutation
    detection (NCBI Reference Gene/Point Mutation databases).
    Gene-level: runs once on the concatenated protein set (all genome
    units from prok_bin_proteins), unlike DefenseFinder.
    Reuses env_annotation.yaml, which already pins ncbi-amrfinderplus.
    """
    input:
        manifest = rules.prok_bin_proteins.output.manifest,
        done     = rules.prok_bin_proteins.output.done,
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/amrfinderplus/done.txt",
        results = f"{OUTDIR}/{{sample}}/bins/amrfinderplus/amrfinder_results.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/amrfinderplus.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/amrfinderplus.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("ncbi-amrfinderplus")
    threads: THREADS
    params:
        outdir  = f"{OUTDIR}/{{sample}}/bins/amrfinderplus",
        enabled = DEFENSE_AMR_ENABLED,
    run:
        import os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        all_faa = os.path.join(params.outdir, "all_genomes.faa")

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.results)).write_text("Protein identifier\tGene symbol\n")
            Path(str(output.done)).touch()

        if (not params.enabled or not os.path.exists(str(input.manifest))
                or os.path.getsize(str(input.manifest)) == 0):
            write_empty("[amrfinderplus] Disabled or no genome units -- skipping")
            return

        if not _concat_proteins(str(input.manifest), all_faa):
            write_empty("[amrfinderplus] No proteins -- skipping")
            return

        shell("amrfinder -u >> {log} 2>&1 || "
              "echo '[amrfinderplus] WARNING: database update failed (may already be current)' >> {log}")
        shell(
            "amrfinder -p {all_faa} --plus -o {output.results} --threads {threads} "
            ">> {log} 2>&1 || echo '[amrfinderplus] WARNING: amrfinder failed' >> {log}"
        )

        if not os.path.exists(str(output.results)) or os.path.getsize(str(output.results)) == 0:
            Path(str(output.results)).write_text("Protein identifier\tGene symbol\n")
        Path(str(output.done)).touch()


rule rgi_card:
    """
    RGI (CARD) -- curated AMR detection via the Comprehensive Antibiotic
    Resistance Database (homology + SNP models). Same concatenated batch
    input as AMRFinderPlus; reported alongside it as "curated" AMR.
    CARD database path is configurable (card_db) -- auto-downloaded and
    loaded into that directory on first use (public CARD download URL),
    so a fresh clone is reproducible without a manual prep step.
    """
    input:
        manifest = rules.prok_bin_proteins.output.manifest,
        done     = rules.prok_bin_proteins.output.done,
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/rgi/done.txt",
        results = f"{OUTDIR}/{{sample}}/bins/rgi/rgi_results.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/rgi.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/rgi.tsv"
    conda: "../envs/env_rgi.yaml"
    container:  CONTAINERS.get("rgi")
    threads: THREADS
    params:
        outdir  = f"{OUTDIR}/{{sample}}/bins/rgi",
        card_db = CARD_DB,
        enabled = DEFENSE_AMR_ENABLED,
    run:
        import os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        all_faa = os.path.join(params.outdir, "all_genomes.faa")

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.results)).write_text("ORF_ID\tBest_Hit_ARO\n")
            Path(str(output.done)).touch()

        if (not params.enabled or not os.path.exists(str(input.manifest))
                or os.path.getsize(str(input.manifest)) == 0):
            write_empty("[rgi] Disabled or no genome units -- skipping")
            return

        if not _concat_proteins(str(input.manifest), all_faa):
            write_empty("[rgi] No proteins -- skipping")
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
            write_empty("[rgi] WARNING: CARD database unavailable -- skipping")
            return

        # --local makes RGI create/read ./localDB relative to the CURRENT
        # working directory, not the --card_json path -- 'load' and 'main'
        # must run from the exact same directory or 'main' can't find the
        # database 'load' just built.
        shell("cd {params.outdir} && rgi load --card_json {card_json} --local >> {log} 2>&1")
        shell(
            "cd {params.outdir} && rgi main -i {all_faa} -t protein --include_loose --local "
            "-o rgi_results -n {threads} >> {log} 2>&1 || "
            "echo '[rgi] WARNING: rgi main failed' >> {log}"
        )

        # 'rgi main -o rgi_results' inside params.outdir already writes
        # directly to output.results (rgi_results.txt) -- same path, no
        # copy needed (shutil.copy onto an identical path raises
        # SameFileError).
        produced = os.path.join(params.outdir, "rgi_results.txt")
        if not os.path.exists(produced) or os.path.getsize(produced) == 0:
            Path(str(output.results)).write_text("ORF_ID\tBest_Hit_ARO\n")
        Path(str(output.done)).touch()


rule deeparg:
    """
    DeepARG -- deep-learning AMR gene prediction (CNN trained on
    CARD/ARDB/UniProt-derived sequences). Exploratory/sensitivity-oriented
    complement to AMRFinderPlus+RGI's curated calls (Serrana et al. 2026,
    iMetaOmics): higher recall on divergent/novel environmental ARGs that
    curated homology search misses, but with no accuracy benchmark --
    kept and reported separately, never merged into the curated AMR count.
    Same concatenated batch input as AMRFinderPlus/RGI, so it automatically
    inherits the low-depth contig fallback from prok_bin_proteins.

    NOTE: bioconda's deeparg=1.0.4 is the classic Python2/Theano codebase,
    not the newer PyTorch/HuggingFace rewrite some docs describe -- there is
    no auto-download and no --threads flag. `deeparg predict` requires an
    explicit --data-path (fetched once via `deeparg download_data`).
    """
    input:
        manifest = rules.prok_bin_proteins.output.manifest,
        done     = rules.prok_bin_proteins.output.done,
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/deeparg/done.txt",
        results = f"{OUTDIR}/{{sample}}/bins/deeparg/deeparg_results.mapping.ARG",
    log:
        f"{OUTDIR}/{{sample}}/logs/deeparg.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/deeparg.tsv"
    conda: "../envs/env_deeparg.yaml"
    container:  CONTAINERS.get("deeparg")
    threads: THREADS
    params:
        outdir   = f"{OUTDIR}/{{sample}}/bins/deeparg",
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
            Path(str(output.results)).write_text("#ARG\tquery-start\n")
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
            Path(str(output.results)).write_text("#ARG\tquery-start\n")
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


rule abricate:
    """
    ABRicate -- BLASTN mass screening of contigs/bins for gene presence.
    Used here only for the two databases AMRFinderPlus/RGI/DeepARG do not
    cover -- virulence factors (VFDB) and plasmid replicons (PlasmidFinder)
    -- NOT for AMR gene calling: ABRicate's bundled AMR databases are
    static flat-file BLASTN screens with no point-mutation detection
    (AMRFinderPlus) or SNP/variant models (RGI/CARD), so they would be a
    strict downgrade if used to replace either tool.
    Runs per genome unit (manifest from prok_bin_proteins) on the
    nucleotide sequence -- ABRicate works on DNA via BLASTN, unlike the
    protein-level AMR/defense tools above.
    """
    input:
        manifest = rules.prok_bin_proteins.output.manifest,
        done     = rules.prok_bin_proteins.output.done,
    output:
        done          = f"{OUTDIR}/{{sample}}/bins/abricate/done.txt",
        vfdb          = f"{OUTDIR}/{{sample}}/bins/abricate/vfdb_results.tsv",
        plasmidfinder = f"{OUTDIR}/{{sample}}/bins/abricate/plasmidfinder_results.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/abricate.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/abricate.tsv"
    conda: "../envs/env_abricate.yaml"
    container:  CONTAINERS.get("abricate")
    threads: THREADS
    params:
        outdir  = f"{OUTDIR}/{{sample}}/bins/abricate",
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


rule argnorm_normalize:
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
        amrfinder      = rules.amrfinderplus.output.results,
        amrfinder_done = rules.amrfinderplus.output.done,
        deeparg        = rules.deeparg.output.results,
        deeparg_done   = rules.deeparg.output.done,
    output:
        done             = f"{OUTDIR}/{{sample}}/bins/argnorm/done.txt",
        amrfinder_normed = f"{OUTDIR}/{{sample}}/bins/argnorm/amrfinderplus_normed.tsv",
        deeparg_normed   = f"{OUTDIR}/{{sample}}/bins/argnorm/deeparg_normed.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/argnorm.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/argnorm.tsv"
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
