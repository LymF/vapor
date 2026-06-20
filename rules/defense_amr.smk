# ══════════════════════════════════════════════════════════════════════
# rules/defense_amr.smk — BLOCK 10.5: Defense Systems + AMR (prok bins)
#
# prok_bin_proteins — shared per-genome Prodigal proteins, feeds all 4
#                      tools below. Falls back to the whole (viral-
#                      filtered) contig set as a single pseudo-genome
#                      when a sample has no bins (low depth/coverage).
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
# hamronize_rgi — converts rgi_card's native output into the standardized
#                 format argnorm_normalize needs (argNorm has no native
#                 RGI parser, only AMRFinderPlus/DeepARG/hAMRonization ones)
# argnorm_normalize — maps AMRFinderPlus/RGI(via hAMRonization)/DeepARG gene
#                 calls onto the Antibiotic Resistance Ontology (ARO), so the
#                 same gene reported under different names by the 3 tools can
#                 be compared. Normalized tables stay separate per tool --
#                 same "never merged" rule as the raw outputs.
#
# All tools soft-fail (header-only TSV + done.txt) when disabled
# (defense_amr_enabled / abricate_enabled / argnorm_enabled: false) or no
# genome units / upstream hits exist.
# Config keys: defense_amr_enabled, defense_amr_contig_fallback, card_db,
# abricate_enabled, argnorm_enabled
# ══════════════════════════════════════════════════════════════════════


rule prok_bin_proteins:
    """
    Per-genome protein prediction (Prodigal) shared by every defense/AMR
    annotation tool below -- runs once regardless of how many downstream
    tools consume it (same reuse pattern as rules.prodigal_viral.output.faa
    for viral taxonomy).

    Normal path: one Prodigal call per Binette final bin (single mode).
    Low-depth fallback (no bins -- e.g. shallow/low-coverage samples where
    binning produced nothing): one Prodigal call (meta mode) on the whole
    viral-filtered contig set used as binning input, treated as a single
    pseudo-genome. Gated by defense_amr_contig_fallback (config.yaml);
    every downstream rule just iterates whatever the manifest contains,
    with no separate low-depth codepath of its own.
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
        bins_dir = lambda wc: f"{OUTDIR}/{wc.sample}/bins/binette/final_bins",
        outdir   = f"{OUTDIR}/{{sample}}/bins/proteins",
        fallback = DEFENSE_AMR_CONTIG_FALLBACK,
        enabled  = DEFENSE_AMR_ENABLED,
    run:
        import glob, os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        manifest_rows = []

        with open(str(log[0]), "w") as lf:
            if not params.enabled:
                lf.write("[prok_bin_proteins] defense_amr_enabled=False -- skipping\n")
            else:
                bins = sorted(glob.glob(os.path.join(params.bins_dir, "*.fa")))
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
                elif (params.fallback and os.path.exists(str(input.contigs))
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
        outdir  = f"{OUTDIR}/{{sample}}/bins/defensefinder",
        enabled = DEFENSE_AMR_ENABLED,
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

        shell("defense-finder update >> {log} 2>&1 || "
              "echo '[defensefinder] WARNING: model update failed (may already be cached)' >> {log}")
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
                "defense-finder run -o {genome_out} --antidefensefinder {faa} "
                ">> {log} 2>&1 || echo '[defensefinder] WARNING: failed on {name}' >> {log}"
            )

        def merge(pattern, exclude_substr=None):
            rows, header = [], None
            for tsv in sorted(glob.glob(pattern)):
                if exclude_substr and exclude_substr in os.path.basename(tsv).lower():
                    continue
                genome = os.path.basename(os.path.dirname(tsv))
                with open(tsv) as f:
                    r = csv.reader(f, delimiter="\t")
                    h = next(r, None)
                    if h is None:
                        continue
                    if header is None:
                        header = h
                    for row in r:
                        rows.append([genome] + row)
            return header, rows

        def write(path, header, rows):
            with open(path, "w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["genome"] + (header or []))
                w.writerows(rows)

        # AntiDefenseFinder (--antidefensefinder) output files are distinguished
        # from defense-system tables by filename; if the installed tool version
        # names them differently, the antidefense table comes back header-only --
        # check {params.outdir}/<genome>/ directly for the actual filenames.
        def_h, def_rows   = merge(os.path.join(params.outdir, "*", "*systems.tsv"),
                                   exclude_substr="anti")
        anti_h, anti_rows = merge(os.path.join(params.outdir, "*", "*anti*systems.tsv"))
        write(str(output.systems), def_h, def_rows)
        write(str(output.antisystems), anti_h, anti_rows)

        with open(str(log[0]), "a") as lf:
            lf.write(f"[defensefinder] Done -- {len(def_rows)} defense, "
                      f"{len(anti_rows)} antidefense system rows\n")
        Path(str(output.done)).touch()


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


rule hamronize_rgi:
    """
    hAMRonization -- converts RGI's native output into the standardized
    format argNorm can read. argNorm has no native RGI parser (only
    AMRFinderPlus/DeepARG/hAMRonization subcommands), so this bridge rule
    lets argnorm_normalize below treat RGI the same way as the other two.
    """
    input:
        results = rules.rgi_card.output.results,
        done    = rules.rgi_card.output.done,
    output:
        done       = f"{OUTDIR}/{{sample}}/bins/argnorm/hamronize_rgi_done.txt",
        hamronized = f"{OUTDIR}/{{sample}}/bins/argnorm/rgi_hamronized.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/hamronize_rgi.log"
    conda: "../envs/env_argnorm.yaml"
    container:  CONTAINERS.get("hamronization")
    threads: 1
    params:
        card_db = CARD_DB,
        enabled = ARGNORM_ENABLED,
    run:
        import json, os
        from pathlib import Path

        os.makedirs(os.path.dirname(str(output.hamronized)), exist_ok=True)

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.hamronized)).write_text("gene_symbol\n")
            Path(str(output.done)).touch()

        if not params.enabled or not _has_data_rows(str(input.results)):
            write_empty("[hamronize_rgi] Disabled or no RGI hits -- skipping")
            return

        # Same default location rgi_card falls back to when card_db is unset.
        card_dir  = params.card_db or os.path.join(OUTDIR, wildcards.sample, "bins", "rgi", "card_db")
        card_json = os.path.join(card_dir, "card.json")
        card_version = "unknown"
        if os.path.exists(card_json):
            try:
                with open(card_json) as f:
                    card_version = json.load(f).get("_version", "unknown")
            except (json.JSONDecodeError, OSError):
                pass

        shell(
            "hamronize rgi --input_file_name rgi_{wildcards.sample} "
            "--analysis_software_version 6.0.5 "
            "--reference_database_version {card_version} "
            "--format tsv --output {output.hamronized} {input.results} "
            ">> {log} 2>&1 || echo '[hamronize_rgi] WARNING: hamronize failed' >> {log}"
        )

        if not os.path.exists(str(output.hamronized)) or os.path.getsize(str(output.hamronized)) == 0:
            Path(str(output.hamronized)).write_text("gene_symbol\n")
        Path(str(output.done)).touch()


rule argnorm_normalize:
    """
    argNorm -- maps AMRFinderPlus/RGI/DeepARG gene calls onto the
    Antibiotic Resistance Ontology (ARO), so the same gene reported under
    different names by different tools can be compared. RGI goes through
    hamronize_rgi first; AMRFinderPlus and DeepARG are read directly via
    their own argnorm subcommands. The 3 normalized tables stay separate
    -- same "never merged" rule as the raw AMR outputs above, this only
    adds a common ARO/drug-class vocabulary to each.
    """
    input:
        amrfinder      = rules.amrfinderplus.output.results,
        amrfinder_done = rules.amrfinderplus.output.done,
        rgi_hamronized = rules.hamronize_rgi.output.hamronized,
        rgi_done       = rules.hamronize_rgi.output.done,
        deeparg        = rules.deeparg.output.results,
        deeparg_done   = rules.deeparg.output.done,
    output:
        done             = f"{OUTDIR}/{{sample}}/bins/argnorm/done.txt",
        amrfinder_normed = f"{OUTDIR}/{{sample}}/bins/argnorm/amrfinderplus_normed.tsv",
        rgi_normed       = f"{OUTDIR}/{{sample}}/bins/argnorm/rgi_normed.tsv",
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
            stub(output.rgi_normed)
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

        if _has_data_rows(str(input.rgi_hamronized)):
            shell(
                "argnorm hamronization -i {input.rgi_hamronized} -o {output.rgi_normed} "
                "--hamronization_skip_unsupported_tool "
                ">> {log} 2>&1 || echo '[argnorm] WARNING: rgi normalization failed' >> {log}"
            )
        if not os.path.exists(str(output.rgi_normed)) or os.path.getsize(str(output.rgi_normed)) == 0:
            stub(output.rgi_normed)

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
