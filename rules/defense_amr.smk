# ══════════════════════════════════════════════════════════════════════
# rules/defense_amr.smk — BLOCK 10.5: Defense Systems + AMR (prok bins)
#
# prok_bin_proteins — shared per-genome Prodigal proteins, feeds all 5
#                      tools below. Falls back to the whole (viral-
#                      filtered) contig set as a single pseudo-genome
#                      when a sample has no bins (low depth/coverage).
#
# defensefinder — MacSyFinder anti-phage defense systems + built-in
#                 AntiDefenseFinder (-a flag); per-genome (architecture-aware)
# padloc        — complementary defense-system detection (own HMM models);
#                 per-genome, reported separately (different nomenclature)
# amrfinderplus — curated AMR genes/point mutations (NCBI); gene-level,
#                 single batch call on the concatenated protein set
# rgi_card      — curated AMR genes (CARD); gene-level, single batch call
# deeparg       — deep-learning AMR genes; exploratory/sensitivity-oriented
#                 complement to amrfinderplus+rgi's curated calls
#                 (Serrana et al. 2026) -- reported separately, never merged
#
# All five tools soft-fail (header-only TSV + done.txt) when disabled
# (defense_amr_enabled: false) or no genome units exist.
# Config keys: defense_amr_enabled, defense_amr_contig_fallback, card_db
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
                            manifest_rows.append((name, "bins", faa, gff))
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
                        manifest_rows.append((name, "contigs", faa, gff))
                else:
                    lf.write("[prok_bin_proteins] No bins and fallback disabled/no contigs "
                             "-- skipping\n")

            with open(str(output.manifest), "w") as mf:
                for name, mode, faa, gff in manifest_rows:
                    mf.write(f"{name}\t{mode}\t{faa}\t{gff}\n")

            lf.write(f"[prok_bin_proteins] {len(manifest_rows)} genome unit(s) in manifest\n")

        Path(str(output.done)).touch()


def _read_manifest(path):
    """Yield (name, mode, faa, gff) tuples from a prok_bin_proteins manifest."""
    import os
    if not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            yield parts[0], parts[1], parts[2], parts[3]


def _concat_proteins(manifest_path, out_faa):
    """Concatenate every genome unit's proteins into one FASTA, prefixing
    each header with its genome/bin name so gene-level hits (AMRFinderPlus,
    RGI, DeepARG) can be attributed back to a bin or the contig fallback."""
    import os
    wrote_any = False
    with open(out_faa, "w") as out_f:
        for name, mode, faa, gff in _read_manifest(manifest_path):
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

        for name, mode, faa, gff in _read_manifest(str(input.manifest)):
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


rule padloc:
    """
    PADLOC -- complementary defense-system detection (own HMM models +
    MacSyFinder-style architecture rules). Reported separately from
    DefenseFinder rather than merged: system nomenclature differs between
    the two tools, and each catches systems the other misses.
    """
    input:
        manifest = rules.prok_bin_proteins.output.manifest,
        done     = rules.prok_bin_proteins.output.done,
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/padloc/done.txt",
        systems = f"{OUTDIR}/{{sample}}/bins/padloc/padloc_systems.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/padloc.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/padloc.tsv"
    conda: "../envs/env_defense.yaml"
    container:  CONTAINERS.get("padloc")
    threads: THREADS
    params:
        outdir  = f"{OUTDIR}/{{sample}}/bins/padloc",
        enabled = DEFENSE_AMR_ENABLED,
    run:
        import csv, glob, os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.systems)).write_text("genome\n")
            Path(str(output.done)).touch()

        if (not params.enabled or not os.path.exists(str(input.manifest))
                or os.path.getsize(str(input.manifest)) == 0):
            write_empty("[padloc] Disabled or no genome units -- skipping")
            return

        shell("padloc --db-update >> {log} 2>&1 || "
              "echo '[padloc] WARNING: db update failed (may already be cached)' >> {log}")

        for name, mode, faa, gff in _read_manifest(str(input.manifest)):
            if not os.path.exists(faa) or os.path.getsize(faa) == 0:
                continue
            genome_out = os.path.join(params.outdir, name)
            os.makedirs(genome_out, exist_ok=True)
            shell(
                "padloc --faa {faa} --gff {gff} --outdir {genome_out} --cpu {threads} "
                ">> {log} 2>&1 || echo '[padloc] WARNING: failed on {name}' >> {log}"
            )

        rows, header = [], None
        for csv_f in sorted(glob.glob(os.path.join(params.outdir, "*", "*_padloc.csv"))):
            genome = os.path.basename(os.path.dirname(csv_f))
            with open(csv_f) as f:
                r = csv.reader(f)
                h = next(r, None)
                if h is None:
                    continue
                if header is None:
                    header = h
                for row in r:
                    rows.append([genome] + row)

        with open(str(output.systems), "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["genome"] + (header or []))
            w.writerows(rows)

        with open(str(log[0]), "a") as lf:
            lf.write(f"[padloc] Done -- {len(rows)} system rows\n")
        Path(str(output.done)).touch()


rule amrfinderplus:
    """
    AMRFinderPlus -- curated, alignment-based AMR gene + point-mutation
    detection (NCBI Reference Gene/Point Mutation databases).
    Gene-level: runs once on the concatenated protein set (all genome
    units from prok_bin_proteins), unlike DefenseFinder/PADLOC.
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
        import os, shutil
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
                ">> {log} 2>&1 && tar -xjf {card_dir}/card_data.tar.bz2 -C {card_dir} card.json "
                ">> {log} 2>&1 || true"
            )

        if not os.path.exists(card_json) or os.path.getsize(card_json) == 0:
            write_empty("[rgi] WARNING: CARD database unavailable -- skipping")
            return

        shell("rgi load --card_json {card_json} --local >> {log} 2>&1")
        shell(
            "cd {params.outdir} && rgi main -i {all_faa} -t protein --include_loose --local "
            "-o rgi_results -n {threads} >> {log} 2>&1 || "
            "echo '[rgi] WARNING: rgi main failed' >> {log}"
        )

        produced = os.path.join(params.outdir, "rgi_results.txt")
        if os.path.exists(produced) and os.path.getsize(produced) > 0:
            shutil.copy(produced, str(output.results))
        else:
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
        outdir  = f"{OUTDIR}/{{sample}}/bins/deeparg",
        enabled = DEFENSE_AMR_ENABLED,
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

        out_prefix = os.path.join(params.outdir, "deeparg_results")
        # First run auto-downloads model + DIAMOND DB from Hugging Face into the
        # env's HF cache -- no config path required.
        shell(
            "deeparg predict --model LS --type prot -i {all_faa} -o {out_prefix} "
            "--threads {threads} >> {log} 2>&1 || "
            "echo '[deeparg] WARNING: deeparg predict failed' >> {log}"
        )

        if not os.path.exists(str(output.results)) or os.path.getsize(str(output.results)) == 0:
            Path(str(output.results)).write_text("#ARG\tquery-start\n")
        Path(str(output.done)).touch()
