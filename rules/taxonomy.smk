# ══════════════════════════════════════════════════════════════════════
# rules/taxonomy.smk — BLOCK 9: Viral Taxonomy
#
# No fixed source priority -- rule viral_taxonomy compares the rank DEPTH
# each source resolves to and takes the deepest, ties broken by trust order.
# Sources:
#   vConTACT3       — protein-sharing network (genus-level, most specific)
#   MMseqs2/INPHARED — real per-query LCA against INPHARED (rule
#                      mmseqs_taxonomy_viral), contig-level LCA-of-LCAs
#                      across that contig's proteins. Replaced an earlier
#                      Diamond/INPHARED best-hit+majority-vote tier (removed
#                      2026-06-22 -- see memory project_viral_taxonomy_merge).
#   MMseqs2/Custom  — optional user-supplied seqTaxDB (e.g. IMG/VR, rule
#                      mmseqs_taxonomy_custom_viral), same LCA approach.
#                      Replaced an earlier diamond_custom_viral best-hit+
#                      majority-vote rule entirely (removed 2026-06-23).
#   GeNomad         — marker-gene lineage string, parsed to deepest level
# See memory project_viral_taxonomy_merge for the full rationale.
#
# Also contains mmseqs_taxonomy_prok for prokaryote bin taxonomy (separate
# from the viral merge above; diamond_custom_prok was removed -- mmseqs
# replaces it outright, not just outranks it).
# ══════════════════════════════════════════════════════════════════════


def _mmseqs_lca_rollup(hits_path, ranks):
    """Shared by mmseqs_taxonomy_viral (INPHARED) and
    mmseqs_taxonomy_custom_viral (e.g. IMG/VR) results -- both produce the
    same (qseqid, taxid, rank, name, lineage) shape over the same 8-level
    ICTV rank scheme (realm..genus), just against different seqTaxDBs.

    mmseqs operates per-PROTEIN; roll up to per-CONTIG by taking the longest
    common prefix (a second, contig-level LCA) across that contig's own
    already-LCA-resolved protein lineages -- avoids one well-conserved
    protein dominating the whole-genome call. Returns
    {contig: {order, family, subfamily, genus, rank, lineage, n_proteins}}.
    """
    import csv, os, re, collections
    _RANK_PREFIX = re.compile(r'^[a-z]+_')
    contig_protein_lineages = collections.defaultdict(list)
    if os.path.exists(hits_path) and os.path.getsize(hits_path) > 0:
        with open(hits_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                protein_id = row.get("qseqid", "")
                rank       = (row.get("rank") or "").strip()
                lineage    = (row.get("lineage") or "").strip()
                if not protein_id or not lineage or rank in ("", "no rank", "root"):
                    continue
                contig = "_".join(protein_id.split("_")[:-1]) or protein_id
                names = [_RANK_PREFIX.sub("", p).strip() for p in lineage.split(";") if p.strip()]
                if names:
                    contig_protein_lineages[contig].append(names)

    result = {}
    for contig, lineages in contig_protein_lineages.items():
        common = []
        for level in zip(*lineages):
            if len(set(level)) == 1:
                common.append(level[0])
            else:
                break
        if not common:
            continue
        depth = len(common)
        result[contig] = {
            "order":      common[4] if depth >= 5 else "",
            "family":     common[5] if depth >= 6 else "",
            "subfamily":  common[6] if depth >= 7 else "",
            "genus":      common[7] if depth >= 8 else "",
            "rank":       ranks[depth - 1] if depth <= len(ranks) else ranks[-1],
            "lineage":    ";".join(common),
            "n_proteins": len(lineages),
        }
    return result


rule prodigal_viral:
    """Predict ORFs from the vOTU MQ+ representatives for taxonomy searches."""
    input:
        viral = rules.votu_catalog_reps.output.mq_fasta,
    output:
        faa  = f"{OUTDIR}/{{sample}}/viral/taxonomy/viral_proteins.faa",
        done = f"{OUTDIR}/{{sample}}/viral/taxonomy/prodigal_done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/prodigal_viral.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/prodigal_viral.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("prodigal")
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.faa})
        if [ ! -s {input.viral} ]; then
            touch {output.faa} {output.done}; exit 0
        fi
        prodigal -i {input.viral} -a {output.faa} -p meta -f gff > {log} 2>&1
        touch {output.done}
        """


rule mmseqs_taxonomy_viral:
    """
    The pipeline's only INPHARED-based taxonomy source (replaced an earlier
    Diamond/INPHARED best-hit+majority-vote rule on 2026-06-22 -- removed
    entirely, see memory project_viral_taxonomy_merge). MMseqs2 `taxonomy`
    computes a real per-query lowest-common-ancestor against an
    INPHARED-derived seqTaxDB, avoiding "spurious specificity" (von
    Meijenfeldt et al. 2019, CAT/BAT) the same way mmseqs_taxonomy_prok does
    for IMG_NR -- relevant here because divergent environmental phages
    (INPHARED's own niche relative to RefSeq) are exactly where best-hit
    identity cutoffs are least reliable.

    Same pattern as mmseqs_taxonomy_prok: the seqTaxDB is NOT built here.
    Pre-build it ONCE with scripts/prepare_mmseqs_taxdb.py --format inphared
    (see INSTALL.md) -- skipped gracefully if missing. Deliberately not
    auto-built inline: the seqTaxDB lives under the shared INPHARED_DB, not
    a per-sample Snakemake output, so auto-building on first sight would
    race across samples under --cores >1.

    viral_taxonomy compares this against vConTACT3/Diamond-custom/GeNomad by
    resolved rank depth, not a fixed priority order, so this can win, lose,
    or tie per contig.
    """
    input:
        faa  = rules.prodigal_viral.output.faa,
        done = rules.prodigal_viral.output.done,
    output:
        hits = f"{OUTDIR}/{{sample}}/viral/taxonomy/mmseqs_vs_inphared.tsv",
        done = f"{OUTDIR}/{{sample}}/viral/taxonomy/mmseqs_inphared_done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/mmseqs_taxonomy_viral.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/mmseqs_taxonomy_viral.tsv"
    conda: "../envs/env_assembly.yaml"
    container:  CONTAINERS.get("mmseqs2")
    threads: THREADS
    params:
        seqtaxdb = f"{INPHARED_DB}/inphared_mmseqs_taxdb/seqTaxDB",
        # derivado do output, nao de {{sample}}: regra herdada por
        # coassembly.smk via `use rule ... as ... with:` (wildcard {group}).
        outdir   = lambda wc, output: os.path.dirname(output.done),
        querydb  = lambda wc, output: os.path.join(os.path.dirname(output.done), "queryDB"),
        result   = lambda wc, output: os.path.join(os.path.dirname(output.done), "result"),
        tmp      = lambda wc, output: os.path.join(os.path.dirname(output.done), "tmp"),
    run:
        import os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        header = "qseqid\ttaxid\trank\tname\tlineage\n"

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.hits)).write_text(header)
            Path(str(output.done)).touch()

        if not os.path.exists(str(params.seqtaxdb) + ".dbtype"):
            write_empty(
                "[mmseqs_taxonomy_viral] No seqTaxDB at " + str(params.seqtaxdb) +
                " -- run scripts/prepare_mmseqs_taxdb.py --format inphared once first (see INSTALL.md). Skipping."
            )
            return

        if not os.path.exists(str(input.faa)) or os.path.getsize(str(input.faa)) == 0:
            write_empty("[mmseqs_taxonomy_viral] No viral proteins -- skipping")
            return

        # mmseqs taxonomy refuses to run if its output DB already exists
        # ("result.dbtype exists already!") -- happens on any retry after a
        # previous attempt got partway through (same fix as mmseqs_taxonomy_prok,
        # confirmed live on litrp4). querydb doesn't need this: createdb's own
        # message confirms it overwrites an existing one on its own.
        shell("rm -rf {params.tmp} {params.result}*; mkdir -p {params.tmp}")
        shell("mmseqs createdb {input.faa} {params.querydb} >> {log} 2>&1")
        shell(
            "mmseqs taxonomy {params.querydb} {params.seqtaxdb} {params.result} {params.tmp} "
            "--threads {threads} --tax-lineage 1 >> {log} 2>&1"
        )
        shell(
            "mmseqs createtsv {params.querydb} {params.result} {output.hits}.raw >> {log} 2>&1"
        )
        Path(str(output.hits)).write_text(header)
        if os.path.exists(str(output.hits) + ".raw"):
            with open(str(output.hits) + ".raw") as f, open(str(output.hits), "a") as out:
                out.writelines(f)
            os.remove(str(output.hits) + ".raw")
        Path(str(output.done)).touch()


rule mmseqs_taxonomy_custom_viral:
    """
    Optional custom viral MMseqs2 seqTaxDB (e.g. IMG/VR) -- real per-query
    LCA, same approach as mmseqs_taxonomy_viral but against a user-supplied
    DB instead of INPHARED. Replaces an earlier diamond_custom_viral
    best-hit+majority-vote rule entirely (removed 2026-06-23, see memory
    project_viral_taxonomy_merge).

    Same pattern as mmseqs_taxonomy_prok/mmseqs_taxonomy_viral: the seqTaxDB
    is NOT built here -- pre-build it once with scripts/prepare_mmseqs_taxdb.py
    (e.g. --format imgvr, see INSTALL.md). Skipped gracefully if
    custom_viral_mmseqs_db isn't configured, same as any other optional DB.
    """
    input:
        faa  = rules.prodigal_viral.output.faa,
        done = rules.prodigal_viral.output.done,
    output:
        hits = f"{OUTDIR}/{{sample}}/viral/taxonomy/mmseqs_vs_custom.tsv",
        done = f"{OUTDIR}/{{sample}}/viral/taxonomy/mmseqs_custom_viral_done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/mmseqs_taxonomy_custom_viral.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/mmseqs_taxonomy_custom_viral.tsv"
    conda: "../envs/env_assembly.yaml"
    container:  CONTAINERS.get("mmseqs2")
    threads: THREADS
    params:
        seqtaxdb = CUSTOM_VIRAL_MMSEQS_DB,
        outdir   = f"{OUTDIR}/{{sample}}/viral/taxonomy/mmseqs_custom",
        querydb  = f"{OUTDIR}/{{sample}}/viral/taxonomy/mmseqs_custom/queryDB",
        result   = f"{OUTDIR}/{{sample}}/viral/taxonomy/mmseqs_custom/result",
        tmp      = f"{OUTDIR}/{{sample}}/viral/taxonomy/mmseqs_custom/tmp",
    run:
        import os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        header = "qseqid\ttaxid\trank\tname\tlineage\n"

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.hits)).write_text(header)
            Path(str(output.done)).touch()

        if not params.seqtaxdb or not os.path.exists(str(params.seqtaxdb) + ".dbtype"):
            write_empty("[mmseqs_taxonomy_custom_viral] No custom_viral_mmseqs_db configured -- skipping")
            return

        if not os.path.exists(str(input.faa)) or os.path.getsize(str(input.faa)) == 0:
            write_empty("[mmseqs_taxonomy_custom_viral] No viral proteins -- skipping")
            return

        # mmseqs taxonomy refuses to run if its output DB already exists --
        # same lesson as mmseqs_taxonomy_viral/mmseqs_taxonomy_prok.
        shell("rm -rf {params.tmp} {params.result}*; mkdir -p {params.tmp}")
        shell("mmseqs createdb {input.faa} {params.querydb} >> {log} 2>&1")
        shell(
            "mmseqs taxonomy {params.querydb} {params.seqtaxdb} {params.result} {params.tmp} "
            "--threads {threads} --tax-lineage 1 >> {log} 2>&1"
        )
        shell(
            "mmseqs createtsv {params.querydb} {params.result} {output.hits}.raw >> {log} 2>&1"
        )
        Path(str(output.hits)).write_text(header)
        if os.path.exists(str(output.hits) + ".raw"):
            with open(str(output.hits) + ".raw") as f, open(str(output.hits), "a") as out:
                out.writelines(f)
            os.remove(str(output.hits) + ".raw")
        Path(str(output.done)).touch()


rule mmseqs_taxonomy_prok:
    """
    Primary source for custom prokaryote taxonomy: MMseqs2 `taxonomy`
    against a custom IMG_NR seqTaxDB (scripts/prepare_mmseqs_taxdb.py)
    instead of DIAMOND blastp + best-hit/majority-vote. Computes a real
    lowest-common-ancestor per query across all its hits -- avoids the
    "spurious specificity" best-hit problem (von Meijenfeldt et al. 2019,
    CAT/BAT), more relevant here than usual since IMG_NR exists specifically
    to cover environmental/divergent genomes standard databases
    under-represent.

    Output (qseqid, taxid, rank, name, lineage) is per-PROTEIN; the report
    loader (load_mmseqs_taxonomy_prok in scripts/report/data_loaders.py)
    aggregates to genome level with a second LCA pass across each genome
    unit's own proteins -- a vote would reintroduce the same problem this
    rule exists to avoid. Under low_depth_mode the genome unit is each
    contig individually, not the whole sample's pooled pseudo-genome (see
    _prok_genome_unit), so taxonomy isn't blended across organisms.

    Replaces diamond_custom_prok entirely (removed) -- custom prokaryote
    taxonomy now runs exclusively through this rule. Skips gracefully if
    custom_prok_mmseqs_db isn't configured, same as any other optional DB.
    """
    input:
        manifest = rules.prok_bin_proteins.output.manifest,
        done     = rules.prok_bin_proteins.output.done,
    output:
        hits = f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/taxonomy.tsv",
        done = f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/mmseqs_taxonomy_prok.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/mmseqs_taxonomy_prok.tsv"
    conda: "../envs/env_assembly.yaml"
    container:  CONTAINERS.get("mmseqs2")
    threads: THREADS
    params:
        seqtaxdb = CUSTOM_PROK_MMSEQS_DB,
        outdir   = f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok",
        prok_faa = f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/all_bins.faa",
        querydb  = f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/queryDB",
        result   = f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/result",
        tmp      = f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/tmp",
    run:
        import os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        header = "qseqid\ttaxid\trank\tname\tlineage\n"

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.hits)).write_text(header)
            Path(str(output.done)).touch()

        if not params.seqtaxdb or not os.path.exists(str(params.seqtaxdb) + ".dbtype"):
            write_empty("[mmseqs_taxonomy_prok] No custom_prok_mmseqs_db configured -- skipping")
            return

        if not _concat_proteins(str(input.manifest), params.prok_faa):
            write_empty("[mmseqs_taxonomy_prok] No genome-unit proteins found")
            return

        shell("rm -rf {params.tmp} {params.result}*; mkdir -p {params.tmp}")
        shell("mmseqs createdb {params.prok_faa} {params.querydb} >> {log} 2>&1")
        shell(
            "mmseqs taxonomy {params.querydb} {params.seqtaxdb} {params.result} {params.tmp} "
            "--threads {threads} --tax-lineage 1 >> {log} 2>&1"
        )
        shell(
            "mmseqs createtsv {params.querydb} {params.result} {output.hits}.raw >> {log} 2>&1"
        )
        Path(str(output.hits)).write_text(header)
        if os.path.exists(str(output.hits) + ".raw"):
            with open(str(output.hits) + ".raw") as f, open(str(output.hits), "a") as out:
                out.writelines(f)
            os.remove(str(output.hits) + ".raw")
        Path(str(output.done)).touch()


rule vcontact3:
    """
    vConTACT3: protein-sharing network clustering, genus-level.
    Generally the most specific source when it resolves (compared by rank
    depth against other sources in viral_taxonomy, not a fixed priority).
    Runs on HQ+/Complete vOTU representatives >= 10 kb only (votu_catalog_reps
    hq_10kb_fasta) — shorter or lower-quality sequences add noise to the
    protein-sharing network without meaningful genus-level signal.
    Options: SqRoot metric, 10 iterations, --reduce-memory, family+genus ranks.
    """
    input:
        viral = rules.votu_catalog_reps.output.hq_10kb_fasta,
    output:
        done    = f"{OUTDIR}/{{sample}}/viral/vcontact3/done.txt",
        network = f"{OUTDIR}/{{sample}}/viral/vcontact3/genome_clusters.tsv",
    log:   f"{OUTDIR}/{{sample}}/logs/vcontact3.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/vcontact3.tsv"
    conda: "../envs/env_vcontact3.yaml"
    container:  CONTAINERS.get("vcontact3")
    threads: THREADS
    params:
        outdir  = f"{OUTDIR}/{{sample}}/viral/vcontact3",
        vc3_db  = VCONTACT3_DB,
        vc3_ver = VCONTACT3_VER,
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        if [ ! -s {input.viral} ]; then
            echo "[vcontact3] No HQ+/>=10kb genomes — skipping" | tee -a {log}
            mkdir -p {params.outdir}/vConTACT3_results
            printf "genome\tVC\tVC_status\tgenus\tfamily\torder\n" > {output.network}
            printf "skipped: no HQ+/>=10kb genomes\n" > {output.done}; exit 0
        fi

        N=$(grep -c "^>" {input.viral} || echo 0)
        echo "[vcontact3] Running on $N genomes (HQ+/>=10kb)" | tee -a {log}

        VC3_EXIT=0
        vcontact3 run \
            --nucleotide      {input.viral} \
            --output          {params.outdir}/vConTACT3_results \
            --threads         {threads} \
            --db-path         {params.vc3_db} \
            --db-version      {params.vc3_ver} \
            --db-domain       prokaryotes \
            --distance-metric SqRoot \
            --max-iterations  10 \
            --reduce-memory \
            --target-rank family genus \
            --exports profiles completeness \
            --no-progress \
            --force-overwrite \
            >> {log} 2>&1 || VC3_EXIT=$?
        if [ $VC3_EXIT -ne 0 ]; then
            echo "[vcontact3] WARNING: vConTACT3 exited with code $VC3_EXIT — check log" | tee -a {log}
        fi

        # vConTACT3 v3 writes final_assignments.csv in exports/
        VC3_OUT="{params.outdir}/vConTACT3_results/exports/final_assignments.csv"

        if [ -f "$VC3_OUT" ]; then
            cp "$VC3_OUT" {output.network}
            echo "[vcontact3] Output: $VC3_OUT ($(wc -l < $VC3_OUT) lines)" | tee -a {log}
        else
            echo "[vcontact3] final_assignments.csv not found — listing:" | tee -a {log}
            find {params.outdir}/vConTACT3_results -type f | tee -a {log}
            printf "Genome,family_prediction,genus_prediction,realm_prediction,order_prediction\n" > {output.network}
        fi
        if [ $VC3_EXIT -ne 0 ]; then
            printf "failed: vcontact3 run exit %s\n" "$VC3_EXIT" > {output.done}
        else
            printf "ok\n" > {output.done}
        fi
        """


rule viral_taxonomy:
    """
    Merge taxonomy from all sources into one table per contig.
    Sources: vConTACT3, MMseqs2/INPHARED (LCA), MMseqs2/Custom (LCA), GeNomad.
    No fixed source priority -- each source proposes a (family, genus, order)
    call, and whichever resolves the DEEPEST rank wins per contig (genus >
    family > order > unclassified). Ties are broken by trust order:
    vConTACT3 > mmseqs_inphared > mmseqs_custom > genomad.
    Why: with a fixed priority, an upstream tier's *shallow* hit would always
    beat a downstream tier's deeper, better-supported call (e.g. GeNomad
    resolving to genus) -- see [[project_viral_taxonomy_merge]] memory for
    the full rationale and how each source's depth is computed. That same
    memory documents why the earlier Diamond/INPHARED (2026-06-22) and
    diamond_custom_viral (2026-06-23) best-hit+majority-vote tiers were
    removed in favour of real per-query LCA (mmseqs_inphared/mmseqs_custom).
    Output: viral_taxonomy_merged.tsv
    """
    input:
        genomad_done    = rules.genomad.output.done,
        mmseqs_hits     = rules.mmseqs_taxonomy_viral.output.hits,
        mmseqs_done     = rules.mmseqs_taxonomy_viral.output.done,
        custom_hits     = rules.mmseqs_taxonomy_custom_viral.output.hits,
        custom_done     = rules.mmseqs_taxonomy_custom_viral.output.done,
        vcontact3_net   = rules.vcontact3.output.network,
        vcontact3_done  = rules.vcontact3.output.done,
        viral           = rules.votu_catalog_reps.output.mq_fasta,
    output:
        tsv  = f"{OUTDIR}/{{sample}}/viral/taxonomy/viral_taxonomy_merged.tsv",
        done = f"{OUTDIR}/{{sample}}/viral/taxonomy/taxonomy_done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/viral_taxonomy.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/viral_taxonomy.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("diamond")
    threads: 1
    run:
        import csv, os, glob, collections
        from pathlib import Path

        lf = open(str(log[0]), "w")

        # All viral contigs
        contigs = []
        if os.path.exists(str(input.viral)):
            with open(str(input.viral)) as f:
                for line in f:
                    if line.startswith(">"): contigs.append(line[1:].split()[0])
        lf.write(f"Total viral contigs: {len(contigs)}\n")

        # ── vConTACT3 ─────────────────────────────────────────────────
        vc3_tax = {}
        vc3_path = str(input.vcontact3_net)
        if os.path.exists(vc3_path) and os.path.getsize(vc3_path) > 100:
            with open(vc3_path) as f:
                first = f.readline(); f.seek(0)
                delim = ',' if first.count(',') > first.count('\t') else '\t'
            with open(vc3_path) as f:
                for row in csv.DictReader(f, delimiter=delim):
                    g = row.get("Genome", row.get("genome","")).strip()
                    if not g: continue
                    if row.get("Reference","").lower() in ("true","1","yes"): continue
                    fam  = row.get("family_prediction", row.get("family", row.get("Family","")))
                    gen  = row.get("genus_prediction",  row.get("genus",  row.get("Genus","")))
                    ord_ = row.get("order_prediction",  row.get("order",  row.get("Order","")))
                    # Reconstruct VC_status: vConTACT3 v3 encodes novelty in prediction names.
                    # "singleton" is output literally for genomes with no network neighbours.
                    _NOVEL_VALS = {"singleton", "unclassified", "nd", "none", ""}
                    if not fam or fam.lower().startswith("novel_") or fam.lower() in _NOVEL_VALS:
                        status = "Novel"
                    elif gen and "|" in gen:
                        status = "Shared"
                    else:
                        status = "Assigned"
                    # Extract best-known parent anchor for Novel genomes
                    # e.g. "novel_family_13_of_novel_order_64_of_Caudoviricetes" → "Caudoviricetes"
                    novel_anchor = ""
                    if status == "Novel" and fam and "of_" in fam:
                        novel_anchor = fam.rsplit("of_", 1)[-1].strip()
                    vc3_tax[g] = {
                        "status":       status,
                        "genus":        gen.split("|")[0].strip() if gen else "",
                        "family":       fam if status != "Novel" else "",
                        "order":        ord_ if ord_ and not ord_.lower().startswith("novel_") and ord_.lower() not in _NOVEL_VALS else "",
                        "novel_anchor": novel_anchor,
                    }
        lf.write(f"vConTACT3: {len(vc3_tax)} genomes\n")

        # ── MMseqs2/INPHARED + MMseqs2/Custom (real per-query LCA) ──────
        # Both rules produce the same (qseqid, taxid, rank, name, lineage)
        # shape over the same 8-level ICTV rank scheme -- _mmseqs_lca_rollup
        # (module-level helper, top of this file) does the per-contig rollup
        # for either one.
        _RANKS = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus']
        mmseqs_tax = _mmseqs_lca_rollup(str(input.mmseqs_hits), _RANKS)
        lf.write(f"MMseqs2/INPHARED: {len(mmseqs_tax)} contigs\n")

        custom_tax = _mmseqs_lca_rollup(str(input.custom_hits), _RANKS)
        lf.write(f"MMseqs2/Custom: {len(custom_tax)} contigs\n")

        # ── GeNomad ───────────────────────────────────────────────────
        # Path derivado de input.genomad_done (idêntico ao Snakefile funcional)
        genomad_tax = {}
        gdir   = os.path.dirname(str(input.genomad_done))
        gfiles = glob.glob(os.path.join(gdir, "**", "*_virus_summary.tsv"), recursive=True)
        gpath  = gfiles[0] if gfiles else ""
        if gpath and os.path.getsize(gpath) > 0:
            with open(gpath) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    name  = row.get("seq_name","")
                    tax   = row.get("taxonomy","")
                    score = row.get("virus_score","0")
                    if not name or not tax: continue
                    parts = [p.strip() for p in tax.split(";")
                             if p.strip() and p.strip() not in ("Viruses", "")]
                    family=""; genus=""; order=""; cls=""; best=""
                    # High-rank names (phylum/kingdom) must not be used as fallback
                    _high = set()
                    for p in parts:
                        if p.endswith("viridae") or p.endswith("virnae"):
                            family = p
                        elif p.endswith("virales"):
                            order  = p
                        elif p.endswith("viricetes"):
                            cls    = p
                        elif p.endswith("virus") or p.endswith("phage"):
                            genus  = p
                        elif any(p.endswith(s) for s in
                                 ("viricota", "virae", "viria", "virites")):
                            _high.add(p)  # phylum/kingdom — skip as taxonomy fallback
                    _low = [p for p in parts if p not in _high]
                    best = genus or family or order or cls or (
                        _low[-1] if _low else (parts[-1] if parts else ""))
                    genomad_tax[name] = {
                        "family": family, "genus": genus, "order": order,
                        "class": cls, "best": best, "lineage": tax,
                        "score": float(score or 0),
                    }
        lf.write(f"GeNomad: {len(genomad_tax)} contigs\n")

        # ── Build final table ─────────────────────────────────────────
        # No fixed tier priority: each source proposes (family, genus, order);
        # whichever resolves the deepest rank wins. Ties broken by trust order
        # (vcontact3 > mmseqs_inphared > mmseqs_custom > genomad).
        _PRIORITY = {"vcontact3": 0, "mmseqs_inphared": 1,
                     "mmseqs_custom": 2, "genomad": 3}

        def _depth(genus, family, order):
            if genus:  return 3
            if family: return 2
            if order:  return 1
            return 0

        rows = []; stats = collections.Counter()
        for contig in contigs:
            vc3 = vc3_tax.get(contig, {})
            mms = mmseqs_tax.get(contig, {})
            gmd = genomad_tax.get(contig, {})
            cms = custom_tax.get(contig, {})

            candidates = []  # (source, ff, fg, fo, lin, conf, best)

            vc3_status = vc3.get("status","").strip()
            if (vc3 and vc3_status in ("Assigned", "Shared") and
                    (vc3.get("family") or vc3.get("genus"))):
                ff, fg = vc3.get("family",""), vc3.get("genus","")
                fo     = vc3.get("order","")
                candidates.append(("vcontact3", ff, fg, fo,
                                    ";".join(filter(None, [fo, ff, fg])),
                                    vc3_status, fg or ff or fo))

            if mms and (mms.get("family") or mms.get("genus") or mms.get("order")):
                ff, fg, fo = mms.get("family",""), mms.get("genus",""), mms.get("order","")
                candidates.append(("mmseqs_inphared", ff, fg, fo, mms.get("lineage",""),
                                    f"{mms.get('rank','')} ({mms.get('n_proteins',0)} proteins)",
                                    fg or ff or fo))

            if cms and (cms.get("family") or cms.get("genus") or cms.get("order")):
                ff, fg, fo = cms.get("family",""), cms.get("genus",""), cms.get("order","")
                candidates.append(("mmseqs_custom", ff, fg, fo, cms.get("lineage",""),
                                    f"{cms.get('rank','')} ({cms.get('n_proteins',0)} proteins)",
                                    fg or ff or fo))

            if gmd and gmd.get("family"):  # only true family; no class/order via this candidate
                ff, fg, fo = gmd.get("family",""), gmd.get("genus",""), gmd.get("order","")
                candidates.append(("genomad", ff, fg, fo, gmd.get("lineage",""),
                                    f"{gmd['score']:.3f}", gmd.get("best","")))

            if candidates:
                candidates.sort(key=lambda c: (-_depth(c[2], c[1], c[3]), _PRIORITY[c[0]]))
                source, ff, fg, fo, lin, conf, best = candidates[0]
            elif gmd:
                # GeNomad's only signal is class/order/higher -- still better than nothing
                source = "genomad"
                ff, fg, fo = "", "", gmd.get("order","")
                best   = gmd.get("best","")
                conf   = f"{gmd['score']:.3f}"
                lin    = gmd.get("lineage","")
            else:
                source = "unclassified"
                ff = fg = fo = conf = lin = best = ""

            stats[source] += 1
            rows.append({
                "seq_name":      contig,
                "final_family":  ff,
                "final_genus":   fg,
                "final_order":   fo,
                "best_taxonomy": best,
                "source":        source,
                "confidence":    conf,
                "lineage":       lin,
                "vc3_status":    vc3.get("status",""),
                "vc3_novel_anchor": vc3.get("novel_anchor",""),
                "genomad_best":  gmd.get("best",""),
                "genomad_class": gmd.get("class",""),
                "genomad_score": gmd.get("score",""),
                "mmseqs_rank":         mms.get("rank",""),
                "mmseqs_lineage":      mms.get("lineage",""),
                "mmseqs_n_proteins":   mms.get("n_proteins",""),
                "custom_rank":         cms.get("rank",""),
                "custom_lineage":      cms.get("lineage",""),
                "custom_n_proteins":   cms.get("n_proteins",""),
            })

        lf.write("\nSummary:\n")
        for k, v in stats.most_common():
            lf.write(f"  {k}: {v}\n")
        total = len(rows)
        unclass = stats.get("unclassified", 0)
        lf.write(f"  Novel (unclassified): {unclass}/{total} = {100*unclass/total:.1f}%\n" if total else "")
        lf.close()

        fields = ["seq_name","final_family","final_genus","final_order","best_taxonomy",
                  "source","confidence","lineage",
                  "vc3_status","vc3_novel_anchor",
                  "genomad_best","genomad_class","genomad_score",
                  "mmseqs_rank","mmseqs_lineage","mmseqs_n_proteins",
                  "custom_rank","custom_lineage","custom_n_proteins"]
        with open(str(output.tsv), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader(); w.writerows(rows)
        Path(str(output.done)).write_text("ok\n")
