# ══════════════════════════════════════════════════════════════════════
# rules/votu_catalog.smk — BLOCK 7.5: Global vOTU catalog
#
# Replaces the former per-sample ANI-clustering chain (removed). A vOTU is
# defined ONCE over the pooled viral sets of every sample and co-assembly
# group, so richness is comparable across samples and against the
# literature. Per-sample presence comes from read recruitment against the
# catalog (see votu_catalog_abundance).
#
#   votu_catalog_pool    — concatenate all viral sets, namespacing IDs
#   votu_catalog_skani   — skani triangle --sparse over the pool
#   votu_catalog_cluster — single-linkage vOTUs (ANI >= 95, AF >= 85)
#   votu_catalog_reps    — all / MQ+ / HQ+>=10kb representative FASTAs
# ══════════════════════════════════════════════════════════════════════

import sys as _sys
_sys.path.insert(0, SCRIPTS_DIR)
from votu_catalog import (
    build_pool, parse_skani_sparse, cluster_votus, write_clusters,
)

CATALOG_DIR = f"{OUTDIR}/votu_catalog"


def _catalog_sources():
    """(source_type, source_id, fasta_path) for every viral set entering the pool."""
    sources = [
        ("sample", s,
         f"{OUTDIR}/{s}/viral/consensus/{s}_viral_nonredundant.fasta")
        for s in SAMPLES
    ]
    if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:
        # coassembly_vrhyme/coassembly_viral_nonredundant only exist for
        # short-read groups (rules/coassembly.smk gates that whole block on
        # `not LONG_READS` — no group-level vRhyme for long reads). Long-read
        # groups keep pointing at the pre-binning trimmed set, same as before
        # item (e); short-read groups now point past vRhyme's bins-first +
        # composite gate (rule coassembly_viral_nonredundant), so group vMAGs
        # actually reach the catalog instead of being computed and discarded.
        sources += [
            ("group", g,
             (f"{OUTDIR}/coassembly/{g}/viral/consensus/{g}_viral_nonredundant.fasta"
              if not LONG_READS else
              f"{OUTDIR}/coassembly/{g}/viral/checkv/{g}_viral_trimmed.fasta"))
            for g in GROUPS
        ]
    return sources


def _catalog_input_fastas(wildcards):
    return [path for _, _, path in _catalog_sources()]


def _catalog_checkv_pairs():
    """(source_id, checkv_summary_path) in the SAME order as _catalog_sources().

    Returned as pairs, not a bare path list, on purpose: re-keying CheckV
    completeness onto the pool's namespaced IDs requires knowing which source
    each summary belongs to. Zipping two independently built lists by position
    would produce a silently mis-keyed catalog -- the same class of quiet
    wrong answer this whole stage exists to eliminate.
    """
    pairs = [(s, f"{OUTDIR}/{s}/viral/checkv/quality_summary.tsv")
             for s in SAMPLES]
    if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:
        pairs += [(g, f"{OUTDIR}/coassembly/{g}/viral/checkv/quality_summary.tsv")
                  for g in GROUPS]
    return pairs


def _catalog_checkv_summaries(wildcards):
    return [path for _, path in _catalog_checkv_pairs()]


def _load_catalog_completeness():
    """CheckV completeness and quality tier, keyed by the pool's namespaced IDs.

    Duas passadas, e a segunda existe por um motivo especifico. O
    quality_summary.tsv do CheckV tem uma linha por contig ORIGINAL
    ("k141_219139"), mas a sequencia que chega ao pool pode ser um fragmento
    de provirus aparado ("k141_219139_1"). Chavear so pelo id original faria
    todo profago cair para completude 0.0 -- reprovando o braco de qualidade
    do portao, sumindo do nivel `mq` de representantes e nunca sendo anotado
    pelo pharokka.

    Ate 2026-08-19 isso nao aparecia porque o aparo do CheckV silenciosamente
    nunca acontecia (o mesmo bug de sufixo, uma camada abaixo -- ver
    scripts/checkv_provirus.py). Consertar aquela camada sem esta apenas
    moveria o descasamento de lugar.

    A segunda passada le o provenance.tsv, que ja lista todo member_id do pool
    com seu source, e copia a evidencia do contig de origem para o id do
    fragmento.
    """
    import csv
    import sys as _sys
    _sys.path.insert(0, SCRIPTS_DIR)
    from checkv_provirus import resolve_original_id

    completeness = {}
    quality = {}
    known_by_source = {}
    for source_id, path in _catalog_checkv_pairs():
        if not os.path.exists(path):
            continue
        known = known_by_source.setdefault(source_id, set())
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                contig = (row.get("contig_id") or "").strip()
                if not contig:
                    continue
                known.add(contig)
                key = f"{source_id}|{contig}"
                quality[key] = (row.get("checkv_quality") or "").strip()
                try:
                    completeness[key] = float(row.get("completeness") or 0)
                except ValueError:
                    completeness[key] = 0.0

    prov_path = f"{CATALOG_DIR}/provenance.tsv"
    if os.path.exists(prov_path):
        with open(prov_path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                member = (row.get("member_id") or "").strip()
                source = (row.get("source_id") or "").strip()
                header = (row.get("original_contig_id") or "").strip()
                if not member or member in completeness:
                    continue
                orig, ok = resolve_original_id(header,
                                               known_by_source.get(source, set()))
                if not ok:
                    continue
                src_key = f"{source}|{orig}"
                if src_key in completeness:
                    completeness[member] = completeness[src_key]
                    quality[member] = quality.get(src_key, "")

    return completeness, quality


def _catalog_genomad_pairs():
    """(source_id, genomad_done_path) for every source that can carry a
    GeNomad virus_summary, in the SAME source enumeration as
    _catalog_checkv_pairs()/_catalog_sources() -- samples then groups.

    Only the directory is needed (the actual *_virus_summary.tsv is found
    by glob, same as the removed per-sample viral_taxonomy rule did), but
    the *_done.txt path is what every genomad rule (sample-level `genomad`
    and its coassembly twin `coassembly_genomad`) actually declares as an
    output, so it is what downstream rules should declare as an input to
    track the dependency correctly.
    """
    pairs = [(s, f"{OUTDIR}/{s}/viral/genomad/done.txt") for s in SAMPLES]
    if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:
        pairs += [(g, f"{OUTDIR}/coassembly/{g}/viral/genomad/done.txt")
                  for g in GROUPS]
    return pairs


def _catalog_genomad_done_files(wildcards):
    return [path for _, path in _catalog_genomad_pairs()]


def _load_catalog_genomad():
    """Raw GeNomad virus_summary rows from every catalog source, keyed by
    the pool's namespaced IDs ({source_id}|seq_name) -- same prefixing
    convention as _load_catalog_completeness() for CheckV.

    Fixes the join that silently dropped GeNomad from viral_taxonomy: the
    old per-sample rule read ONE sample's genomad_dir and matched its bare
    seq_name against namespaced contig IDs from the global representative
    FASTA, which never matched (see docs/ROADMAP_SIMPLIFICACAO.md, "O bug
    de namespace de ID"). Aggregating and prefixing every source's GeNomad
    output here, once, makes the keys land in the same namespace as the
    representative FASTA headers votu_taxonomy iterates over.
    """
    import csv, glob
    genomad_rows = {}
    n_sources_with_output = 0
    for source_id, done_path in _catalog_genomad_pairs():
        gdir = os.path.dirname(done_path)
        gfiles = glob.glob(os.path.join(gdir, "**", "*_virus_summary.tsv"), recursive=True)
        if not gfiles:
            continue
        n_sources_with_output += 1
        with open(gfiles[0]) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                name = (row.get("seq_name") or "").strip()
                if not name:
                    continue
                genomad_rows[f"{source_id}|{name}"] = row
    return genomad_rows, n_sources_with_output


def _mmseqs_lca_rollup(hits_path, ranks):
    """Shared by votu_mmseqs_taxonomy (INPHARED) and
    votu_mmseqs_taxonomy_custom (e.g. IMG/VR) results -- both produce the
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


rule votu_catalog_pool:
    """Concatenate every viral non-redundant set into one namespaced pool.

    Contig IDs are only unique within an assembly, so they are prefixed with
    their source. Without this the pool silently merges unrelated contigs.
    """
    input:
        fastas = _catalog_input_fastas,
    output:
        pool       = f"{CATALOG_DIR}/pool.fasta",
        provenance = f"{CATALOG_DIR}/provenance.tsv",
    log:
        f"{OUTDIR}/logs/votu_catalog_pool.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_pool.tsv"
    run:
        stats = build_pool(_catalog_sources(),
                           str(output.pool), str(output.provenance))
        with open(str(log[0]), "w") as lf:
            lf.write(f"[votu_catalog_pool] sequences: {stats['n_sequences']}\n")
            lf.write(f"[votu_catalog_pool] sources used: {stats['n_sources']}\n")
            lf.write(f"[votu_catalog_pool] sources skipped (missing/empty): "
                     f"{stats['n_skipped']}\n")


rule votu_catalog_skani:
    """skani triangle over the pooled catalog.

    --sparse is REQUIRED, not an optimisation: the default dense matrix
    emits ANI only, with no aligned fraction, so the ICTV AF >= 85 criterion
    cannot be evaluated from it at all.
    """
    input:
        pool = rules.votu_catalog_pool.output.pool,
    output:
        ani = f"{CATALOG_DIR}/skani_ani.tsv",
    log:
        f"{OUTDIR}/logs/votu_catalog_skani.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_skani.tsv"
    conda: "../envs/env_derep.yaml"
    container: CONTAINERS.get("skani")
    threads: THREADS
    shell:
        """
        mkdir -p $(dirname {output.ani})
        N_SEQ=$(grep -c '^>' {input.pool} 2>/dev/null || echo 0)
        if [ "$N_SEQ" -eq 0 ]; then
            echo "[votu_catalog_skani] Empty pool" | tee {log}
            printf "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n" > {output.ani}
            exit 0
        fi
        echo "[votu_catalog_skani] $N_SEQ genomes in pool" | tee {log}
        skani triangle \
            -i {input.pool} \
            -o {output.ani} \
            -t {threads} \
            --slow \
            --sparse \
            >> {log} 2>&1
        """


rule votu_catalog_cluster:
    """Single-linkage vOTU clustering over the global pool.

    Fails loudly when nothing collapses -- N sequences producing N clusters
    is the exact signature of the format/parser mismatch this stage replaces.
    """
    input:
        ani    = rules.votu_catalog_skani.output.ani,
        pool   = rules.votu_catalog_pool.output.pool,
        checkv = _catalog_checkv_summaries,
    output:
        clusters = f"{CATALOG_DIR}/vOTU_clusters.tsv",
    log:
        f"{OUTDIR}/logs/votu_catalog_cluster.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_cluster.tsv"
    params:
        ani_min = VOTU_ANI,
        af_min  = VOTU_AF,
    run:
        import csv

        ids = []
        with open(str(input.pool)) as fh:
            for line in fh:
                if line.startswith(">"):
                    ids.append(line[1:].strip().split()[0])

        completeness, _quality = _load_catalog_completeness()

        edges = parse_skani_sparse(str(input.ani), params.ani_min,
                                   params.af_min, set(ids))
        clusters = cluster_votus(ids, edges, completeness)
        n_clusters = write_clusters(clusters, len(ids), str(output.clusters),
                                    completeness)

        with open(str(log[0]), "w") as lf:
            lf.write(f"[votu_catalog_cluster] genomes={len(ids)} "
                     f"edges={len(edges)} clusters={n_clusters} "
                     f"ani>={params.ani_min} af>={params.af_min}\n")
            if len(ids):
                lf.write(f"[votu_catalog_cluster] reduction="
                         f"{100 * (1 - n_clusters / len(ids)):.1f}%\n")


rule votu_catalog_reps:
    """Extract the two representative tiers used downstream.

    Same quality gates as the removed per-sample representative-extraction
    rule, applied once over the global catalog:
      all      — one per vOTU; recruitment reference and report base
      mq       — MQ+ (Complete/HQ/MQ or completeness >= 50%); taxonomy,
                 PHIST, annotation
    """
    input:
        pool     = rules.votu_catalog_pool.output.pool,
        clusters = rules.votu_catalog_cluster.output.clusters,
        checkv   = _catalog_checkv_summaries,
    output:
        all_fasta     = f"{CATALOG_DIR}/catalog_all_reps.fasta",
        mq_fasta      = f"{CATALOG_DIR}/catalog_mq_reps.fasta",
        done          = f"{CATALOG_DIR}/done.txt",
    log:
        f"{OUTDIR}/logs/votu_catalog_reps.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_reps.tsv"
    params:
        keep_tiers = VIRAL_KEEP_TIERS,
    run:
        import csv

        reps = set()
        with open(str(input.clusters)) as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                reps.add(row["representative"])

        completeness, quality = _load_catalog_completeness()

        def is_mq(rid):
            return (quality.get(rid, "") in params.keep_tiers
                    or completeness.get(rid, 0.0) >= 50.0)

        seqs = {}
        cur = None
        with open(str(input.pool)) as fh:
            for line in fh:
                if line.startswith(">"):
                    cur = line[1:].strip().split()[0]
                    if cur in reps:
                        seqs[cur] = []
                    else:
                        cur = None
                elif cur is not None:
                    seqs[cur].append(line.strip())

        n_all = n_mq = 0
        with open(str(output.all_fasta), "w") as fa, \
             open(str(output.mq_fasta), "w") as fm:
            for rid, chunks in seqs.items():
                seq = "".join(chunks)
                record = f">{rid}\n{seq}\n"
                fa.write(record)
                n_all += 1
                if is_mq(rid):
                    fm.write(record)
                    n_mq += 1

        with open(str(log[0]), "w") as lf:
            lf.write(f"[votu_catalog_reps] all: {n_all}\n")
            lf.write(f"[votu_catalog_reps] MQ+ (taxonomy/PHIST/annotation): {n_mq}\n")
        if n_all == 0:
            write_status(str(output.done), "skipped: empty catalog")
        else:
            write_status(str(output.done), "ok")


rule votu_prodigal:
    """Predict ORFs from the vOTU MQ+ representatives, once for the whole
    catalog.

    Moved from the per-sample `prodigal_viral` (rules/taxonomy.smk) on
    2026-08-18: its input was ALREADY exclusively the global representative
    FASTA (rules.votu_catalog_reps.output.mq_fasta), so every sample ran an
    identical prodigal pass over byte-identical input and produced
    byte-identical output -- see docs/ROADMAP_SIMPLIFICACAO.md "(h)".
    """
    input:
        viral = rules.votu_catalog_reps.output.mq_fasta,
    output:
        faa  = f"{CATALOG_DIR}/taxonomy/viral_proteins.faa",
        done = f"{CATALOG_DIR}/taxonomy/prodigal_done.txt",
    log:   f"{OUTDIR}/logs/votu_prodigal.log"
    benchmark: f"{OUTDIR}/benchmarks/votu_prodigal.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("prodigal")
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.faa})
        if [ ! -s {input.viral} ]; then
            touch {output.faa}
            echo "skipped: no viral contigs" > {output.done}; exit 0
        fi
        prodigal -i {input.viral} -a {output.faa} -p meta -f gff > {log} 2>&1
        echo "ok" > {output.done}
        """


rule votu_mmseqs_taxonomy:
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

    votu_taxonomy compares this against votu_mmseqs_taxonomy_custom/GeNomad
    by resolved rank depth, not a fixed priority order, so this can win,
    lose, or tie per contig.

    Moved from the per-sample `mmseqs_taxonomy_viral` (rules/taxonomy.smk)
    on 2026-08-18 -- its input (prodigal_viral's .faa) was already global,
    see rule votu_prodigal above.
    """
    input:
        faa  = rules.votu_prodigal.output.faa,
        done = rules.votu_prodigal.output.done,
    output:
        hits = f"{CATALOG_DIR}/taxonomy/mmseqs_vs_inphared.tsv",
        done = f"{CATALOG_DIR}/taxonomy/mmseqs_inphared_done.txt",
    log:   f"{OUTDIR}/logs/votu_mmseqs_taxonomy.log"
    benchmark: f"{OUTDIR}/benchmarks/votu_mmseqs_taxonomy.tsv"
    conda: "../envs/env_assembly.yaml"
    container:  CONTAINERS.get("mmseqs2")
    threads: THREADS
    params:
        seqtaxdb = f"{INPHARED_DB}/inphared_mmseqs_taxdb/seqTaxDB",
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
            write_status(output.done, "skipped: " + msg.split("] ", 1)[-1])

        if not os.path.exists(str(params.seqtaxdb) + ".dbtype"):
            write_empty(
                "[votu_mmseqs_taxonomy] No seqTaxDB at " + str(params.seqtaxdb) +
                " -- run scripts/prepare_mmseqs_taxdb.py --format inphared once first (see INSTALL.md). Skipping."
            )
            return

        if not os.path.exists(str(input.faa)) or os.path.getsize(str(input.faa)) == 0:
            write_empty("[votu_mmseqs_taxonomy] No viral proteins -- skipping")
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
        write_status(output.done, "ok")


rule votu_mmseqs_taxonomy_custom:
    """
    Optional custom viral MMseqs2 seqTaxDB (e.g. IMG/VR) -- real per-query
    LCA, same approach as votu_mmseqs_taxonomy but against a user-supplied
    DB instead of INPHARED. Replaces an earlier diamond_custom_viral
    best-hit+majority-vote rule entirely (removed 2026-06-23, see memory
    project_viral_taxonomy_merge).

    Same pattern as mmseqs_taxonomy_prok/votu_mmseqs_taxonomy: the seqTaxDB
    is NOT built here -- pre-build it once with scripts/prepare_mmseqs_taxdb.py
    (e.g. --format imgvr, see INSTALL.md). Skipped gracefully if
    custom_viral_mmseqs_db isn't configured, same as any other optional DB.

    Moved from the per-sample `mmseqs_taxonomy_custom_viral`
    (rules/taxonomy.smk) on 2026-08-18, same reasoning as votu_mmseqs_taxonomy.
    """
    input:
        faa  = rules.votu_prodigal.output.faa,
        done = rules.votu_prodigal.output.done,
    output:
        hits = f"{CATALOG_DIR}/taxonomy/mmseqs_vs_custom.tsv",
        done = f"{CATALOG_DIR}/taxonomy/mmseqs_custom_viral_done.txt",
    log:   f"{OUTDIR}/logs/votu_mmseqs_taxonomy_custom.log"
    benchmark: f"{OUTDIR}/benchmarks/votu_mmseqs_taxonomy_custom.tsv"
    conda: "../envs/env_assembly.yaml"
    container:  CONTAINERS.get("mmseqs2")
    threads: THREADS
    params:
        seqtaxdb = CUSTOM_VIRAL_MMSEQS_DB,
        outdir   = lambda wc, output: os.path.join(os.path.dirname(output.done), "mmseqs_custom"),
        querydb  = lambda wc, output: os.path.join(os.path.dirname(output.done), "mmseqs_custom", "queryDB"),
        result   = lambda wc, output: os.path.join(os.path.dirname(output.done), "mmseqs_custom", "result"),
        tmp      = lambda wc, output: os.path.join(os.path.dirname(output.done), "mmseqs_custom", "tmp"),
    run:
        import os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        header = "qseqid\ttaxid\trank\tname\tlineage\n"

        def write_empty(msg):
            with open(str(log[0]), "a") as lf:
                lf.write(msg + "\n")
            Path(str(output.hits)).write_text(header)
            write_status(output.done, "skipped: " + msg.split("] ", 1)[-1])

        if not params.seqtaxdb or not os.path.exists(str(params.seqtaxdb) + ".dbtype"):
            write_empty("[votu_mmseqs_taxonomy_custom] No custom_viral_mmseqs_db configured -- skipping")
            return

        if not os.path.exists(str(input.faa)) or os.path.getsize(str(input.faa)) == 0:
            write_empty("[votu_mmseqs_taxonomy_custom] No viral proteins -- skipping")
            return

        # mmseqs taxonomy refuses to run if its output DB already exists --
        # same lesson as votu_mmseqs_taxonomy/mmseqs_taxonomy_prok.
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
        write_status(output.done, "ok")


rule votu_taxonomy:
    """
    Merge taxonomy from all sources into one table per vOTU representative,
    once for the whole catalog.
    Sources: MMseqs2/INPHARED (LCA), MMseqs2/Custom (LCA), GeNomad.
    No fixed source priority -- each source proposes a (family, genus, order)
    call, and whichever resolves the DEEPEST rank wins per contig (genus >
    family > order > unclassified). Ties are broken by trust order:
    mmseqs_inphared > mmseqs_custom > genomad.
    Why: with a fixed priority, an upstream tier's *shallow* hit would always
    beat a downstream tier's deeper, better-supported call (e.g. GeNomad
    resolving to genus) -- see [[project_viral_taxonomy_merge]] memory for
    the full rationale and how each source's depth is computed. That same
    memory documents why the earlier Diamond/INPHARED (2026-06-22) and
    diamond_custom_viral (2026-06-23) best-hit+majority-vote tiers were
    removed in favour of real per-query LCA (mmseqs_inphared/mmseqs_custom).

    Moved from the per-sample `viral_taxonomy` (rules/taxonomy.smk) on
    2026-08-18, fixing a silent GeNomad dropout in the process: the old
    per-sample rule globbed ONE sample's genomad_dir for its
    *_virus_summary.tsv and matched `seq_name` (bare contig ID) against
    contigs read from the GLOBAL representative FASTA (namespaced
    "{source}|{contig}") -- never matched, so GeNomad silently contributed
    zero rows to every sample's merge (see docs/ROADMAP_SIMPLIFICACAO.md,
    "O bug de namespace de ID"). `_load_catalog_genomad()` (this file)
    aggregates GeNomad output from EVERY source (sample and group) and
    prefixes each `seq_name` with its own source_id, landing it in the same
    namespace as the representative FASTA headers this rule iterates over.

    The MMseqs2 tiers were never affected by that bug -- they already ran
    on the global prodigal .faa (votu_prodigal), so their contig keys were
    already namespaced correctly.

    Output: {CATALOG_DIR}/taxonomy/viral_taxonomy_merged.tsv, one row per
    vOTU representative (from votu_catalog_reps.output.mq_fasta), seq_name
    NAMESPACED ("{source_id}|{original_contig_id}"). Per-sample/per-group
    consumers read this through the `viral_taxonomy_view` rule
    (rules/taxonomy.smk), which filters to that source and strips the
    prefix back to a bare ID.
    """
    input:
        genomad_done = _catalog_genomad_done_files,
        mmseqs_hits  = rules.votu_mmseqs_taxonomy.output.hits,
        mmseqs_done  = rules.votu_mmseqs_taxonomy.output.done,
        custom_hits  = rules.votu_mmseqs_taxonomy_custom.output.hits,
        custom_done  = rules.votu_mmseqs_taxonomy_custom.output.done,
        viral        = rules.votu_catalog_reps.output.mq_fasta,
    output:
        tsv  = f"{CATALOG_DIR}/taxonomy/viral_taxonomy_merged.tsv",
        done = f"{CATALOG_DIR}/taxonomy/taxonomy_done.txt",
    log:   f"{OUTDIR}/logs/votu_taxonomy.log"
    benchmark: f"{OUTDIR}/benchmarks/votu_taxonomy.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("diamond")
    threads: 1
    run:
        import csv, os, collections
        from pathlib import Path

        lf = open(str(log[0]), "w")

        # All vOTU representatives (namespaced: "{source_id}|{contig}")
        contigs = []
        if os.path.exists(str(input.viral)):
            with open(str(input.viral)) as f:
                for line in f:
                    if line.startswith(">"): contigs.append(line[1:].split()[0])
        lf.write(f"Total vOTU representatives (namespaced seq_name): {len(contigs)}\n")

        # ── MMseqs2/INPHARED + MMseqs2/Custom (real per-query LCA) ──────
        # Both rules produce the same (qseqid, taxid, rank, name, lineage)
        # shape over the same 8-level ICTV rank scheme -- _mmseqs_lca_rollup
        # (module-level helper, top of this file) does the per-contig rollup
        # for either one. Keys here are already namespaced (both ran on the
        # global votu_prodigal .faa), so they match `contigs` directly.
        _RANKS = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus']
        mmseqs_tax = _mmseqs_lca_rollup(str(input.mmseqs_hits), _RANKS)
        lf.write(f"MMseqs2/INPHARED: {len(mmseqs_tax)} contigs (namespaced keys)\n")

        custom_tax = _mmseqs_lca_rollup(str(input.custom_hits), _RANKS)
        lf.write(f"MMseqs2/Custom: {len(custom_tax)} contigs (namespaced keys)\n")

        # ── GeNomad ───────────────────────────────────────────────────
        # Aggregated across every catalog source (sample + group) and
        # namespaced by _load_catalog_genomad() -- see docstring above for
        # why this replaces a single-sample glob.
        genomad_raw, n_genomad_sources = _load_catalog_genomad()
        lf.write(f"GeNomad: {n_genomad_sources}/{len(_catalog_genomad_pairs())} "
                 f"sources had a virus_summary; {len(genomad_raw)} raw rows "
                 f"(namespaced keys)\n")

        genomad_tax = {}
        for name, row in genomad_raw.items():
            tax   = row.get("taxonomy","")
            score = row.get("virus_score","0")
            if not tax: continue
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
        lf.write(f"GeNomad parsed: {len(genomad_tax)} contigs (namespaced keys)\n")
        n_genomad_matched = sum(1 for c in contigs if c in genomad_tax)
        lf.write(f"GeNomad matched against vOTU representatives: "
                 f"{n_genomad_matched}/{len(contigs)}\n")

        # ── Build final table ─────────────────────────────────────────
        # No fixed tier priority: each source proposes (family, genus, order);
        # whichever resolves the deepest rank wins. Ties broken by trust order
        # (mmseqs_inphared > mmseqs_custom > genomad).
        _PRIORITY = {"mmseqs_inphared": 1, "mmseqs_custom": 2, "genomad": 3}

        def _depth(genus, family, order):
            if genus:  return 3
            if family: return 2
            if order:  return 1
            return 0

        rows = []; stats = collections.Counter()
        for contig in contigs:
            mms = mmseqs_tax.get(contig, {})
            gmd = genomad_tax.get(contig, {})
            cms = custom_tax.get(contig, {})

            candidates = []  # (source, ff, fg, fo, lin, conf, best)


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
                  "genomad_best","genomad_class","genomad_score",
                  "mmseqs_rank","mmseqs_lineage","mmseqs_n_proteins",
                  "custom_rank","custom_lineage","custom_n_proteins"]
        with open(str(output.tsv), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader(); w.writerows(rows)
        Path(str(output.done)).write_text("ok\n")


rule bacphlip_votu:
    """
    BACPHLIP lifestyle (lytic/lysogenic) call, once per vOTU representative.

    Default gate is Complete/High-quality only, deliberately conservative:
    BACPHLIP decides lysogeny from the presence of integrase/recombinase/
    repressor domains, and an incomplete genome can simply be missing the
    fragment that would have carried one. The error is NOT symmetric —
    fragmentation biases calls towards "lytic" and inflates the apparent
    lytic fraction — so genomes below `bacphlip_min_quality` are excluded
    rather than called anyway.
    """
    input:
        fasta  = rules.votu_catalog_reps.output.all_fasta,
        checkv = _catalog_checkv_summaries,
    output:
        lifestyle = f"{CATALOG_DIR}/bacphlip/votu_lifestyle.tsv",
        done      = f"{CATALOG_DIR}/bacphlip/done.txt",
    log:
        f"{OUTDIR}/logs/bacphlip_votu.log"
    benchmark:
        f"{OUTDIR}/benchmarks/bacphlip_votu.tsv"
    conda: "../envs/env_reads_classify.yaml"
    container: CONTAINERS.get("bacphlip")
    threads: 1
    params:
        min_quality         = BACPHLIP_MIN_QUALITY,
        virulence_threshold = BACPHLIP_VIRULENCE_THRESHOLD,
        outdir              = lambda wc, output: os.path.dirname(output.lifestyle),
    run:
        import subprocess

        HEADER = ["votu_id", "lifestyle", "virulent_score", "temperate_score",
                  "checkv_quality", "quality_tier_used"]
        os.makedirs(params.outdir, exist_ok=True)

        def write_header_only():
            with open(str(output.lifestyle), "w") as fh:
                fh.write("\t".join(HEADER) + "\n")

        completeness, quality = _load_catalog_completeness()

        tmp_fasta = os.path.join(params.outdir, "votu_input.fasta")
        kept = {}
        cur_id = None
        cur_keep = False
        n_kept = 0
        with open(str(input.fasta)) as fin, open(tmp_fasta, "w") as fout:
            for line in fin:
                if line.startswith(">"):
                    cur_id = line[1:].strip().split()[0]
                    cur_keep = quality.get(cur_id, "") in params.min_quality
                    if cur_keep:
                        kept[cur_id] = quality.get(cur_id, "")
                        n_kept += 1
                        fout.write(line)
                elif cur_keep:
                    fout.write(line)

        with open(str(log[0]), "w") as lf:
            lf.write(f"[bacphlip_votu] quality tiers accepted: {params.min_quality}\n")
            lf.write(f"[bacphlip_votu] vOTU representatives kept: {n_kept}\n")

        if n_kept == 0:
            write_header_only()
            write_status(str(output.done), "skipped: no vOTU at required CheckV quality")
            return

        rc = subprocess.run(
            ["bacphlip", "-i", tmp_fasta, "--multi_fasta"],
            stdout=open(str(log[0]), "a"), stderr=subprocess.STDOUT,
        ).returncode

        bacphlip_out = tmp_fasta + ".bacphlip"
        if rc != 0 or not os.path.exists(bacphlip_out):
            write_header_only()
            write_status(str(output.done), f"failed: bacphlip exit {rc}")
            return

        import csv as _csv
        with open(bacphlip_out) as fh:
            rdr = _csv.reader(fh, delimiter="\t")
            bp_header = next(rdr, None)
            rows = list(rdr)

        with open(str(output.lifestyle), "w") as fh:
            fh.write("\t".join(HEADER) + "\n")
            for row in rows:
                if len(row) < 3:
                    continue
                votu_id = row[0]
                try:
                    virulent_score  = float(row[1])
                    temperate_score = float(row[2])
                except ValueError:
                    continue
                lifestyle = "lytic" if virulent_score >= params.virulence_threshold else "lysogenic"
                fh.write("\t".join([
                    votu_id, lifestyle, f"{virulent_score:.6g}", f"{temperate_score:.6g}",
                    kept.get(votu_id, ""), "/".join(params.min_quality),
                ]) + "\n")

        write_status(str(output.done), "ok")


rule eggnog_viral:
    """
    eggNOG-mapper over ORF predictions from the MQ+ vOTU representatives,
    once for the whole catalog (prodigal is cheap, eggNOG is the expensive
    step, so this runs one pass instead of per-sample).

    Downstream filtering (filter_putative_amgs.py) flags candidates only —
    a KEGG_ko hit landing on a metabolism KEGG pathway map. These are
    "putative AMGs": AMG calling from annotation alone is known to be
    error-prone and requires genomic-context inspection to confirm, so the
    pipeline never claims a confirmed AMG anywhere.
    """
    input:
        fasta = rules.votu_catalog_reps.output.mq_fasta,
    output:
        annot = f"{CATALOG_DIR}/eggnog_viral/eggnog_annotations.tsv",
        amg   = f"{CATALOG_DIR}/eggnog_viral/putative_amgs.tsv",
        done  = f"{CATALOG_DIR}/eggnog_viral/done.txt",
    log:
        f"{OUTDIR}/logs/eggnog_viral.log"
    benchmark:
        f"{OUTDIR}/benchmarks/eggnog_viral.tsv"
    conda: "../envs/env_annotation.yaml"
    container: CONTAINERS.get("eggnog_mapper")
    threads: THREADS
    params:
        outdir   = lambda wc, output: os.path.dirname(output.annot),
        proteins = lambda wc, output: os.path.join(os.path.dirname(output.annot), "votu_proteins.faa"),
        amg_script = os.path.join(SCRIPTS_DIR, "filter_putative_amgs.py"),
    shell:
        """
        mkdir -p {params.outdir}

        AMG_HEADER="votu_id\\tprotein_id\\tKEGG_ko\\tKEGG_Pathway\\tCOG_category\\tDescription"

        if [ "{USE_EGGNOG_VIRAL}" != "True" ]; then
            echo "[eggnog_viral] use_eggnog_viral: false — skipping" | tee {log}
            touch {output.annot}
            printf "$AMG_HEADER\\n" > {output.amg}
            echo "skipped: use_eggnog_viral is false" > {output.done}
            exit 0
        fi

        if [ -z "{EGGNOG_DB}" ] || [ ! -d "{EGGNOG_DB}" ]; then
            echo "[eggnog_viral] EGGNOG_DB not configured — skipping" | tee {log}
            touch {output.annot}
            printf "$AMG_HEADER\\n" > {output.amg}
            echo "skipped: EGGNOG_DB not configured" > {output.done}
            exit 0
        fi

        if [ ! -s {input.fasta} ]; then
            echo "[eggnog_viral] empty input FASTA (no MQ+ vOTU representatives) — skipping" | tee -a {log}
            touch {output.annot}
            printf "$AMG_HEADER\\n" > {output.amg}
            echo "skipped: empty MQ+ representative FASTA" > {output.done}
            exit 0
        fi

        prodigal -i {input.fasta} -a {params.proteins} -p meta >> {log} 2>&1 && RC=0 || RC=$?
        if [ $RC -ne 0 ]; then
            echo "[eggnog_viral] prodigal failed (exit $RC)" | tee -a {log}
            touch {output.annot}
            printf "$AMG_HEADER\\n" > {output.amg}
            echo "failed: prodigal exit $RC" > {output.done}
            exit 0
        fi

        emapper.py \
            -m diamond \
            --itype proteins \
            -i {params.proteins} \
            -o eggnog_viral \
            --output_dir {params.outdir} \
            --cpu {threads} \
            --data_dir {EGGNOG_DB} \
            --override \
            >> {log} 2>&1 && RC=0 || RC=$?

        if [ $RC -ne 0 ] || [ ! -f {params.outdir}/eggnog_viral.emapper.annotations ]; then
            echo "[eggnog_viral] emapper.py failed (exit $RC)" | tee -a {log}
            touch {output.annot}
            printf "$AMG_HEADER\\n" > {output.amg}
            echo "failed: emapper.py exit $RC" > {output.done}
            exit 0
        fi

        cp {params.outdir}/eggnog_viral.emapper.annotations {output.annot}

        python3 {params.amg_script} {output.annot} {output.amg} && RC=0 || RC=$?
        if [ $RC -ne 0 ]; then
            echo "[eggnog_viral] filter_putative_amgs.py failed (exit $RC)" | tee -a {log}
            printf "$AMG_HEADER\\n" > {output.amg}
            echo "failed: filter_putative_amgs.py exit $RC" > {output.done}
            exit 0
        fi

        echo "ok" > {output.done}
        """


# ══════════════════════════════════════════════════════════════════════
# Read recruitment against the catalog + presence/abundance matrices
#
# Two-branch pattern. Both variants produce the exact same temp SAM path
# (`{CATALOG_DIR}/mapping/{sample}.catalog.sam`), so votu_catalog_sort,
# votu_catalog_coverm and votu_catalog_matrices run without knowing which
# technology was used -- same strategy as rules/mapping.smk with
# bwa_mem / minimap2_lr.
#
# Alignment and sorting are separate rules, mandatorily. The containers
# are distinct images (bwa-mem2, minimap2, samtools in containers.yaml),
# so a `| samtools sort` inside the aligner's rule would fail because
# samtools does not exist in that image. It would work under --use-conda
# (env_mapping.yaml carries all three), which is exactly why the defect
# would go unnoticed until the first container run. rules/mapping.smk
# already splits for this reason; the catalog does the same.
# ══════════════════════════════════════════════════════════════════════

if LONG_READS:

    rule votu_catalog_map:
        """Long-read recruitment of one sample against the whole catalog.

        minimap2 indexes the reference on the fly, so there is no separate
        index rule on this branch. ONT: -ax map-ont; HiFi: -ax map-hifi.
        """
        input:
            reads = _clean_lr,
            fasta = rules.votu_catalog_reps.output.all_fasta,
        output:
            sam = temp(f"{CATALOG_DIR}/mapping/{{sample}}.catalog.sam"),
        log:
            f"{OUTDIR}/{{sample}}/logs/votu_catalog_map.log"
        benchmark:
            f"{OUTDIR}/{{sample}}/benchmarks/votu_catalog_map.tsv"
        conda: "../envs/env_mapping.yaml"
        container: CONTAINERS.get("minimap2")
        threads: THREADS
        shell:
            """
            mkdir -p $(dirname {output.sam})
            if [ ! -s {input.fasta} ]; then
                echo "[votu_catalog_map] Empty catalog -- header-only SAM" | tee {log}
                printf "@HD\tVN:1.6\tSO:unsorted\n" > {output.sam}
                exit 0
            fi
            if [ "{LR_TECH}" = "hifi" ]; then
                PRESET="map-hifi"
            else
                PRESET="map-ont"
            fi
            minimap2 -ax $PRESET \
                -t {threads} \
                {input.fasta} {input.reads} \
                > {output.sam} 2> {log}
            """

else:

    rule votu_catalog_index:
        """Index the catalog once; every sample maps against the same reference."""
        input:
            fasta = rules.votu_catalog_reps.output.all_fasta,
        output:
            idx = f"{CATALOG_DIR}/mapping/catalog_index.bwt.2bit.64",
        log:
            f"{OUTDIR}/logs/votu_catalog_index.log"
        conda: "../envs/env_mapping.yaml"
        container: CONTAINERS.get("bwa_mem2")
        params:
            prefix = f"{CATALOG_DIR}/mapping/catalog_index",
        shell:
            """
            mkdir -p $(dirname {output.idx})
            if [ ! -s {input.fasta} ]; then
                echo "[votu_catalog_index] Empty catalog -- skipping index" | tee {log}
                touch {output.idx}
                exit 0
            fi
            bwa-mem2 index -p {params.prefix} {input.fasta} > {log} 2>&1
            """

    rule votu_catalog_map:
        """Competitive short-read recruitment of one sample against the catalog.

        This is what makes per-sample presence assembly-independent: a virus
        present but too low-coverage to assemble in a sample is still detected
        here, which assembly-derived presence cannot do.
        """
        input:
            tr1   = _clean_r1,
            tr2   = _clean_r2,
            idx   = rules.votu_catalog_index.output.idx,
            fasta = rules.votu_catalog_reps.output.all_fasta,
        output:
            sam = temp(f"{CATALOG_DIR}/mapping/{{sample}}.catalog.sam"),
        log:
            f"{OUTDIR}/{{sample}}/logs/votu_catalog_map.log"
        benchmark:
            f"{OUTDIR}/{{sample}}/benchmarks/votu_catalog_map.tsv"
        conda: "../envs/env_mapping.yaml"
        container: CONTAINERS.get("bwa_mem2")
        threads: THREADS
        params:
            prefix     = f"{CATALOG_DIR}/mapping/catalog_index",
            single_end = SINGLE_END,
        shell:
            """
            mkdir -p $(dirname {output.sam})
            if [ ! -s {input.fasta} ]; then
                echo "[votu_catalog_map] Empty catalog -- header-only SAM" | tee {log}
                printf "@HD\tVN:1.6\tSO:unsorted\n" > {output.sam}
                exit 0
            fi
            if [ "{params.single_end}" = "True" ]; then
                bwa-mem2 mem -t {threads} {params.prefix} {input.tr1} \
                    > {output.sam} 2> {log}
            else
                bwa-mem2 mem -t {threads} {params.prefix} {input.tr1} {input.tr2} \
                    > {output.sam} 2> {log}
            fi
            """


rule votu_catalog_sort:
    """Sort and index the catalog alignment. Separate rule because samtools
    lives in its own container image, not in the aligners'."""
    input:
        sam = rules.votu_catalog_map.output.sam,
    output:
        bam = f"{CATALOG_DIR}/mapping/{{sample}}.catalog.sorted.bam",
        bai = f"{CATALOG_DIR}/mapping/{{sample}}.catalog.sorted.bam.bai",
    log:
        f"{OUTDIR}/{{sample}}/logs/votu_catalog_sort.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/votu_catalog_sort.tsv"
    conda: "../envs/env_mapping.yaml"
    container: CONTAINERS.get("samtools")
    threads: THREADS
    shell:
        """
        samtools sort -@ {threads} -o {output.bam} {input.sam} 2> {log}
        samtools index {output.bam} 2>> {log}
        """


rule votu_catalog_coverm:
    """CoverM over the catalog BAM, with the same filters as coverm_viral."""
    input:
        bam   = rules.votu_catalog_sort.output.bam,
        fasta = rules.votu_catalog_reps.output.all_fasta,
    output:
        tsv = f"{CATALOG_DIR}/coverm/{{sample}}.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/votu_catalog_coverm.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/votu_catalog_coverm.tsv"
    conda: "../envs/env_coverm.yaml"
    container: CONTAINERS.get("coverm")
    threads: THREADS
    params:
        method   = COVERM_METHOD,
        min_id   = VOTU_RECRUIT_MIN_ID,
        long_reads = LONG_READS,
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        if [ ! -s {input.fasta} ]; then
            echo "[votu_catalog_coverm] Empty catalog" | tee {log}
            printf "Contig\t{params.method}\tcovered_fraction\n" > {output.tsv}
            exit 0
        fi
        # --min-read-aligned-length is a short-read filter: a 45 bp floor is
        # meaningless for reads averaging kilobases, and on ONT it would only
        # discard genuinely short alignments that carry real signal.
        if [ "{params.long_reads}" = "True" ]; then
            LEN_FILTER=""
        else
            LEN_FILTER="--min-read-aligned-length 45"
        fi
        echo "[votu_catalog_coverm] min identity: {params.min_id}%" | tee {log}
        coverm contig \
            --bam-files {input.bam} \
            --min-read-percent-identity {params.min_id} \
            $LEN_FILTER \
            --contig-end-exclusion 75 \
            --methods {params.method} covered_fraction \
            --threads {threads} \
            --output-file {output.tsv} \
            >> {log} 2>&1
        """


rule votu_catalog_matrices:
    """Build the vOTU x sample presence and abundance matrices.

    Two independent presence signals, reported side by side and never merged:
      assembled — the vOTU has a member contig from that sample (provenance)
      recruited — >= votu_presence_min_coverage of the representative is
                  covered by that sample's reads (Roux et al. 2017)
    """
    input:
        clusters   = rules.votu_catalog_cluster.output.clusters,
        provenance = rules.votu_catalog_pool.output.provenance,
        coverm     = expand(f"{CATALOG_DIR}/coverm/{{sample}}.tsv", sample=SAMPLES),
    output:
        presence  = f"{CATALOG_DIR}/presence_matrix.tsv",
        abundance = f"{CATALOG_DIR}/votu_abundance_matrix.tsv",
        done      = f"{CATALOG_DIR}/matrices_done.txt",
    log:
        f"{OUTDIR}/logs/votu_catalog_matrices.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_matrices.tsv"
    params:
        min_cov = VOTU_PRESENCE_MIN_COV,
        method  = COVERM_METHOD,
    run:
        import csv
        from collections import defaultdict

        member_to_votu = {}
        votu_rep = {}
        with open(str(input.clusters)) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                member_to_votu[row["member"]] = row["votu_id"]
                votu_rep[row["votu_id"]] = row["representative"]

        # assembled presence: straight from provenance x clusters, free.
        assembled = defaultdict(set)
        with open(str(input.provenance)) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                votu = member_to_votu.get(row["member_id"])
                if votu and row["source_type"] == "sample":
                    assembled[votu].add(row["source_id"])

        rep_to_votu = {rep: votu for votu, rep in votu_rep.items()}
        recruited = defaultdict(set)
        abundance = defaultdict(dict)
        for sample, path in zip(SAMPLES, list(input.coverm)):
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                rdr = csv.reader(fh, delimiter="\t")
                header = next(rdr, None)
                if not header:
                    continue
                # CoverM names columns "<bam> <method>"; positions are stable:
                # 0 = Contig, 1 = configured method, 2 = covered_fraction.
                for parts in rdr:
                    if len(parts) < 3:
                        continue
                    votu = rep_to_votu.get(parts[0])
                    if not votu:
                        continue
                    try:
                        value = float(parts[1])
                        covered = float(parts[2])
                    except ValueError:
                        continue
                    # CoverM reports covered_fraction in [0,1]; the config
                    # threshold is a percentage.
                    if covered * 100.0 >= params.min_cov:
                        recruited[votu].add(sample)
                        abundance[votu][sample] = value

        votus = sorted(votu_rep)
        with open(str(output.presence), "w") as fh:
            fh.write("votu_id\trepresentative\t" + "\t".join(SAMPLES) + "\n")
            for votu in votus:
                cells = []
                for s in SAMPLES:
                    a = s in assembled[votu]
                    r = s in recruited[votu]
                    cells.append("both" if (a and r) else
                                 "assembled" if a else
                                 "recruited" if r else "absent")
                fh.write(f"{votu}\t{votu_rep[votu]}\t" + "\t".join(cells) + "\n")

        with open(str(output.abundance), "w") as fh:
            fh.write("votu_id\t" + "\t".join(SAMPLES) + "\n")
            for votu in votus:
                vals = [f"{abundance[votu].get(s, 0.0):.6g}" for s in SAMPLES]
                fh.write(f"{votu}\t" + "\t".join(vals) + "\n")

        with open(str(log[0]), "w") as lf:
            lf.write(f"[votu_catalog_matrices] vOTUs: {len(votus)}\n")
            lf.write(f"[votu_catalog_matrices] metric: {params.method}\n")
            lf.write(f"[votu_catalog_matrices] presence cutoff: "
                     f"{params.min_cov}% of representative covered\n")
            for s in SAMPLES:
                n_a = sum(1 for v in votus if s in assembled[v])
                n_r = sum(1 for v in votus if s in recruited[v])
                lf.write(f"  {s}: assembled={n_a} recruited={n_r}\n")
        write_status(str(output.done), "ok")




rule votu_pharokka:
    """
    Bacteriophage genome annotation with Pharokka (PHROGS database), once
    for the whole vOTU catalog. Moved from the per-sample `pharokka`
    (rules/annotation.smk) on 2026-08-18 -- see module docstring above for
    the namespace bug this fixes.

    Selection criteria (unchanged from the per-sample rule):
      - completeness >= PHAROKKA_MIN_COMPLETENESS  (quality tier NOT required —
        novel phages often receive "Not-determined" from CheckV despite high
        completeness because no close reference cluster exists)
      - No cap on genome count: every qualifying vOTU representative is
        annotated.

    Completeness gate now uses `_load_catalog_completeness()` (namespaced,
    keyed exactly like `input.viral_nr`'s headers) instead of any single
    sample's CheckV summary -- this IS the fix, not just the move.

    Output GBK/TSV are later used by votu_phold.
    Skipped if PHAROKKA_DB is not configured (empty string).
    """
    input:
        viral_nr = rules.votu_catalog_reps.output.mq_fasta,
        checkv   = _catalog_checkv_summaries,
    output:
        done = f"{CATALOG_DIR}/annotation/pharokka/done.txt",
        gbk  = f"{CATALOG_DIR}/annotation/pharokka/pharokka.gbk",
        tsv  = f"{CATALOG_DIR}/annotation/pharokka/pharokka_cds_final_merged_output.tsv",
    log:
        f"{OUTDIR}/logs/votu_pharokka.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_pharokka.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("pharokka")
    threads: THREADS
    params:
        outdir   = lambda wc, output: os.path.dirname(output.done),
        db       = PHAROKKA_DB,
        min_comp = PHAROKKA_MIN_COMPLETENESS,
        # NOT inside outdir: pharokka.py --force deletes/recreates its own
        # -o directory on startup, which would delete this -i input too if
        # it lived underneath it.
        hq_fa    = lambda wc, output: os.path.join(
            os.path.dirname(os.path.dirname(output.done)), "pharokka_hq_phages.fasta"),
    run:
        import os
        from pathlib import Path

        os.makedirs(params.outdir, exist_ok=True)
        log_path = str(log[0])

        def touch_empty(path):
            Path(path).touch()

        if not params.db or not os.path.isdir(str(params.db)):
            with open(log_path, "w") as lf:
                lf.write("[votu_pharokka] PHAROKKA_DB not configured — skipping\n")
            write_status(output.done, "skipped: PHAROKKA_DB not configured")
            touch_empty(output.gbk)
            touch_empty(output.tsv)
            return

        # Filter by the CATALOG's own completeness (namespaced), not a
        # sample's CheckV -- see rule docstring for why. Quality tier is
        # NOT required: novel phages (Caudovirales, etc.) often receive
        # "Not-determined" quality from CheckV even at high completeness,
        # because CheckV cannot assign a tier without a close reference
        # cluster; filtering by tier would silently drop bona-fide phages.
        completeness, _quality = _load_catalog_completeness()
        hq_set = {cid for cid, comp in completeness.items()
                  if comp >= float(params.min_comp)}

        # Extract sequences
        with open(str(params.hq_fa), "w") as out_fa, \
             open(str(input.viral_nr)) as in_fa, \
             open(log_path, "w") as lf:
            lf.write(f"[votu_pharokka] catalog namespaced IDs with completeness: "
                     f"{len(completeness)}\n")
            lf.write(f"[votu_pharokka] vOTU representatives selected "
                     f"(completeness >= {params.min_comp}%): {len(hq_set)}\n")
            write = False
            n_written = 0
            for line in in_fa:
                if line.startswith(">"):
                    name = line[1:].split()[0]
                    write = name in hq_set
                    if write:
                        n_written += 1
                if write:
                    out_fa.write(line)
            lf.write(f"[votu_pharokka] sequences extracted from mq_fasta: {n_written}\n")

        if not hq_set or os.path.getsize(str(params.hq_fa)) == 0:
            with open(log_path, "a") as lf:
                lf.write("[votu_pharokka] No HQ phages found — skipping\n")
            write_status(output.done,
                         "skipped: no vOTU representative reached "
                         "%s%% completeness" % params.min_comp)
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
            lf.write(f"[votu_pharokka] {n_seqs} sequence(s) in input — "
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

        write_status(output.done, "ok")


rule votu_phold:
    """
    Phold — structure-based annotation of viral hypothetical proteins, once
    for the whole vOTU catalog. Moved from the per-sample `phold`
    (rules/annotation.smk) on 2026-08-18; pure sequence-dependent
    annotation of votu_pharokka's GBK, no per-sample input at all.
    """
    input:
        pharokka_done = rules.votu_pharokka.output.done,
        pharokka_gbk  = rules.votu_pharokka.output.gbk,
    output:
        done = f"{CATALOG_DIR}/annotation/phold/done.txt",
        gbk  = f"{CATALOG_DIR}/annotation/phold/phold.gbk",
    log:
        f"{OUTDIR}/logs/votu_phold.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_phold.tsv"
    conda: "../envs/env_annotation.yaml"
    container:  CONTAINERS.get("phold")
    threads: THREADS
    params:
        outdir = lambda wc, output: os.path.dirname(output.done),
        db     = PHOLD_DB,
    shell:
        """
        mkdir -p {params.outdir}
        if [ ! -s {input.pharokka_gbk} ]; then
            echo "[votu_phold] Pharokka GBK empty — skipping" | tee {log}
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
            > {log} 2>&1 && RC=0 || RC=$?
        [ -f {params.outdir}/phold.gbk ] && \
            cp {params.outdir}/phold.gbk {output.gbk} || \
            touch {output.gbk}
        if [ "$RC" -ne 0 ]; then
            echo "failed: phold exit $RC" > {output.done}
        else
            echo "ok" > {output.done}
        fi
        """


rule votu_defensefinder_viral:
    """
    Anti-defense systems on VIRAL proteins (Han et al. 2026, Nat Commun
    cold seep defensome paper, Fig 6b/6c) -- DefenseFinder's
    --antidefensefinder side of the SAME tool/models/container as
    `defensefinder` (rules/defense_amr.smk, host bins), just pointed at the
    whole vOTU catalog's viral ORFs instead of bin proteins. Today
    DefenseFinder/AntiDefenseFinder only runs on the host (bin) side, so
    every Host<->Virus defense/anti-defense cross-link is one-sided (we
    know which bins defend, never which phages counter-defend).

    Moved from the per-sample `defensefinder_viral` (rules/defense_amr.smk)
    and its per-group twin `coassembly_defensefinder_viral`
    (rules/coassembly.smk) on 2026-08-18 (second half of "(h)",
    docs/ROADMAP_SIMPLIFICACAO.md). Both already read the GLOBAL prodigal
    output (rules.votu_prodigal, this file) since the first half of "(h)"
    earlier the same day -- input was already byte-identical across every
    sample/group, so the per-sample/per-group fan-out (N+G jobs) was pure
    waste on top of that. Worse: each per-group run still processed the
    WHOLE catalog (same global .faa), so
    coassembly/{group}/viral/defensefinder/viral_defense_systems.tsv
    silently held every source's systems, not just that group's, while
    rules/coassembly.smk's finalize step copied it into
    final/viral/defense_amr/ as if it were group-scoped -- a silent
    meaning change, not just wasted compute. One global run fixes both.
    """
    input:
        faa  = rules.votu_prodigal.output.faa,
        done = rules.votu_prodigal.output.done,
    output:
        done        = f"{CATALOG_DIR}/defensefinder/done.txt",
        systems     = f"{CATALOG_DIR}/defensefinder/viral_defense_systems.tsv",
        antisystems = f"{CATALOG_DIR}/defensefinder/viral_antidefense_systems.tsv",
    log:
        f"{OUTDIR}/logs/votu_defensefinder_viral.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_defensefinder_viral.tsv"
    conda: "../envs/env_defense.yaml"
    container:  CONTAINERS.get("defense_finder")
    threads: THREADS
    params:
        # No {{sample}}/{{group}} wildcard any more -- this runs once for
        # the whole catalog, so these no longer need to be derived from
        # `output` the way the removed per-sample rule did.
        outdir     = f"{CATALOG_DIR}/defensefinder",
        unit_label = "votu_catalog",
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
            write_status(output.done, "skipped: " + msg.split("] ", 1)[-1])

        if (not params.enabled or not os.path.exists(str(input.faa))
                or os.path.getsize(str(input.faa)) == 0):
            write_empty("[votu_defensefinder_viral] Disabled or no viral proteins -- skipping")
            return

        models_dir = params.models_dir or os.path.join(params.outdir, "models")
        os.makedirs(models_dir, exist_ok=True)
        if os.listdir(models_dir):
            with open(str(log[0]), "a") as lf:
                lf.write(f"[votu_defensefinder_viral] Models already cached in {models_dir}\n")
        else:
            shell("defense-finder update --models-dir {models_dir} >> {log} 2>&1 || "
                  "echo '[votu_defensefinder_viral] WARNING: model update failed' >> {log}")

        run_out = os.path.join(params.outdir, "run")
        os.makedirs(run_out, exist_ok=True)
        # O `|| echo WARNING` que estava aqui engolia a falha: as tabelas
        # saiam vazias e o done.txt vinha tocado em branco, ou seja, um
        # DefenseFinder quebrado era lido como "nenhum sistema anti-defesa
        # neste virome" -- exatamente o modo de falha silenciosa que o
        # write_status existe para impedir.
        _df_error = None
        try:
            shell(
                "defense-finder run -o {run_out} --models-dir {models_dir} "
                "--antidefensefinder {input.faa} >> {log} 2>&1"
            )
        except Exception as exc:
            _df_error = exc
            with open(str(log[0]), "a") as lf:
                lf.write("[votu_defensefinder_viral] WARNING: run failed: %s\n" % exc)

        # Same split-by-"activity" logic as the host-side rule (3.0.0
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
                    # "genome" column here is the catalog-wide protein set
                    # label (one DefenseFinder call across the whole
                    # catalog), not a per-bin name -- the protein_in_syst
                    # column (already part of `header`) carries the actual
                    # viral contig/ORF IDs needed to attribute hits
                    # downstream, and those IDs are namespaced
                    # "{source_id}|{contig}" since votu_prodigal ran on the
                    # namespaced catalog FASTA.
                    (anti_rows if is_anti else def_rows).append([params.unit_label] + row)

        def write(path, rows):
            with open(path, "w", newline="") as f:
                w = csv.writer(f, delimiter="\t")
                w.writerow(["genome"] + (header or []))
                w.writerows(rows)

        write(str(output.systems), def_rows)
        write(str(output.antisystems), anti_rows)

        with open(str(log[0]), "a") as lf:
            lf.write(f"[votu_defensefinder_viral] Done -- {len(def_rows)} defense, "
                      f"{len(anti_rows)} antidefense system rows\n")
        if _df_error is not None:
            write_status(output.done, "failed: defense-finder run %s" % _df_error)
        else:
            write_status(output.done, "ok")


rule votu_dbapis_viral:
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
    merged with votu_defensefinder_viral (same "never merge tiers" rule as
    AMR curated/exploratory, see rules/defense_amr.smk module docstring).
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

    Moved from the per-sample `dbapis_viral` (rules/defense_amr.smk) and its
    per-group twin `coassembly_dbapis_viral` (rules/coassembly.smk) on
    2026-08-18 (second half of "(h)", docs/ROADMAP_SIMPLIFICACAO.md) -- same
    move/rationale as votu_defensefinder_viral above (both already read the
    GLOBAL rules.votu_prodigal output, so this eliminates the N+G redundant
    runs AND the silent per-group meaning change).
    """
    input:
        faa  = rules.votu_prodigal.output.faa,
        done = rules.votu_prodigal.output.done,
    output:
        done = f"{CATALOG_DIR}/dbapis/done.txt",
        hits = f"{CATALOG_DIR}/dbapis/dbapis_hits.tsv",
    log:
        f"{OUTDIR}/logs/votu_dbapis_viral.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_dbapis_viral.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("diamond")
    threads: THREADS
    params:
        # No {{sample}}/{{group}} wildcard any more -- runs once globally.
        outdir   = f"{CATALOG_DIR}/dbapis",
        apis_dir = APIS_DB or f"{OUTDIR}/dbapis_db",
        enabled  = DEFENSE_AMR_VIRAL_ENABLED,
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}
        if [ "{params.enabled}" != "True" ] || [ ! -s {input.faa} ]; then
            echo "[votu_dbapis_viral] Disabled or no viral proteins -- skipping" | tee {log}
            printf "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\n" > {output.hits}
            echo "skipped: disabled or no viral proteins" > {output.done}
            exit 0
        fi

        APIS_DIR="{params.apis_dir}"
        mkdir -p "$APIS_DIR"

        if [ ! -s "$APIS_DIR/APIS_db.dmnd" ]; then
            echo "[votu_dbapis_viral] Building dbAPIS Diamond DB in $APIS_DIR" | tee -a {log}
            wget -q -O "$APIS_DIR/anti_defense.pep" \
                "https://pro.unl.edu/dbAPIS/download_file.php?file=anti_defense.pep" \
                >> {log} 2>&1 || echo "[votu_dbapis_viral] WARNING: wget anti_defense.pep failed" >> {log}
            # Family -> (gene name, inhibited defense-type) mapping, one row
            # per APIS family -- confirmed against a real download 2026-06-23
            # (columns: "APIS families", "APIS genes", "Defense systems", ...).
            wget -q -O "$APIS_DIR/seed_and_familyrep_all_infor.tsv" \
                "https://pro.unl.edu/dbAPIS/download_file.php?file=seed_and_familyrep_all_infor.tsv" \
                >> {log} 2>&1 || echo "[votu_dbapis_viral] WARNING: wget seed_and_familyrep_all_infor.tsv failed (mapping disabled, family IDs still reported)" >> {log}
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
            echo "ok" > {output.done}
        else
            echo "[votu_dbapis_viral] No dbAPIS DB available -- writing empty hits" | tee -a {log}
            printf "qseqid\\tsseqid\\tpident\\tlength\\tmismatch\\tgapopen\\tqstart\\tqend\\tsstart\\tsend\\tevalue\\tbitscore\\tqlen\\tslen\\n" > {output.hits}
            echo "skipped: dbAPIS DB not configured" > {output.done}
        fi
        echo "[votu_dbapis_viral] Done -- $(( $(wc -l < {output.hits}) - 1 )) hits" | tee -a {log}
        """
