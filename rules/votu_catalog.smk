# ══════════════════════════════════════════════════════════════════════
# rules/votu_catalog.smk — BLOCK 7.5: Global vOTU catalog
#
# Replaces the former per-sample skani chain (skani_votu / skani_cluster /
# viral_votu_reps). A vOTU is defined ONCE over the pooled viral sets of
# every sample and co-assembly group, so richness is comparable across
# samples and against the literature. Per-sample presence comes from
# read recruitment against the catalog (see votu_catalog_abundance).
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
        sources += [
            ("group", g,
             f"{OUTDIR}/coassembly/{g}/viral/checkv/{g}_viral_trimmed.fasta")
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
    """CheckV completeness and quality tier, keyed by the pool's namespaced IDs."""
    import csv
    completeness = {}
    quality = {}
    for source_id, path in _catalog_checkv_pairs():
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                contig = (row.get("contig_id") or "").strip()
                if not contig:
                    continue
                key = f"{source_id}|{contig}"
                quality[key] = (row.get("checkv_quality") or "").strip()
                try:
                    completeness[key] = float(row.get("completeness") or 0)
                except ValueError:
                    completeness[key] = 0.0
    return completeness, quality


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
    """Extract the three representative tiers used downstream.

    Same quality gates as the removed per-sample viral_votu_reps, applied
    once over the global catalog:
      all      — one per vOTU; recruitment reference and report base
      mq       — MQ+ (Complete/HQ/MQ or completeness >= 50%); taxonomy,
                 PHIST, annotation
      hq_10kb  — HQ+/Complete and >= 10 kb; vConTACT3
    """
    input:
        pool     = rules.votu_catalog_pool.output.pool,
        clusters = rules.votu_catalog_cluster.output.clusters,
        checkv   = _catalog_checkv_summaries,
    output:
        all_fasta     = f"{CATALOG_DIR}/catalog_all_reps.fasta",
        mq_fasta      = f"{CATALOG_DIR}/catalog_mq_reps.fasta",
        hq_10kb_fasta = f"{CATALOG_DIR}/catalog_hq_10kb_reps.fasta",
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

        def is_hq(rid):
            return quality.get(rid, "") in ("Complete", "High-quality")

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

        n_all = n_mq = n_hq = 0
        with open(str(output.all_fasta), "w") as fa, \
             open(str(output.mq_fasta), "w") as fm, \
             open(str(output.hq_10kb_fasta), "w") as fh10:
            for rid, chunks in seqs.items():
                seq = "".join(chunks)
                record = f">{rid}\n{seq}\n"
                fa.write(record)
                n_all += 1
                if is_mq(rid):
                    fm.write(record)
                    n_mq += 1
                if is_hq(rid) and len(seq) >= 10000:
                    fh10.write(record)
                    n_hq += 1

        with open(str(log[0]), "w") as lf:
            lf.write(f"[votu_catalog_reps] all: {n_all}\n")
            lf.write(f"[votu_catalog_reps] MQ+ (taxonomy/PHIST/annotation): {n_mq}\n")
            lf.write(f"[votu_catalog_reps] HQ+/>=10kb (vConTACT3): {n_hq}\n")
        write_status(str(output.done), "ok")
