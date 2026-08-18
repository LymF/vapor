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
