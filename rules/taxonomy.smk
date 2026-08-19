# ══════════════════════════════════════════════════════════════════════
# rules/taxonomy.smk — BLOCK 9: Viral Taxonomy (per-sample view) +
#                       Prokaryote custom taxonomy
#
# The viral taxonomy computation itself (prodigal_viral -> mmseqs_taxonomy_
# viral/custom -> the GeNomad/MMseqs merge) moved to rules/votu_catalog.smk
# on 2026-08-18 as `votu_prodigal` / `votu_mmseqs_taxonomy` /
# `votu_mmseqs_taxonomy_custom` / `votu_taxonomy` -- see
# docs/ROADMAP_SIMPLIFICACAO.md "(h) Princípio: computar no representante,
# herdar no membro". A vOTU is the same biological entity in every sample
# that carries a member of it, so its taxonomy is computed ONCE on the
# global vOTU representative, not once per sample on byte-identical input.
#
# `viral_taxonomy` below is what is left on the per-sample side (kept its
# old name so rules.viral_taxonomy references, e.g. finalize.smk, don't
# need to change): it
# reads the global table and re-exposes it at the old per-sample path with
# the old per-sample schema (bare contig IDs), so make_votu_table.py,
# rules/report.smk and rules/finalize.smk keep working unmodified.
#
# Also contains mmseqs_taxonomy_prok for prokaryote bin taxonomy -- this one
# stays per-sample on purpose: prokaryote bins are not deduplicated into a
# global catalog the way vOTUs are, so there is no single global
# representative to compute this against.
# ══════════════════════════════════════════════════════════════════════


rule viral_taxonomy:
    """
    Per-sample/group view of the global vOTU taxonomy table (rule
    votu_taxonomy, rules/votu_catalog.smk). Implements INHERITANCE, not
    identity: every member of a vOTU that has a representative annotated in
    the global table gets THAT representative's annotation, not just
    members that happen to be their own vOTU's representative. This is the
    premise the whole (h) architecture rests on -- a vOTU is an ANI>=95% /
    AF>=85% cluster (ICTV species level), so every member is by
    construction the same biological entity as its representative, which is
    exactly what authorizes computing the analysis once and reusing it.

    (2026-08-18 correction: an earlier version of this rule only filtered
    the global table by source_id prefix, which is IDENTITY, not
    inheritance -- it only recovered a row for a member that IS the
    representative of its own vOTU. With 32 samples most members' vOTU
    representative belongs to a different sample, so that version left the
    taxonomy column silently empty for most contigs. See
    docs/ROADMAP_SIMPLIFICACAO.md "(h)" for the full writeup.)

    JOIN, spelled out (three tables, two join steps):
      1. clusters (votu_catalog_cluster.output.clusters): `member`
         (NAMESPACED) -> `representative` (NAMESPACED). Built into
         member_to_rep.
      2. provenance (votu_catalog_pool.output.provenance): `member_id`
         (NAMESPACED) -> `original_contig_id` (BARE), filtered to rows
         whose `source_id` == this sample/group -- this enumerates every
         member contig belonging to this source, not just the ones chosen
         as a representative.
      3. global_tsv (votu_taxonomy.output.tsv): `seq_name` (NAMESPACED,
         representative IDs only -- one row per MQ+ vOTU representative)
         -> annotation row.
      For every member of this source (step 2): look up its representative
      (step 1), then look up that representative's row in global_tsv (step
      3). A vOTU whose representative isn't MQ+ has no row in global_tsv --
      the member then gets NO output row. There is no fallback to another
      member of the same vOTU; that gap is intentional (see "quality gate"
      note below), not a bug.

    Output side: `seq_name` is written as the MEMBER's own BARE
    `original_contig_id` (never the representative's), so every other
    per-sample table (CheckV, PHIST) that make_votu_table.py joins against
    keeps working. Two columns are appended at the end, additive only --
    confirmed both make_votu_table.py (load_taxonomy, csv.DictReader +
    row.get by name) and the report loader
    (scripts/report/data_loaders.py:load_viral_taxonomy, same pattern)
    read named columns and ignore trailing ones, so this does not break
    either contract:
      votu_representative -- namespaced ID of the representative the row's
                              annotation was inherited from
      is_representative    -- "True" if this member IS its own vOTU's
                              representative, else "False"

    Quality-gate note: global_tsv only has rows for MQ+ representatives
    (votu_catalog_reps.output.mq_fasta gate). A vOTU whose representative
    falls below that gate has no annotation ANYWHERE in the catalog -- not
    a per-sample limitation, the same vOTU would be unannotated in every
    sample that carries one of its members. Logged below as "representative
    has no row in global table (not MQ+)".
    """
    input:
        global_tsv = rules.votu_taxonomy.output.tsv,
        provenance = rules.votu_catalog_pool.output.provenance,
        clusters   = rules.votu_catalog_cluster.output.clusters,
    output:
        tsv  = f"{OUTDIR}/{{sample}}/viral/taxonomy/viral_taxonomy_merged.tsv",
        done = f"{OUTDIR}/{{sample}}/viral/taxonomy/taxonomy_done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/viral_taxonomy_view.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/viral_taxonomy_view.tsv"
    threads: 1
    params:
        # derivado do output, nao de {{sample}}: regra herdada por
        # coassembly.smk via `use rule ... as ... with:` (wildcard {group}).
        source_id = lambda wc, output: os.path.basename(
            os.path.dirname(os.path.dirname(os.path.dirname(output.tsv)))),
    run:
        import csv, os
        from pathlib import Path

        os.makedirs(os.path.dirname(str(output.tsv)), exist_ok=True)

        # 1. clusters: member (namespaced) -> representative (namespaced)
        member_to_rep = {}
        with open(str(input.clusters)) as fh:
            for row in csv.DictReader(fh, delimiter="	"):
                member_to_rep[row["member"]] = row["representative"]

        # 3. global annotation table: representative seq_name (namespaced,
        #    MQ+ only) -> full annotation row
        rep_annotation = {}
        fields = None
        with open(str(input.global_tsv)) as fh:
            rdr = csv.DictReader(fh, delimiter="	")
            fields = rdr.fieldnames or []
            for row in rdr:
                rep_annotation[row.get("seq_name", "")] = row

        extra_fields = ["votu_representative", "is_representative"]
        out_fields = list(fields) + extra_fields

        # 2. provenance, filtered to this source's members, drives the
        #    output: one row attempted per member, not per representative.
        n_members = 0
        n_no_rep = 0
        n_rep_not_mq = 0
        n_written = 0
        n_is_rep = 0
        rows_out = []
        with open(str(input.provenance)) as fh:
            for row in csv.DictReader(fh, delimiter="	"):
                if row.get("source_id", "") != params.source_id:
                    continue
                n_members += 1
                member_id = row["member_id"]
                bare_id   = row["original_contig_id"]

                rep = member_to_rep.get(member_id)
                if not rep:
                    n_no_rep += 1
                    continue

                annot = rep_annotation.get(rep)
                if annot is None:
                    n_rep_not_mq += 1
                    continue

                out_row = dict(annot)
                out_row["seq_name"] = bare_id
                out_row["votu_representative"] = rep
                is_rep = (member_id == rep)
                out_row["is_representative"] = "True" if is_rep else "False"
                if is_rep:
                    n_is_rep += 1
                rows_out.append(out_row)
                n_written += 1

        with open(str(output.tsv), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=out_fields, delimiter="	")
            w.writeheader()
            w.writerows(rows_out)

        with open(str(log[0]), "w") as lf:
            lf.write(f"[viral_taxonomy_view] source: {params.source_id}\n")
            lf.write(f"[viral_taxonomy_view] members of this source (provenance): {n_members}\n")
            lf.write(f"[viral_taxonomy_view] members with NO representative found in clusters.tsv "
                     f"(should be 0 -- pool/clusters mismatch if not): {n_no_rep}\n")
            lf.write(f"[viral_taxonomy_view] members whose representative has NO row in the "
                     f"global taxonomy table (representative not MQ+): {n_rep_not_mq}\n")
            lf.write(f"[viral_taxonomy_view] rows written (member -> inherited annotation): "
                     f"{n_written}/{n_members}\n")
            lf.write(f"[viral_taxonomy_view]   of which this source's own vOTU representative: "
                     f"{n_is_rep}\n")
            if n_members and n_written == 0:
                lf.write("[viral_taxonomy_view] WARNING: 0 rows written for a non-empty source "
                          "-- check clusters.tsv / global_tsv namespace alignment\n")

        Path(str(output.done)).write_text("ok\n")


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
    rule exists to avoid.

    Replaces diamond_custom_prok entirely (removed) -- custom prokaryote
    taxonomy now runs exclusively through this rule. Skips gracefully if
    custom_prok_mmseqs_db isn't configured, same as any other optional DB.

    Stays per-sample: prokaryote bins are not deduplicated into a global
    catalog the way vOTUs are (see rules/votu_catalog.smk), so there is no
    single global representative to move this computation to.
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
