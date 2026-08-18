# ══════════════════════════════════════════════════════════════════════
# rules/finalize.smk — BLOCK 11: Organize Final Outputs
#
# Cria estrutura limpa em {sample}/final/:
#
#   final/assembly/
#     quast_report.tsv              — assembly quality metrics (QUAST)
#
#   final/viral/
#     viral_consensus.fasta         — contigs virais confirmados (multi-tool)
#     viral_nonredundant.fasta      — deduplicados (95% ANI)
#     checkv_quality.tsv            — qualidade/completude CheckV
#     viral_bins/                   — vMAGs do vRhyme
#     taxonomy/
#       viral_taxonomy_merged.tsv   — taxonomia 3-tier (MMseqs2/GeNomad/INPHARED)
#     host_prediction/
#       phist_results.tsv           — predição hospedeiro (PHIST)
#     (defense_amr/ viral não existe mais aqui -- DefenseFinder/dbAPIS em
#      ORFs virais rodam uma única vez no catálogo global desde 2026-08-18,
#      ver {OUTDIR}/votu_catalog/{defensefinder,dbapis}/)
#
#   final/bins/
#     bacteria/                     — MAGs Bacteria (GTDB-Tk)
#     archaea/                      — MAGs Archaea (GTDB-Tk)
#     unclassified/                 — MAGs sem domínio definido
#     all_bins_checkm2.tsv          — qualidade CheckM2 todos os bins
#     binette_quality.tsv           — sumário Binette (dereplication)
#     taxonomy/
#       gtdbtk_bacteria.tsv
#       gtdbtk_archaea.tsv
#       mmseqs_taxonomy_prok.tsv    — taxonomia LCA (MMseqs2)
#     defense_amr/
#       defensefinder_systems.tsv   — sistemas de defesa anti-fago (DefenseFinder)
#       antidefensefinder_systems.tsv — sistemas anti-defesa em bins
#       amrfinder_results.tsv       — AMR curado (AMRFinderPlus/NCBI)
#       rgi_results.tsv             — AMR curado (RGI/CARD)
#       deeparg_results.tsv         — AMR exploratório (DeepARG)
#       vfdb_results.tsv            — fatores de virulência (ABRicate/VFDB)
#       plasmidfinder_results.tsv   — replicons plasmidiais (ABRicate/PlasmidFinder)
#       amrfinderplus_normed.tsv    — AMRFinderPlus normalizado ARO (argNorm)
#       deeparg_normed.tsv          — DeepARG normalizado ARO (argNorm)
#
#   final/reads_classify/           — opcional (reads_classify: true em config)
#     merged_relative_abundance_filtered.tsv — abundância relativa sylph-tax (filtrada)
#     merged_relative_abundance.tsv
#     merged_sequence_abundance.tsv
#     otu_table.tsv
#     viral_abundance_by_host.tsv   — abundância viral por hospedeiro predito
# ══════════════════════════════════════════════════════════════════════


rule aggregate_benchmarks:
    """
    Agrega todos os TSVs de benchmark (Snakemake benchmark: directive) num
    único summary: pipeline_timing_summary.tsv.
    Colunas: sample, rule, wall_s, cpu_s, max_rss_mb, mean_load, io_in_mb, io_out_mb.
    Depende de organize_outputs para garantir que todos os rules já rodaram.
    """
    input:
        org = expand(f"{OUTDIR}/{{sample}}/final/done.txt", sample=SAMPLES),
    output:
        summary = f"{OUTDIR}/benchmarks/pipeline_timing_summary.tsv",
    run:
        import csv, glob, os
        from pathlib import Path

        os.makedirs(f"{OUTDIR}/benchmarks", exist_ok=True)
        rows = []
        fields = ["sample", "rule", "wall_s", "cpu_s", "max_rss_mb",
                  "mean_load", "io_in_mb", "io_out_mb"]

        def _read_bench(tsv_path, sample_name):
            try:
                with open(tsv_path) as f:
                    for row in csv.DictReader(f, delimiter="\t"):
                        return {
                            "sample":      sample_name,
                            "rule":        Path(tsv_path).stem,
                            "wall_s":      row.get("s", ""),
                            "cpu_s":       row.get("cpu_time", ""),
                            "max_rss_mb":  row.get("max_rss", ""),
                            "mean_load":   row.get("mean_load", ""),
                            "io_in_mb":    row.get("io_in", ""),
                            "io_out_mb":   row.get("io_out", ""),
                        }
            except Exception:
                return None

        for sample in SAMPLES:
            for tsv in sorted(glob.glob(f"{OUTDIR}/{sample}/benchmarks/*.tsv")):
                r = _read_bench(tsv, sample)
                if r:
                    rows.append(r)

        for tsv in sorted(glob.glob(f"{OUTDIR}/benchmarks/*.tsv")):
            if Path(tsv).name == "pipeline_timing_summary.tsv":
                continue
            r = _read_bench(tsv, "global")
            if r:
                rows.append(r)

        with open(output.summary, "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader()
            w.writerows(rows)


rule organize_outputs:
    """
    Classifica bins por domínio usando GTDB-Tk (fallback: CheckM2 lineage).
    Copia outputs virais, metagenômicos, defense/AMR e assembly para final/
    com estrutura limpa.
    """
    input:
        # Assembly
        **({
            "quast": rules.quast.output.report,
        } if (TRACK_VIRAL or TRACK_PROK) else {}),
        # Viral core + taxonomy/host. Viral defense/anti-defense (vdef/
        # vantidef/dbapis) removed 2026-08-18 (second half of "(h)",
        # docs/ROADMAP_SIMPLIFICACAO.md): defensefinder_viral/dbapis_viral
        # now run once globally (votu_defensefinder_viral/votu_dbapis_viral,
        # rules/votu_catalog.smk), so there is no per-sample file left to
        # copy here -- copying the (now nonexistent) per-sample path would
        # either break the DAG or, before this rule moved, silently copy
        # the WHOLE catalog into a sample-scoped final/ path. The single
        # correct copy is the global one at
        # {OUTDIR}/votu_catalog/{defensefinder,dbapis}/.
        **({
            "checkv":     rules.checkv.output.summary,
            "viral":      rules.viral_consensus.output.fasta,
            "viral_nr":   rules.viral_nonredundant.output.fasta,
            "vrhyme":     rules.vrhyme.output.done,
            "taxonomy":   rules.viral_taxonomy.output.tsv,
        } if TRACK_VIRAL else {}),
        # Prok bins + taxonomy + defense/anti-defense + AMR
        **({
            "checkm2":      rules.checkm2.output.report,
            "gtdbtk_b":     rules.gtdbtk.output.bac_tsv,
            "gtdbtk_a":     rules.gtdbtk.output.ar_tsv,
            "binette":      rules.binette.output.summary,
            "mmseqs_prok":  rules.mmseqs_taxonomy_prok.output.hits,
            "pdef":         rules.defensefinder.output.systems,
            "pantidef":     rules.defensefinder.output.antisystems,
            "amrfinder":    rules.amrfinderplus.output.results,
            "rgi":          rules.rgi_card.output.results,
            "deeparg":      rules.deeparg.output.results,
            "vfdb":         rules.abricate.output.vfdb,
            "plasmidfind":  rules.abricate.output.plasmidfinder,
            "amf_normed":   rules.argnorm_normalize.output.amrfinder_normed,
            "deeparg_norm": rules.argnorm_normalize.output.deeparg_normed,
        } if TRACK_PROK else {}),
        # Host prediction / integration
        **({
            "phist": rules.phist.output.results,
        } if INTEGRATION_ENABLED else {}),
    output:
        done = f"{OUTDIR}/{{sample}}/final/done.txt",
    log:
        f"{OUTDIR}/{{sample}}/logs/organize_outputs.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/organize_outputs.tsv"
    params:
        s = f"{OUTDIR}/{{sample}}",
    run:
        import os, glob, shutil

        s     = params.s
        final = f"{s}/final"

        for d in [
            f"{final}/assembly",
            f"{final}/viral/viral_bins",
            f"{final}/viral/taxonomy",
            f"{final}/viral/host_prediction",
            f"{final}/bins/bacteria",
            f"{final}/bins/archaea",
            f"{final}/bins/unclassified",
            f"{final}/bins/taxonomy",
            f"{final}/bins/defense_amr",
        ]:
            os.makedirs(d, exist_ok=True)

        def cp(src, dst):
            if src and os.path.exists(str(src)):
                shutil.copy(str(src), dst)

        g = lambda name: getattr(input, name, None)

        with open(log[0], "w") as lf:

            # ── Assembly ───────────────────────────────────────────────
            cp(g('quast'), f"{final}/assembly/quast_report.tsv")

            # ── Viral core ────────────────────────────────────────────
            cp(g('viral'),    f"{final}/viral/viral_consensus.fasta")
            cp(g('viral_nr'), f"{final}/viral/viral_nonredundant.fasta")
            cp(g('checkv'),   f"{final}/viral/checkv_quality.tsv")

            if TRACK_VIRAL:
                vrhyme_bins = glob.glob(f"{s}/bins/vrhyme/vRhyme_best_bins.*.fasta")
                for bf in vrhyme_bins:
                    shutil.copy(bf, f"{final}/viral/viral_bins/")
                lf.write(f"vRhyme bins: {len(vrhyme_bins)}\n")

            # ── Viral taxonomy + host ─────────────────────────────────
            cp(g('taxonomy'),  f"{final}/viral/taxonomy/viral_taxonomy_merged.tsv")
            cp(g('phist'),     f"{final}/viral/host_prediction/phist_results.tsv")

            # Viral defense / anti-defense: no longer copied here (see the
            # removal note above `input:`) -- only the global vOTU catalog
            # copy at {OUTDIR}/votu_catalog/{defensefinder,dbapis}/ exists.

            # ── Prokaryotic bins — classify with GTDB-Tk ──────────────
            archaea_bins, bacteria_bins = set(), set()

            if TRACK_PROK:
                gtdbtk_b, gtdbtk_a = g('gtdbtk_b'), g('gtdbtk_a')
                for tsv_path in [str(gtdbtk_b) if gtdbtk_b else '', str(gtdbtk_a) if gtdbtk_a else '']:
                    is_arc = 'ar53' in tsv_path or 'ar_tsv' in tsv_path
                    if not tsv_path or not os.path.exists(tsv_path):
                        continue
                    with open(tsv_path) as f:
                        hdr = None
                        for line in f:
                            parts = line.strip().split('\t')
                            if hdr is None:
                                hdr = [h.lower() for h in parts]; continue
                            if not parts or len(parts) < 2: continue
                            bin_name = parts[0]
                            if is_arc:
                                archaea_bins.add(bin_name)
                            else:
                                bacteria_bins.add(bin_name)

                # Fallback to CheckM2 lineage if GTDB-Tk empty
                checkm2_path = g('checkm2')
                if not archaea_bins and not bacteria_bins and checkm2_path and os.path.exists(str(checkm2_path)):
                    with open(str(checkm2_path)) as f:
                        hdr = None
                        for line in f:
                            parts = line.strip().split('\t')
                            if hdr is None:
                                hdr = [h.lower() for h in parts]; continue
                            if not parts or len(parts) < 2: continue
                            bin_name = parts[0]
                            tax = ''
                            for col in ['taxonomic_lineage', 'lineage', 'taxonomy']:
                                if col in hdr and hdr.index(col) < len(parts):
                                    tax = parts[hdr.index(col)]; break
                            if 'Archaea' in tax: archaea_bins.add(bin_name)
                            elif 'Bacteria' in tax: bacteria_bins.add(bin_name)

                lf.write(f"GTDB-Tk — Bacteria: {len(bacteria_bins)}, Archaea: {len(archaea_bins)}\n")

                copied = {'bacteria': 0, 'archaea': 0, 'unclassified': 0}
                bins_dir = f"{s}/bins/binette/final_bins"
                for bf in glob.glob(f"{bins_dir}/*.fa"):
                    bin_name = os.path.basename(bf).replace('.fa', '')
                    if bin_name in archaea_bins:
                        shutil.copy(bf, f"{final}/bins/archaea/"); copied['archaea'] += 1
                    elif bin_name in bacteria_bins:
                        shutil.copy(bf, f"{final}/bins/bacteria/"); copied['bacteria'] += 1
                    else:
                        shutil.copy(bf, f"{final}/bins/unclassified/"); copied['unclassified'] += 1

                lf.write(f"Bins copied: {copied}\n")

            # ── Prok QC + taxonomy ────────────────────────────────────
            cp(g('checkm2'),     f"{final}/bins/all_bins_checkm2.tsv")
            cp(g('binette'),     f"{final}/bins/binette_quality.tsv")
            cp(g('gtdbtk_b'),    f"{final}/bins/taxonomy/gtdbtk_bacteria.tsv")
            cp(g('gtdbtk_a'),    f"{final}/bins/taxonomy/gtdbtk_archaea.tsv")
            cp(g('mmseqs_prok'), f"{final}/bins/taxonomy/mmseqs_taxonomy_prok.tsv")

            # ── Prok defense / anti-defense ───────────────────────────
            cp(g('pdef'),     f"{final}/bins/defense_amr/defensefinder_systems.tsv")
            cp(g('pantidef'), f"{final}/bins/defense_amr/antidefensefinder_systems.tsv")

            # ── Prok AMR ──────────────────────────────────────────────
            cp(g('amrfinder'),    f"{final}/bins/defense_amr/amrfinder_results.tsv")
            cp(g('rgi'),          f"{final}/bins/defense_amr/rgi_results.tsv")
            cp(g('deeparg'),      f"{final}/bins/defense_amr/deeparg_results.tsv")
            cp(g('vfdb'),         f"{final}/bins/defense_amr/vfdb_results.tsv")
            cp(g('plasmidfind'),  f"{final}/bins/defense_amr/plasmidfinder_results.tsv")
            cp(g('amf_normed'),   f"{final}/bins/defense_amr/amrfinderplus_normed.tsv")
            cp(g('deeparg_norm'), f"{final}/bins/defense_amr/deeparg_normed.tsv")

        with open(output.done, 'w') as f:
            f.write('ok\n')


if READS_CLASSIFY_ENABLED:
    rule finalize_reads_classify:
        """
        Copia outputs do reads_classify (sylph) para final/reads_classify/.
        Executa uma vez por pipeline (global, não por sample).
        """
        input:
            done     = f"{OUTDIR}/reads_classify/reads_classify_done.txt",
        output:
            sentinel = f"{OUTDIR}/final/reads_classify/done.txt",
        run:
            import os, shutil
            src = f"{OUTDIR}/reads_classify"
            dst = f"{OUTDIR}/final/reads_classify"
            os.makedirs(dst, exist_ok=True)
            for fname in [
                "merged_relative_abundance_filtered.tsv",
                "merged_relative_abundance.tsv",
                "merged_sequence_abundance.tsv",
                "otu_table.tsv",
                "viral_abundance_by_host.tsv",
            ]:
                src_f = os.path.join(src, fname)
                if os.path.exists(src_f):
                    shutil.copy(src_f, os.path.join(dst, fname))
            with open(output.sentinel, 'w') as f:
                f.write('ok\n')


rule finalize_votu_catalog:
    """
    Copia os outputs globais do catalogo de vOTU para final/votu_catalog/.
    Roda uma vez por pipeline, como o finalize_reads_classify.

    Existe porque o item (h) (docs/ROADMAP_SIMPLIFICACAO.md) globalizou
    pharokka, phold, genome maps, defensefinder viral e dbAPIS: eles deixaram
    de ter caminho por amostra, e o organize_outputs por amostra/grupo passou
    a criar final/viral/defense_amr/ sem que nenhuma regra escrevesse dentro.
    O dado nao sumiu -- vive em {outdir}/votu_catalog/ -- mas parou de ser
    entregue em final/. Esta regra reata a entrega, agora num lugar so, que e
    o que a natureza global do dado pede.
    """
    input:
        catalog_done = f"{OUTDIR}/votu_catalog/done.txt",
        taxonomy     = f"{OUTDIR}/votu_catalog/taxonomy/viral_taxonomy_merged.tsv",
        pharokka     = f"{OUTDIR}/votu_catalog/annotation/pharokka/done.txt",
        phold        = f"{OUTDIR}/votu_catalog/annotation/phold/done.txt",
        **({
            "defense": f"{OUTDIR}/votu_catalog/defensefinder/done.txt",
            "dbapis":  f"{OUTDIR}/votu_catalog/dbapis/done.txt",
        } if DEFENSE_AMR_VIRAL_ENABLED else {}),
    output:
        sentinel = f"{OUTDIR}/final/votu_catalog/done.txt",
    log:
        f"{OUTDIR}/logs/finalize_votu_catalog.log"
    run:
        import os, shutil

        src = f"{OUTDIR}/votu_catalog"
        dst = f"{OUTDIR}/final/votu_catalog"
        os.makedirs(dst, exist_ok=True)

        # (origem relativa a votu_catalog/, destino relativo a final/votu_catalog/)
        wanted = [
            ("vOTU_clusters.tsv",                              "vOTU_clusters.tsv"),
            ("catalog_all_reps.fasta",                         "catalog_all_reps.fasta"),
            ("catalog_mq_reps.fasta",                          "catalog_mq_reps.fasta"),
            ("provenance.tsv",                                 "provenance.tsv"),
            ("presence_matrix.tsv",                            "presence_matrix.tsv"),
            ("votu_abundance_matrix.tsv",                      "votu_abundance_matrix.tsv"),
            ("taxonomy/viral_taxonomy_merged.tsv",             "viral_taxonomy_merged.tsv"),
            ("bacphlip/votu_lifestyle.tsv",                    "votu_lifestyle.tsv"),
            ("eggnog_viral/putative_amgs.tsv",                 "putative_amgs.tsv"),
            ("annotation/pharokka/pharokka_cds_final_merged_output.tsv",
                                                               "pharokka_cds.tsv"),
            ("annotation/phold/phold.gbk",                     "phold.gbk"),
            ("defensefinder/viral_defense_systems.tsv",        "viral_defense_systems.tsv"),
            ("defensefinder/viral_antidefense_systems.tsv",    "viral_antidefense_systems.tsv"),
            ("dbapis/dbapis_hits.tsv",                         "dbapis_hits.tsv"),
        ]

        copied, missing = [], []
        for rel_src, rel_dst in wanted:
            src_f = os.path.join(src, rel_src)
            # Ausente e legitimo: a fonte pode ter sido pulada (DB nao
            # configurado) ou nao ter produzido nada. Fica registrado no log
            # em vez de virar erro -- mas fica registrado, para "sumiu" e
            # "nao rodou" nao se confundirem.
            if os.path.exists(src_f) and os.path.getsize(src_f) > 0:
                shutil.copy(src_f, os.path.join(dst, rel_dst))
                copied.append(rel_dst)
            else:
                missing.append(rel_src)

        with open(str(log[0]), "w") as lf:
            lf.write(f"[finalize_votu_catalog] copiados: {len(copied)}\n")
            for name in copied:
                lf.write(f"  ok      {name}\n")
            for name in missing:
                lf.write(f"  ausente {name}\n")

        with open(output.sentinel, "w") as f:
            f.write("ok\n")
