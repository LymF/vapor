# ══════════════════════════════════════════════════════════════════════
# rules/finalize.smk — BLOCK 11: Organize Final Outputs
#
# Cria estrutura limpa em {sample}/final/:
#
#   final/viral/
#     viral_consensus.fasta       — contigs virais confirmados
#     checkv_quality.tsv          — qualidade CheckV
#     viral_bins/                 — vMAGs do vRhyme
#     taxonomy/
#       viral_taxonomy_merged.tsv — taxonomia 3-tier
#       vcontact3_clusters.tsv    — clusters vConTACT3
#     host_prediction/
#       phist_results.tsv
#
#   final/bins/
#     bacteria/                   — MAGs classificados como Bactéria
#     archaea/                    — MAGs classificados como Archaea
#     unclassified/               — MAGs sem domínio definido
#     taxonomy/
#       gtdbtk_bacteria.tsv
#       gtdbtk_archaea.tsv
#     all_bins_checkm2.tsv
#     binette_quality.tsv
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
    Copia outputs virais e metagenômicos para final/ com estrutura limpa.
    """
    input:
        checkm2  = rules.checkm2.output.report,
        checkv   = rules.checkv.output.summary,
        viral    = rules.viral_consensus.output.fasta,
        viral_nr = rules.viral_nonredundant.output.fasta,
        vrhyme   = rules.vrhyme.output.done,
        gtdbtk_b = rules.gtdbtk.output.bac_tsv,
        gtdbtk_a = rules.gtdbtk.output.ar_tsv,
        phist    = rules.phist.output.results,
        taxonomy  = rules.viral_taxonomy.output.tsv,
        vcontact3 = rules.vcontact3.output.network,
        binette  = rules.binette.output.summary,
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

        for d in [f"{final}/viral/viral_bins",
                  f"{final}/viral/taxonomy",
                  f"{final}/viral/host_prediction",
                  f"{final}/bins/bacteria",
                  f"{final}/bins/archaea",
                  f"{final}/bins/unclassified",
                  f"{final}/bins/taxonomy"]:
            os.makedirs(d, exist_ok=True)

        with open(log[0], "w") as lf:

            # ── Viral outputs ──────────────────────────────────────────
            shutil.copy(input.viral,    f"{final}/viral/viral_consensus.fasta")
            shutil.copy(input.viral_nr, f"{final}/viral/viral_nonredundant.fasta")
            shutil.copy(input.checkv,   f"{final}/viral/checkv_quality.tsv")

            # vRhyme bins
            vrhyme_bins = glob.glob(f"{s}/bins/vrhyme/vRhyme_best_bins.*.fasta")
            for bf in vrhyme_bins:
                shutil.copy(bf, f"{final}/viral/viral_bins/")
            lf.write(f"vRhyme bins: {len(vrhyme_bins)}\n")

            # Viral taxonomy (3-tier: vConTACT3 > INPHARED > GeNomad)
            if os.path.exists(str(input.taxonomy)):
                shutil.copy(input.taxonomy, f"{final}/viral/taxonomy/viral_taxonomy_merged.tsv")
            if os.path.exists(str(input.vcontact3)):
                shutil.copy(input.vcontact3, f"{final}/viral/taxonomy/vcontact3_clusters.tsv")

            # Host predictions
            if os.path.exists(str(input.phist)):
                shutil.copy(input.phist, f"{final}/viral/host_prediction/phist_results.tsv")

            # ── Prokaryotic bins — classify with GTDB-Tk ──────────────
            archaea_bins, bacteria_bins = set(), set()

            # Parse GTDB-Tk bacteria summary
            for tsv_path in [str(input.gtdbtk_b), str(input.gtdbtk_a)]:
                is_arc = 'ar53' in tsv_path or 'ar_tsv' in tsv_path
                if not os.path.exists(tsv_path):
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
            if not archaea_bins and not bacteria_bins:
                with open(str(input.checkm2)) as f:
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

            # Copy Binette bins to final/
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

            # Copy QC/taxonomy tables
            shutil.copy(str(input.checkm2),  f"{final}/bins/all_bins_checkm2.tsv")
            shutil.copy(str(input.binette),  f"{final}/bins/binette_quality.tsv")
            if os.path.exists(str(input.gtdbtk_b)):
                shutil.copy(str(input.gtdbtk_b), f"{final}/bins/taxonomy/gtdbtk_bacteria.tsv")
            if os.path.exists(str(input.gtdbtk_a)):
                shutil.copy(str(input.gtdbtk_a), f"{final}/bins/taxonomy/gtdbtk_archaea.tsv")
            lf.write(f"Bins copied: {copied}\n")

        with open(output.done, 'w') as f:
            f.write('ok\n')
