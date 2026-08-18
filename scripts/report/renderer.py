"""renderer.py — Assembles the VAPOR HTML report from template + data.

Reads shell.html, inlines CSS/JS/assets via {{PLACEHOLDER}} replacement.
No f-strings in this file — all substitution is explicit string.replace().
"""
import os
import glob
import json
import traceback

from .data_loaders import (
    safe_float, safe_int, parse_tsv, load_tsv, load_csv,
    parse_quast_all, parse_support, parse_support_combos, parse_fastp_json, parse_mapping_rate,
    parse_total_reads, build_host_collapse,
    collect_depth_data, parse_fasta_lengths,
    collect_viral_tool_counts, collect_vrhyme_stats, collect_binner_counts,
    parse_checkm2_phyla,
    load_viral_taxonomy, load_viral_source_distribution, load_gtdbtk,
    load_mmseqs_taxonomy_prok, load_phist,
    load_defensefinder, load_antidefensefinder,
    load_antidefensefinder_viral, load_dbapis_viral, compute_defense_islands,
    load_amr_consensus,
    build_host_defense_links, build_bin_annotation_summary,
    enrich_taxonomy_with_checkv, collapse_taxonomy_to_votu, merge_prok_taxonomy,
    load_alpha_diversity, load_pcoord, load_eggnog, load_phrogs,
    load_genome_maps, load_reads_classify, load_coassembly, load_coassembly_rich,
    load_votu_accumulation,
    load_votu_catalog, load_votu_presence,
    load_votu_lifestyle, load_putative_amgs,
    load_tool_status, summarize_tool_status,
    path_dict, collect_tool_versions,
)

_HERE = os.path.dirname(__file__)
_COMP = os.path.join(_HERE, "components")
_ASSETS = os.path.join(_HERE, "assets")


def _read(path):
    with open(path, encoding='utf-8') as f:
        return f.read()


def _jsstr(obj):
    """JSON-encode and escape closing script tags."""
    return json.dumps(obj, ensure_ascii=False).replace("</", "<\\/")


def _drop_field(records, field):
    """Strip a key from every record before embedding in the HTML -- for
    fields only needed server-side (e.g. 'Proteins', the raw per-system
    protein-ID list defensefinder loaders carry purely for
    compute_defense_islands; no JS component reads it, so embedding it
    again in DEFENSE_DATA/ANTIDEFENSE_DATA is pure dead weight in the
    browser bundle)."""
    return [{k: v for k, v in r.items() if k != field} for r in records]


def build_report(snakemake):
    """Entry point called from generate_report.py."""
    try:
        _build(snakemake)
    except Exception:
        out = snakemake.output.html
        os.makedirs(os.path.dirname(out) or '.', exist_ok=True)
        tb = traceback.format_exc()
        with open(out, 'w', encoding='utf-8') as f:
            f.write(f"<pre>VAPOR report generation failed:\n\n{tb}</pre>")
        raise


def _build(snakemake):
    samples        = snakemake.params.samples
    outdir         = snakemake.params.outdir
    low_depth_mode = bool(getattr(snakemake.params, 'low_depth_mode', False))
    out_html       = snakemake.output.html

    tracks_param = dict(getattr(snakemake.params, "tracks", {}) or {})
    tracks_param.setdefault("reads", True)
    tracks_param.setdefault("viral", True)
    tracks_param.setdefault("prok", True)
    tracks_param.setdefault("integration", True)
    tracks_param.setdefault("coassembly", False)

    def _inp(name):
        return list(getattr(snakemake.input, name, []) or [])

    # ── Path dicts ────────────────────────────────────────────────────────────
    quast_paths      = path_dict(_inp('quast'),         samples)
    checkv_paths     = path_dict(_inp('checkv'),        samples)
    checkv_vrh_paths = path_dict(_inp('checkv_vrhyme'), samples)
    checkm2_paths    = path_dict(_inp('checkm2'),       samples)
    support_paths    = path_dict(_inp('support'),       samples)
    depth_paths      = path_dict(_inp('depth'),         samples)
    taxonomy_paths   = path_dict(_inp('taxonomy'),      samples)
    mmseqs_prok_paths= path_dict(_inp('mmseqs_prok'),   samples)
    gtdbtk_bac_l     = _inp('gtdbtk_bac')
    gtdbtk_arc_l     = _inp('gtdbtk_arc')
    phist_paths_l    = [path_dict(_inp('phist'), samples).get(s, '') for s in samples]
    defensefinder_paths_l     = [path_dict(_inp('defensefinder'), samples).get(s, '') for s in samples]
    antidefensefinder_paths_l = [path_dict(_inp('antidefensefinder'), samples).get(s, '') for s in samples]
    antidefense_viral_paths_l = [path_dict(_inp('antidefense_viral'), samples).get(s, '') for s in samples]
    dbapis_viral_paths_l      = [path_dict(_inp('dbapis_viral'), samples).get(s, '') for s in samples]
    prok_protein_manifest_l   = [path_dict(_inp('prok_protein_manifest'), samples).get(s, '') for s in samples]
    amr_consensus_paths_l     = [path_dict(_inp('amr_consensus'), samples).get(s, '') for s in samples]

    # ── Load per-sample data ──────────────────────────────────────────────────
    checkm2_data   = {s: parse_tsv(checkm2_paths[s])     for s in samples}
    quast_data     = {s: parse_quast_all(quast_paths[s]) for s in samples}
    checkv_data    = {s: parse_tsv(checkv_paths[s])      for s in samples}
    checkv_vrh_data= {s: parse_tsv(checkv_vrh_paths[s])  for s in samples}
    support_data   = {s: parse_support(support_paths[s]) for s in samples}
    support_combos = {s: parse_support_combos(support_paths[s]) for s in samples}
    fastp_data     = {s: parse_fastp_json(outdir, s)     for s in samples}
    mapping_data   = {s: parse_mapping_rate(outdir, s)   for s in samples}
    viral_tool_counts = {s: collect_viral_tool_counts(outdir, s) for s in samples}
    checkm2_tax_data  = {s: parse_checkm2_phyla(checkm2_data[s]) for s in samples}
    binner_counts  = {s: collect_binner_counts(outdir, s, checkm2_data[s]) for s in samples}
    vrhyme_data    = {s: collect_vrhyme_stats(outdir, s) for s in samples}
    viral_contig_lengths = {}
    for s in samples:
        vcfa = glob.glob(os.path.join(outdir, s, "viral", "consensus", "*_viral_consensus.fasta"))
        viral_contig_lengths[s] = parse_fasta_lengths(vcfa[0]) if vcfa else []

    # ── Load taxonomy + host prediction ───────────────────────────────────────
    tax_data     = load_viral_taxonomy([taxonomy_paths.get(s, '') for s in samples], samples)
    viral_source_dist = load_viral_source_distribution(
        [taxonomy_paths.get(s, '') for s in samples], samples)
    gtdb_data    = load_gtdbtk(gtdbtk_bac_l, gtdbtk_arc_l, samples)
    phist_data   = load_phist(phist_paths_l, samples)
    defensefinder_data     = load_defensefinder(defensefinder_paths_l, samples)
    antidefensefinder_data = load_antidefensefinder(antidefensefinder_paths_l, samples)
    antidefense_viral_df_data    = load_antidefensefinder_viral(antidefense_viral_paths_l, samples)
    antidefense_viral_dbapis_data = load_dbapis_viral(
        dbapis_viral_paths_l, samples, getattr(snakemake.params, 'apis_db_dir', ''))
    defense_islands = compute_defense_islands(prok_protein_manifest_l, samples, defensefinder_data)
    amr_consensus_data = load_amr_consensus(amr_consensus_paths_l, samples,
                                             low_depth_mode=low_depth_mode)
    host_defense_links = build_host_defense_links(phist_data, gtdb_data)
    bin_annotations = build_bin_annotation_summary(
        defensefinder_data, antidefensefinder_data, amr_consensus_data)
    mmseqs_prok_data = load_mmseqs_taxonomy_prok(mmseqs_prok_paths, samples)
    reads_classify_data = load_reads_classify(
        getattr(snakemake.input, 'reads_classify_abundance', None) or '',
        getattr(snakemake.input, 'reads_classify_host', None) or '',
        samples,
    )
    _coas_groups = list(getattr(snakemake.params, "coassembly_groups", []) or [])
    coassembly_data = load_coassembly(outdir, _coas_groups)
    coassembly_rich = load_coassembly_rich(outdir, _coas_groups)
    votu_accum      = load_votu_accumulation(outdir, _coas_groups)
    # Real per-rule outcome, so a crashed tool renders as a gap, not a zero.
    tool_status     = load_tool_status(outdir, samples, _coas_groups)
    tool_status_issues = summarize_tool_status(tool_status)
    votu_catalog  = load_votu_catalog(outdir)
    votu_presence = load_votu_presence(outdir, samples)
    votu_lifestyle = load_votu_lifestyle(outdir)
    putative_amgs  = load_putative_amgs(outdir)

    tax_data = enrich_taxonomy_with_checkv(tax_data, checkv_data)
    # Collapse rep_seq-level (MMseqs2, 95% identity) rows down to one per
    # vOTU representative (skani, 95% ANI + 85% AF) so the same viral
    # population isn't counted more than once in taxonomy charts/tables.
    tax_data = collapse_taxonomy_to_votu(tax_data, outdir, samples)

    # ── vOTU table + lifestyle ────────────────────────────────────────────────
    votu_data = {}
    for s in samples:
        p = os.path.join(outdir, s, "viral", "votu", f"{s}_vOTU_table.tsv")
        votu_data[s] = load_tsv(p)

    lifestyle_data = {}
    for s in samples:
        lytic = lysogenic = 0
        reps = [r for r in votu_data[s] if str(r.get('is_rep', 'True')).lower() in ('true', '1', 'yes')]
        if not reps:
            reps = votu_data[s]  # fallback: old table without is_rep column
        for row in reps:
            ls = (row.get('lifestyle', '') or '').lower()
            if 'lytic' in ls or 'virulent' in ls: lytic += 1
            elif 'lysogenic' in ls or 'temperate' in ls: lysogenic += 1
        total = len(reps)
        lifestyle_data[s] = {
            'lytic': lytic, 'lysogenic': lysogenic,
            'unknown': total - lytic - lysogenic, 'total': total,
        }

    # ── Abundance ─────────────────────────────────────────────────────────────
    viral_abund = {s: load_tsv(os.path.join(outdir, s, "abundance", "viral_abundance.tsv")) for s in samples}
    prok_abund  = {s: load_tsv(os.path.join(outdir, s, "abundance", "prok_abundance.tsv"))  for s in samples}
    # vOTU-level abundance (skani cluster representatives) — raw read counts
    # summed across cluster members, RPKM/TPM/mean recomputed for the
    # representative; see rules/abundance.smk votu_abundance for the maths.
    votu_abund  = {s: load_tsv(os.path.join(outdir, s, "viral", "votu", f"{s}_vOTU_abundance.tsv"))
                   for s in samples}

    # Enrich votu_data rows with covered_fraction (breadth) and RPMPM.
    # RPMPM = (raw_reads / length_mb) / (total_reads_sample / 1e6)
    #       = raw_reads * 1e12 / (length_bp * total_reads_sample)
    total_reads_per_sample = {s: parse_total_reads(outdir, s) for s in samples}
    for s in samples:
        abund_lookup = {r.get('representative', ''): r for r in votu_abund[s]}
        tot = total_reads_per_sample[s] or 1
        for row in votu_data[s]:
            rep = row.get('representative', row.get('member', ''))
            ab  = abund_lookup.get(rep, {})
            breadth    = safe_float(ab.get('covered_fraction', 0))
            raw_reads  = safe_float(ab.get('total_reads', 0))
            length_bp  = safe_float(row.get('rep_length_bp', 0)) or safe_float(ab.get('length', 0)) or 1
            rpmpm = raw_reads * 1e12 / (length_bp * tot) if raw_reads > 0 else 0.0
            row['breadth']  = f"{breadth:.3f}"
            row['rpmpm']    = f"{rpmpm:.4f}"

    # Host collapse: viral RPKM aggregated by predicted host genus per sample.
    host_collapse_data = build_host_collapse(phist_data, votu_abund, samples)

    # ── Diversity ─────────────────────────────────────────────────────────────
    div_base = os.path.join(outdir, "diversity")
    alpha_rows      = load_alpha_diversity(os.path.join(div_base, "alpha_diversity.tsv"))
    pcoa_viral      = load_pcoord(os.path.join(div_base, "beta_pcoord_viral.tsv"))
    pcoa_prok       = load_pcoord(os.path.join(div_base, "beta_pcoord_prok.tsv"))
    pcoa_combined   = load_pcoord(os.path.join(div_base, "beta_pcoord_combined.tsv"))
    procrustes_rows = load_tsv(os.path.join(div_base, "procrustes_coords.tsv"))

    # ── Annotation ───────────────────────────────────────────────────────────
    eggnog_data  = load_eggnog(outdir, samples)
    phrogs_data  = load_phrogs(outdir, samples)
    genome_maps  = load_genome_maps(outdir, samples)

    # ── Prokaryotic merged taxonomy ───────────────────────────────────────────
    merged_prok = merge_prok_taxonomy(gtdb_data, mmseqs_prok_data, checkm2_data,
                                       low_depth_mode=low_depth_mode)

    # ── Build overview dict ───────────────────────────────────────────────────
    overview = {}
    for s in samples:
        qd  = quast_data[s].get("deduplicated", {})
        cm  = checkm2_data[s]
        cv  = checkv_data[s]
        sp  = support_data[s]
        fq  = fastp_data[s]["reads"]
        raw_r1 = next((r for r in fq if r["stage"] == "raw"     and r["read"] == "R1"), None)
        tr_r1  = next((r for r in fq if r["stage"] == "trimmed" and r["read"] == "R1"), None)
        das = binner_counts[s].get("Binette (final)", {})
        # viral consensus count
        vcfa = os.path.join(outdir, s, "viral", "consensus", f"{s}_viral_consensus.fasta")
        vc_count = sum(1 for l in open(vcfa) if l.startswith('>')) if os.path.exists(vcfa) else \
                   sum(cnt for n, cnt in sp.items() if n >= 2)
        overview[s] = {
            "total_raw_reads":    raw_r1["total_sequences"] if raw_r1 else fastp_data[s]["trim"].get("reads_in", 0),
            "total_trimmed_reads":tr_r1["total_sequences"]  if tr_r1  else 0,
            "mean_qual":          f"{tr_r1['mean_quality']:.1f}" if tr_r1 else "N/A",
            "gc_pct":             f"{tr_r1['gc_percent']:.1f}%"  if tr_r1 else "N/A",
            "mapping_rate":       f"{mapping_data[s]:.1f}%"      if mapping_data[s] > 0 else "N/A",
            "n_contigs":          qd.get("# contigs", "N/A"),
            "n50":                qd.get("N50", "N/A"),
            "viral_consensus":    vc_count,
            # Use checkv_data (all consensus contigs) rather than votu_data (mmseqs
            # cluster representatives only) — the latter often excludes the HQ/Complete
            # contigs entirely if they weren't picked as cluster representatives,
            # which made this KPI read 0% even when the CheckV scatter/donuts show HQ hits.
            "complete_viral":     sum(1 for r in cv if r.get("checkv_quality", "") in ("Complete", "High-quality")),
            "vmags":              vrhyme_data[s]["n_bins"],
            "total_bins":         das.get("total", 0),
            "hq_bins":            sum(1 for r in cm if safe_float(r.get("Completeness", 0)) >= 90
                                      and safe_float(r.get("Contamination", 100)) <= 5),
            "bacteria_bins":      das.get("bacteria", 0),
            "archaea_bins":       das.get("archaea", 0),
            "taxonomy_classified":sum(1 for r in tax_data if r.get("sample") == s),
            "host_pred_total":    len({r["Virus"] for r in phist_data if r.get("sample") == s and r.get("Virus")}),
            "lytic_count":        lifestyle_data[s]["lytic"],
            "lysogenic_count":    lifestyle_data[s]["lysogenic"],
            "lytic_ratio":        round(lifestyle_data[s]["lytic"] / lifestyle_data[s]["total"], 3)
                                  if lifestyle_data[s]["total"] > 0 else 0.0,
            "gtdb_classified":    len({r["Bin"] for r in gtdb_data
                                       if r.get("sample") == s and r.get("Bin")}),
            "total_defense":      sum(1 for r in defensefinder_data if r.get("sample") == s),
            "total_amr_hq":       sum(1 for r in amr_consensus_data
                                      if r.get("sample") == s
                                      and safe_int(r.get("n_tools", 0)) == 3),
        }

    # Back-fill Binette domain from GTDB-Tk
    for s in samples:
        gb = {r["Bin"] for r in gtdb_data if r.get("sample") == s and r.get("Domain") == "Bacteria"}
        ga = {r["Bin"] for r in gtdb_data if r.get("sample") == s and r.get("Domain") == "Archaea"}
        bd = os.path.join(outdir, s, "bins", "binette", "final_bins")
        das = {"total": 0, "bacteria": 0, "archaea": 0, "unknown": 0}
        for bf in glob.glob(os.path.join(bd, "*.fa")):
            name = os.path.basename(bf).replace(".fa", "")
            das["total"] += 1
            if name in ga:   das["archaea"]  += 1
            elif name in gb: das["bacteria"] += 1
            else:            das["unknown"]  += 1
        binner_counts[s]["Binette (final)"] = das
        overview[s].update({"bacteria_bins": das["bacteria"],
                             "archaea_bins":  das["archaea"],
                             "total_bins":    das["total"]})

    # ── Novelty + MIMAG ───────────────────────────────────────────────────────
    novelty_data = {}
    for s in samples:
        total_viral  = overview[s].get("viral_consensus", 0)
        classified   = overview[s].get("taxonomy_classified", 0)
        unclassified = max(0, total_viral - classified)
        pct_novel = round(100.0 * unclassified / total_viral, 1) if total_viral > 0 else 0.0
        novelty_data[s] = {
            'total': total_viral, 'classified': classified,
            'unclassified': unclassified,
            'pct_novel': pct_novel,
        }
        overview[s]['pct_novel'] = pct_novel

    mimag_data = {}
    for s in samples:
        hq = mq = lq = 0
        for row in checkm2_data[s]:
            comp = safe_float(row.get('Completeness', 0))
            cont = safe_float(row.get('Contamination', 100))
            if comp >= 90 and cont <= 5:    hq += 1
            elif comp >= 50 and cont <= 10: mq += 1
            else:                           lq += 1
        mimag_data[s] = {'HQ': hq, 'MQ': mq, 'LQ': lq, 'total': hq + mq + lq}

    # ── Pipeline config (for About tab) ──────────────────────────────────────
    _p = snakemake.params
    cfg_params = {
        "Threads":                     str(getattr(_p, "threads",         "?")),
        "Min contig length":           f"{getattr(_p, 'min_contig',       '?')} bp",
        "SPAdes memory":               f"{getattr(_p, 'spades_mem',       '?')} GB",
        "MEGAHIT memory":              f"{getattr(_p, 'megahit_mem',      '?')} GB",
        "MMseqs2 min identity":        str(getattr(_p, "min_seq_id",      "?")),
        "SemiBin2 environment":        str(getattr(_p, "semibin_env",     "?")),
        "Viral consensus (min tools)": f"{getattr(_p, 'min_viral_tools', '?')} / 3",
    }

    # ── Benchmark data ────────────────────────────────────────────────────────
    bench_path = getattr(snakemake.input, 'benchmark_summary', None)
    bench_data = load_tsv(bench_path) if bench_path and os.path.exists(bench_path) else []

    # ── Tool versions ─────────────────────────────────────────────────────────
    tool_versions = collect_tool_versions()

    # ── Viral depth per-contig ────────────────────────────────────────────────
    viral_depth_data = {}
    for s in samples:
        viral_names = {r.get('contig_id', r.get('contig', '')) for r in checkv_data[s]}
        depth_path  = os.path.join(outdir, s, "mapping", f"{s}_depth.txt")
        vd = []
        if os.path.exists(depth_path):
            for row in parse_tsv(depth_path):
                name = row.get('contigName', row.get('Contig', ''))
                if name in viral_names:
                    d = safe_float(row.get('totalAvgDepth', row.get('avgDepth', 0)))
                    if 0 < d < 2000: vd.append(round(d, 2))
        viral_depth_data[s] = vd

    # ── Serialize all data to JSON ────────────────────────────────────────────
    data_script = _build_data_script({
        "SAMPLES":      samples,
        "OVERVIEW":     overview,
        "CHECKV":       {s: [dict(r) for r in checkv_data[s]] for s in samples},
        "CHECKV_VRH":   {s: [dict(r) for r in checkv_vrh_data[s]] for s in samples},
        "CHECKM2":      {s: [dict(r) for r in checkm2_data[s]] for s in samples},
        "SUPPORT":      support_data,
        "SUPPORT_COMBOS": support_combos,
        "FASTP":        {s: {
            "reads": fastp_data[s]["reads"],
            "trim":  fastp_data[s]["trim"],
        } for s in samples},
        "MAPPING":      mapping_data,
        "QUAST":        {s: quast_data[s] for s in samples},
        "BINNER":       binner_counts,
        "VRHYME":       {s: {
            "n_bins": vrhyme_data[s]["n_bins"],
            "total_members": vrhyme_data[s]["total_members"],
            "rows": vrhyme_data[s]["rows"],
        } for s in samples},
        "VIRAL_TOOLS":  viral_tool_counts,
        "VIRAL_LENGTHS":{s: viral_contig_lengths[s] for s in samples},
        "VIRAL_DEPTH":  viral_depth_data,
        "VIRAL_ABUND":  viral_abund,
        "VOTU_ABUND":   votu_abund,
        "PROK_ABUND":   prok_abund,
        "TAX_DATA":     tax_data,
        "VIRAL_SOURCE_DIST": viral_source_dist,
        "GTDB_DATA":    gtdb_data,
        "MERGED_PROK":  merged_prok,
        "PHIST_DATA":   phist_data,
        "HOST_COLLAPSE":host_collapse_data,
        "DEFENSE_DATA": _drop_field(defensefinder_data, 'Proteins'),
        "ANTIDEFENSE_DATA": _drop_field(antidefensefinder_data, 'Proteins'),
        "ANTIDEFENSE_VIRAL_DF":     antidefense_viral_df_data,
        "ANTIDEFENSE_VIRAL_DBAPIS": antidefense_viral_dbapis_data,
        "DEFENSE_ISLANDS": defense_islands,
        "AMR_CONSENSUS": amr_consensus_data,
        "HOST_DEFENSE_LINKS": host_defense_links,
        "BIN_ANNOTATIONS": bin_annotations,
        "VOTU_DATA":    votu_data,
        "LIFESTYLE":    lifestyle_data,
        "NOVELTY":      novelty_data,
        "MIMAG":        mimag_data,
        "EGGNOG_DATA":  eggnog_data,
        "PHROGS_DATA":  phrogs_data,
        "GENOME_MAPS":  genome_maps,
        "ALPHA_DATA":   alpha_rows,
        "PCOA_VIRAL":   pcoa_viral,
        "PCOA_PROK":    pcoa_prok,
        "PCOA_COMBINED":pcoa_combined,
        "PROCRUSTES":   procrustes_rows,
        "CFG_PARAMS":   cfg_params,
        "TOOL_VERSIONS":tool_versions,
        "BENCH_DATA":   bench_data,
        "READS_CLASSIFY": reads_classify_data,
        "TRACKS":       tracks_param,
        "COASSEMBLY_DATA": coassembly_data,
        "COAS_RICH":       coassembly_rich,
        "VOTU_ACCUM":      votu_accum,
        "TOOL_STATUS":        tool_status,
        "TOOL_STATUS_ISSUES": tool_status_issues,
        "VOTU_CATALOG":  votu_catalog,
        "VOTU_PRESENCE": votu_presence,
        "VOTU_LIFESTYLE": votu_lifestyle,
        "PUTATIVE_AMGS":  putative_amgs,
    })

    # ── Assemble HTML ─────────────────────────────────────────────────────────
    css       = _read(os.path.join(_COMP, "base.css"))
    shell     = _read(os.path.join(_COMP, "shell.html"))
    echarts   = _read(os.path.join(_ASSETS, "echarts.min.js"))
    d3_js     = _read(os.path.join(_ASSETS, "d3.min.js"))

    js_parts  = []
    for js_file in ["app.js", "export.js", "overview.js", "sequencing.js", "viral.js",
                    "prokaryotic.js", "hostdefense.js", "diversity.js", "annotation.js",
                    "reads_classify.js", "coassembly.js", "about.js"]:
        js_parts.append(_read(os.path.join(_COMP, js_file)))
    app_js = "\n".join(js_parts)

    html = (shell
            .replace("{{CSS}}", css)
            .replace("{{ECHARTS_JS}}", echarts)
            .replace("{{D3_JS}}", d3_js)
            .replace("{{DATA_JSON}}", data_script)
            .replace("{{APP_JS}}", app_js))

    os.makedirs(os.path.dirname(out_html) or '.', exist_ok=True)
    with open(out_html, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f"[VAPOR] Report written to {out_html}")


def _build_data_script(data_dict):
    """Build a <script> block defining all JS constants."""
    lines = ["<script>"]
    for name, obj in data_dict.items():
        lines.append(f"const {name} = {_jsstr(obj)};")
    lines.append("</script>")
    return "\n".join(lines)
