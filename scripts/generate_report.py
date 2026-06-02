#!/usr/bin/env python3
"""
generate_report.py — MITE (Metagenomic Integrated Taxonomic Engine) HTML report.
Tabs: Overview | Read QC | Assembly QC | Viral Analysis | Taxonomy | Host Prediction | Bin Quality | Abundance | About
"""
import json, os, re, glob, sys, zipfile, gzip, traceback
from collections import defaultdict, Counter

try:
    import plotly.graph_objects as go
    import plotly.express as px
    import plotly.io as pio
    from plotly.subplots import make_subplots
except ImportError:
    sys.exit("ERROR: pip install plotly")

try:
    samples  = snakemake.params.samples
    outdir   = snakemake.params.outdir
    out_html = snakemake.output.html

    # New module inputs — handle gracefully if not present
    def _inp(name):
        return list(getattr(snakemake.input, name, []) or [])
    taxonomy_paths  = _inp('taxonomy')
    gtdbtk_bac      = _inp('gtdbtk_bac')
    gtdbtk_arc      = _inp('gtdbtk_arc')
    phist_paths     = _inp('phist')
    binette_paths   = _inp('binette')
    benchmark_path  = getattr(snakemake.input, 'benchmark_summary', None)

    # ── Pipeline configuration (from Snakefile) ───────────────────────────────
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
    
    # ══════════════════════════════════════════════════════════════════════════════
    #  HELPERS
    # ══════════════════════════════════════════════════════════════════════════════
    
    def safe_float(v, d=0.0):
        try: return float(str(v).replace(",","").replace(">","").strip())
        except: return d
    
    def safe_int(v, d=0):
        try: return int(str(v).replace(",","").replace(">","").strip())
        except: return d
    
    def parse_tsv(path):
        rows = []
        if not path or not os.path.exists(path): return rows
        with open(path) as f:
            hdr = None
            for line in f:
                line = line.strip()
                if not line: continue
                parts = line.split("\t")
                if hdr is None: hdr = parts; continue
                if len(parts) == len(hdr): rows.append(dict(zip(hdr, parts)))
        return rows
    
    def parse_quast_all(path):
        data = {}
        if not path or not os.path.exists(path): return data
        labels = None
        with open(path) as f:
            for line in f:
                parts = line.strip().split("\t")
                if not parts: continue
                if parts[0] == "Assembly":
                    labels = parts[1:]
                    for lbl in labels: data[lbl] = {}
                    continue
                if labels and len(parts) >= 2:
                    metric = parts[0].strip()
                    for i, lbl in enumerate(labels):
                        data[lbl][metric] = parts[i+1].strip() if i+1 < len(parts) else ""
        return data
    
    def parse_support(path):
        counts = {1:0, 2:0, 3:0, 4:0}
        for row in parse_tsv(path):
            try:
                n = int(row.get("n_tools", 0))
                if n in counts: counts[n] += 1
            except: pass
        return counts
    
    def compute_bin_abundance(depth_path, s2b_dir):
        cov, clen = {}, {}
        if os.path.exists(depth_path):
            with open(depth_path) as f:
                first = True
                for line in f:
                    if first: first = False; continue
                    p = line.strip().split("\t")
                    if len(p) >= 3:
                        cov[p[0]]  = safe_float(p[2])
                        clen[p[0]] = safe_int(p[1])
        bin_totals = defaultdict(lambda: {"cov":0.0,"len":0})
        for s2b in glob.glob(os.path.join(s2b_dir, "*_s2b.tsv")):
            binner = os.path.basename(s2b).replace("_s2b.tsv","")
            with open(s2b) as f:
                for line in f:
                    p = line.strip().split("\t")
                    if len(p) == 2:
                        contig, bn = p; key = f"{binner}::{bn}"
                        l = clen.get(contig, 0)
                        bin_totals[key]["cov"] += cov.get(contig,0.0)*l
                        bin_totals[key]["len"] += l
        return {k: v["cov"]/v["len"] for k,v in bin_totals.items() if v["len"]>0}
    
    
    # ── fastp ─────────────────────────────────────────────────────────────────────
    def parse_fastp_json(outdir, sample):
        """Parse fastp JSON → {"reads": [...], "trim": {...}}
        reads: list of dicts with total_sequences, gc_percent, mean_quality, stage, read, r_label
        trim:  dict with reads_in, reads_written, adapter_r1_pct, adapter_r2_pct, bp_removed_pct
        """
        import json as _json
        _empty = {"reads": [], "trim": {
            "reads_in": 0, "reads_written": 0,
            "adapter_r1_pct": 0.0, "adapter_r2_pct": 0.0, "bp_removed_pct": 0.0}}
        json_path = os.path.join(outdir, sample, "qc_raw", f"{sample}_fastp.json")
        if not os.path.exists(json_path): return _empty
        try:
            with open(json_path) as f: j = _json.load(f)
        except Exception: return _empty

        reads = []
        for stage_key, stage_lbl in [("before_filtering", "raw"), ("after_filtering", "trimmed")]:
            for read_key, read_lbl in [("read1", "R1"), ("read2", "R2")]:
                rd = j.get(f"{read_key}_{stage_key}") or j.get("summary", {}).get(stage_key, {})
                q30 = safe_float(rd.get("q30_rate", 0) or 0)
                reads.append({
                    "stage":           stage_lbl,
                    "read":            read_lbl,
                    "r_label":         f"{stage_lbl} {read_lbl}",
                    "total_sequences": safe_int(rd.get("total_reads", 0)),
                    "gc_percent":      safe_float(rd.get("gc_content", 0)) * 100,
                    "mean_quality":    q30 * 35 + (1 - q30) * 10,
                })

        summ   = j.get("summary", {})
        bf     = summ.get("before_filtering", {})
        af     = summ.get("after_filtering",  {})
        reads_in   = safe_int(bf.get("total_reads", 0)) or 1
        bases_in   = safe_int(bf.get("total_bases", 0)) or 1
        bases_out  = safe_int(af.get("total_bases", 0))
        reads_out  = safe_int(af.get("total_reads", 0))
        adp        = j.get("adapter_cutting", {})
        adp_reads  = safe_int(adp.get("adapter_trimmed_reads", 0))
        adp_pct    = adp_reads / reads_in * 100

        trim = {
            "reads_in":        reads_in,
            "reads_written":   reads_out,
            "adapter_r1_pct":  adp_pct,
            "adapter_r2_pct":  adp_pct,
            "bp_removed_pct":  (bases_in - bases_out) / bases_in * 100,
        }
        return {"reads": reads, "trim": trim}

    def parse_mapping_rate(outdir, sample):
        flagstat_path = os.path.join(outdir, sample, "mapping", "flagstat.txt")
        if not os.path.exists(flagstat_path): return 0.0
        with open(flagstat_path) as f: content = f.read()
        m = re.search(r"mapped \(([\d.]+)%", content)
        if m: return safe_float(m.group(1))
        return 0.0
    
    
    # ── Depth / contig data ───────────────────────────────────────────────────────
    def collect_depth_data(depth_path):
        lengths, depths = [], []
        if not os.path.exists(depth_path): return lengths, depths
        with open(depth_path) as f:
            hdr = None
            for line in f:
                # depth file can be tab-separated or space-separated
                raw = line.rstrip("\n")
                p = raw.split("\t") if "\t" in raw else raw.split()
                if len(p) < 3: continue
                if hdr is None:
                    hdr = p; continue  # skip header
                try:
                    lengths.append(int(p[1]))
                    depths.append(float(p[2]))
                except (ValueError, IndexError):
                    continue
        return lengths, depths

    def parse_fasta_lengths(fasta_path):
        lengths = []; curr = 0
        if not fasta_path or not os.path.exists(fasta_path): return lengths
        with open(fasta_path) as f:
            for line in f:
                line = line.strip()
                if line.startswith(">"):
                    if curr > 0: lengths.append(curr)
                    curr = 0
                else: curr += len(line)
        if curr > 0: lengths.append(curr)
        return lengths
    
    
    # ── Viral data ────────────────────────────────────────────────────────────────
    def collect_viral_tool_counts(outdir, sample):
        counts = {}
        vs2 = os.path.join(outdir, sample, "viral", "virsorter2", "final-viral-combined.fa")
        counts["VirSorter2"] = sum(1 for l in open(vs2) if l.startswith(">")) if os.path.exists(vs2) else 0
        counts["GeNomad"] = sum(
            len(parse_tsv(p))
            for p in glob.glob(os.path.join(outdir, sample, "viral", "genomad", "*_summary", "*_virus_summary.tsv"))
        )
        # VIBRANT: count scaffolds in VIBRANT_phages_* directory
        vib_fastas = glob.glob(os.path.join(outdir, sample, "viral", "vibrant",
                                            "**", "VIBRANT_phages_*", "*.fna"), recursive=True)
        vib_combined = glob.glob(os.path.join(outdir, sample, "viral", "vibrant",
                                              "**", "*phages_combined*"), recursive=True)
        vib_seqs = set()
        for fa in vib_fastas + vib_combined:
            try:
                for l in open(fa):
                    if l.startswith(">"): vib_seqs.add(l[1:].split()[0])
            except: pass
        counts["VIBRANT"] = len(vib_seqs)
        return counts
    
    def collect_viral_scores(outdir, sample):
        scores = {"GeNomad":[], "VirSorter2":[]}
        for p in glob.glob(os.path.join(outdir, sample, "viral", "genomad", "*_summary", "*_virus_summary.tsv")):
            for row in parse_tsv(p):
                v = safe_float(row.get("virus_score", row.get("Virus_score",-1)))
                if 0 <= v <= 1: scores["GeNomad"].append(v)
        vs2_score = os.path.join(outdir, sample, "viral", "virsorter2", "final-viral-score.tsv")
        for row in parse_tsv(vs2_score):
            v = safe_float(row.get("max_score", row.get("score",-1)))
            if 0 <= v <= 1: scores["VirSorter2"].append(v)
        return scores
    
    
    # ── vRhyme ────────────────────────────────────────────────────────────────────
    def collect_vrhyme_stats(outdir, sample):
        vdir      = os.path.join(outdir, sample, "bins", "vrhyme")
        fasta_dir = os.path.join(vdir, "vRhyme_best_bins_fasta")
        n_bins = len([f for f in os.listdir(fasta_dir)
                      if f.endswith(".fasta")]) if os.path.isdir(fasta_dir) else 0
        summary_rows = []
        total_members = 0
        for sf in glob.glob(os.path.join(vdir, "vRhyme_best_bins.*.summary.tsv")):
            rows = parse_tsv(sf)
            summary_rows.extend(rows)
            total_members += sum(safe_int(r.get("members",0)) for r in rows)
        if n_bins == 0 and summary_rows:
            n_bins = len(summary_rows)
        return {"n_bins": n_bins, "total_members": total_members, "rows": summary_rows}
    
    
    # ── Binning ───────────────────────────────────────────────────────────────────
    def collect_binner_counts(outdir, sample, checkm2_rows):
        """
        Count bins per tool using confirmed output structures:
          MetaBAT2 : bin.*.fa
          VAMB     : clusters.tsv (unique bin IDs)
          SemiBin2 : SemiBin_*.fa.gz in output_bins/
                     fallback: unique bins in contig_bins.tsv
          Binette  : binette/final_bins/*.fa
        """
        b = os.path.join(outdir, sample, "bins")

        # MetaBAT2
        mb_count = len(glob.glob(os.path.join(b, "metabat2", "bin.*.fa")))
    
        # VAMB - count unique bin IDs in clusters.tsv
        vamb_count = 0
        vamb_clusters = os.path.join(b, "vamb", "clusters.tsv")
        if os.path.exists(vamb_clusters):
            bins = set()
            with open(vamb_clusters) as f:
                for line in f:
                    parts = line.strip().split("\t")
                    if len(parts) == 2 and not parts[0].startswith("clusterid"):
                        bins.add(parts[0])
            vamb_count = len(bins)
        if vamb_count == 0:
            vamb_count = len(glob.glob(os.path.join(b, "vamb", "bins", "*.fna")))
    
        # SemiBin2 — count SemiBin_*.fa.gz files; fallback: unique bins in contig_bins.tsv
        sb_gz = glob.glob(os.path.join(b, "semibin2", "output_bins", "SemiBin_*.fa.gz"))
        if sb_gz:
            sb_count = len(sb_gz)
        else:
            sb_tsv = os.path.join(b, "semibin2", "contig_bins.tsv")
            if os.path.exists(sb_tsv):
                bins = set()
                with open(sb_tsv) as f:
                    header = True
                    for line in f:
                        if header: header = False; continue
                        parts = line.strip().split("\t")
                        if len(parts) == 2: bins.add(parts[1])
                sb_count = len(bins)
            else:
                sb_count = 0
    
        # Binette final bins
        das_dir  = os.path.join(b, "binette", "final_bins")
        das_bins = glob.glob(os.path.join(das_dir, "*.fa"))
        das = {"total": len(das_bins), "bacteria":0, "archaea":0, "unknown":len(das_bins)}
    
        return {
            "MetaBAT2":        {"total": mb_count,   "bacteria":0,"archaea":0,"unknown":0},
            "VAMB":            {"total": vamb_count, "bacteria":0,"archaea":0,"unknown":0},
            "SemiBin2":        {"total": sb_count,   "bacteria":0,"archaea":0,"unknown":0},
            "Binette (final)":das,
        }
    
    
    # ── Taxonomy ──────────────────────────────────────────────────────────────────
    def parse_genomad_taxonomy(outdir, sample):
        records = []
        for path in glob.glob(os.path.join(outdir, sample, "viral", "genomad", "*_summary", "*_virus_summary.tsv")):
            for row in parse_tsv(path):
                tax = row.get("taxonomy", row.get("Taxonomy","")).strip()
                if not tax or tax in ("NA","n/a",""): continue
                parts = [re.sub(r"^[-_ ]+","",p).strip() for p in tax.split(";")]
                while len(parts) < 4: parts.append("Unclassified")
                records.append({"realm":parts[0] or "Unclassified","kingdom":parts[1] or "Unclassified",
                                "phylum":parts[2] or "Unclassified","cls":parts[3] or "Unclassified"})
        return records
    
    def parse_checkm2_phyla(checkm2_rows):
        phyla = Counter()
        for row in checkm2_rows:
            tax = row.get("Taxonomic_lineage",row.get("taxonomic_lineage","")).strip()
            if not tax: phyla["Unclassified"] += 1; continue
            phylum = "Unclassified"
            for p in tax.split(";"):
                if p.strip().startswith("p__"): phylum = p.strip()[3:] or "Unclassified"; break
            phyla[phylum] += 1
        return phyla
    
    
    # ══════════════════════════════════════════════════════════════════════════════

    # ── Tool versions (collected at runtime via conda run) ────────────────────
    import subprocess as _sp

    def _get_ver(cmd):
        try:
            r = _sp.run(cmd, shell=True, capture_output=True, text=True, timeout=20)
            out = (r.stdout + r.stderr).strip()
            lines = [l for l in out.splitlines() if l.strip()]
            if not lines: return "n/a"
            import re as _re
            for l in lines:
                m = _re.search(r"[\d]+\.[\d]+[\d.]*", l)
                if m: return m.group(0)
            return lines[0][:40]
        except Exception:
            return "n/a"

    _C = "conda run --no-capture-output -n"

    # CheckV: version appears in first line of bare `checkv` call (no --version flag)
    def _get_ver2(cmd, line_idx=0):
        """Like _get_ver but choose which output line to parse."""
        try:
            r = _sp.run(cmd, shell=True, capture_output=True, text=True, timeout=20)
            out = (r.stdout + r.stderr).strip()
            lines = [l for l in out.splitlines() if l.strip()]
            if not lines: return "n/a"
            import re as _re
            target = lines[line_idx] if line_idx < len(lines) else lines[0]
            m = _re.search(r"[\d]+\.[\d]+[\d.]*", target)
            return m.group(0) if m else target[:40]
        except Exception:
            return "n/a"

    tool_versions = {
        # QC
        "fastp":         _get_ver(f"{_C} env_qc fastp --version"),
        "MultiQC":       _get_ver(f"{_C} env_qc multiqc --version"),
        # Assembly
        "MEGAHIT":       _get_ver(f"{_C} env_assembly megahit --version"),
        "metaSPAdes":    _get_ver(f"{_C} env_assembly spades.py --version"),
        "MMseqs2":       _get_ver(f"{_C} env_mmseqs mmseqs version"),
        "QUAST":         _get_ver(f"{_C} env_qc quast.py --version"),
        # Mapping
        "BWA-MEM2":      _get_ver(f"{_C} env_mapping bwa-mem2 version"),
        "minimap2":      _get_ver(f"{_C} env_mapping minimap2 --version"),
        "samtools":      _get_ver(f"{_C} env_mapping samtools --version"),
        "CoverM":        _get_ver(f"{_C} env_coverm coverm --version"),
        # Viral detection
        "VirSorter2":    _get_ver(f"{_C} env_viral virsorter --version"),
        "GeNomad":       _get_ver(f"{_C} env_genomad genomad --version"),
        "VIBRANT":       _get_ver(f"{_C} env_viral VIBRANT_run.py 2>&1 | head -3"),
        "DeepVirFinder": _get_ver(f"{_C} env_dvf python 2>&1 | head -1"),
        "CenoteTaker3":  _get_ver(f"{_C} env_cenote cenote-taker 2>&1 | head -2"),
        # Viral QC & taxonomy
        "CheckV":        _get_ver2(f"{_C} env_viral checkv 2>&1 | head -1", line_idx=0),
        "vRhyme":        _get_ver(f"{_C} env_vrhyme vRhyme --version"),
        "vConTACT3":     _get_ver(f"{_C} env_vcontact3 vcontact3 --version"),
        "Diamond":       _get_ver(f"{_C} env_viral diamond version"),
        "Prodigal":      _get_ver(f"{_C} env_viral prodigal -v 2>&1 | head -2"),
        # Prokaryote binning
        "MetaBAT2":      _get_ver(f"{_C} env_binning metabat2 2>&1 | head -3"),
        "MaxBin2":       _get_ver(f"{_C} env_binning run_MaxBin.pl 2>&1 | head -2"),
        "VAMB":          _get_ver(f"{_C} env_binning vamb --version"),
        "SemiBin2":      _get_ver(f"{_C} env_binning SemiBin2 --version"),
        "Binette":       _get_ver(f"{_C} env_binette binette --version"),
        "CheckM2":       _get_ver(f"{_C} env_checkm2 checkm2 --version"),
        "GTDB-Tk":       _get_ver(f"{_C} env_gtdbtk gtdbtk --version"),
        # Host prediction
        "PHIST":         _get_ver(f"{_C} env_phist phist 2>&1 | head -1"),
        # Annotation
        "Pharokka":      _get_ver(f"{_C} env_annotation pharokka --version"),
        "EggNOG-mapper": _get_ver(f"{_C} env_annotation emapper.py --version"),
    }
    print("[generate_report] versions:", {k:v for k,v in tool_versions.items() if v != "n/a"}, flush=True)

    #  LOAD DATA
    # ══════════════════════════════════════════════════════════════════════════════
    
    # Build path dicts by matching sample name in path (robust to Snakemake expand ordering)
    def _path_dict(paths, samples):
        d = {}
        for p in paths:
            for s in samples:
                if f"/{s}/" in p or p.endswith(f"/{s}"):
                    d[s] = p; break
        # fallback: zip order if pattern matching fails
        if len(d) < len(samples):
            d = {s:p for s,p in zip(samples, paths)}
        return d

    quast_paths      = _path_dict(snakemake.input.quast,         samples)
    checkv_paths     = _path_dict(snakemake.input.checkv,        samples)
    checkv_vrh_paths = _path_dict(snakemake.input.checkv_vrhyme, samples)
    checkm2_paths    = _path_dict(snakemake.input.checkm2,       samples)
    support_paths    = _path_dict(snakemake.input.support,       samples)
    depth_paths      = _path_dict(snakemake.input.depth,         samples)
    
    quast_data={}; checkv_data={}; checkv_vrh_data={}; checkm2_data={}; support_data={}
    abundance_data={}; fastp_data={}; mapping_data={}
    viral_tool_counts={}; viral_scores={}; viral_contig_lengths={}
    final_contig_lengths={}; coverage_data={}
    checkm2_tax_data={}
    binner_counts={}; vrhyme_data={}
    
    for sample in samples:
        checkm2_data[sample]      = parse_tsv(checkm2_paths[sample])
        quast_data[sample]        = parse_quast_all(quast_paths[sample])
        checkv_data[sample]       = parse_tsv(checkv_paths[sample])
        checkv_vrh_data[sample]   = parse_tsv(checkv_vrh_paths[sample])
        support_data[sample]      = parse_support(support_paths[sample])
        s2b                       = os.path.join(outdir, sample, "bins", "scaffold2bin")
        abundance_data[sample]    = compute_bin_abundance(depth_paths[sample], s2b)
        fastp_data[sample]        = parse_fastp_json(outdir, sample)
        mapping_data[sample]      = parse_mapping_rate(outdir, sample)
        viral_tool_counts[sample] = collect_viral_tool_counts(outdir, sample)
        viral_scores[sample]      = collect_viral_scores(outdir, sample)
        checkm2_tax_data[sample]  = parse_checkm2_phyla(checkm2_data[sample])
        binner_counts[sample]     = collect_binner_counts(outdir, sample, checkm2_data[sample])
        vrhyme_data[sample]       = collect_vrhyme_stats(outdir, sample)
        final_contig_lengths[sample], coverage_data[sample] = collect_depth_data(depth_paths[sample])
        vcfa = glob.glob(os.path.join(outdir, sample, "viral", "consensus", "*_viral_consensus.fasta"))
        viral_contig_lengths[sample] = parse_fasta_lengths(vcfa[0]) if vcfa else []
    
    
    # ══════════════════════════════════════════════════════════════════════════════
    #  FIGURES
    # ══════════════════════════════════════════════════════════════════════════════
    
    # ── Design system colours ──────────────────────────────────────────────────
    TEAL   = "#0d9488"   # virology primary
    AMBER  = "#d97706"   # prokaryotes
    GREEN  = "#16a34a"
    RED    = "#ef4444"
    GRAY   = "#6b7280"
    PURPLE = "#7c3aed"
    CYAN   = "#0891b2"
    BLUE   = "#2563eb"
    # keep legacy aliases for code that still references BLUE/ORANGE
    ORANGE = AMBER

    PALETTE = [TEAL, AMBER, PURPLE, CYAN, GREEN, "#f59e0b", "#9333ea",
               RED, "#0e7490", "#b45309"]

    # ── Register custom Plotly template ───────────────────────────────────────
    _light_layout = dict(
        font=dict(family="Inter, 'Source Sans Pro', system-ui, sans-serif", size=12, color="#1e293b"),
        paper_bgcolor="rgba(0,0,0,0)",
        plot_bgcolor="#f8fafc",
        colorway=PALETTE,
        xaxis=dict(gridcolor="#e2e8f0", linecolor="#cbd5e1", zerolinecolor="#e2e8f0",
                   tickfont=dict(size=11)),
        yaxis=dict(gridcolor="#e2e8f0", linecolor="#cbd5e1", zerolinecolor="#e2e8f0",
                   tickfont=dict(size=11)),
        legend=dict(bgcolor="rgba(255,255,255,0.85)", bordercolor="#e2e8f0",
                    borderwidth=1, font=dict(size=11)),
        margin=dict(t=40, b=60, l=60, r=20),
        hoverlabel=dict(bgcolor="#1e293b", font=dict(color="#f1f5f9", size=12)),
        title=dict(font=dict(size=14, color="#0f172a"), x=0.02),
    )
    pio.templates["mite_light"] = go.layout.Template(layout=_light_layout)
    T = "mite_light"
    ncols = max(len(samples), 1)
    
    # ─ fastp ──────────────────────────────────────────────────────────────────────
    stage_colors={"raw R1":CYAN,"raw R2":"#67e8f9","trimmed R1":TEAL,"trimmed R2":"#5eead4"}
    all_fq_labels=[]
    for s in samples:
        for rec in fastp_data[s]["reads"]:
            if rec["r_label"] not in all_fq_labels: all_fq_labels.append(rec["r_label"])
    fig_fq_reads=go.Figure(); fig_fq_qual=go.Figure(); fig_fq_gc=go.Figure()
    for lbl in all_fq_labels:
        rv,qv,gv=[],[],[]
        for s in samples:
            rec=next((r for r in fastp_data[s]["reads"] if r["r_label"]==lbl),None)
            rv.append(rec["total_sequences"] if rec else 0)
            qv.append(rec["mean_quality"]    if rec else 0)
            gv.append(rec["gc_percent"]      if rec else 0)
        color=stage_colors.get(lbl,GRAY); kw=dict(name=lbl,x=samples,marker_color=color)
        fig_fq_reads.add_trace(go.Bar(**kw,y=rv,hovertemplate="<b>%{x}</b><br>"+lbl+": %{y:,}<extra></extra>"))
        fig_fq_qual.add_trace( go.Bar(**kw,y=qv,hovertemplate="<b>%{x}</b><br>"+lbl+": %{y:.1f}<extra></extra>"))
        fig_fq_gc.add_trace(   go.Bar(**kw,y=gv,hovertemplate="<b>%{x}</b><br>"+lbl+": %{y:.1f}%<extra></extra>"))
    for fig,title,yl in [(fig_fq_reads,"Total Reads — Raw vs Trimmed","Read count"),
                         (fig_fq_qual,"Mean Per-Base Quality Score","Phred score"),
                         (fig_fq_gc,"GC Content — Raw vs Trimmed","% GC")]:
        fig.update_layout(barmode="group",title=title,yaxis_title=yl,legend_title="Read set",template=T,height=420)
    fig_fq_qual.add_hline(y=30,line_dash="dash",line_color=GREEN,annotation_text="Q30",annotation_position="right")
    fig_fq_qual.add_hline(y=20,line_dash="dot", line_color=AMBER,annotation_text="Q20",annotation_position="right")
    fig_fq_gc.add_hline(y=60,line_dash="dot",line_color=GRAY,annotation_text="60%",annotation_position="right")
    fig_fq_gc.add_hline(y=40,line_dash="dot",line_color=GRAY,annotation_text="40%",annotation_position="right")

    fig_trim=make_subplots(rows=1,cols=3,horizontal_spacing=0.10,
        subplot_titles=["Adapter reads — R1 (%)","Adapter reads — R2 (%)","Bases quality-trimmed (%)"])
    for col,(key,color) in enumerate([("adapter_r1_pct",AMBER),("adapter_r2_pct","#fbbf24"),("bp_removed_pct",RED)],1):
        vals=[fastp_data[s]["trim"].get(key,0) for s in samples]
        fig_trim.add_trace(go.Bar(x=samples,y=vals,marker_color=color,showlegend=False,
            hovertemplate=f"<b>%{{x}}</b><br>%{{y:.1f}}%<extra></extra>"),row=1,col=col)
    fig_trim.update_layout(title="Fastp Trimming Efficiency",height=400,template=T)
    
    fig_mapping=go.Figure()
    mvals=[mapping_data[s] for s in samples]
    mcolors=[GREEN if v>=70 else ORANGE if v>=50 else RED if v>0 else GRAY for v in mvals]
    fig_mapping.add_trace(go.Bar(x=samples,y=mvals,marker_color=mcolors,showlegend=False,
        text=[f"{v:.1f}%" if v>0 else "N/A" for v in mvals],textposition="outside",
        hovertemplate="<b>%{x}</b><br>Mapped: %{y:.1f}%<extra></extra>"))
    fig_mapping.add_hline(y=70,line_dash="dash",line_color=GREEN,annotation_text="70% (good)",annotation_position="right")
    fig_mapping.add_hline(y=50,line_dash="dot", line_color=AMBER,annotation_text="50% (acceptable)",annotation_position="right")
    fig_mapping.update_layout(title="Read Mapping Rate (BWA-MEM2 → Assembly)",
        yaxis_title="% reads mapped",yaxis=dict(range=[0,115]),template=T,height=380)
    
    # ─ Assembly ───────────────────────────────────────────────────────────────────
    asm_stages=["MEGAHIT","metaSPAdes","metaviralSPAdes","merged_filtered","deduplicated"]
    asm_colors=[CYAN,"#67e8f9",PURPLE,AMBER,TEAL]
    fig_asm_prog=make_subplots(rows=1,cols=3,horizontal_spacing=0.10,
        subplot_titles=["N50 (bp)","Total Length (bp)","Number of Contigs"])
    for col,(metric,) in enumerate([("N50",),("Total length",),("# contigs",)],1):
        for si,sample in enumerate(samples):
            for stage,color in zip(asm_stages,asm_colors):
                val=safe_float(quast_data[sample].get(stage,{}).get(metric,0))
                fig_asm_prog.add_trace(go.Bar(name=stage,x=[sample],y=[val],marker_color=color,
                    legendgroup=stage,showlegend=(col==1 and si==0),
                    hovertemplate=f"<b>%{{x}} — {stage}</b><br>{metric}: %{{y:,.0f}}<extra></extra>"),row=1,col=col)
    fig_asm_prog.update_layout(barmode="group",title="Assembly Progression",height=450,template=T,legend_title="Stage")
    
    fig_assembly=make_subplots(rows=1,cols=3,horizontal_spacing=0.10,
        subplot_titles=["N50 (bp)","Total Length (bp)","Number of Contigs"])
    for col,(metric,color) in enumerate([("N50",TEAL),("Total length",GREEN),("# contigs",AMBER)],1):
        vals=[safe_float(quast_data[s].get("deduplicated",{}).get(metric,0)) for s in samples]
        fig_assembly.add_trace(go.Bar(x=samples,y=vals,marker_color=color,showlegend=False,
            hovertemplate=f"<b>%{{x}}</b><br>{metric}: %{{y:,.0f}}<extra></extra>"),row=1,col=col)
    fig_assembly.update_layout(title="Final Assembly Quality (QUAST — deduplicated)",height=420,template=T)
    
    # Ensure contig lengths: depth.txt primary, fasta fallback
    for sample in samples:
        if not final_contig_lengths[sample]:
            for fpath in [
                os.path.join(outdir, sample, "mmseqs", f"{sample}_rep_seq.fasta"),
                os.path.join(outdir, sample, "assembly", "merged", f"{sample}_merged_contigs.fasta"),
            ]:
                if os.path.exists(fpath):
                    final_contig_lengths[sample] = parse_fasta_lengths(fpath)
                    break

    # Ensure contig lengths — fallback to mmseqs fasta if depth.txt missing/empty
    for sample in samples:
        if not final_contig_lengths[sample]:
            for fpath in [
                os.path.join(outdir, sample, "mmseqs", f"{sample}_rep_seq.fasta"),
                os.path.join(outdir, sample, "assembly", "merged", f"{sample}_merged_contigs.fasta"),
            ]:
                if os.path.exists(fpath):
                    final_contig_lengths[sample] = parse_fasta_lengths(fpath); break

    fig_contig_len=go.Figure()
    for i,sample in enumerate(samples):
        lengths=[l for l in final_contig_lengths[sample] if l>0]
        if lengths:
            fig_contig_len.add_trace(go.Histogram(x=lengths,name=sample,opacity=0.7,nbinsx=60,
                marker_color=PALETTE[i%len(PALETTE)],
                hovertemplate="Length: %{x:,} bp<br>Count: %{y}<extra>"+sample+"</extra>"))
    fig_contig_len.update_layout(barmode="overlay",xaxis_type="log",
        title="Contig Length Distribution (log x)",xaxis_title="Length (bp)",yaxis_title="Count",template=T,height=420)
    
    fig_cov=go.Figure()
    for i,sample in enumerate(samples):
        depths=[d for d in coverage_data[sample] if 0<d<500]
        if depths:
            fig_cov.add_trace(go.Histogram(x=depths,name=sample,opacity=0.7,nbinsx=60,
                marker_color=PALETTE[i%len(PALETTE)],
                hovertemplate="Coverage: %{x:.1f}×<br>Count: %{y}<extra>"+sample+"</extra>"))
    fig_cov.update_layout(barmode="overlay",title="Coverage Distribution (clipped 500×)",
        xaxis_title="Mean depth (×)",yaxis_title="Count",template=T,height=420)
    
    # ─ Viral ──────────────────────────────────────────────────────────────────────
    tools=["VirSorter2","GeNomad","VIBRANT"]
    tcols={"VirSorter2":CYAN,"GeNomad":TEAL,"VIBRANT":"#9333ea"}
    fig_tool_counts=go.Figure()
    for tool in tools:
        fig_tool_counts.add_trace(go.Bar(name=tool,x=samples,
            y=[viral_tool_counts[s].get(tool,0) for s in samples],marker_color=tcols[tool],
            hovertemplate="<b>%{x}</b><br>"+tool+": %{y}<extra></extra>"))
    fig_tool_counts.update_layout(barmode="group",title="Viral Contigs per Tool (Before Consensus Filter)",
        yaxis_title="Contigs",legend_title="Tool",template=T,height=420)
    
    fig_viral=go.Figure()
    for n,lbl,color in zip([1,2,3],["1 tool","2 tools","3 tools (consensus)"],[GRAY,AMBER,TEAL]):
        fig_viral.add_trace(go.Bar(name=lbl,x=samples,
            y=[support_data[s].get(n,0) for s in samples],marker_color=color,
            hovertemplate="<b>%{x}</b><br>"+lbl+": %{y}<extra></extra>"))
    fig_viral.update_layout(barmode="stack",title="Tool Agreement — Consensus Filter",
        yaxis_title="Contigs",legend_title="Agreement",template=T,height=420)
    
    fig_scores=go.Figure()
    for sample in samples:
        for tool,color in [("VirSorter2",CYAN),("GeNomad",TEAL)]:
            vals=viral_scores[sample].get(tool,[])
            if vals:
                fig_scores.add_trace(go.Violin(y=vals,name=f"{sample} — {tool}",
                    box_visible=True,meanline_visible=True,marker_color=color,opacity=0.7))
    fig_scores.update_layout(title="Viral Score Distributions",yaxis_title="Score (0–1)",template=T,height=480)
    
    checkv_cats=["Complete","High-quality","Medium-quality","Low-quality","Not-determined"]
    cv_colors=[GREEN,"#4ade80","#fbbf24",AMBER,RED]
    fig_checkv=make_subplots(rows=1,cols=ncols,specs=[[{"type":"pie"}]*ncols],subplot_titles=samples)
    for i,sample in enumerate(samples,1):
        cc=Counter(r.get("checkv_quality","Not-determined") for r in checkv_data[sample])
        fig_checkv.add_trace(go.Pie(labels=checkv_cats,values=[cc.get(c,0) for c in checkv_cats],
            marker_colors=cv_colors,showlegend=(i==1),
            hovertemplate="<b>%{label}</b><br>%{value} (%{percent})<extra></extra>"),row=1,col=i)
    fig_checkv.update_layout(title="CheckV — Consensus Contigs Quality",height=440,template=T)
    
    fig_checkv_vrh=make_subplots(rows=1,cols=ncols,specs=[[{"type":"pie"}]*ncols],subplot_titles=samples)
    for i,sample in enumerate(samples,1):
        cc=Counter(r.get("checkv_quality","Not-determined") for r in checkv_vrh_data[sample])
        fig_checkv_vrh.add_trace(go.Pie(labels=checkv_cats,values=[cc.get(c,0) for c in checkv_cats],
            marker_colors=cv_colors,showlegend=(i==1),
            hovertemplate="<b>%{label}</b><br>%{value} (%{percent})<extra></extra>"),row=1,col=i)
    fig_checkv_vrh.update_layout(title="CheckV — vRhyme vMAGs Quality",height=440,template=T)
    
    fig_checkv_scatter=go.Figure()
    qual_cmap={"Complete":GREEN,"High-quality":"#4ade80","Medium-quality":"#fbbf24","Low-quality":AMBER,"Not-determined":GRAY}
    shown_cv=set()
    for sample in samples:
        for source,rows,sym in [("consensus",checkv_data[sample],"circle"),
                                 ("vRhyme",checkv_vrh_data[sample],"diamond")]:
            for row in rows:
                try:
                    length=safe_int(row.get("contig_length",0)); comp=safe_float(row.get("completeness",0))
                    qual=row.get("checkv_quality","Not-determined"); name=row.get("contig_id",row.get("contig",""))
                    color=qual_cmap.get(qual,GRAY); key=f"{qual}_{source}"; show=key not in shown_cv
                    if show: shown_cv.add(key)
                    fig_checkv_scatter.add_trace(go.Scatter(x=[length],y=[comp],mode="markers",
                        name=f"{qual} ({source})",legendgroup=qual,showlegend=show,
                        marker=dict(size=9,symbol=sym,color=color,line=dict(width=1,color="white")),
                        text=[f"{name} ({sample}, {source})"],
                        hovertemplate="<b>%{text}</b><br>Length: %{x:,} bp<br>Completeness: %{y:.1f}%<extra></extra>"))
                except: continue
    fig_checkv_scatter.update_layout(title="CheckV — Length vs Completeness (● consensus  ◆ vRhyme)",
        xaxis_title="Contig length (bp)",yaxis_title="Completeness (%)",
        xaxis_type="log",legend_title="Quality",template=T,height=500)
    
    fig_viral_len=go.Figure()
    for i,sample in enumerate(samples):
        vl=viral_contig_lengths[sample]
        if vl:
            fig_viral_len.add_trace(go.Histogram(x=vl,name=sample,opacity=0.75,nbinsx=30,
                marker_color=PALETTE[i%len(PALETTE)],
                hovertemplate="Length: %{x:,} bp<br>Count: %{y}<extra>"+sample+"</extra>"))
    fig_viral_len.update_layout(barmode="overlay",title="Viral Consensus Contig Length Distribution",
        xaxis_title="Length (bp)",yaxis_title="Count",template=T,height=420)
    
    fig_vrhyme=go.Figure()
    fig_vrhyme.add_trace(go.Bar(name="Consensus contigs (input)",x=samples,
        y=[support_data[s].get(3,0)+support_data[s].get(4,0) for s in samples],marker_color=CYAN))
    fig_vrhyme.add_trace(go.Bar(name="vMAGs formed",x=samples,
        y=[vrhyme_data[s]["n_bins"] for s in samples],marker_color=TEAL))
    fig_vrhyme.add_trace(go.Bar(name="Contigs binned",x=samples,
        y=[vrhyme_data[s]["total_members"] for s in samples],marker_color=AMBER))
    fig_vrhyme.update_layout(barmode="group",title="vRhyme — vMAG Summary",
        yaxis_title="Count",template=T,height=400)
    
    fig_vrhyme_detail=go.Figure()
    for i,sample in enumerate(samples):
        rows=vrhyme_data[sample].get("rows",[])
        if rows:
            fig_vrhyme_detail.add_trace(go.Scatter(
                x=[safe_int(r.get("members",0)) for r in rows],
                y=[safe_int(r.get("proteins",0)) for r in rows],
                mode="markers",name=sample,
                marker=dict(size=[max(10,safe_int(r.get("redundancy",0))/2+10) for r in rows],
                            color=PALETTE[i%len(PALETTE)],opacity=0.8,line=dict(width=1,color="white")),
                text=[f"bin {r.get('bin','')} ({sample})" for r in rows],
                hovertemplate="<b>%{text}</b><br>Members: %{x}<br>Proteins: %{y}<extra></extra>"))
    fig_vrhyme_detail.update_layout(title="vRhyme — Per-Bin Detail (size = redundancy %)",
        xaxis_title="Member contigs",yaxis_title="Proteins",template=T,height=420)
    
    # ─ Bins ───────────────────────────────────────────────────────────────────────
    all_binners=["MetaBAT2","VAMB","SemiBin2","Binette (final)"]
    fig_binner_total=go.Figure()
    for sample in samples:
        fig_binner_total.add_trace(go.Bar(name=sample,x=all_binners,
            y=[binner_counts[sample].get(b,{}).get("total",0) for b in all_binners],
            hovertemplate="<b>%{x}</b><br>"+sample+": %{y} bins<extra></extra>"))
    fig_binner_total.update_layout(barmode="group",title="Total Bins per Tool",
        yaxis_title="Number of bins",legend_title="Sample",template=T,height=420)
    
    fig_das_tax=go.Figure()
    for cat,color in [("bacteria",TEAL),("archaea",AMBER),("unknown",GRAY)]:
        fig_das_tax.add_trace(go.Bar(name=cat.capitalize(),x=samples,
            y=[binner_counts[s].get("Binette (final)",{}).get(cat,0) for s in samples],
            marker_color=color,hovertemplate="<b>%{x}</b><br>"+cat+": %{y}<extra></extra>"))
    fig_das_tax.update_layout(barmode="stack",title="Binette Final Bins — Domain Classification (GTDB-Tk)",
        yaxis_title="Bins",legend_title="Domain",template=T,height=400)
    
    fig_checkm2=go.Figure()
    dom_colors={"Bacteria":TEAL,"Archaea":AMBER,"Unknown":GRAY}; added=set()
    for sample in samples:
        for row in checkm2_data[sample]:
            try:
                comp=safe_float(row.get("Completeness",0)); cont=safe_float(row.get("Contamination",0))
                tax=row.get("Taxonomic_lineage",row.get("taxonomic_lineage",""))
                name=row.get("Name",row.get("name",""))
                dom=("Archaea" if "Archaea" in tax else "Bacteria" if "Bacteria" in tax or "bacteria" in tax else "Unknown")
                show=dom not in added
                if show: added.add(dom)
                fig_checkm2.add_trace(go.Scatter(x=[cont],y=[comp],mode="markers",
                    name=dom,legendgroup=dom,showlegend=show,
                    marker=dict(size=10,color=dom_colors[dom],line=dict(width=1,color="white")),
                    text=[f"{name} ({sample})"],
                    hovertemplate="<b>%{text}</b><br>Comp: %{y:.1f}%<br>Cont: %{x:.1f}%<extra></extra>"))
            except: continue
    fig_checkm2.add_hline(y=90,line_dash="dash",line_color=GREEN,annotation_text="≥90% (HQ)",annotation_position="right")
    fig_checkm2.add_hline(y=50,line_dash="dot",line_color=AMBER,annotation_text="≥50% (MQ)",annotation_position="right")
    fig_checkm2.add_vline(x=5,line_dash="dash",line_color=RED,annotation_text="≤5% contam.")
    fig_checkm2.add_vline(x=10,line_dash="dot",line_color=AMBER,annotation_text="≤10% (MQ)",annotation_position="top right")
    fig_checkm2.update_layout(
        title="Bin Quality — Completeness vs Contamination (CheckM2)",
        xaxis_title="Contamination (%)",yaxis_title="Completeness (%)",
        xaxis=dict(range=[-1,15]),yaxis=dict(range=[-5,105]),
        legend_title="Domain",template=T,height=500,
        shapes=[
            dict(type="rect",x0=-1,y0=90,x1=5,y1=106,
                 fillcolor="rgba(16,185,129,0.10)",line=dict(width=0),layer="below"),
            dict(type="rect",x0=-1,y0=50,x1=10,y1=90,
                 fillcolor="rgba(245,158,11,0.07)",line=dict(width=0),layer="below"),
            dict(type="rect",x0=-1,y0=-6,x1=15,y1=50,
                 fillcolor="rgba(239,68,68,0.05)",line=dict(width=0),layer="below"),
        ],
        annotations=[
            dict(x=2,y=102,text="HQ zone",showarrow=False,
                 font=dict(size=9,color="rgba(16,185,129,0.7)"),xanchor="center"),
            dict(x=4,y=88,text="MQ zone",showarrow=False,
                 font=dict(size=9,color="rgba(245,158,11,0.7)"),xanchor="center"),
        ]
    )
    
    fig_cm2_hist=make_subplots(rows=1,cols=2,horizontal_spacing=0.12,
        subplot_titles=["Completeness","Contamination"])
    for i,sample in enumerate(samples):
        comps=[safe_float(r.get("Completeness",0)) for r in checkm2_data[sample]]
        conts=[safe_float(r.get("Contamination",0)) for r in checkm2_data[sample]]
        c=PALETTE[i%len(PALETTE)]
        if comps:
            fig_cm2_hist.add_trace(go.Histogram(x=comps,name=sample,opacity=0.7,nbinsx=20,marker_color=c),row=1,col=1)
        if conts:
            fig_cm2_hist.add_trace(go.Histogram(x=conts,name=sample,opacity=0.7,nbinsx=20,marker_color=c,showlegend=False),row=1,col=2)
    fig_cm2_hist.update_xaxes(title_text="Completeness (%)",row=1,col=1)
    fig_cm2_hist.update_xaxes(title_text="Contamination (%)",row=1,col=2)
    fig_cm2_hist.update_layout(barmode="overlay",title="Bin Quality Distributions",height=400,template=T)
    
    fig_bin_size=go.Figure(); any_size=False
    for i,sample in enumerate(samples):
        sizes=[safe_int(row.get("Genome_Size","0"))/1e6 for row in checkm2_data[sample] if row.get("Genome_Size","")]
        if sizes:
            any_size=True
            fig_bin_size.add_trace(go.Histogram(x=sizes,name=sample,opacity=0.75,nbinsx=30,
                marker_color=PALETTE[i%len(PALETTE)],
                hovertemplate="Size: %{x:.1f} Mb<br>Count: %{y}<extra>"+sample+"</extra>"))
    if not any_size:
        fig_bin_size.add_annotation(text="Genome_Size not in CheckM2 output",xref="paper",yref="paper",x=0.5,y=0.5,showarrow=False)
    fig_bin_size.update_layout(barmode="overlay",title="Bin Genome Size Distribution",
        xaxis_title="Estimated size (Mb)",yaxis_title="Count",template=T,height=400)
    
    # ─ Abundance ──────────────────────────────────────────────────────────────────
    all_bin_ids=[]
    for s in samples:
        for bid in sorted(abundance_data[s]):
            lbl=f"{s} | {bid}"
            if lbl not in all_bin_ids: all_bin_ids.append(lbl)
    z_matrix=[]; bin_labels=[]
    for lbl in all_bin_ids:
        src=lbl.split(" | ")[0]; bid=lbl.split(" | ",1)[1]
        z_matrix.append([abundance_data[s].get(bid,0.0) if s==src else 0.0 for s in samples])
        bin_labels.append(lbl)
    if z_matrix:
        fig_abundance=go.Figure(go.Heatmap(z=z_matrix,x=samples,y=bin_labels,colorscale="Viridis",
            colorbar=dict(title="Mean depth (×)"),
            hovertemplate="<b>%{y}</b><br>%{x}: %{z:.2f}×<extra></extra>"))
        fig_abundance.update_layout(title="Bin Abundance",height=max(450,len(bin_labels)*14+100),template=T,
            yaxis=dict(tickfont=dict(size=8),autorange="reversed"))
    else:
        fig_abundance=go.Figure()
        fig_abundance.add_annotation(text="No bin abundance data.",xref="paper",yref="paper",x=0.5,y=0.5,showarrow=False)
        fig_abundance.update_layout(height=300,template=T)
    
    # ─ Taxonomy ───────────────────────────────────────────────────────────────────

    all_phyla=Counter(); phylum_counts={}
    for s in samples:
        phylum_counts[s]=checkm2_tax_data[s]; all_phyla.update(checkm2_tax_data[s])
    top_phyla=[p for p,_ in all_phyla.most_common(20) if p!="Unclassified"]
    fig_tax_bac=go.Figure()
    for s in samples:
        fig_tax_bac.add_trace(go.Bar(name=s,x=top_phyla,y=[phylum_counts[s].get(p,0) for p in top_phyla]))
    fig_tax_bac.update_layout(barmode="group",title="⚠️ Preliminary Prokaryote Taxonomy — CheckM2 lineage (phylum, top 20)",
        xaxis_tickangle=-35,template=T,height=450)
    
    # ══════════════════════════════════════════════════════════════════════════════
    #  OVERVIEW
    # ══════════════════════════════════════════════════════════════════════════════
    overview={}
    for sample in samples:
        qd=quast_data[sample].get("deduplicated",{}); cm=checkm2_data[sample]
        cv=checkv_data[sample]; sp=support_data[sample]; fq=fastp_data[sample]["reads"]
        raw_r1=next((r for r in fq if r["stage"]=="raw"     and r["read"]=="R1"),None)
        tr_r1 =next((r for r in fq if r["stage"]=="trimmed" and r["read"]=="R1"),None)
        das=binner_counts[sample].get("Binette (final)",{})
        overview[sample]={
            "total_raw_reads": raw_r1["total_sequences"] if raw_r1 else fastp_data[sample]["trim"].get("reads_in","N/A"),
            "mean_qual":  f"{tr_r1['mean_quality']:.1f}" if tr_r1 else "N/A",
            "gc_pct":     f"{tr_r1['gc_percent']:.1f}%"  if tr_r1 else "N/A",
            "mapping_rate": f"{mapping_data[sample]:.1f}%" if mapping_data[sample]>0 else "N/A",
            "n_contigs":  qd.get("# contigs","N/A"),
            "n50":        qd.get("N50","N/A"),
            "viral_consensus": (lambda p: sum(1 for l in open(p) if l.startswith('>')) if os.path.exists(p) else
                                sum(cnt for n, cnt in sp.items() if n >= 2))(
                                os.path.join(outdir, sample, "viral", "consensus",
                                             f"{sample}_viral_consensus.fasta")),
            "complete_viral":  sum(1 for r in cv if r.get("checkv_quality","")=="Complete"),
            "vmags":      vrhyme_data[sample]["n_bins"],
            "unbinned_viral": (lambda p: sum(1 for l in open(p) if l.startswith('>')) - vrhyme_data[sample]["total_members"]
                               if os.path.exists(p) else 0)(
                               os.path.join(outdir, sample, "viral", "consensus",
                                            f"{sample}_viral_nonredundant.fasta")),
            "total_bins": das.get("total",0),
            "hq_bins":    sum(1 for r in cm if safe_float(r.get("Completeness",0))>=90
                              and safe_float(r.get("Contamination",100))<=5),
            "bacteria_bins":     das.get("bacteria",0),
            "archaea_bins":      das.get("archaea",0),
            "taxonomy_classified": 0,  # filled after tax_data load below
            "host_pred_total":   0,  # filled after phist_data load below
        }
    
    # ══════════════════════════════════════════════════════════════════════════════
    #  SERIALIZE
    # ══════════════════════════════════════════════════════════════════════════════
    figs_json={
        "fq_reads":       fig_fq_reads.to_json(),
        "fq_qual":        fig_fq_qual.to_json(),
        "fq_gc":          fig_fq_gc.to_json(),
        "trim":           fig_trim.to_json(),
        "mapping":        fig_mapping.to_json(),
        "asm_prog":       fig_asm_prog.to_json(),
        "assembly":       fig_assembly.to_json(),
        "contig_len":     fig_contig_len.to_json(),
        "coverage":       fig_cov.to_json(),
        "tool_counts":    fig_tool_counts.to_json(),
        "viral":          fig_viral.to_json(),
        "scores":         fig_scores.to_json(),
        "checkv":         fig_checkv.to_json(),
        "checkv_vrh":     fig_checkv_vrh.to_json(),
        "checkv_scatter": fig_checkv_scatter.to_json(),
        "viral_len":      fig_viral_len.to_json(),
        "vrhyme":         fig_vrhyme.to_json(),
        "vrhyme_detail":  fig_vrhyme_detail.to_json(),
        "checkm2":        fig_checkm2.to_json(),
        "cm2_hist":       fig_cm2_hist.to_json(),
        "bin_size":       fig_bin_size.to_json(),
        "binner_total":   fig_binner_total.to_json(),
        "das_tax":        fig_das_tax.to_json(),
        "abundance":      fig_abundance.to_json(),
        "tax_bac":        fig_tax_bac.to_json(),
    }
    samples_json     = json.dumps(samples).replace("</", "<\\/")
    overview_json    = json.dumps(overview).replace("</", "<\\/")
    # figs_json_str serialized later after Task-5 figures are created
    sample_list      = ", ".join(samples)
    consensus_mode   = str(getattr(snakemake.params, 'consensus_mode',    'count'))
    min_viral_tools  = str(getattr(snakemake.params, 'min_viral_tools',   '2'))
    
    # ══════════════════════════════════════════════════════════════════════════════
    #  LOAD TAXONOMY + HOST PREDICTION DATA
    # ══════════════════════════════════════════════════════════════════════════════

    import csv as _csv_mod

    def load_tsv(path, skip_empty=True):
        rows = []
        try:
            with open(path) as f:
                rdr = _csv_mod.DictReader(f, delimiter='\t')
                for row in rdr:
                    if skip_empty and not any(row.values()): continue
                    rows.append(row)
        except Exception:
            pass
        return rows

    def load_csv(path):
        rows = []
        try:
            with open(path) as f:
                rdr = _csv_mod.DictReader(f)
                for row in rdr:
                    rows.append(row)
        except Exception:
            pass
        return rows

    def load_vibrant(outdir_base, samples):
        """Parse VIBRANT summary + AMG results for each sample."""
        scaffold_records = []
        amg_records = []
        for s in samples:
            vdir = os.path.join(outdir_base, s, 'viral', 'vibrant')
            # Scaffold summary
            for tsv in glob.glob(os.path.join(vdir, '**', 'VIBRANT_summary_results_*.tsv'), recursive=True):
                for row in load_tsv(tsv):
                    scaffold = row.get('scaffold','')
                    if not scaffold: continue
                    scaffold_records.append({'sample': s,
                        'Scaffold':    scaffold,
                        'Total_genes': row.get('total genes',''),
                        'VOG_score':   row.get('VOG v-score',''),
                        'Pfam_score':  row.get('Pfam v-score',''),
                        'KEGG_score':  row.get('KEGG v-score',''),
                    })
            # AMG pathways — VIBRANT_AMG_pathways_*.tsv (pathway-level summary)
            for tsv in glob.glob(os.path.join(vdir, '**', 'VIBRANT_AMG_pathways_*.tsv'), recursive=True):
                for row in load_tsv(tsv):
                    pathway = row.get('Pathway','') or row.get('pathway','')
                    if not pathway: continue
                    n_amgs = row.get('Total AMGs','') or row.get('total amgs','1')
                    kos    = row.get('Present AMG KOs','') or row.get('present amg kos','')
                    amg_records.append({'sample':    s,
                        'Pathway':    pathway,
                        'Metabolism': row.get('Metabolism','') or row.get('metabolism',''),
                        'KEGG_map':   row.get('KEGG Entry','') or row.get('kegg entry',''),
                        'Total_AMGs': n_amgs,
                        'KOs':        kos,
                    })
            # AMG individuals — VIBRANT_AMG_individuals_*.tsv (per-protein detail)
            for tsv in glob.glob(os.path.join(vdir, '**', 'VIBRANT_AMG_individuals_*.tsv'), recursive=True):
                for row in load_tsv(tsv):
                    prot = row.get('protein','') or row.get('gene','')
                    if not prot: continue
                    amg_records.append({'sample':    s,
                        'Protein':    prot,
                        'Scaffold':   row.get('scaffold',''),
                        'AMG_KO':     row.get('AMG KO','') or row.get('KO',''),
                        'AMG_gene':   row.get('AMG gene','') or row.get('gene name',''),
                        'Pathway':    row.get('metabolic pathway','') or row.get('pathway',''),
                        'Metabolism': row.get('AMG category','') or row.get('category',''),
                        'Total_AMGs': '1',
                        'KOs':        row.get('AMG KO',''),
                    })
        return scaffold_records, amg_records

    def load_vcontact3(paths_d, samples):
        """Load vConTACT3 final_assignments.csv (v3 output format)."""
        records = []
        for s in samples:
            p = paths_d.get(s,'')
            # v3 uses CSV with columns: Genome, family_prediction, genus_prediction, etc.
            rows = load_csv(p) if p.endswith('.csv') else load_tsv(p)
            for row in rows:
                genome = row.get('Genome', row.get('genome',''))
                if not genome: continue
                # Skip reference genomes (Reference == True)
                if row.get('Reference','').lower() in ('true','1','yes'): continue
                fam = row.get('family_prediction', row.get('Family',''))
                gen = row.get('genus_prediction',  row.get('Genus',''))
                rlm = row.get('realm_prediction',  row.get('realm_reference',''))
                cls = row.get('class_prediction',  '')
                ord_ = row.get('order_prediction', '')
                # Best taxonomy = deepest assigned level
                best = gen or fam or ord_ or cls or rlm or ''
                # Novel flag
                is_novel = any('novel' in str(v).lower()
                               for v in [fam, gen, ord_] if v)
                records.append({'sample':   s,
                    'Genome':       genome,
                    'VC':           row.get('VC', ''),
                    'VC_status':    'novel' if is_novel else ('classified' if best else 'singleton'),
                    'Family':       fam,
                    'Genus':        gen,
                    'Order':        ord_,
                    'Realm':        rlm,
                    'Best_taxonomy': best,
                    'Is_novel':     str(is_novel),
                })
        return records

    def load_viral_taxonomy(paths, samples):
        RANKS = ['realm','kingdom','phylum','class','order','family','genus','species']
        def deepest_level(lineage):
            """Extract deepest non-empty taxonomic level from GeNomad lineage string."""
            if not lineage: return '', '', ''
            parts = [p.strip() for p in lineage.split(';')]
            # GeNomad lineage: Viruses;Realm;Kingdom;Phylum;Class;Order;Family;Genus
            # Find deepest non-empty, non-generic part
            family = ''; genus = ''; order = ''
            if len(parts) >= 8 and parts[7]: genus  = parts[7]
            if len(parts) >= 7 and parts[6]: family = parts[6]
            if len(parts) >= 6 and parts[5]: order  = parts[5]
            # If family/genus empty, return deepest available level as "best_level"
            return family, genus, order

        records = []
        for p, s in zip(paths, samples):
            for row in load_tsv(p):
                name = row.get('seq_name','')
                if not name: continue
                # Get stored values
                final_family = row.get('final_family','')
                final_genus  = row.get('final_genus','')
                final_order  = row.get('final_order','')
                lineage      = row.get('lineage','')
                source       = row.get('source','')
                # If family/genus empty, try to extract from lineage
                if not final_family or not final_genus:
                    lf, lg, lo = deepest_level(lineage)
                    if not final_family: final_family = lf
                    if not final_genus:  final_genus  = lg
                    if not final_order:  final_order  = lo
                # Build best_taxonomy: deepest assigned level for display
                best_tax = final_genus or final_family or final_order
                if not best_tax and lineage:
                    parts = [p.strip() for p in lineage.split(';') if p.strip() and p.strip() != 'Viruses']
                    best_tax = parts[-1] if parts else ''
                records.append({'sample': s,
                    'Genome':        name,
                    'final_family':  final_family,
                    'final_genus':   final_genus,
                    'final_order':   final_order,
                    'Family':        final_family,
                    'Genus':         final_genus,
                    'Order':         final_order,
                    'Best_taxonomy': row.get('best_taxonomy', best_tax),
                    'GeNomad_best':  row.get('genomad_best', ''),
                    'GeNomad_class': row.get('genomad_class', ''),
                    'Source':        source,
                    'Confidence':    row.get('confidence',''),
                    'Lineage':       lineage,
                    'INPHARED_name': row.get('inphared_name',''),
                })
        return records

    def load_gtdbtk(bac_paths, arc_paths, samples):
        records = []
        for bac_p, arc_p, s in zip(bac_paths, arc_paths, samples):
            for fpath in [bac_p, arc_p]:
                for row in load_tsv(fpath):
                    classif = row.get('classification','')
                    if not classif or classif in ('N/A','NA',''): continue
                    def _g(prefix, c=classif):
                        if f';{prefix}__' in c:
                            v = c.split(f';{prefix}__')[-1].split(';')[0]
                            return v.strip() if v.strip() else ''
                        return ''
                    domain = ('Archaea' if 'd__Archaea' in classif
                              else 'Bacteria' if 'd__Bacteria' in classif
                              else 'Unknown')
                    bin_name = row.get('user_genome','').strip()
                    if not bin_name: continue
                    records.append({'sample':s, 'Bin':bin_name,
                        'Domain':domain, 'Phylum':_g('p'), 'Class':_g('c'),
                        'Order':_g('o'), 'Family':_g('f'), 'Genus':_g('g'),
                        'Species':_g('s'), 'Full_classification':classif,
                        'RED_value':row.get('red_value',''), 'Note':row.get('note',''),
                    })
        return records

    def load_custom_prok(paths_d, samples, meta_path=''):
        """Load custom Diamond prokaryote hits joined with metadata TSV.
        Uses sscinames/sskingdoms from Diamond output as fallback when metadata unavailable."""
        from collections import Counter as _Counter, defaultdict as _dd
        meta = {}
        if meta_path and os.path.exists(str(meta_path)):
            for row in load_tsv(str(meta_path)):
                acc = row.get('accession','').strip()
                if acc: meta[acc] = row
        records = []
        for s in samples:
            p = paths_d.get(s,'')
            if not p or not os.path.exists(p): continue
            votes      = _dd(_Counter)
            pidents    = _dd(list)
            sscinames  = _dd(list)
            sskingdoms = {}
            for row in load_tsv(p):
                q = row.get('qseqid',''); h = row.get('sseqid','')
                if not q or not h: continue
                bin_name = '_'.join(q.split('_')[:-1]) or q
                acc      = '_'.join(h.split('_')[:-1]) or h
                votes[bin_name][acc] += 1
                try: pidents[bin_name].append(float(row.get('pident',0)))
                except: pass
                ss = row.get('sscinames','')
                if ss and ss not in ('N/A', '0', ''):
                    sscinames[bin_name].append(ss)
                sk = row.get('sskingdoms','')
                if sk and sk not in ('N/A', '0', ''):
                    sskingdoms[bin_name] = sk
            for bin_name, v in votes.items():
                top_acc, _ = v.most_common(1)[0]
                m   = meta.get(top_acc, {})
                avg = sum(pidents[bin_name])/len(pidents[bin_name]) if pidents[bin_name] else 0
                sci_fb = sscinames[bin_name][0] if sscinames[bin_name] else ''
                records.append({'sample': s, 'Bin': bin_name,
                    'Phylum':   m.get('phylum',''),  'Class':    m.get('class',''),
                    'Order':    m.get('order',''),   'Family':   m.get('family',''),
                    'Genus':    m.get('genus',''),
                    'Organism': m.get('organism','') or sci_fb,
                    'Domain':   m.get('domain','') or sskingdoms.get(bin_name,''),
                    'Top_acc':  top_acc,
                    'Avg_pident': f"{avg:.1f}",      'Source':   'diamond_custom',
                    'Sci_name': sci_fb,
                })
        return records

    def load_phist(paths, samples):
        records = []
        for p, s in zip(paths, samples):
            # PHIST output: phage,host,#common-kmers,pvalue,adj-pvalue
            # phage column = path to fasta file (e.g. viral_fastas/contig_MEGA_k141_xxx.fasta)
            for row in load_csv(p):
                virus_raw = row.get('phage', row.get('Phage', row.get('virus', row.get('Virus', ''))))
                host_raw  = row.get('host',  row.get('Host', ''))
                score     = row.get('#common-kmers', row.get('Score', row.get('score', '')))
                pval      = row.get('adj-pvalue', row.get('pvalue', ''))
                if not virus_raw or not host_raw: continue
                # Clean names: strip paths, extensions, prefixes added during splitting
                virus_clean = os.path.basename(virus_raw).replace('.fasta','').replace('.fa','')
                virus_clean = virus_clean.replace('contig_','',1)   # remove split prefix
                host_clean  = os.path.basename(host_raw).replace('.fa','').replace('.fasta','')
                records.append({'sample':s, 'Virus':virus_clean, 'Host':host_clean,
                                'Score':str(score), 'P_value':str(pval)})
        return records

    # ── Load taxonomy + host prediction data ──────────────────────────────────
    taxonomy_paths_d = _path_dict(list(snakemake.input.taxonomy), samples)
    # vcontat3 input points to genome_clusters.tsv but real data is in exports/
    vcontact3_raw    = _path_dict(list(snakemake.input.vcontact3), samples)
    vcontact3_paths  = {}
    for s in samples:
        raw = vcontact3_raw.get(s,'')
        # Try exports/final_assignments.csv first (vConTACT3 v3 actual output)
        if raw:
            exports_csv = os.path.join(os.path.dirname(raw), 
                                       'vConTACT3_results', 'exports', 'final_assignments.csv')
            vcontact3_paths[s] = exports_csv if os.path.exists(exports_csv) else raw
        else:
            vcontact3_paths[s] = raw
    gtdbtk_bac_l     = list(snakemake.input.gtdbtk_bac)
    gtdbtk_arc_l     = list(snakemake.input.gtdbtk_arc)
    phist_paths_l    = _path_dict(list(snakemake.input.phist), samples)
    phist_paths_l    = [phist_paths_l.get(s,"") for s in samples]

    tax_data      = load_viral_taxonomy([taxonomy_paths_d.get(s,'') for s in samples], samples)
    custom_prok_paths = _path_dict(list(snakemake.input.custom_prok), samples)
    custom_prok_data  = load_custom_prok(
        custom_prok_paths, samples,
        getattr(snakemake.params, 'custom_prok_meta', ''))
    vcontact3_data = load_vcontact3(vcontact3_paths, samples)
    vibrant_data, vibrant_amg = load_vibrant(snakemake.params.outdir, samples)
    gtdb_data  = load_gtdbtk(gtdbtk_bac_l, gtdbtk_arc_l, samples)
    phist_data = load_phist(phist_paths_l, samples)

    print(f"[generate_report] Viral taxonomy: {len(tax_data)} classified viruses")
    print(f"[generate_report] vConTACT3: {len(vcontact3_data)} genomes clustered")
    print(f"[generate_report] VIBRANT: {len(vibrant_data)} scaffolds, {len(vibrant_amg)} AMGs")
    print(f"[generate_report] GTDB-Tk: {len(gtdb_data)} MAGs")
    print(f"[generate_report] Custom prok: {len(custom_prok_data)} bins")
    print(f"[generate_report] PHIST: {len(phist_data)} predictions")

    # Novelty metrics per sample
    novelty_data = {}
    for sample in samples:
        total_viral = sum(1 for r in tax_data if r.get('sample')==sample)
        unclassified = sum(1 for r in tax_data if r.get('sample')==sample
                          and r.get('Source','') == 'unclassified')
        classified = total_viral - unclassified
        pct_novel = round(100.0 * unclassified / total_viral, 1) if total_viral > 0 else 0.0
        novelty_data[sample] = {
            'total': total_viral, 'classified': classified,
            'unclassified': unclassified, 'pct_novel': pct_novel
        }

    # MIMAG quality tiers for MAGs
    mimag_data = {}
    for sample in samples:
        hq = mq = lq = 0
        for row in checkm2_data[sample]:
            comp = safe_float(row.get('Completeness',0))
            cont = safe_float(row.get('Contamination',100))
            if comp >= 90 and cont <= 5:   hq += 1
            elif comp >= 50 and cont <= 10: mq += 1
            else:                           lq += 1
        mimag_data[sample] = {'HQ': hq, 'MQ': mq, 'LQ': lq, 'total': hq+mq+lq}

    # Back-fill host predictions (PHIST unique viruses per sample)
    for sample in samples:
        phist_viruses = {r["Virus"] for r in phist_data if r.get("sample")==sample and r.get("Virus")}
        overview[sample]["host_pred_total"] = len(phist_viruses)

    # Back-fill Binette domain from GTDB-Tk
    for sample in samples:
        gb = set(r["Bin"] for r in gtdb_data if r.get("sample")==sample and r.get("Domain")=="Bacteria")
        ga = set(r["Bin"] for r in gtdb_data if r.get("sample")==sample and r.get("Domain")=="Archaea")
        bd = os.path.join(outdir, sample, "bins", "binette", "final_bins")
        das = {"total":0,"bacteria":0,"archaea":0,"unknown":0}
        for bf in glob.glob(os.path.join(bd,"*.fa")):
            name = os.path.basename(bf).replace(".fa","")
            das["total"] += 1
            if name in ga:   das["archaea"]  += 1
            elif name in gb: das["bacteria"] += 1
            else:            das["unknown"]  += 1
        binner_counts[sample]["Binette (final)"] = das
        overview[sample].update({"bacteria_bins":das["bacteria"],
                                  "archaea_bins":das["archaea"],
                                  "total_bins":das["total"]})

    # ── Merge vConTACT3 records into unified tax_data ─────────────────────────
    tax_genome_keys = {(r.get('sample',''), r.get('Genome','')) for r in tax_data}
    for vc_row in vcontact3_data:
        key = (vc_row.get('sample',''), vc_row.get('Genome',''))
        if key not in tax_genome_keys:
            tax_data.append({
                'sample':        vc_row['sample'],
                'Genome':        vc_row['Genome'],
                'final_family':  vc_row.get('Family',''),
                'final_genus':   vc_row.get('Genus',''),
                'Family':        vc_row.get('Family',''),
                'Genus':         vc_row.get('Genus',''),
                'Order':         vc_row.get('Order',''),
                'Best_taxonomy': vc_row.get('Best_taxonomy',''),
                'Source':        'vcontact3',
                'Confidence':    '',
                'Lineage':       '',
                'Completeness':  '',
                'Genome_length': '',
                'CheckV_quality':'',
            })
            tax_genome_keys.add(key)

    # Back-fill taxonomy_classified: sequences with at least family- or genus-level assignment.
    # Uses Best_taxonomy (deepest assigned rank) so sequences with empty classification
    # (Source set but no actual rank assigned) are not counted.
    for sample in samples:
        overview[sample]["taxonomy_classified"] = sum(
            1 for r in tax_data
            if r.get("sample") == sample
            and (r.get("Best_taxonomy","") or r.get("final_family","") or r.get("final_genus","")))

    # ── Task 2: Enrich viral taxonomy with CheckV data ────────────────────────
    def enrich_taxonomy_with_checkv(tax_records, checkv_dict):
        cv_lookup = {}
        for s, rows in checkv_dict.items():
            cv_lookup[s] = {r.get('contig_id', r.get('contig','')): r for r in rows}
        for rec in tax_records:
            cv_row = cv_lookup.get(rec.get('sample',''), {}).get(rec.get('Genome',''), {})
            rec['Completeness']   = cv_row.get('completeness','')
            rec['Genome_length']  = cv_row.get('contig_length','')
            rec['CheckV_quality'] = cv_row.get('checkv_quality','')
        return tax_records

    tax_data = enrich_taxonomy_with_checkv(tax_data, checkv_data)

    # ── Task 3: Merge prokaryotic taxonomy (GTDB-Tk priority + Diamond fallback)
    def merge_prok_taxonomy(gtdb_records, custom_prok_records, checkm2_dict):
        gtdb_bins   = {(r.get('sample',''), r.get('Bin','').replace('.fa','')): r
                       for r in gtdb_records}
        custom_bins = {(r.get('sample',''), r.get('Bin','').replace('.fa','')): r
                       for r in custom_prok_records}
        merged = []
        for sample, rows in checkm2_dict.items():
            for cm_row in rows:
                bin_name = cm_row.get('Name', cm_row.get('name','')).replace('.fa','')
                if not bin_name: continue
                key = (sample, bin_name)
                if key in gtdb_bins:
                    base = dict(gtdb_bins[key])
                    base['Source_tax'] = 'GTDB-Tk'
                elif key in custom_bins:
                    base = dict(custom_bins[key])
                    base['Bin'] = bin_name
                    base['sample'] = sample
                    base['Source_tax'] = 'Diamond-Custom'
                else:
                    base = {'sample': sample, 'Bin': bin_name,
                            'Domain':'','Phylum':'','Class':'','Order':'',
                            'Family':'','Genus':'','Species':'',
                            'Source_tax':'Unclassified'}
                base['Completeness']  = cm_row.get('Completeness','')
                base['Contamination'] = cm_row.get('Contamination','')
                base['Genome_size']   = cm_row.get('Genome_Size', cm_row.get('genome_size',''))
                merged.append(base)
        return merged

    merged_prok_data = merge_prok_taxonomy(gtdb_data, custom_prok_data, checkm2_data)
    print(f"[generate_report] Merged prok taxonomy: {len(merged_prok_data)} MAGs")

    # ── Task 5: Tool support matrix for heatmap ───────────────────────────────
    def load_tool_support_matrix(outdir_base, samples_list):
        tools_order = ["VirSorter2","GeNomad","VIBRANT"]
        rows = []
        for s in samples_list:
            for p in glob.glob(os.path.join(outdir_base, s, "viral", "consensus",
                                            "*_tool_support.tsv")):
                for row in parse_tsv(p):
                    contig    = row.get('contig', row.get('contig_id',''))
                    tool_str  = row.get('tools','')
                    n_tools   = safe_int(row.get('n_tools', 0))
                    if not contig: continue
                    rec = {'sample': s, 'contig': contig, 'n_tools': n_tools}
                    for t in tools_order:
                        rec[t] = 1 if t in tool_str else 0
                    rows.append(rec)
        return rows

    tool_matrix_data = load_tool_support_matrix(outdir, samples)

    # ── Task 5: Read utilization funnel ──────────────────────────────────────
    fig_read_funnel = go.Figure()
    funnel_stages = ["Raw reads", "Post-trimming", "Mapped to assembly",
                     "Viral consensus (≥3 tools)", "Prokaryotic MAG bins"]
    for i, sample in enumerate(samples):
        raw     = safe_int(fastp_data[sample]["trim"].get('reads_in', 0))
        trimmed = safe_int(fastp_data[sample]["trim"].get('reads_written', 0))
        mapped  = int(raw * mapping_data[sample] / 100) if mapping_data[sample] > 0 and raw > 0 else 0
        viral_c = support_data[sample].get(3,0) + support_data[sample].get(4,0)
        bin_c   = sum(1 for r in checkm2_data[sample])
        # estimate reads in viral contigs / bins as proxy (per-contig depth × len / read_len)
        viral_reads_proxy = viral_c * 5000 // 150 if viral_c else 0
        bin_reads_proxy   = bin_c   * 20000 // 150 if bin_c   else 0
        funnel_vals = [raw, trimmed, mapped, viral_reads_proxy, bin_reads_proxy]
        if any(v > 0 for v in funnel_vals):
            fig_read_funnel.add_trace(go.Funnel(
                name=sample,
                y=funnel_stages,
                x=funnel_vals,
                textinfo="value+percent previous",
                marker=dict(color=PALETTE[i % len(PALETTE)]),
                opacity=0.85,
            ))
    fig_read_funnel.update_layout(
        title="Read Utilization — From Raw to Assigned",
        height=480, template=T,
        annotations=[dict(
            text="Note: viral/bin read counts are approximate (contig-based proxy)",
            xref="paper", yref="paper", x=0.5, y=-0.08,
            showarrow=False, font=dict(size=10, color=GRAY)
        )]
    )

    # ── Task 5: Tool agreement heatmap ────────────────────────────────────────
    fig_tool_heatmap = go.Figure()
    if tool_matrix_data:
        _tools_ord = ["VirSorter2","GeNomad","VIBRANT"]
        _z = [[r.get(t, 0) for t in _tools_ord] for r in tool_matrix_data]
        _contigs = [r['contig'] for r in tool_matrix_data]
        # show at most 100 contigs for readability
        _max = min(100, len(_z))
        fig_tool_heatmap = go.Figure(go.Heatmap(
            z=_z[:_max],
            x=_tools_ord,
            y=_contigs[:_max],
            colorscale=[[0,"#f1f5f9"],[1,TEAL]],
            showscale=False,
            hovertemplate="<b>%{y}</b><br>%{x}: %{z}<extra></extra>",
        ))
        fig_tool_heatmap.update_layout(
            title=f"Tool Agreement Heatmap (first {_max} of {len(_z)} viral contigs)",
            height=max(400, _max * 8 + 80),
            yaxis=dict(tickfont=dict(size=7), autorange="reversed"),
            template=T,
        )
    else:
        fig_tool_heatmap.add_annotation(
            text="No tool support data available",
            xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False)
        fig_tool_heatmap.update_layout(height=300, template=T)

    # ── Task 5: Viral contig depth distribution ────────────────────────────
    fig_viral_depth = go.Figure()
    for i, sample in enumerate(samples):
        # Get depth values for viral consensus contigs only
        viral_names = {r.get('contig_id', r.get('contig',''))
                       for r in checkv_data[sample]}
        depth_path = os.path.join(outdir, sample, "mapping", f"{sample}_depth.txt")
        vd = []
        if os.path.exists(depth_path):
            for row in parse_tsv(depth_path):
                name = row.get('contigName', row.get('Contig',''))
                if name in viral_names:
                    d = safe_float(row.get('totalAvgDepth', row.get('avgDepth', 0)))
                    if 0 < d < 2000:
                        vd.append(d)
        if vd:
            fig_viral_depth.add_trace(go.Histogram(
                x=vd, name=sample, opacity=0.75, nbinsx=50,
                marker_color=PALETTE[i % len(PALETTE)],
                hovertemplate="Depth: %{x:.1f}×<br>Count: %{y}<extra>" + sample + "</extra>",
            ))
    fig_viral_depth.update_layout(
        barmode="overlay",
        title="Viral Consensus Contigs — Coverage Depth Distribution",
        xaxis_title="Mean read depth (×)", yaxis_title="Count",
        template=T, height=400,
    )

    # ════════════════════════════════════════════════════════════════════════
    #  BLOCO 8 — New module data loading + figures
    # ════════════════════════════════════════════════════════════════════════

    # ── vOTU table ────────────────────────────────────────────────────────
    votu_data = {}
    for _s in samples:
        _p = os.path.join(outdir, _s, "viral", "votu", f"{_s}_vOTU_table.tsv")
        votu_data[_s] = load_tsv(_p)

    # ── Lifestyle from vOTU table ─────────────────────────────────────────
    lifestyle_data = {}
    for _s in samples:
        lytic = lysogenic = 0
        for row in votu_data[_s]:
            ls = (row.get('lifestyle', '') or '').lower()
            if 'lytic' in ls or 'virulent' in ls:
                lytic += 1
            elif 'lysogenic' in ls or 'temperate' in ls:
                lysogenic += 1
        _total = len(votu_data[_s])
        lifestyle_data[_s] = {
            'lytic': lytic, 'lysogenic': lysogenic,
            'unknown': _total - lytic - lysogenic, 'total': _total,
        }

    # ── PHROGS categories (Pharokka TSV) ──────────────────────────────────
    phrogs_data = {}
    for _s in samples:
        _p = os.path.join(outdir, _s, "annotation", "pharokka",
                          "pharokka_cds_final_merged_output.tsv")
        _cnt = Counter()
        for row in load_tsv(_p):
            cat = (row.get('phrog_category', '') or row.get('category', '') or
                   'unknown function').lower().strip()
            if not cat:
                cat = 'unknown function'
            _cnt[cat] += 1
        phrogs_data[_s] = dict(_cnt)

    # ── EggNOG COG categories ─────────────────────────────────────────────
    _COG_LABEL = {
        'J': 'Translation', 'K': 'Transcription', 'L': 'Replication/Repair',
        'D': 'Cell cycle', 'M': 'Cell membrane', 'C': 'Energy production',
        'E': 'Amino acid met.', 'G': 'Carbohydrate met.', 'P': 'Ion transport',
        'T': 'Signal transd.', 'V': 'Defense', 'O': 'Protein modif.',
        'U': 'Secretion', 'N': 'Cell motility', 'S': 'Unknown function',
    }
    eggnog_data = {}
    for _s in samples:
        _p = os.path.join(outdir, _s, "annotation", "eggnog", "eggnog_annotations.tsv")
        _cnt = Counter()
        for row in load_tsv(_p):
            cog = str(row.get('COG_category', '') or row.get('cog_category', '') or 'S').strip()
            for c in cog:
                if c in _COG_LABEL:
                    _cnt[_COG_LABEL[c]] += 1
        eggnog_data[_s] = dict(_cnt)

    # ── CoverM abundance ──────────────────────────────────────────────────
    viral_abund_data = {}
    prok_abund_data  = {}
    for _s in samples:
        viral_abund_data[_s] = load_tsv(os.path.join(outdir, _s, "abundance", "viral_abundance.tsv"))
        prok_abund_data[_s]  = load_tsv(os.path.join(outdir, _s, "abundance", "prok_abundance.tsv"))

    # ── Alpha/Beta diversity + Procrustes ─────────────────────────────────
    alpha_rows       = load_tsv(os.path.join(outdir, "diversity", "alpha_diversity.tsv"))
    pcoa_viral_rows  = load_tsv(os.path.join(outdir, "diversity", "beta_pcoord_viral.tsv"))
    pcoa_prok_rows   = load_tsv(os.path.join(outdir, "diversity", "beta_pcoord_prok.tsv"))
    pcoa_comb_rows   = load_tsv(os.path.join(outdir, "diversity", "beta_pcoord_combined.tsv"))
    procrustes_rows  = load_tsv(os.path.join(outdir, "diversity", "procrustes_coords.tsv"))

    # ── Genome map SVGs ───────────────────────────────────────────────────
    def load_svg(svg_path):
        if not os.path.exists(svg_path):
            return ""
        try:
            with open(svg_path, encoding='utf-8') as _f:
                return _f.read()
        except Exception:
            return ""

    genome_maps_dict = {}
    for _s in samples:
        _gm_base = os.path.join(outdir, _s, "annotation", "genome_maps")
        genome_maps_dict[_s] = {"phage": [], "virus": [], "prok": []}
        for _mode in ("phage", "virus", "prok"):
            _mdir = os.path.join(_gm_base, _mode)
            for _svg_f in sorted(glob.glob(os.path.join(_mdir, "*.svg")))[:5]:
                _gid = os.path.basename(_svg_f).replace("_map.svg", "")
                _svg = load_svg(_svg_f)
                if _svg:
                    genome_maps_dict[_s][_mode].append({"id": _gid, "svg": _svg})

    # ════════════════════════════════════════════════════════════════════════
    #  BLOCO 8 Plotly figures
    # ════════════════════════════════════════════════════════════════════════

    _B8_PAL = [TEAL, AMBER, PURPLE, GREEN, CYAN, RED, BLUE, ORANGE,
               "#f97316", "#a78bfa", "#34d399", "#fb923c", "#e879f9", "#38bdf8"]

    # V3: Lifestyle donut per sample
    _ncols_ls = max(len(samples), 1)
    fig_lifestyle = go.Figure()
    for _i, _s in enumerate(samples):
        _ls = lifestyle_data[_s]
        fig_lifestyle.add_trace(go.Pie(
            name=_s,
            labels=['Lytic', 'Lysogenic', 'Unknown'],
            values=[_ls['lytic'], _ls['lysogenic'], _ls['unknown']],
            hole=0.42,
            domain=dict(column=_i),
            marker_colors=[RED, PURPLE, GRAY],
            hovertemplate="<b>%{label}</b><br>Count: %{value}<br>%{percent}<extra>" + _s + "</extra>",
        ))
    fig_lifestyle.update_layout(
        title="Viral Lifestyle Prediction (VIBRANT) — Lytic vs Lysogenic",
        grid=dict(rows=1, columns=_ncols_ls),
        annotations=[dict(text=_s, x=(_i + 0.5) / _ncols_ls,
                          y=-0.12, font_size=11, showarrow=False)
                     for _i, _s in enumerate(samples)],
        height=380, template=T, showlegend=True,
    )

    # V5: PHROGS category stacked bar
    _phrogs_cats = sorted({cat for d in phrogs_data.values() for cat in d})
    fig_phrogs = go.Figure()
    for _ci, _cat in enumerate(_phrogs_cats):
        fig_phrogs.add_trace(go.Bar(
            name=_cat,
            x=samples,
            y=[phrogs_data[_s].get(_cat, 0) for _s in samples],
            marker_color=_B8_PAL[_ci % len(_B8_PAL)],
            hovertemplate=f"<b>%{{x}}</b><br>{_cat}: %{{y}}<extra></extra>",
        ))
    fig_phrogs.update_layout(
        barmode='stack',
        title="PHROGS Functional Categories (Pharokka)",
        yaxis_title="Gene count", height=420, template=T,
        legend=dict(font=dict(size=10)),
    )

    # P1: COG category stacked bar
    _cog_cats = sorted({cat for d in eggnog_data.values() for cat in d})
    fig_cog = go.Figure()
    for _ci, _cat in enumerate(_cog_cats):
        fig_cog.add_trace(go.Bar(
            name=_cat,
            x=samples,
            y=[eggnog_data[_s].get(_cat, 0) for _s in samples],
            marker_color=_B8_PAL[_ci % len(_B8_PAL)],
            hovertemplate=f"<b>%{{x}}</b><br>{_cat}: %{{y}}<extra></extra>",
        ))
    fig_cog.update_layout(
        barmode='stack',
        title="COG Functional Categories (EggNOG-mapper)",
        yaxis_title="Gene count", height=420, template=T,
        legend=dict(font=dict(size=10)),
    )

    # E1-E3: Alpha diversity box plots
    _alpha_by_dom = {}
    for _row in alpha_rows:
        _dom = _row.get('domain', 'combined')
        _idx = _row.get('index', 'shannon')
        _sam = _row.get('sample', '')
        _val = safe_float(_row.get('value', 0))
        _alpha_by_dom.setdefault(_dom, {}).setdefault(_idx, []).append((_sam, _val))

    def _make_alpha_fig(domain):
        _fig = go.Figure()
        for _idx_n in ['shannon', 'simpson', 'observed', 'chao1']:
            _pts = _alpha_by_dom.get(domain, {}).get(_idx_n, [])
            if not _pts:
                continue
            _xs, _ys = zip(*_pts)
            _fig.add_trace(go.Box(
                name=_idx_n.capitalize(), x=list(_xs), y=list(_ys),
                boxpoints='all', jitter=0.3, pointpos=-1.6,
                marker=dict(size=7), line=dict(width=1.5),
            ))
        _fig.update_layout(
            title=f"Alpha Diversity — {domain.capitalize()}",
            yaxis_title="Index value", height=400, template=T, boxmode='group',
        )
        return _fig

    fig_alpha_viral    = _make_alpha_fig('viral')
    fig_alpha_prok     = _make_alpha_fig('prok')
    fig_alpha_combined = _make_alpha_fig('combined')

    # E4-E6: PCoA scatter
    def _make_pcoa_fig(pcoa_rows_in, title):
        _fig = go.Figure()
        if not pcoa_rows_in:
            _fig.add_annotation(text="No PCoA data available — run with ≥2 samples",
                                xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
                                font=dict(color=GRAY, size=13))
            _fig.update_layout(height=400, template=T)
            return _fig
        for _i, _row in enumerate(pcoa_rows_in):
            _sam = _row.get('sample', f'S{_i+1}')
            _x1 = safe_float(_row.get('PC1', _row.get('pc1', 0)))
            _x2 = safe_float(_row.get('PC2', _row.get('pc2', 0)))
            _fig.add_trace(go.Scatter(
                x=[_x1], y=[_x2], mode='markers+text',
                text=[_sam], textposition='top center',
                name=_sam, marker=dict(size=13, color=PALETTE[_i % len(PALETTE)]),
                hovertemplate=f"<b>{_sam}</b><br>PC1: {_x1:.4f}<br>PC2: {_x2:.4f}<extra></extra>",
            ))
        _fig.update_layout(
            title=title,
            xaxis_title="PC1 (Bray-Curtis)", yaxis_title="PC2",
            height=440, template=T, showlegend=False,
        )
        return _fig

    fig_pcoa_viral    = _make_pcoa_fig(pcoa_viral_rows,  "Beta Diversity PCoA — Viral Community")
    fig_pcoa_prok     = _make_pcoa_fig(pcoa_prok_rows,   "Beta Diversity PCoA — Prokaryotic Community")
    fig_pcoa_combined = _make_pcoa_fig(pcoa_comb_rows,   "Beta Diversity PCoA — Combined Community")

    # E7: Procrustes overlay
    fig_procrustes = go.Figure()
    if procrustes_rows:
        _disp = safe_float(procrustes_rows[0].get('disparity', '') if procrustes_rows else '')
        for _i, _row in enumerate(procrustes_rows):
            _sam = _row.get('sample', f'S{_i+1}')
            _vx = safe_float(_row.get('viral_PC1', _row.get('PC1_viral', 0)))
            _vy = safe_float(_row.get('viral_PC2', _row.get('PC2_viral', 0)))
            _px = safe_float(_row.get('prok_PC1',  _row.get('PC1_prok',  0)))
            _py = safe_float(_row.get('prok_PC2',  _row.get('PC2_prok',  0)))
            _col = PALETTE[_i % len(PALETTE)]
            fig_procrustes.add_trace(go.Scatter(
                x=[_vx], y=[_vy], mode='markers',
                name=f"{_sam} (viral)",
                marker=dict(size=11, symbol='circle', color=_col),
                hovertemplate=f"<b>{_sam} — viral</b><br>D1: {_vx:.4f}<br>D2: {_vy:.4f}<extra></extra>",
            ))
            fig_procrustes.add_trace(go.Scatter(
                x=[_px], y=[_py], mode='markers',
                name=f"{_sam} (prok)",
                marker=dict(size=11, symbol='square', color=_col),
                hovertemplate=f"<b>{_sam} — prok</b><br>D1: {_px:.4f}<br>D2: {_py:.4f}<extra></extra>",
            ))
            fig_procrustes.add_shape(type='line',
                x0=_vx, y0=_vy, x1=_px, y1=_py,
                line=dict(color=_col, width=1.5, dash='dot'))
        _disp_lbl = f"disparity = {_disp:.4f}" if _disp else ""
        _interp = ("strong co-variation" if _disp < 0.2 else
                   "moderate" if _disp < 0.5 else "weak co-variation")
        fig_procrustes.update_layout(
            title=f"Procrustes — Viral vs Prokaryotic Structure ({_disp_lbl} → {_interp})",
            xaxis_title="Dimension 1", yaxis_title="Dimension 2",
            height=480, template=T,
        )
    else:
        fig_procrustes.add_annotation(
            text="No Procrustes data — requires ≥2 samples with abundance data",
            xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
            font=dict(color=GRAY, size=13))
        fig_procrustes.update_layout(height=300, template=T)

    # E8-E9: Relative abundance stacked bars (top 20 features)
    def _make_abund_fig(abund_dict, label):
        _fig = go.Figure()
        _all_ids = Counter()
        _by_feat = {}
        for _s, _rows in abund_dict.items():
            for _row in _rows:
                _keys = list(_row.keys())
                if not _keys:
                    continue
                _feat_col = _keys[0]
                _feat = _row[_feat_col]
                _rpkm_col = next((k for k in _keys[1:] if 'rpkm' in k.lower()), None)
                if not _rpkm_col and len(_keys) > 1:
                    _rpkm_col = _keys[1]
                _val = safe_float(_row.get(_rpkm_col, 0)) if _rpkm_col else 0
                _all_ids[_feat] += _val
                _by_feat.setdefault(_feat, {})[_s] = _val
        _top20 = [k for k, _ in _all_ids.most_common(20)]
        if not _top20:
            _fig.add_annotation(text="No abundance data available",
                                xref="paper", yref="paper", x=0.5, y=0.5, showarrow=False,
                                font=dict(color=GRAY, size=13))
            _fig.update_layout(height=400, template=T)
            return _fig
        for _ci, _feat in enumerate(_top20):
            _fig.add_trace(go.Bar(
                name=_feat,
                x=samples,
                y=[_by_feat.get(_feat, {}).get(_s, 0) for _s in samples],
                marker_color=_B8_PAL[_ci % len(_B8_PAL)],
                hovertemplate=f"<b>%{{x}}</b><br>{_feat}: %{{y:.2f}} RPKM<extra></extra>",
            ))
        _fig.update_layout(
            barmode='stack',
            title=f"Relative Abundance — Top 20 {label} (RPKM)",
            yaxis_title="RPKM", height=460, template=T,
            legend=dict(font=dict(size=9)),
        )
        return _fig

    fig_abund_viral = _make_abund_fig(viral_abund_data, "vOTUs")
    fig_abund_prok  = _make_abund_fig(prok_abund_data,  "MAGs")

    # ── Benchmark / timing data ────────────────────────────────────────────────
    bench_rows = load_tsv(benchmark_path) if benchmark_path and os.path.exists(str(benchmark_path)) else []

    # Build matrix: rule → sample → wall_minutes
    _bench_rules   = []
    _bench_samples = []
    _bench_by      = {}   # {rule: {sample: wall_min}}
    _bench_total   = {}   # {rule: total_wall_min across all samples}
    for _row in bench_rows:
        _rule   = _row.get('rule', '')
        _sample = _row.get('sample', 'global')
        _wall_s = safe_float(_row.get('wall_s', 0))
        _wall_m = round(_wall_s / 60, 2)
        if not _rule:
            continue
        _bench_by.setdefault(_rule, {})[_sample] = _wall_m
        _bench_total[_rule] = _bench_total.get(_rule, 0.0) + _wall_m
        if _rule not in _bench_rules:
            _bench_rules.append(_rule)
        if _sample not in _bench_samples:
            _bench_samples.append(_sample)
    # Sort rules by total time descending
    _bench_rules = sorted(_bench_rules, key=lambda r: _bench_total.get(r, 0), reverse=True)

    # Heatmap: rules (y-axis) × samples (x-axis), colour = wall_minutes
    _hm_z    = [[_bench_by.get(r, {}).get(s, None) for s in _bench_samples] for r in _bench_rules]
    _hm_text = [[(f"{_bench_by.get(r,{}).get(s,0):.1f} min" if _bench_by.get(r,{}).get(s) is not None else "") for s in _bench_samples] for r in _bench_rules]
    fig_bench_heatmap = go.Figure(go.Heatmap(
        z=_hm_z, x=_bench_samples, y=_bench_rules,
        text=_hm_text, texttemplate="%{text}",
        colorscale="YlOrRd", colorbar=dict(title="min"),
        hovertemplate="<b>%{y}</b> — %{x}<br>%{text}<extra></extra>",
    ))
    fig_bench_heatmap.update_layout(
        title="Rule Runtime Heatmap (wall-clock minutes)",
        xaxis_title="Sample", yaxis_title="Rule",
        height=max(400, len(_bench_rules) * 22 + 120),
        margin=dict(l=220), template=T,
    )

    # Bar chart: total time per rule, sorted descending
    _bar_rules = _bench_rules[:30]  # cap at 30 for readability
    fig_bench_bar = go.Figure(go.Bar(
        x=[_bench_total.get(r, 0) for r in _bar_rules],
        y=_bar_rules,
        orientation='h',
        marker_color=TEAL,
        hovertemplate="<b>%{y}</b><br>Total: %{x:.1f} min<extra></extra>",
        text=[f"{_bench_total.get(r,0):.1f} min" for r in _bar_rules],
        textposition='outside',
    ))
    fig_bench_bar.update_layout(
        title="Total Runtime per Rule (all samples, wall-clock minutes)",
        xaxis_title="Minutes", yaxis_title="",
        height=max(400, len(_bar_rules) * 22 + 120),
        margin=dict(l=220), template=T,
    )

    # Table rows for JS
    bench_table_rows = []
    for _row in bench_rows:
        _wall_s = safe_float(_row.get('wall_s', 0))
        bench_table_rows.append({
            'rule':      _row.get('rule', ''),
            'sample':    _row.get('sample', ''),
            'wall_min':  round(_wall_s / 60, 2),
            'cpu_min':   round(safe_float(_row.get('cpu_s', 0)) / 60, 2),
            'max_rss_mb': round(safe_float(_row.get('max_rss_mb', 0)), 1),
            'mean_load': round(safe_float(_row.get('mean_load', 0)), 2),
        })

    # Serialize — replace </ with <\/ so the browser HTML parser never sees
    # "</script>" inside a <script> block and truncates the JS early.
    def _js(obj):
        return json.dumps(obj).replace("</", "<\\/")

    tax_json        = _js(tax_data)
    vcontact3_json  = _js(vcontact3_data)
    vibrant_json    = _js(vibrant_data)
    vibrant_amg_json = _js(vibrant_amg)
    novelty_json    = _js(novelty_data)
    mimag_json      = _js(mimag_data)
    gtdb_json        = _js(gtdb_data)
    custom_prok_json = _js(custom_prok_data)
    merged_prok_json = _js(merged_prok_data)
    tool_matrix_json = _js(tool_matrix_data)
    phist_json = _js(phist_data)
    overview_json = _js(overview)

    # BLOCO 8 serializations
    votu_json      = _js(votu_data)
    lifestyle_json = _js(lifestyle_data)
    eggnog_json    = _js(eggnog_data)
    alpha_json     = _js(alpha_rows)
    genome_maps_json = _js(genome_maps_dict)
    bench_json     = _js(bench_table_rows)

    # Add Task-5 figures to figs_json and finalize the string
    figs_json["read_funnel"]   = fig_read_funnel.to_json()
    figs_json["tool_heatmap"]  = fig_tool_heatmap.to_json()
    figs_json["viral_depth"]   = fig_viral_depth.to_json()

    # BLOCO 8 figures
    figs_json["lifestyle"]      = fig_lifestyle.to_json()
    figs_json["phrogs"]         = fig_phrogs.to_json()
    figs_json["cog"]            = fig_cog.to_json()
    figs_json["alpha_viral"]    = fig_alpha_viral.to_json()
    figs_json["alpha_prok"]     = fig_alpha_prok.to_json()
    figs_json["alpha_combined"] = fig_alpha_combined.to_json()
    figs_json["pcoa_viral"]     = fig_pcoa_viral.to_json()
    figs_json["pcoa_prok"]      = fig_pcoa_prok.to_json()
    figs_json["pcoa_combined"]  = fig_pcoa_combined.to_json()
    figs_json["procrustes"]     = fig_procrustes.to_json()
    figs_json["abund_viral"]    = fig_abund_viral.to_json()
    figs_json["abund_prok"]     = fig_abund_prok.to_json()
    figs_json["bench_heatmap"]  = fig_bench_heatmap.to_json()
    figs_json["bench_bar"]      = fig_bench_bar.to_json()

    figs_json_str = _js(figs_json)


    def _tv(t): return tool_versions.get(t, "n/a")
    def _cp(k): return cfg_params.get(k, "?")

    params_table_html = (
        '<table class="tool-table"><thead><tr><th>Parameter</th><th>Value</th></tr></thead><tbody>'
        + "".join(f'<tr><td>{k}</td><td><b>{v}</b></td></tr>' for k,v in cfg_params.items())
        + "</tbody></table>"
    )

    versions_table_html = (
        '<table class="tool-table"><thead><tr><th>Tool</th><th>Version</th></tr></thead><tbody>'
        + "".join(
            f'<tr><td>{t}</td><td><span class="ver-badge">{v if v not in ("n/a","") else "—"}</span></td></tr>'
            for t, v in sorted(tool_versions.items())
        )
        + "</tbody></table>"
    )

    tool_table_html = (
        '<table class="tool-table"><thead><tr><th>Step</th><th>Tool</th><th>Version</th><th>Key parameters</th></tr></thead><tbody>'
        + '<tr><td class="cat" colspan="4">QC &amp; Trimming</td></tr>'
        + f'<tr><td>QC &amp; Trimming</td><td>fastp</td><td>{_tv("fastp")}</td><td>—</td></tr>'
        + f'<tr><td>MultiQC</td><td>MultiQC</td><td>{_tv("MultiQC")}</td><td>—</td></tr>'
        + '<tr><td class="cat" colspan="4">Assembly</td></tr>'
        + f'<tr><td>Assembly 1</td><td>MEGAHIT</td><td>{_tv("MEGAHIT")}</td><td>--min-contig-len {_cp("Min contig length")} / mem {_cp("MEGAHIT memory")}</td></tr>'
        + f'<tr><td>Assembly 2</td><td>metaSPAdes</td><td>{_tv("metaSPAdes")}</td><td>--meta / mem {_cp("SPAdes memory")}</td></tr>'
        + f'<tr><td>Deduplication</td><td>MMseqs2</td><td>{_tv("MMseqs2")}</td><td>easy-linclust --min-seq-id {_cp("MMseqs2 min identity")} -c 0.85</td></tr>'
        + f'<tr><td>Assembly QC</td><td>QUAST</td><td>{_tv("QUAST")}</td><td>4 assemblies compared</td></tr>'
        + '<tr><td class="cat" colspan="4">Read Mapping</td></tr>'
        + f'<tr><td>Alignment</td><td>BWA-MEM2</td><td>{_tv("BWA-MEM2")}</td><td>mem -t {_cp("Threads")}</td></tr>'
        + f'<tr><td>Sorting/indexing</td><td>samtools</td><td>{_tv("samtools")}</td><td>sort + index + flagstat</td></tr>'
        + '<tr><td>Coverage</td><td>jgi_summarize_bam_contig_depths</td><td>—</td><td>MetaBAT2 utility</td></tr>'
        + '<tr><td class="cat" colspan="4">Viral Detection</td></tr>'
        + f'<tr><td>Viral tool 1</td><td>VirSorter2</td><td>{_tv("VirSorter2")}</td><td>all classifier groups</td></tr>'
        + f'<tr><td>Viral tool 2</td><td>GeNomad</td><td>{_tv("GeNomad")}</td><td>--disable-nn-classification</td></tr>'
        + f'<tr><td>Viral tool 3</td><td>VIBRANT</td><td>{_tv("VIBRANT")}</td><td>lytic/lysogenic binary call</td></tr>'
        + f'<tr><td>Consensus</td><td>custom script</td><td>—</td><td>{_cp("Viral consensus (min tools)")} agreement required</td></tr>'
        + '<tr><td class="cat" colspan="4">Viral QC &amp; Binning</td></tr>'
        + f'<tr><td>Viral completeness</td><td>CheckV</td><td>{_tv("CheckV")}</td><td>end_to_end — consensus + vMAGs</td></tr>'
        + f'<tr><td>Viral binning</td><td>vRhyme</td><td>{_tv("vRhyme")}</td><td>coverage + protein homology</td></tr>'
        + '<tr><td class="cat" colspan="4">Prokaryote Binning</td></tr>'
        + f'<tr><td>Binner 1</td><td>MetaBAT2</td><td>{_tv("MetaBAT2")}</td><td>--minContig 1500</td></tr>'
        + f'<tr><td>Binner 2</td><td>VAMB</td><td>{_tv("VAMB")}</td><td>--minfasta 200000</td></tr>'
        + f'<tr><td>Binner 3</td><td>SemiBin2</td><td>{_tv("SemiBin2")}</td><td>single_easy_bin --environment {_cp("SemiBin2 environment")}</td></tr>'
        + f'<tr><td>Bin integration</td><td>Binette</td><td>{_tv("Binette")}</td><td>Successor to DAS Tool — same scaffold2bin input</td></tr>'
        + f'<tr><td>Bin taxonomy</td><td>GTDB-Tk</td><td>{_tv("GTDB-Tk")}</td><td>classify_wf — GTDB-Tk r220</td></tr>'
        + '<tr><td class="cat" colspan="4">Viral Taxonomy &amp; Host Prediction</td></tr>'
        + f'<tr><td>Viral taxonomy 1</td><td>GeNomad</td><td>{_tv("GeNomad")}</td><td>ML taxonomy — family/genus from virus_summary</td></tr>'
        + f'<tr><td>Viral taxonomy 2</td><td>Diamond+INPHARED</td><td>{_tv("Diamond")}</td><td>BLASTp vs INPHARED proteins — majority vote</td></tr>'
        + f'<tr><td>Host prediction</td><td>PHIST</td><td>{_tv("PHIST")}</td><td>k-mer similarity — sample MAGs as hosts</td></tr>'
        + f'<tr><td>Bin QC</td><td>CheckM2</td><td>{_tv("CheckM2")}</td><td>predict --universal</td></tr>'
        + "</tbody></table>"
    )


    # ── New CSS design system (plain string — no f-string escaping needed) ──────
    _CSS = """
*{box-sizing:border-box;margin:0;padding:0}
:root{
  --bg0:#f8f9fa;--bg1:#ffffff;--bg2:#f1f5f9;--bg3:#e2e8f0;
  --sidebar-bg:#0f172a;--sidebar-brd:#1e293b;--sidebar-acc:#2dd4bf;
  --acc:#0d9488;--acc2:#14b8a6;--amb:#d97706;--grn:#16a34a;
  --pur:#7c3aed;--red:#ef4444;--txt:#1e293b;--txt2:#475569;
  --txt3:#94a3b8;--brd:#e2e8f0;--card-bg:#ffffff;
  --shadow:0 1px 3px rgba(0,0,0,0.08),0 1px 2px rgba(0,0,0,0.04);
  --font:'Inter','Source Sans Pro',system-ui,-apple-system,sans-serif;
}
[data-theme="dark"]{
  --bg0:#080e1a;--bg1:#0f1827;--bg2:#141f30;--bg3:#1a2a3f;
  --acc:#2dd4bf;--acc2:#5eead4;--amb:#f59e0b;--grn:#10b981;
  --txt:#e8f0f8;--txt2:#a8bdd4;--txt3:#7a9ab8;
  --brd:#182336;--card-bg:#141f30;
  --shadow:0 1px 3px rgba(0,0,0,0.4);
}
body{font-family:var(--font);background:var(--bg0);color:var(--txt);display:flex;min-height:100vh;font-size:14px}
/* ── Sidebar ── */
#sidebar{width:220px;min-height:100vh;background:var(--sidebar-bg);border-right:1px solid var(--sidebar-brd);display:flex;flex-direction:column;flex-shrink:0;position:sticky;top:0;height:100vh;overflow-y:auto}
#sidebar .brand{padding:20px 18px 16px;border-bottom:1px solid var(--sidebar-brd)}
#sidebar .brand-logo{display:flex;align-items:center;gap:9px;margin-bottom:6px}
#sidebar .brand-title{font-size:.88rem;font-weight:700;color:#f1f5f9;line-height:1.3;letter-spacing:-.3px}
#sidebar .brand-name{font-size:.65rem;font-weight:600;color:var(--sidebar-acc);text-transform:uppercase;letter-spacing:1px}
#sidebar .brand-sub{font-size:.65rem;color:#64748b;margin-top:4px;line-height:1.4}
#sidebar .nav-section{padding:16px 14px 4px;font-size:.60rem;font-weight:600;color:#475569;text-transform:uppercase;letter-spacing:.9px}
.tab{display:flex;align-items:center;gap:9px;padding:8px 14px;cursor:pointer;font-size:.78rem;color:#94a3b8;border-radius:6px;margin:1px 7px;transition:all .15s;border:none;background:none;text-align:left;width:calc(100% - 14px)}
.tab:hover{background:rgba(255,255,255,0.06);color:#e2e8f0}
.tab.active{background:rgba(45,212,191,0.12);color:var(--sidebar-acc);border-left:2px solid var(--sidebar-acc)}
.tab .ticon{font-size:.88rem;width:16px;text-align:center;flex-shrink:0}
.tab .tbadge{margin-left:auto;background:rgba(255,255,255,0.08);color:#64748b;font-size:.58rem;padding:1px 5px;border-radius:99px}
.tab.active .tbadge{background:rgba(45,212,191,0.2);color:var(--sidebar-acc)}
#theme-btn{margin:12px 14px 16px;padding:6px 12px;font-size:.68rem;border:1px solid #1e293b;border-radius:5px;background:rgba(255,255,255,0.04);color:#64748b;cursor:pointer;text-align:center;transition:all .15s}
#theme-btn:hover{background:rgba(255,255,255,0.08);color:#94a3b8}
/* ── Main ── */
#main{flex:1;min-width:0;display:flex;flex-direction:column}
#topbar{background:var(--bg1);border-bottom:1px solid var(--brd);padding:13px 28px;display:flex;align-items:center;justify-content:space-between;flex-shrink:0;box-shadow:0 1px 0 var(--brd)}
#topbar .page-title{font-size:.95rem;font-weight:600;color:var(--txt)}
#topbar .chips{display:flex;gap:6px;align-items:center}
#topbar .chip{font-size:.65rem;padding:3px 9px;border-radius:99px;border:1px solid var(--brd);color:var(--txt3);background:var(--bg2)}
.content{padding:22px 28px;flex:1}
.panel{display:none}.panel.active{display:block}
/* ── Sections ── */
.sec{font-size:.68rem;font-weight:700;color:var(--txt3);text-transform:uppercase;letter-spacing:.7px;margin:24px 0 10px;padding-left:2px;display:flex;align-items:center;gap:8px}
.sec::after{content:'';flex:1;height:1px;background:var(--brd);margin-left:6px}
/* ── Cards ── */
.card{background:var(--card-bg);border:1px solid var(--brd);border-radius:10px;padding:18px 20px;margin-bottom:14px;box-shadow:var(--shadow)}
.card h2{font-size:.84rem;font-weight:600;color:var(--txt);margin-bottom:3px}
.sub{font-size:.72rem;color:var(--txt3);margin-bottom:12px;line-height:1.5}
.warn{font-size:.74rem;color:var(--amb);padding:5px 11px;background:rgba(217,119,6,.08);border:1px solid rgba(217,119,6,.2);border-radius:5px;margin-bottom:10px;display:inline-block}
.note{font-size:.72rem;color:var(--txt3);padding:8px 0;border-top:1px solid var(--brd);margin-top:10px}
/* ── Stats grid / metric cards ── */
.stats-grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(148px,1fr));gap:9px;margin-bottom:18px}
.s-hdr{grid-column:1/-1;background:var(--bg2);border:1px solid var(--brd);color:var(--acc);border-radius:7px;padding:9px 15px;font-weight:600;font-size:.78rem;margin-top:5px}
.sc{background:var(--bg1);border:1px solid var(--brd);border-top:3px solid var(--acc);border-radius:8px;padding:13px 15px}
.sc .sv{display:flex;align-items:center;gap:7px}
.sc .v{font-size:1.5rem;font-weight:700;color:var(--acc);line-height:1.1}
.sc .l{font-size:.66rem;color:var(--txt3);margin-top:5px;line-height:1.35}
.sc-ok{border-top-color:var(--grn)}.sc-ok .v{color:var(--grn)}
.sc-warn{border-top-color:var(--amb)}.sc-warn .v{color:var(--amb)}
.sc-bad{border-top-color:var(--red)}.sc-bad .v{color:var(--red)}
.sc-neutral{border-top-color:var(--acc)}.sc-neutral .v{color:var(--acc)}
/* ── Layout helpers ── */
.r2{display:grid;grid-template-columns:1fr 1fr;gap:14px}
@media(max-width:900px){.r2{grid-template-columns:1fr}}
.pw{width:100%}
/* ── Sample dropdown ── */
.sample-ctrl{margin-bottom:12px;display:flex;align-items:center;gap:10px}
.sample-ctrl label{font-size:.72rem;font-weight:600;color:var(--txt2)}
.sample-sel{font-size:.75rem;padding:5px 10px;border:1px solid var(--brd);border-radius:5px;background:var(--bg2);color:var(--txt);cursor:pointer;min-width:160px}
.sample-sel:focus{outline:none;border-color:var(--acc)}
/* ── Quality badges ── */
.qbadge{display:inline-block;font-size:.62rem;font-weight:600;padding:2px 7px;border-radius:3px;white-space:nowrap;vertical-align:middle}
.qbadge-complete,.qbadge-hq{background:rgba(22,163,74,.12);color:var(--grn);border:1px solid rgba(22,163,74,.25)}
.qbadge-mq{background:rgba(217,119,6,.12);color:var(--amb);border:1px solid rgba(217,119,6,.25)}
.qbadge-lq,.qbadge-nd{background:rgba(107,114,128,.10);color:var(--txt3);border:1px solid rgba(107,114,128,.2)}
/* ── Source badges ── */
.src-badge{display:inline-block;font-size:.60rem;font-weight:600;padding:1px 6px;border-radius:3px;white-space:nowrap}
/* ── Tables ── */
.tbl-search{width:100%;padding:6px 10px;font-size:.75rem;border:1px solid var(--brd);border-radius:5px;background:var(--bg2);color:var(--txt);margin-bottom:8px}
.tbl-search:focus{outline:none;border-color:var(--acc)}
.tbl-wrap{clear:both;overflow-x:auto;overflow-y:auto;max-height:inherit}
.dl-btn{font-size:.66rem;padding:3px 10px;border-radius:4px;border:1px solid var(--brd);color:var(--txt3);background:var(--bg2);cursor:pointer;float:right;margin-bottom:6px;transition:all .15s}
.dl-btn:hover{color:var(--acc);border-color:var(--acc)}
table.dtbl{width:100%;border-collapse:collapse;font-size:.76rem}
table.dtbl th{background:var(--bg3);color:var(--txt2);padding:7px 12px;text-align:left;white-space:nowrap;font-weight:600;border-bottom:2px solid var(--brd);position:sticky;top:0}
table.dtbl td{padding:6px 12px;border-bottom:1px solid var(--brd);white-space:nowrap;color:var(--txt)}
table.dtbl tbody tr:nth-child(even){background:var(--bg2)}
table.dtbl tbody tr:hover{background:rgba(13,148,136,.06)}
/* ── Pipeline flow ── */
.flow-wrap{overflow-x:auto;padding:8px 0}
.flow{display:flex;align-items:center;flex-wrap:nowrap;min-width:900px}
.fstep{background:var(--bg2);border:1px solid var(--acc);border-radius:8px;padding:11px 15px;min-width:130px;text-align:center}
.fstep .ftitle{font-weight:600;font-size:.76rem;color:var(--acc)}
.fstep .ftools{font-size:.64rem;color:var(--txt3);margin-top:4px;line-height:1.4}
.fstep .ftag{display:inline-block;background:var(--bg3);color:var(--acc);border-radius:3px;font-size:.60rem;padding:1px 5px;margin:2px 1px}
.farrow{color:var(--txt3);font-size:1.1rem;padding:0 5px;flex-shrink:0}
.fstep-viral{border-color:var(--pur)}.fstep-viral .ftitle{color:var(--pur)}
.fstep-viral .ftag{color:var(--pur)}
.fstep-bin{border-color:var(--grn)}.fstep-bin .ftitle{color:var(--grn)}
.fstep-bin .ftag{color:var(--grn)}
.fstep-qc{border-color:var(--amb)}.fstep-qc .ftitle{color:var(--amb)}
.fstep-qc .ftag{color:var(--amb)}
/* ── Tool table (About) ── */
.tool-table{width:100%;border-collapse:collapse;font-size:.78rem}
.tool-table th{background:var(--bg3);color:var(--txt2);padding:8px 13px;text-align:left;font-weight:600;border-bottom:1px solid var(--brd)}
.tool-table td{padding:7px 13px;border-bottom:1px solid var(--brd);vertical-align:top;color:var(--txt)}
.tool-table tr:hover td{background:var(--bg2)}
.tool-table .cat{font-weight:700;color:var(--acc);background:rgba(13,148,136,.05)!important}
/* ── Code blocks ── */
.cmd-block{background:#0f172a;color:#cdd6f4;border:1px solid var(--brd);border-radius:8px;padding:15px 18px;font-family:'Fira Code','Cascadia Code',monospace;font-size:.76rem;line-height:1.7;overflow-x:auto;margin-top:8px}
.cmd-block .kw{color:#89dceb}.cmd-block .flag{color:#a6e3a1}.cmd-block .val{color:#f9e2af}.cmd-block .comment{color:#475569}
/* ── Misc ── */
::-webkit-scrollbar{width:5px;height:5px}
::-webkit-scrollbar-track{background:var(--bg2)}
::-webkit-scrollbar-thumb{background:var(--brd);border-radius:3px}
::-webkit-scrollbar-thumb:hover{background:var(--acc)}
.ov-group-hdr{margin:18px 0 8px;font-size:.68rem;font-weight:700;text-transform:uppercase;letter-spacing:.7px;display:flex;align-items:center;gap:8px}
.ov-group-hdr::after{content:'';flex:1;height:1px;background:var(--brd);margin-left:6px}
.run-summary{background:var(--bg1);border:1px solid var(--brd);border-left:3px solid var(--acc);border-radius:8px;padding:13px 18px;margin-bottom:18px;font-size:.82rem;color:var(--txt);line-height:1.7;box-shadow:var(--shadow)}
.run-summary b{color:var(--acc)}
.run-summary .rs-label{font-size:.66rem;color:var(--txt3);text-transform:uppercase;letter-spacing:.5px;margin-right:6px}
.about-grid{display:grid;grid-template-columns:1fr 1fr;gap:14px;margin-bottom:14px}
@media(max-width:900px){.about-grid{grid-template-columns:1fr}}
.ver-badge{display:inline-block;background:var(--bg3);color:var(--txt3);border-radius:3px;font-size:.60rem;padding:1px 6px;margin-left:4px;font-family:monospace}
.sunburst-row{display:grid;gap:14px;margin-bottom:14px}
.sunburst-row.cols-2{grid-template-columns:1fr 1fr}
@media(max-width:900px){.sunburst-row.cols-2{grid-template-columns:1fr}}
@media(min-width:1440px){body{font-size:15px}#sidebar{width:240px}.card{padding:22px 24px}}
"""

    html = (
        "<!DOCTYPE html>\n<html lang='en' data-theme='light'>\n<head>\n"
        "<meta charset='UTF-8'><meta name='viewport' content='width=device-width,initial-scale=1.0'>\n"
        "<title>MITE — Metagenomic Integrated Taxonomic Engine</title>\n"
        "<link rel='preconnect' href='https://fonts.googleapis.com'>\n"
        "<link rel='preconnect' href='https://fonts.gstatic.com' crossorigin>\n"
        "<link href='https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap' rel='stylesheet'>\n"
        "<script src='https://cdn.plot.ly/plotly-2.27.0.min.js'></script>\n"
        "<style>\n" + _CSS + "\n</style>\n"
        "</head>\n<body>\n"
    ) + f"""<nav id="sidebar">
  <div class="brand">
    <div class="brand-logo">
      <svg width="28" height="28" viewBox="0 0 28 28" fill="none" xmlns="http://www.w3.org/2000/svg">
        <circle cx="14" cy="14" r="13" stroke="#2dd4bf" stroke-width="1.5" fill="none"/>
        <circle cx="14" cy="14" r="5" fill="#2dd4bf" opacity=".8"/>
        <line x1="14" y1="1" x2="14" y2="6" stroke="#2dd4bf" stroke-width="1.5"/>
        <line x1="14" y1="22" x2="14" y2="27" stroke="#2dd4bf" stroke-width="1.5"/>
        <line x1="1" y1="14" x2="6" y2="14" stroke="#2dd4bf" stroke-width="1.5"/>
        <line x1="22" y1="14" x2="27" y2="14" stroke="#2dd4bf" stroke-width="1.5"/>
        <line x1="4.5" y1="4.5" x2="8" y2="8" stroke="#2dd4bf" stroke-width="1" opacity=".5"/>
        <line x1="20" y1="20" x2="23.5" y2="23.5" stroke="#2dd4bf" stroke-width="1" opacity=".5"/>
        <line x1="23.5" y1="4.5" x2="20" y2="8" stroke="#2dd4bf" stroke-width="1" opacity=".5"/>
        <line x1="8" y1="20" x2="4.5" y2="23.5" stroke="#2dd4bf" stroke-width="1" opacity=".5"/>
      </svg>
      <div>
        <div class="brand-name">MITE</div>
        <div class="brand-title">Metagenomics<br>Pipeline Report</div>
      </div>
    </div>
    <div class="brand-sub">{sample_list}</div>
  </div>
  <div class="nav-section">Analysis</div>
  <button class="tab active" onclick="showTab('overview', this)"><span class="ticon">▦</span>Overview<span class="tbadge">5</span></button>
  <button class="tab" onclick="showTab('readqc',   this)"><span class="ticon">🔎</span>Read QC</button>
  <button class="tab" onclick="showTab('assembly', this)"><span class="ticon">🔬</span>Assembly</button>
  <button class="tab" onclick="showTab('viral',    this)"><span class="ticon">🦠</span>Viral Analysis</button>
  <button class="tab" onclick="showTab('taxonomy', this)"><span class="ticon">🌳</span>Taxonomy</button>
  <button class="tab" onclick="showTab('hostpred', this)"><span class="ticon">🔗</span>Host Prediction</button>
  <button class="tab" onclick="showTab('bins',     this)"><span class="ticon">🧫</span>Bin Quality</button>
  <button class="tab" onclick="showTab('abundance', this)"><span class="ticon">📈</span>Abundance</button>
  <button class="tab" onclick="showTab('ecology',  this)"><span class="ticon">🌿</span>Ecology</button>
  <button class="tab" onclick="showTab('annotation',this)"><span class="ticon">🗺</span>Annotation</button>
  <div class="nav-section">Pipeline</div>
  <button class="tab" onclick="showTab('performance',this)"><span class="ticon">⏱</span>Performance</button>
  <button class="tab" onclick="showTab('about',    this)"><span class="ticon">ℹ</span>About</button>
  <div style="flex:1"></div>
  <button id="theme-btn" onclick="toggleTheme()">☀ Light / 🌙 Dark</button>
</nav>
<div id="main">
<div id="topbar">
  <span class="page-title" id="topbar-title">Overview</span>
  <div class="chips">
    <span class="chip">MITE pipeline</span>
    <span class="chip" id="chip-samples">{len(samples)} sample{"s" if len(samples)!=1 else ""}</span>
    <span class="chip" id="chip-date"></span>
  </div>
</div>
<div class="content">

<div id="panel-overview" class="panel active">
  <div id="run-summary" class="run-summary"></div>
  <div class="stats-grid" id="ovg"></div>
</div>

<div id="panel-readqc" class="panel">
  <div id="sample-ctrl-readqc" class="sample-ctrl"></div>
  <div class="sec">Read Quality</div>
  <div class="card"><h2>Total Reads — Raw vs Trimmed</h2><p class="sub">R1 and R2 before/after fastp.</p><div id="p-fq-reads" class="pw"></div></div>
  <div class="r2">
    <div class="card"><h2>Mean Per-Base Quality</h2><p class="sub">Q30 = 1 error per 1,000 bases.</p><div id="p-fq-qual" class="pw"></div></div>
    <div class="card"><h2>GC Content</h2><p class="sub">Soil metagenomes typically 55–65%.</p><div id="p-fq-gc" class="pw"></div></div>
  </div>
  <div class="sec">Trimming &amp; Mapping</div>
  <div class="card"><h2>Fastp Trimming Efficiency</h2><p class="sub">From fastp JSON report.</p><div id="p-trim" class="pw"></div></div>
  <div class="card"><h2>Read Mapping Rate (BWA-MEM2 → Assembly)</h2><p class="sub">samtools flagstat. &lt;50% may indicate assembly quality issues.</p><div id="p-mapping" class="pw"></div></div>
  <div class="card"><h2>Read Utilization Funnel</h2><p class="sub">Reads progressing from raw sequencing to assigned bins/viruses.</p><div id="p-read-funnel" class="pw"></div></div>
</div>

<div id="panel-assembly" class="panel">
  <div id="sample-ctrl-assembly" class="sample-ctrl"></div>
  <div class="card"><h2>Assembly Progression — All 4 Stages</h2><p class="sub">MEGAHIT → metaSPAdes → merged+filtered → MMseqs2 deduplicated.</p><div id="p-asm-prog" class="pw"></div></div>
  <div class="card"><h2>Final Assembly Quality (QUAST)</h2><div id="p-assembly" class="pw"></div></div>
  <div class="r2">
    <div class="card"><h2>Contig Length Distribution</h2><p class="sub">Log x-axis.</p><div id="p-contig-len" class="pw"></div></div>
    <div class="card"><h2>Contig Coverage Distribution</h2><p class="sub">Clipped at 500×.</p><div id="p-coverage" class="pw"></div></div>
  </div>
</div>

<div id="panel-viral" class="panel">
  <div id="sample-ctrl-viral" class="sample-ctrl"></div>
  <div class="sec">Detection</div>
  <div class="card"><h2>Per-Tool Counts (Before Consensus Filter)</h2><div id="p-tool-counts" class="pw"></div></div>
  <div class="card"><h2>Tool Agreement — Stacked Bar</h2><p class="sub">Pipeline threshold: ≥{min_viral_tools} tools.</p><div id="p-viral" class="pw"></div></div>
  <div class="card"><h2>Score Distributions</h2><div id="p-scores" class="pw"></div></div>
  <div class="card"><h2>Tool Agreement Heatmap</h2><p class="sub">Each row = one viral contig. Columns = detection tools. Teal = detected.</p><div id="p-tool-heatmap" class="pw"></div></div>
  <div class="sec">Quality Assessment</div>
  <div class="card"><h2>Viral Quality Pyramid</h2><p class="sub">CheckV quality tiers across all samples (MIUViG standard).</p><div id="p-viral-pyramid" style="height:340px"></div></div>
  <div class="card"><h2>CheckV — Consensus Contigs</h2><div id="p-checkv" class="pw"></div></div>
  <div class="card"><h2>CheckV — vRhyme vMAGs</h2><div id="p-checkv-vrh" class="pw"></div></div>
  <div class="card"><h2>CheckV — Length vs Completeness (● consensus  ◆ vRhyme)</h2><div id="p-checkv-scatter" class="pw"></div></div>
  <div class="r2">
    <div class="card"><h2>Viral Consensus Contig Length Distribution</h2><div id="p-viral-len" class="pw"></div></div>
    <div class="card"><h2>Viral Contig Coverage Depth</h2><p class="sub">Read depth distribution for viral consensus contigs.</p><div id="p-viral-depth" class="pw"></div></div>
  </div>
  <div class="r2">
    <div class="card"><h2>VIBRANT — AMGs per Sample</h2><p class="sub">Auxiliary Metabolic Genes detected by VIBRANT.</p><div id="p-vibrant-amg" class="pw"></div></div>
    <div class="card"><h2>VIBRANT — Top AMG Pathways</h2><p class="sub">Most frequent metabolic pathways targeted by viral AMGs.</p><div id="p-vibrant-pathways" class="pw"></div></div>
  </div>
  <div class="card"><h2>VIBRANT — AMG Table</h2><p class="sub">All AMG proteins with KO, gene name, metabolic pathway and category.</p>
    <div id="tbl-vibrant-amg" style="max-height:340px;overflow:auto"></div></div>
  <div class="sec">vRhyme</div>
  <div class="card"><h2>vMAG Summary</h2><div id="p-vrhyme" class="pw"></div></div>
  <div class="card"><h2>Per-Bin Detail (size = redundancy %)</h2><div id="p-vrhyme-detail" class="pw"></div></div>
  <div class="sec">vOTU Summary (MIUViG)</div>
  <div id="votu-summary-cards" style="margin-bottom:12px"></div>
  <div class="card"><h2>vOTU Master Table</h2><p class="sub">All viral OTUs with CheckV quality, taxonomy, lifestyle, host, and AMG count.</p>
    <div id="tbl-votu" style="max-height:400px;overflow:auto"></div></div>
  <div class="sec">Lifestyle &amp; Functional Annotation</div>
  <div class="card"><h2>Lifestyle Prediction — Lytic vs Lysogenic (VIBRANT)</h2><div id="p-lifestyle" class="pw"></div></div>
  <div class="card"><h2>PHROGS Functional Categories (Pharokka)</h2><p class="sub">Gene-level functional annotation of bacteriophages (top-N HQ per sample).</p><div id="p-phrogs" class="pw"></div></div>
</div>

<div id="panel-taxonomy" class="panel">
  <div id="sample-ctrl-taxonomy" class="sample-ctrl"></div>
  <div class="sec">Viral Taxonomy — Unified (vConTACT3 + Diamond/INPHARED + Diamond/Custom + GeNomad)</div>
  <div class="r2">
    <div class="card"><h2>Classification Source Distribution</h2><p class="sub">Proportion of viral contigs classified by each method tier.</p><div id="p-tax-source" class="pw"></div></div>
    <div class="card"><h2>Top 15 Viral Families</h2><p class="sub">Most frequent families — use dropdown above to filter by sample.</p><div id="p-tax-family" class="pw"></div></div>
  </div>
  <div class="card"><h2>Unified Viral Taxonomy — Hierarchical Sunburst</h2><p class="sub">Family → Genus breakdown (click to drill down). All classification tiers merged.</p>
    <div id="p-tax-sunburst" style="height:520px"></div></div>
  <div class="card"><h2>Viral Taxonomy — Master Table</h2><p class="sub">All viral contigs with classification source badge, completeness (CheckV) and genome length.</p>
    <div id="tbl-vcontact" style="max-height:420px;overflow:auto"></div></div>
  <div class="sec">Viral Novelty — Unclassified Fraction</div>
  <div class="r2">
    <div class="card"><h2>% Novel Viruses per Sample</h2><p class="sub">Fraction without family-level classification.</p><div id="p-novelty" class="pw"></div></div>
    <div class="card"><h2>Classification Depth</h2><p class="sub">Stack: classified vs unclassified.</p><div id="p-tax-depth" class="pw"></div></div>
  </div>
  <div class="sec">Prokaryotic Taxonomy — GTDB-Tk r220 + Diamond Custom (merged)</div>
  <div class="card"><h2>MAG Taxonomy — Master Table</h2><p class="sub">GTDB-Tk has priority; Diamond/Custom fills unclassified MAGs. Source shown as badge.</p>
    <div id="tbl-prok-master" style="max-height:380px;overflow:auto"></div></div>
  <div class="r2">
    <div class="card"><h2>Domain Summary per Sample</h2><p class="sub">Bacteria / Archaea / Unclassified (stacked bar).</p><div id="p-prok-domain" class="pw"></div></div>
    <div class="card"><h2>Top Bacterial Phyla</h2><p class="sub">Based on GTDB-Tk + Diamond/Custom merged taxonomy.</p><div id="p-prok-phyla" class="pw"></div></div>
  </div>
  <div class="card"><h2>Prokaryotic Taxonomy — Hierarchical Sunburst</h2><p class="sub">Domain → Phylum → Class (click to drill down).</p>
    <div id="p-gtdb-sunburst" style="height:500px"></div></div>
  <div class="card"><h2>GTDB-Tk Full Classification Table</h2>
    <div id="tbl-gtdb" style="max-height:340px;overflow:auto"></div></div>
  <div id="custom-prok-section" style="display:none">
    <div class="sec">Diamond Custom DB — Prokaryote Taxonomy</div>
    <div class="r2">
      <div class="card"><h2>Top Genera</h2><div id="p-custom-prok-genus" class="pw"></div></div>
      <div class="card"><h2>Top Phyla</h2><div id="p-custom-prok-phylum" class="pw"></div></div>
    </div>
    <div class="card"><h2>Custom DB — Classification Table</h2>
      <div id="tbl-custom-prok" style="max-height:340px;overflow:auto"></div></div>
  </div>
</div>

<div id="panel-hostpred" class="panel">
  <div id="sample-ctrl-hostpred" class="sample-ctrl"></div>
  <div class="sec">Host Prediction — PHIST</div>
  <div class="r2">
    <div class="card"><h2>Viruses with Host Prediction</h2><p class="sub">PHIST k-mer similarity predictions per sample.</p><div id="p-hostpred-cov" class="pw"></div></div>
    <div class="card"><h2>Top Predicted Host MAGs</h2><p class="sub">Most frequent host MAGs across PHIST predictions.</p><div id="p-hostpred-hosts" class="pw"></div></div>
  </div>
  <div class="card"><h2>Host Prediction — Results Table</h2><p class="sub">PHIST predictions with host MAG and significance score.</p>
    <div id="tbl-hostpred" style="max-height:400px;overflow:auto"></div></div>
  <div class="sec">Virus-Host Interaction Network</div>
  <div class="card"><h2>Virus ↔ Host Network</h2><p class="sub">Virus (purple) ↔ Host MAG (teal). Edge = PHIST prediction.</p><div id="p-vhost-network" style="height:520px"></div></div>
</div>

<div id="panel-bins" class="panel">
  <div id="sample-ctrl-bins" class="sample-ctrl"></div>
  <div class="r2">
    <div class="card"><h2>Bins per Tool</h2><div id="p-binner-total" class="pw"></div></div>
    <div class="card"><h2>Binette Final Bins — Domain Classification</h2><p class="sub">Bacterial / Archaeal / Unknown from GTDB-Tk.</p><div id="p-das-tax" class="pw"></div></div>
  </div>
  <div class="r2">
    <div class="card"><h2>MAG Quality Tiers (MIMAG)</h2><p class="sub">HQ: ≥90% / ≤5%. MQ: ≥50% / ≤10%. Published standard.</p><div id="p-mimag" class="pw"></div></div>
    <div class="card"><h2>MIMAG Summary Table</h2><div id="tbl-mimag" style="max-height:300px;overflow:auto"></div></div>
  </div>
  <div class="card"><h2>Completeness vs Contamination (CheckM2)</h2><p class="sub">HQ zone (green) / MQ zone (amber). Hover for bin name and taxonomy.</p><div id="p-checkm2" class="pw"></div></div>
  <div class="r2">
    <div class="card"><h2>Quality Distributions</h2><div id="p-cm2-hist" class="pw"></div></div>
    <div class="card"><h2>Bin Genome Size Distribution</h2><p class="sub">Bacteria ~2–6 Mb; Archaea often larger.</p><div id="p-bin-size" class="pw"></div></div>
  </div>
</div>

<div id="panel-bins-extra" style="display:none">
  <div class="sec">Functional Annotation (COG)</div>
  <div class="card"><h2>COG Functional Categories (EggNOG-mapper)</h2><p class="sub">Broad functional annotation of MAG genes from EggNOG-mapper v2.</p><div id="p-cog" class="pw"></div></div>
</div>

<div id="panel-abundance" class="panel">
  <div class="card"><h2>Bin Abundance Across Samples</h2><p class="sub">Weighted mean read depth per bin.</p>
  <div id="p-abundance" class="pw"></div>
  <p class="note">⚠️ Within-sample coverage only. Cross-sample comparison requires dRep + re-mapping.</p></div>
</div>

<!-- ══════════════════════ ECOLOGY TAB ══════════════════════ -->
<div id="panel-ecology" class="panel">
  <div class="sec">Alpha Diversity</div>
  <div class="card"><h2>Alpha Diversity — Viral Community</h2><p class="sub">Shannon, Simpson, Observed richness and Chao1 per sample (vOTU-based RPKM ≥ 0).</p><div id="p-alpha-viral" class="pw"></div></div>
  <div class="card"><h2>Alpha Diversity — Prokaryotic Community</h2><p class="sub">MAG-based. Requires ≥ 2 MAGs with RPKM > 0.</p><div id="p-alpha-prok" class="pw"></div></div>
  <div class="card"><h2>Alpha Diversity — Combined (Viral + Prokaryotic)</h2><div id="p-alpha-combined" class="pw"></div></div>
  <div class="sec">Beta Diversity — PCoA (Bray-Curtis)</div>
  <div class="r2">
    <div class="card"><h2>Viral Community PCoA</h2><p class="sub">Principal Coordinates Analysis. Requires ≥ 2 samples with viral abundance.</p><div id="p-pcoa-viral" class="pw"></div></div>
    <div class="card"><h2>Prokaryotic Community PCoA</h2><p class="sub">MAG-based Bray-Curtis dissimilarity.</p><div id="p-pcoa-prok" class="pw"></div></div>
  </div>
  <div class="r2">
    <div class="card"><h2>Combined Community PCoA</h2><div id="p-pcoa-combined" class="pw"></div></div>
    <div class="card"><h2>Procrustes Overlay — Viral vs Prokaryotic</h2><p class="sub">● = viral position, ■ = prokaryotic position (same sample). Disparity → 0 = strong co-variation.</p><div id="p-procrustes" class="pw"></div></div>
  </div>
  <div class="sec">Relative Abundance (CoverM RPKM)</div>
  <div class="card"><h2>Viral vOTU Relative Abundance — Top 20</h2><p class="sub">RPKM-normalized. Stacked bar across samples. Top 20 by total RPKM.</p><div id="p-abund-viral" class="pw"></div></div>
  <div class="card"><h2>Prokaryotic MAG Relative Abundance — Top 20</h2><div id="p-abund-prok" class="pw"></div></div>
</div>

<!-- ══════════════════════ ANNOTATION TAB ══════════════════════ -->
<div id="panel-annotation" class="panel">
  <div class="sec">Circular Genome Maps</div>
  <div class="card">
    <h2>Genome Map Viewer</h2>
    <p class="sub">Circular maps generated with pycirclize. Phages: PHROGS categories. Viruses: VOGDB categories. Prokaryotes: COG categories.</p>
    <div style="display:flex;gap:10px;flex-wrap:wrap;margin-bottom:14px">
      <select id="gm-sample-sel" style="padding:6px 10px;border-radius:6px;border:1px solid var(--brd);background:var(--bg1);color:var(--txt);font-size:.82rem"></select>
      <button class="dl-btn" onclick="setGmMode('phage')" id="gm-btn-phage" style="border-radius:6px;padding:6px 14px">🦠 Phages</button>
      <button class="dl-btn" onclick="setGmMode('virus')" id="gm-btn-virus" style="border-radius:6px;padding:6px 14px">🔵 Other Viruses</button>
      <button class="dl-btn" onclick="setGmMode('prok')"  id="gm-btn-prok"  style="border-radius:6px;padding:6px 14px">🧫 Prokaryotes</button>
      <select id="gm-genome-sel" style="padding:6px 10px;border-radius:6px;border:1px solid var(--brd);background:var(--bg1);color:var(--txt);font-size:.82rem"></select>
    </div>
    <div id="genome-map-viewer" style="width:100%;overflow:auto;text-align:center;min-height:200px"></div>
  </div>
</div>

<!-- ══════════════════════ ABOUT TAB ══════════════════════ -->
<div id="panel-about" class="panel">

  <div class="sec">Pipeline Overview</div>
  <div class="card">
    <h2>Workflow Diagram</h2>
    <p class="sub">Each box represents a Snakemake rule. Click through the tabs for detailed results at each step.</p>
    <div class="flow-wrap">
      <div class="flow">
        <div class="fstep fstep-qc">
          <div class="ftitle">📥 Raw Reads</div>
          <div class="ftools">Illumina paired-end<br><span class="ftag">FASTQ R1+R2</span></div>
        </div>
        <div class="farrow">→</div>
        <div class="fstep fstep-qc">
          <div class="ftitle">🔎 QC &amp; Trimming</div>
          <div class="ftools"><span class="ftag">fastp</span></div>
        </div>
        <div class="farrow">→</div>
        <div class="fstep">
          <div class="ftitle">🔬 Assembly</div>
          <div class="ftools"><span class="ftag">MEGAHIT</span><span class="ftag">metaSPAdes</span><br><span class="ftag">MMseqs2 dedup</span></div>
        </div>
        <div class="farrow">→</div>
        <div class="fstep">
          <div class="ftitle">🗺️ Mapping</div>
          <div class="ftools"><span class="ftag">BWA-MEM2</span><span class="ftag">CoverM</span></div>
        </div>
        <div class="farrow">→</div>
        <div class="fstep fstep-viral">
          <div class="ftitle">🦠 Viral Detection</div>
          <div class="ftools"><span class="ftag">VirSorter2</span><span class="ftag">GeNomad</span><br><span class="ftag">VIBRANT</span><br>≥{min_viral_tools} tools consensus</div>
        </div>
        <div class="farrow">→</div>
        <div class="fstep fstep-viral">
          <div class="ftitle">✅ Viral QC &amp; Binning</div>
          <div class="ftools"><span class="ftag">CheckV</span><span class="ftag">vRhyme</span></div>
        </div>
        <div class="farrow">→</div>
        <div class="fstep fstep-bin">
          <div class="ftitle">🧫 Prokaryote Binning</div>
          <div class="ftools"><span class="ftag">MetaBAT2</span><span class="ftag">VAMB</span><br><span class="ftag">SemiBin2</span><span class="ftag">Binette</span></div>
        </div>
        <div class="farrow">→</div>
        <div class="fstep fstep-bin">
          <div class="ftitle">📊 Bin QC</div>
          <div class="ftools"><span class="ftag">CheckM2</span><span class="ftag">GTDB-Tk</span></div>
        </div>
        <div class="farrow">→</div>
        <div class="fstep fstep-qc">
          <div class="ftitle">📈 This Report</div>
          <div class="ftools"><span class="ftag">Plotly HTML</span></div>
        </div>
      </div>
    </div>
  </div>

  <div class="sec">Run Parameters</div>
  <div class="about-grid">
    <div class="card">
      <h2>Pipeline Configuration</h2>
      <p class="sub">Values set in the Snakefile for this run.</p>
      {params_table_html}
    </div>
    <div class="card">
      <h2>Tool Versions</h2>
      <p class="sub">Versions collected at runtime via conda run.</p>
      {versions_table_html}
    </div>
  </div>

</div>
<!-- end about -->

<div id="panel-performance" class="panel">
  <div class="sec">Runtime Overview</div>
  <div class="card"><h2>Total Runtime per Rule</h2><p class="sub">Wall-clock minutes summed across all samples, sorted by duration. Top 30 shown.</p><div id="p-bench-bar" class="pw"></div></div>
  <div class="sec">Per-Sample Heatmap</div>
  <div class="card"><h2>Rule × Sample Runtime Heatmap</h2><p class="sub">Wall-clock minutes. Empty cells = rule did not run for that sample.</p><div id="p-bench-heatmap" class="pw"></div></div>
  <div class="sec">Detailed Table</div>
  <div class="card"><h2>All Benchmark Records</h2><p class="sub">Wall time (min), CPU time (min), peak RAM (MB), mean CPU load.</p><div id="bench-tbl"></div></div>
</div>
<!-- end performance -->

</div><!-- end content -->
<script>
const SAMPLES={samples_json};const OVERVIEW={overview_json};const FIGS={figs_json_str};
const TAX_DATA={tax_json};const GTDB_DATA={gtdb_json};const CUSTOM_PROK={custom_prok_json};
const MERGED_PROK={merged_prok_json};const TOOL_MATRIX={tool_matrix_json};
const VC3_DATA={vcontact3_json};const NOVELTY={novelty_json};const MIMAG={mimag_json};
const VIBRANT_DATA={vibrant_json};const VIBRANT_AMG={vibrant_amg_json};
const PHIST_DATA={phist_json};
const VOTU_DATA={votu_json};
const LIFESTYLE={lifestyle_json};
const ALPHA_DATA={alpha_json};
const EGGNOG_DATA={eggnog_json};
const GENOME_MAPS={genome_maps_json};
const BENCH_DATA={bench_json};
const TAB_TITLES={{'overview':'Overview','readqc':'Read QC','assembly':'Assembly QC',
  'viral':'Viral Analysis','taxonomy':'Taxonomy','hostpred':'Host Prediction',
  'bins':'Bin Quality','abundance':'Abundance',
  'ecology':'Ecology','annotation':'Annotation',
  'performance':'Performance','about':'About'}};
let ACTIVE_SAMPLE='__all__';

// ── Theme toggle ────────────────────────────────────────────────────────────
function toggleTheme(){{
  const html=document.documentElement;
  const isDark=html.getAttribute('data-theme')==='dark';
  html.setAttribute('data-theme',isDark?'light':'dark');
  // Re-render all Plotly charts with matching bgcolor
  const bg=isDark?'rgba(0,0,0,0)':'rgba(0,0,0,0)';
  document.querySelectorAll('.js-plotly-plot').forEach(el=>{{
    try{{Plotly.relayout(el,{{paper_bgcolor:bg,plot_bgcolor:isDark?'#f8fafc':'#141f30'}});}}catch(e){{}}
  }});
}}

// ── Tab navigation ──────────────────────────────────────────────────────────
function showTab(name,el){{
  document.querySelectorAll(".panel").forEach(p=>p.classList.remove("active"));
  document.querySelectorAll(".tab").forEach(t=>t.classList.remove("active"));
  document.getElementById("panel-"+name).classList.add("active");
  el.classList.add("active");
  document.getElementById("topbar-title").textContent=TAB_TITLES[name]||name;
  window.dispatchEvent(new Event("resize"));
}}

document.getElementById("chip-date").textContent=new Date().toISOString().slice(0,10);
const cfg={{responsive:true,displayModeBar:true,scrollZoom:false,
  modeBarButtonsToRemove:['sendDataToCloud','select2d','lasso2d']}};

// ── Sample dropdown factory ──────────────────────────────────────────────────
function makeSampleDropdown(containerId, onChangeFn){{
  const wrap=document.getElementById(containerId);
  if(!wrap||SAMPLES.length<2)return;
  const label=document.createElement('label');label.textContent='Sample:';
  const sel=document.createElement('select');sel.className='sample-sel';
  sel.innerHTML='<option value="__all__">All samples — aggregated</option>';
  SAMPLES.forEach(s=>{{sel.innerHTML+=`<option value="${{s}}">${{s}}</option>`;}});
  sel.addEventListener('change',()=>{{ACTIVE_SAMPLE=sel.value;onChangeFn(sel.value);}});
  wrap.appendChild(label);wrap.appendChild(sel);
}}

// ── Helper: render pre-built Plotly figure from FIGS dict ───────────────────
function rf(id,key){{const f=JSON.parse(FIGS[key]);Plotly.newPlot(id,f.data,f.layout,cfg);}}

// ── Helper: filter+render a Plotly figure by sample ────────────────────────
function rfFiltered(id,key,sample){{
  const f=JSON.parse(FIGS[key]);
  if(sample==='__all__'){{Plotly.newPlot(id,f.data,f.layout,cfg);return;}}
  const filtered=f.data.filter(t=>
    (t.name&&t.name.includes(sample))||
    (Array.isArray(t.x)&&t.x.some&&t.x.some(v=>String(v)===sample)));
  Plotly.newPlot(id,filtered.length?filtered:f.data,f.layout,cfg);
}}

// ── Source badge helper ──────────────────────────────────────────────────────
const SOURCE_STYLE={{
  'vcontact3':       {{label:'vConTACT3',       bg:'rgba(124,58,237,.12)',color:'#7c3aed',brd:'rgba(124,58,237,.3)'}},
  'diamond_inphared':{{label:'Diamond/INPHARED', bg:'rgba(37,99,235,.12)', color:'#2563eb',brd:'rgba(37,99,235,.3)'}},
  'diamond_custom':  {{label:'Diamond/Custom',   bg:'rgba(217,119,6,.12)', color:'#d97706',brd:'rgba(217,119,6,.3)'}},
  'genomad':         {{label:'GeNomad',          bg:'rgba(22,163,74,.12)', color:'#16a34a',brd:'rgba(22,163,74,.3)'}},
  'unclassified':    {{label:'Unclassified',     bg:'rgba(107,114,128,.1)',color:'#6b7280',brd:'rgba(107,114,128,.2)'}},
}};
function sourceBadge(src){{
  const k=(src||'').toLowerCase();
  const m=SOURCE_STYLE[k]||{{label:src||'Unknown',bg:'rgba(107,114,128,.1)',color:'#6b7280',brd:'rgba(107,114,128,.2)'}};
  return `<span class="src-badge" style="background:${{m.bg}};color:${{m.color}};border:1px solid ${{m.brd}}">${{m.label}}</span>`;
}}

// ── Quality badge helper ──────────────────────────────────────────────────────
function qualBadge(q){{
  const s=(q||'').toLowerCase();
  if(s==='complete')return`<span class="qbadge qbadge-complete">Complete</span>`;
  if(s.includes('high'))return`<span class="qbadge qbadge-hq">HQ</span>`;
  if(s.includes('medium'))return`<span class="qbadge qbadge-mq">MQ</span>`;
  if(s.includes('low'))return`<span class="qbadge qbadge-lq">LQ</span>`;
  return`<span class="qbadge qbadge-nd">${{q||'ND'}}</span>`;
}}

// ── Enhanced makeTable with search, renderers, sticky header ────────────────
function makeTable(id,rows,cols,renderers,labels){{
  renderers=renderers||{{}};labels=labels||{{}};
  const el=document.getElementById(id);
  if(!el)return;
  if(!rows||!rows.length){{el.innerHTML='<p style="color:var(--txt3);padding:10px;font-size:.78rem">No data available</p>';return;}}
  const btnId='dl-'+id;
  let h=`<input type="text" class="tbl-search" placeholder="Search table…" oninput="filterTable('${{id}}',this.value)">`;
  h+=`<button class="dl-btn" id="${{btnId}}">↓ TSV</button>`;
  h+='<div class="tbl-wrap"><table class="dtbl" id="t-'+id+'"><thead><tr>';
  cols.forEach(c=>h+=`<th>${{labels[c]||c}}</th>`);
  h+='</tr></thead><tbody>';
  rows.forEach(r=>{{
    h+='<tr>';
    cols.forEach(c=>{{
      const raw=r[c]??'';
      const cell=renderers[c]?renderers[c](raw,r):String(raw).replace(/</g,'&lt;');
      h+=`<td>${{cell}}</td>`;
    }});
    h+='</tr>';
  }});
  h+='</tbody></table></div>';
  el.innerHTML=h;
  document.getElementById(btnId).addEventListener('click',()=>{{
    const tsv=[cols.map(c=>labels[c]||c).join('\\t'),...rows.map(r=>cols.map(c=>String(r[c]??'')).join('\\t'))].join('\\n');
    const a=document.createElement('a');
    a.href=URL.createObjectURL(new Blob([tsv],{{type:'text/tab-separated-values'}}));
    a.download=(id.replace(/^tbl-/,'')||id)+'.tsv';a.click();
  }});
}}

function filterTable(id,query){{
  const tbl=document.getElementById('t-'+id);
  if(!tbl)return;
  const q=query.toLowerCase();
  tbl.querySelectorAll('tbody tr').forEach(tr=>{{
    tr.style.display=tr.textContent.toLowerCase().includes(q)?'':'none';
  }});
}}

// ── Render all pre-built figures ────────────────────────────────────────────
rf("p-fq-reads","fq_reads");rf("p-fq-qual","fq_qual");rf("p-fq-gc","fq_gc");
rf("p-trim","trim");rf("p-mapping","mapping");
rf("p-asm-prog","asm_prog");rf("p-assembly","assembly");
rf("p-contig-len","contig_len");rf("p-coverage","coverage");
rf("p-tool-counts","tool_counts");rf("p-viral","viral");rf("p-scores","scores");
rf("p-checkv","checkv");rf("p-checkv-vrh","checkv_vrh");rf("p-checkv-scatter","checkv_scatter");
rf("p-viral-len","viral_len");rf("p-vrhyme","vrhyme");rf("p-vrhyme-detail","vrhyme_detail");
rf("p-checkm2","checkm2");rf("p-cm2-hist","cm2_hist");rf("p-bin-size","bin_size");
rf("p-binner-total","binner_total");rf("p-das-tax","das_tax");rf("p-abundance","abundance");
rf("p-read-funnel","read_funnel");rf("p-tool-heatmap","tool_heatmap");rf("p-viral-depth","viral_depth");
// BLOCO 8 pre-built figures
rf("p-lifestyle","lifestyle");rf("p-phrogs","phrogs");rf("p-cog","cog");
rf("p-alpha-viral","alpha_viral");rf("p-alpha-prok","alpha_prok");rf("p-alpha-combined","alpha_combined");
rf("p-pcoa-viral","pcoa_viral");rf("p-pcoa-prok","pcoa_prok");rf("p-pcoa-combined","pcoa_combined");
rf("p-procrustes","procrustes");rf("p-abund-viral","abund_viral");rf("p-abund-prok","abund_prok");

// ── Per-sample dropdowns ─────────────────────────────────────────────────────
makeSampleDropdown('sample-ctrl-readqc', s=>{{
  rfFiltered('p-fq-reads','fq_reads',s);rfFiltered('p-fq-qual','fq_qual',s);
  rfFiltered('p-fq-gc','fq_gc',s);rfFiltered('p-trim','trim',s);
  rfFiltered('p-mapping','mapping',s);rfFiltered('p-read-funnel','read_funnel',s);
}});
makeSampleDropdown('sample-ctrl-assembly', s=>{{
  rfFiltered('p-asm-prog','asm_prog',s);rfFiltered('p-assembly','assembly',s);
  rfFiltered('p-contig-len','contig_len',s);rfFiltered('p-coverage','coverage',s);
}});
makeSampleDropdown('sample-ctrl-viral', s=>{{
  rfFiltered('p-tool-counts','tool_counts',s);rfFiltered('p-viral','viral',s);
  rfFiltered('p-scores','scores',s);rfFiltered('p-tool-heatmap','tool_heatmap',s);
  rfFiltered('p-checkv','checkv',s);rfFiltered('p-checkv-vrh','checkv_vrh',s);
  rfFiltered('p-checkv-scatter','checkv_scatter',s);
  rfFiltered('p-viral-len','viral_len',s);rfFiltered('p-viral-depth','viral_depth',s);
  rfFiltered('p-vrhyme','vrhyme',s);rfFiltered('p-vrhyme-detail','vrhyme_detail',s);
  renderVibrantSection(s);
  renderVotuSection(s);
}});
makeSampleDropdown('sample-ctrl-taxonomy', s=>{{
  renderTaxFamilyBar(s);renderTaxSunburst(s);renderTaxSourcePie(s);
  renderProkDomainBar(s);renderProkTopPhyla(s);renderProkMasterTable(s);
  makeTable('tbl-vcontact',
    (s==='__all__'?TAX_DATA:TAX_DATA.filter(r=>r.sample===s)),
    ['Source','Genome','Best_taxonomy','final_family','Confidence','CheckV_quality','Completeness','sample'],
    {{Source:v=>sourceBadge(v),CheckV_quality:v=>qualBadge(v),
      Genome:v=>`<span title="${{v}}" style="display:inline-block;max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;vertical-align:bottom">${{v}}</span>`}},
    {{Best_taxonomy:'Best Taxonomy',final_family:'Family',CheckV_quality:'CheckV Quality',Completeness:'Complete %',sample:'Sample'}});
}});
makeSampleDropdown('sample-ctrl-hostpred', s=>{{
  renderHostPredSection(s);
}});
makeSampleDropdown('sample-ctrl-bins', s=>{{
  rfFiltered('p-binner-total','binner_total',s);rfFiltered('p-das-tax','das_tax',s);
  rfFiltered('p-checkm2','checkm2',s);rfFiltered('p-cm2-hist','cm2_hist',s);
  rfFiltered('p-bin-size','bin_size',s);
  renderMimagSection(s);
}});

// ── VIBRANT AMG results ───────────────────────────────────────────────────────
function renderVibrantSection(sample){{
  if(!VIBRANT_AMG.length)return;
  const amgRows=sample==='__all__'?VIBRANT_AMG:VIBRANT_AMG.filter(r=>r.sample===sample);
  const sm=sample==='__all__'?SAMPLES:[sample];
  const bys={{}};SAMPLES.forEach(s=>bys[s]=0);
  amgRows.forEach(r=>{{if(r.sample)bys[r.sample]=(bys[r.sample]||0)+1;}});
  Plotly.newPlot('p-vibrant-amg',[{{type:'bar',x:sm,
    y:sm.map(s=>bys[s]||0),text:sm.map(s=>bys[s]||0),textposition:'auto',
    marker:{{color:'#7c3aed'}}}}],
    {{yaxis:{{title:'AMG count'}},margin:{{t:10,b:80}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
  const pc={{}};
  amgRows.forEach(r=>{{
    const p=r.Pathway||'Unknown';const n=parseInt(r.Total_AMGs)||1;pc[p]=(pc[p]||0)+n;
  }});
  const topP=Object.entries(pc).sort((a,b)=>b[1]-a[1]).slice(0,15);
  if(topP.length){{
    Plotly.newPlot('p-vibrant-pathways',[{{type:'bar',
      x:topP.map(x=>x[0]),y:topP.map(x=>x[1]),
      marker:{{color:'#7c3aed'}}}}],
      {{margin:{{t:10,b:160,l:50,r:10}},
        xaxis:{{tickangle:-45,title:'Pathway'}},
        yaxis:{{title:'AMG count'}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
  }}
  makeTable('tbl-vibrant-amg',amgRows,
    ['sample','Pathway','Metabolism','KEGG_map','Total_AMGs','KOs']);
}}
renderVibrantSection('__all__');

// ── BLOCO 8: vOTU summary section ─────────────────────────────────────────────
function renderVotuSection(sample){{
  const allRows=[];
  if(sample==='__all__'){{SAMPLES.forEach(s=>{{(VOTU_DATA[s]||[]).forEach(r=>{{r._s=s;allRows.push(r);}});}});}}
  else{{(VOTU_DATA[sample]||[]).forEach(r=>{{r._s=sample;allRows.push(r);}});}}
  const total=allRows.length;
  const hq=allRows.filter(r=>(r.checkv_quality||'').match(/High|Complete/i)).length;
  const classified=allRows.filter(r=>r.taxonomy_family&&r.taxonomy_family!=='NA'&&r.taxonomy_family!=='').length;
  const withHost=allRows.filter(r=>r.host&&r.host!=='NA'&&r.host!=='').length;
  const lytic=allRows.filter(r=>(r.lifestyle||'').toLowerCase().includes('lytic')).length;
  const el=document.getElementById('votu-summary-cards');
  if(el){{
    const sc=(v,lbl)=>`<div style="background:var(--bg2);border:1px solid var(--brd);border-radius:8px;padding:12px 18px;text-align:center;min-width:100px"><div style="font-size:1.5rem;font-weight:700;color:var(--acc)">${{v}}</div><div style="font-size:.72rem;color:var(--txt2)">${{lbl}}</div></div>`;
    el.innerHTML=`<div style="display:flex;gap:12px;flex-wrap:wrap;padding:6px 0">
      ${{sc(total,'Total vOTUs')}}${{sc(hq,'High-quality (CheckV)')}}${{sc(classified,'Classified (family)')}}${{sc(withHost,'With host')}}${{sc(lytic,'Lytic (VIBRANT)')}}
    </div>`;
  }}
  const cols=['vOTU_id','n_members','checkv_quality','checkv_completeness','genome_type','lifestyle','taxonomy_family','taxonomy_source','n_AMGs','host'];
  makeTable('tbl-votu',allRows,cols,{{checkv_quality:v=>qualBadge(v)}});
}}
renderVotuSection('__all__');

// ── BLOCO 8: Genome map viewer ────────────────────────────────────────────────
let _gmSample=SAMPLES[0]||'';let _gmMode='phage';
function setGmMode(mode){{
  _gmMode=mode;
  ['phage','virus','prok'].forEach(m=>{{
    const b=document.getElementById('gm-btn-'+m);
    if(b)b.style.background=m===mode?'var(--acc)':'';
  }});
  renderGenomeMapViewer();
}}
function renderGenomeMapViewer(){{
  const maps=(GENOME_MAPS[_gmSample]||{{}})[_gmMode]||[];
  const sel=document.getElementById('gm-genome-sel');
  const viewer=document.getElementById('genome-map-viewer');
  if(!sel||!viewer)return;
  sel.innerHTML='';
  if(!maps.length){{
    viewer.innerHTML='<p style="padding:30px;color:var(--txt3);text-align:center">No genome maps available for this selection.<br>Run the pipeline with pharokka_db / phold_db / bakta_db configured.</p>';
    return;
  }}
  maps.forEach((m,i)=>{{
    const o=document.createElement('option');o.value=i;o.text=m.id;sel.appendChild(o);
  }});
  function showMap(idx){{
    const item=maps[idx];
    if(!item||!item.svg){{viewer.innerHTML='<p style="padding:20px;color:var(--txt3)">No SVG available</p>';return;}}
    viewer.innerHTML=item.svg;
    const svgEl=viewer.querySelector('svg');
    if(svgEl){{svgEl.style.maxWidth='100%';svgEl.style.height='auto';}}
  }}
  sel.onchange=()=>showMap(parseInt(sel.value));
  showMap(0);
}}
(function(){{
  const sampleSel=document.getElementById('gm-sample-sel');
  if(sampleSel){{
    SAMPLES.forEach(s=>{{const o=document.createElement('option');o.value=s;o.text=s;sampleSel.appendChild(o);}});
    sampleSel.onchange=()=>{{_gmSample=sampleSel.value;renderGenomeMapViewer();}};
  }}
  setGmMode('phage');
}})();

// ── Task 5: Viral quality pyramid ─────────────────────────────────────────────
(function(){{
  const tiers=['Complete','High-quality','Medium-quality','Low-quality','Not-determined'];
  const colors=['#16a34a','#4ade80','#fbbf24','#d97706','#6b7280'];
  const allCV=[...Object.values(FIGS)].filter(()=>false);  // fallback — count from TAX_DATA
  const tierCounts={{}};tiers.forEach(t=>tierCounts[t]=0);
  // Use MIMAG data as proxy for quality counts
  SAMPLES.forEach(s=>{{
    const m=MIMAG[s]||{{}};
    tierCounts['High-quality']=(tierCounts['High-quality']||0)+(m.HQ||0);
    tierCounts['Medium-quality']=(tierCounts['Medium-quality']||0)+(m.MQ||0);
    tierCounts['Low-quality']=(tierCounts['Low-quality']||0)+(m.LQ||0);
  }});
  const vals=tiers.map(t=>tierCounts[t]||0);
  const total=vals.reduce((a,b)=>a+b,0)||1;
  if(total<2){{document.getElementById('p-viral-pyramid').innerHTML='<p style="padding:20px;color:var(--txt3);font-size:.8rem">Quality tier data not available</p>';return;}}
  Plotly.newPlot('p-viral-pyramid',[{{
    type:'funnel',y:tiers,x:vals,
    textinfo:'value+percent initial',
    marker:{{color:colors}},
    connector:{{line:{{color:'rgba(0,0,0,0)'}}}},
  }}],{{
    margin:{{t:20,b:20,l:160,r:20}},
    paper_bgcolor:'rgba(0,0,0,0)',
    plot_bgcolor:'rgba(0,0,0,0)',
  }},cfg);
}})();

// ── Novelty metrics ─────────────────────────────────────────────────────────
(function(){{
  const pcts=SAMPLES.map(s=>NOVELTY[s]&&NOVELTY[s].pct_novel||0);
  Plotly.newPlot('p-novelty',[{{type:'bar',x:SAMPLES,y:pcts,text:pcts.map(v=>v+'%'),
    textposition:'auto',marker:{{color:pcts.map(v=>v>70?'#dc2626':v>40?'#d97706':'#16a34a')}}}}],
    {{yaxis:{{title:'% Novel (unclassified)',range:[0,100]}},margin:{{t:10,b:80}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
  Plotly.newPlot('p-tax-depth',[
    {{name:'Classified',type:'bar',x:SAMPLES,y:SAMPLES.map(s=>NOVELTY[s]&&NOVELTY[s].classified||0),marker:{{color:'#0d9488'}}}},
    {{name:'Unclassified',type:'bar',x:SAMPLES,y:SAMPLES.map(s=>NOVELTY[s]&&NOVELTY[s].unclassified||0),marker:{{color:'#dc2626'}}}}],
    {{barmode:'stack',margin:{{t:10,b:80}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
}})();

// ── MIMAG quality tiers ──────────────────────────────────────────────────────
function renderMimagSection(sample){{
  const sm=sample==='__all__'?SAMPLES:[sample];
  Plotly.newPlot('p-mimag',[
    {{name:'High-quality (≥90%/≤5%)', type:'bar',x:sm,y:sm.map(s=>MIMAG[s]&&MIMAG[s].HQ||0),marker:{{color:'#16a34a'}}}},
    {{name:'Medium-quality (≥50%/≤10%)',type:'bar',x:sm,y:sm.map(s=>MIMAG[s]&&MIMAG[s].MQ||0),marker:{{color:'#d97706'}}}},
    {{name:'Low-quality',type:'bar',x:sm,y:sm.map(s=>MIMAG[s]&&MIMAG[s].LQ||0),marker:{{color:'#ef4444'}}}}],
    {{barmode:'stack',yaxis:{{title:'MAGs'}},margin:{{t:10,b:80}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
  const rows=sm.map(s=>{{
    const m=MIMAG[s]||{{HQ:0,MQ:0,LQ:0,total:0}};
    return {{Sample:s,'High-quality':m.HQ,'Medium-quality':m.MQ,'Low-quality':m.LQ,'Total':m.total,
      'HQ %':m.total>0?Math.round(100*m.HQ/m.total)+'%':'—'}};
  }});
  makeTable('tbl-mimag',rows,['Sample','High-quality','Medium-quality','Low-quality','Total','HQ %']);
}}
renderMimagSection('__all__');

// ── Task 2: Unified Viral Taxonomy functions ──────────────────────────────────
function renderTaxSourcePie(sample){{
  if(!TAX_DATA.length)return;
  const rows=sample==='__all__'?TAX_DATA:TAX_DATA.filter(r=>r.sample===sample);
  const sc={{}};
  rows.forEach(r=>{{const s=(r.Source||r.source||'unclassified').toLowerCase();sc[s]=(sc[s]||0)+1;}});
  const labels=Object.keys(sc).map(k=>SOURCE_STYLE[k]?SOURCE_STYLE[k].label:k);
  const colors=Object.keys(sc).map(k=>SOURCE_STYLE[k]?SOURCE_STYLE[k].color:'#6b7280');
  Plotly.newPlot('p-tax-source',[{{type:'pie',labels:labels,values:Object.values(sc),
    hole:0.38,marker:{{colors:colors}},
    hovertemplate:'<b>%{{label}}</b><br>%{{value}} contigs (%{{percent}})<extra></extra>'}}],
    {{margin:{{t:20,b:20,l:20,r:80}},showlegend:true,paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
}}

function renderTaxFamilyBar(sample){{
  if(!TAX_DATA.length)return;
  const rows=sample==='__all__'?TAX_DATA:TAX_DATA.filter(r=>r.sample===sample);
  const fc={{}};
  rows.forEach(r=>{{
    // Use best available taxonomic rank: family > genus (labelled) > order (labelled) > Unclassified
    const f = r.final_family
      || (r.final_genus  ? r.final_genus  + ' (genus)'  : '')
      || (r.final_order  ? r.final_order  + ' (order)'  : '')
      || 'Unclassified';
    fc[f]=(fc[f]||0)+1;
  }});
  const top=Object.entries(fc).sort((a,b)=>b[1]-a[1]).slice(0,15);
  // Colour: teal for family-level, amber for genus/order fallbacks, grey for Unclassified
  const colours=top.map(x=>x[0]==='Unclassified'?'#9ca3af':x[0].includes('(genus)')||x[0].includes('(order)')?'#f59e0b':'#0d9488');
  Plotly.newPlot('p-tax-family',[{{type:'bar',orientation:'h',
    x:top.map(x=>x[1]).reverse(),y:top.map(x=>x[0]).reverse(),
    marker:{{color:colours.slice().reverse()}}}}],
    {{margin:{{t:10,b:40,l:180,r:20}},xaxis:{{title:'Count'}},
      paper_bgcolor:'rgba(0,0,0,0)',
      annotations:[{{xref:'paper',yref:'paper',x:1,y:1.02,xanchor:'right',yanchor:'bottom',
        text:'teal=family · amber=genus/order · grey=unclassified',
        showarrow:false,font:{{size:10,color:'#6b7280'}}}}]}},cfg);
}}

function renderTaxSunburst(sample){{
  try{{
  if(!TAX_DATA||!TAX_DATA.length)return;
  const rows=sample==='__all__'?TAX_DATA:TAX_DATA.filter(r=>r.sample===sample);
  const ROOT='Viruses';const nodeP={{}};const leafC={{}};nodeP[ROOT]='';
  rows.forEach(r=>{{
    const fam=r.final_family||(r.final_genus?r.final_genus+' (genus)':'')
             ||(r.final_order?r.final_order+' (order)':'')||'Unclassified';
    const gen=r.final_genus||r.Genus||'';
    const famId='f::'+fam;
    if(nodeP[famId]===undefined)nodeP[famId]=ROOT;
    if(gen){{const genId='g::'+fam+'::'+gen;if(nodeP[genId]===undefined)nodeP[genId]=famId;leafC[genId]=(leafC[genId]||0)+1;}}
    else{{leafC[famId]=(leafC[famId]||0)+1;}}
  }});
  const allNodes=Object.keys(nodeP);
  const vals={{}};allNodes.forEach(n=>vals[n]=leafC[n]||0);
  allNodes.filter(n=>n.startsWith('g::')).forEach(n=>{{const p=nodeP[n];vals[p]=(vals[p]||0)+(leafC[n]||0);}});
  allNodes.filter(n=>n.startsWith('f::')).forEach(n=>{{vals[ROOT]=(vals[ROOT]||0)+(vals[n]||0);}});
  const ids=allNodes;
  const labels=ids.map(n=>n.startsWith('g::')?n.split('::')[2]:n.startsWith('f::')?n.split('::')[1]:n);
  Plotly.newPlot('p-tax-sunburst',[{{
    type:'sunburst',ids:ids,labels:labels,parents:ids.map(n=>nodeP[n]),
    values:ids.map(n=>vals[n]||0),branchvalues:'total',
    insidetextorientation:'radial',maxdepth:3,
    marker:{{colorscale:'Tealgrn'}},
    hovertemplate:'<b>%{{label}}</b><br>%{{value}} virus sequences<extra></extra>'
  }}],{{margin:{{t:10,b:10,l:10,r:10}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
  }}catch(e){{console.warn('viral sunburst error:',e);}}
}}

// ── Task 3: Prokaryotic taxonomy functions ───────────────────────────────────
function renderProkDomainBar(sample){{
  const rows=sample==='__all__'?MERGED_PROK:MERGED_PROK.filter(r=>r.sample===sample);
  const byS={{}};SAMPLES.forEach(s=>byS[s]={{Bacteria:0,Archaea:0,Unclassified:0}});
  rows.forEach(r=>{{
    const s=r.sample;if(!byS[s])return;
    const d=r.Domain||'Unclassified';
    byS[s][d in byS[s]?d:'Unclassified']++;
  }});
  const sm=sample==='__all__'?SAMPLES:[sample];
  Plotly.newPlot('p-prok-domain',[
    {{name:'Bacteria',   type:'bar',x:sm,y:sm.map(s=>byS[s]&&byS[s].Bacteria||0),    marker:{{color:'#0d9488'}}}},
    {{name:'Archaea',    type:'bar',x:sm,y:sm.map(s=>byS[s]&&byS[s].Archaea||0),     marker:{{color:'#d97706'}}}},
    {{name:'Unclassified',type:'bar',x:sm,y:sm.map(s=>byS[s]&&byS[s].Unclassified||0),marker:{{color:'#6b7280'}}}}],
    {{barmode:'stack',margin:{{t:10,b:80,l:50,r:10}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
}}

function renderProkTopPhyla(sample){{
  const rows=sample==='__all__'?MERGED_PROK:MERGED_PROK.filter(r=>r.sample===sample);
  const pc={{}};
  rows.filter(r=>r.Domain==='Bacteria').forEach(r=>{{const p=r.Phylum||'Unknown';pc[p]=(pc[p]||0)+1;}});
  const top=Object.entries(pc).sort((a,b)=>b[1]-a[1]).slice(0,10);
  Plotly.newPlot('p-prok-phyla',[{{type:'bar',orientation:'h',
    x:top.map(x=>x[1]).reverse(),y:top.map(x=>x[0]).reverse(),
    marker:{{color:'#d97706'}}}}],
    {{margin:{{t:10,b:40,l:180,r:20}},xaxis:{{title:'MAG count'}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
}}

function renderProkMasterTable(sample){{
  const rows=sample==='__all__'?MERGED_PROK:MERGED_PROK.filter(r=>r.sample===sample);
  makeTable('tbl-prok-master',rows,
    ['sample','Bin','Domain','Phylum','Class','Order','Family','Genus','Source_tax','Completeness','Contamination','Genome_size'],
    {{Source_tax:v=>sourceBadge(v)}});
}}

// ── GTDB-Tk sunburst (unchanged logic, new colors) ───────────────────────────
(function(){{
  try{{
  if(!GTDB_DATA||!GTDB_DATA.length)return;
  const ROOT='MAGs';const nodeP={{}};const leafC={{}};nodeP[ROOT]='';
  GTDB_DATA.forEach(r=>{{
    const dom=r.Domain||'Unknown',phy=r.Phylum||'Unknown',cls=r.Class||'Unknown';
    const domId='d::'+dom,phyId='p::'+dom+'::'+phy,clsId='c::'+dom+'::'+phy+'::'+cls;
    if(nodeP[domId]===undefined)nodeP[domId]=ROOT;
    if(nodeP[phyId]===undefined)nodeP[phyId]=domId;
    if(nodeP[clsId]===undefined)nodeP[clsId]=phyId;
    leafC[clsId]=(leafC[clsId]||0)+1;
  }});
  const allNodes=Object.keys(nodeP);const vals={{}};
  allNodes.forEach(n=>vals[n]=leafC[n]||0);
  allNodes.filter(n=>n.startsWith('c::')).forEach(n=>{{const p=nodeP[n];vals[p]=(vals[p]||0)+(leafC[n]||0);}});
  allNodes.filter(n=>n.startsWith('p::')).forEach(n=>{{const p=nodeP[n];vals[p]=(vals[p]||0)+(vals[n]||0);}});
  allNodes.filter(n=>n.startsWith('d::')).forEach(n=>{{vals[ROOT]=(vals[ROOT]||0)+(vals[n]||0);}});
  const ids=allNodes;
  const labels=ids.map(n=>{{const p=n.split('::');return p[p.length-1]||n;}});
  Plotly.newPlot('p-gtdb-sunburst',[{{
    type:'sunburst',ids:ids,labels:labels,
    parents:ids.map(n=>nodeP[n]),values:ids.map(n=>vals[n]||0),
    branchvalues:'total',insidetextorientation:'radial',maxdepth:3,
    marker:{{colorscale:'Tealrose'}},
    hovertemplate:'<b>%{{label}}</b><br>%{{value}} MAGs<extra></extra>'
  }}],{{margin:{{t:10,b:10,l:10,r:10}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
  }}catch(e){{console.warn('gtdb sunburst error:',e);}}
}})();

// ── Custom prok Diamond ──────────────────────────────────────────────────────
(function(){{
  if(!CUSTOM_PROK||!CUSTOM_PROK.length)return;
  document.getElementById('custom-prok-section').style.display='';
  const gc={{}};
  CUSTOM_PROK.forEach(r=>{{const g=r.Genus||'Unknown';gc[g]=(gc[g]||0)+1;}});
  const gn=Object.entries(gc).sort((a,b)=>b[1]-a[1]).slice(0,20);
  Plotly.newPlot('p-custom-prok-genus',[{{type:'bar',orientation:'h',
    x:gn.map(x=>x[1]).reverse(),y:gn.map(x=>x[0]).reverse(),
    marker:{{color:'#0d9488'}}}}],{{margin:{{t:10,b:40,l:160,r:20}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
  const pc={{}};
  CUSTOM_PROK.forEach(r=>{{const p=r.Phylum||'Unknown';pc[p]=(pc[p]||0)+1;}});
  const ph=Object.entries(pc).sort((a,b)=>b[1]-a[1]).slice(0,15);
  Plotly.newPlot('p-custom-prok-phylum',[{{type:'bar',orientation:'h',
    x:ph.map(x=>x[1]).reverse(),y:ph.map(x=>x[0]).reverse(),
    marker:{{color:'#d97706'}}}}],{{margin:{{t:10,b:40,l:180,r:20}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
  makeTable('tbl-custom-prok',CUSTOM_PROK,
    ['sample','Bin','Domain','Phylum','Class','Order','Family','Genus','Organism','Avg_pident']);
}})();

// ── GTDB-Tk table ────────────────────────────────────────────────────────────
makeTable('tbl-gtdb',GTDB_DATA,['sample','Bin','Domain','Phylum','Class','Order','Family','Genus','Species','RED_value']);

// ── Host Prediction — PHIST ───────────────────────────────────────────────────
function renderHostPredSection(sample){{
  const allRows=[];
  const pRows=sample==='__all__'?PHIST_DATA:PHIST_DATA.filter(r=>r.sample===sample);
  pRows.forEach(r=>allRows.push({{
    sample:r.sample, Virus:r.Virus, Host:r.Host||'', Method:'PHIST',
    Score:r.Score||'', Confidence:r.P_value||'',
  }}));

  if(!allRows.length)return;

  const sm=sample==='__all__'?SAMPLES:[sample];
  const pBys={{}};SAMPLES.forEach(s=>pBys[s]=0);pRows.forEach(r=>{{if(r.sample)pBys[r.sample]=(pBys[r.sample]||0)+1;}});
  Plotly.newPlot('p-hostpred-cov',[
    {{name:'PHIST',type:'bar',x:sm,y:sm.map(s=>pBys[s]||0),marker:{{color:'#7c3aed'}}}}],
    {{yaxis:{{title:'Predictions'}},margin:{{t:10,b:80,l:50,r:10}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);

  const hc={{}};
  allRows.forEach(r=>{{const h=r.Host||'Unknown';hc[h]=(hc[h]||0)+1;}});
  const topH=Object.entries(hc).sort((a,b)=>b[1]-a[1]).slice(0,15);
  Plotly.newPlot('p-hostpred-hosts',[{{type:'bar',orientation:'h',
    x:topH.map(x=>x[1]).reverse(),y:topH.map(x=>x[0]).reverse(),
    marker:{{color:'#0d9488'}}}}],
    {{margin:{{t:10,b:40,l:180,r:20}},xaxis:{{title:'Predictions'}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);

  makeTable('tbl-hostpred',allRows,['sample','Virus','Host','Method','Score','Confidence'],
    {{Method:_=>
      '<span class="src-badge" style="background:rgba(124,58,237,.12);color:#7c3aed;border:1px solid rgba(124,58,237,.3)">PHIST</span>'
    }});
}}

// ── Virus-Host Network ────────────────────────────────────────────────────────
(function(){{
  if(!PHIST_DATA.length)return;
  const nodes={{}};const edges=[];
  const addNode=(id,type)=>{{if(!nodes[id])nodes[id]={{id,type,size:1}};else nodes[id].size+=0.5;}};
  PHIST_DATA.forEach(r=>{{if(!r.Virus||!r.Host)return;addNode(r.Virus,'virus');addNode(r.Host,'host');edges.push({{from:r.Virus,to:r.Host,method:'PHIST'}});}});
  const vList=[...new Set(edges.map(e=>e.from))].slice(0,50);
  const hList=[...new Set(edges.map(e=>e.to))].slice(0,30);
  const pos={{}};
  vList.forEach((v,i)=>pos[v]={{x:0.1,y:i/Math.max(vList.length-1,1)}});
  hList.forEach((h,i)=>pos[h]={{x:0.9,y:i/Math.max(hList.length-1,1)}});
  const edgeX=[],edgeY=[];
  edges.filter(e=>pos[e.from]&&pos[e.to]).forEach(e=>{{edgeX.push(pos[e.from].x,pos[e.to].x,null);edgeY.push(pos[e.from].y,pos[e.to].y,null);}});
  Plotly.newPlot('p-vhost-network',[
    {{type:'scatter',x:edgeX,y:edgeY,mode:'lines',line:{{color:'rgba(150,150,150,0.35)',width:1}},hoverinfo:'none',showlegend:false}},
    {{type:'scatter',x:vList.filter(v=>pos[v]).map(v=>pos[v].x),y:vList.filter(v=>pos[v]).map(v=>pos[v].y),
      mode:'markers+text',name:'Virus',marker:{{size:10,color:'#7c3aed'}},
      text:vList.filter(v=>pos[v]).map(v=>v.length>22?v.slice(0,20)+'…':v),
      textposition:'middle left',textfont:{{size:8}},hovertext:vList.filter(v=>pos[v]),hoverinfo:'text'}},
    {{type:'scatter',x:hList.filter(h=>pos[h]).map(h=>pos[h].x),y:hList.filter(h=>pos[h]).map(h=>pos[h].y),
      mode:'markers+text',name:'Host MAG',marker:{{size:12,color:'#0d9488',symbol:'diamond'}},
      text:hList.filter(h=>pos[h]).map(h=>h.length>22?h.slice(0,20)+'…':h),
      textposition:'middle right',textfont:{{size:8}},hovertext:hList.filter(h=>pos[h]),hoverinfo:'text'}}],
    {{showlegend:true,xaxis:{{visible:false,range:[-0.15,1.15]}},yaxis:{{visible:false}},
      margin:{{t:10,b:10,l:130,r:130}},paper_bgcolor:'rgba(0,0,0,0)'}},cfg);
}})();

// ── Performance panel ─────────────────────────────────────────────────────────
(function(){{
  if(!BENCH_DATA||!BENCH_DATA.length){{
    ['p-bench-bar','p-bench-heatmap','bench-tbl'].forEach(id=>{{
      const el=document.getElementById(id);
      if(el)el.innerHTML='<p style="padding:20px;color:var(--txt3);font-size:.8rem">No benchmark data available yet. Run the full pipeline to generate timing data.</p>';
    }});
    return;
  }}
  rf('p-bench-bar','bench_bar');
  rf('p-bench-heatmap','bench_heatmap');
  makeTable('bench-tbl', BENCH_DATA,
    ['rule','sample','wall_min','cpu_min','max_rss_mb','mean_load'],
    {{
      wall_min: v=>{{const n=parseFloat(v);return isNaN(n)?v:`<span style="color:${{n>60?'#ef4444':n>10?'#d97706':'inherit'}}">${{n.toFixed(1)}}</span>`;}},
      max_rss_mb: v=>{{const n=parseFloat(v);return isNaN(n)?v:`<span style="color:${{n>50000?'#ef4444':n>10000?'#d97706':'inherit'}}">${{n.toFixed(0)}}</span>`;}},
    }});
}})();

// ── Initial render of all JS-driven sections ─────────────────────────────────
renderTaxSourcePie('__all__');
renderTaxFamilyBar('__all__');
renderTaxSunburst('__all__');
makeTable('tbl-vcontact',TAX_DATA,
  ['Source','Genome','Best_taxonomy','final_family','Confidence','CheckV_quality','Completeness','sample'],
  {{Source:v=>sourceBadge(v),CheckV_quality:v=>qualBadge(v),
    Genome:v=>`<span title="${{v}}" style="display:inline-block;max-width:220px;overflow:hidden;text-overflow:ellipsis;white-space:nowrap;vertical-align:bottom">${{v}}</span>`}},
  {{Best_taxonomy:'Best Taxonomy',final_family:'Family',CheckV_quality:'CheckV Quality',Completeness:'Complete %',sample:'Sample'}});
renderProkDomainBar('__all__');
renderProkTopPhyla('__all__');
renderProkMasterTable('__all__');
renderHostPredSection('__all__');

// ── Traffic-light helper ──────────────────────────────────────────────────────
function tlClass(val,hi,lo,higher){{
  if(hi===null)return'sc-neutral';
  const n=parseFloat(String(val).replace(/[^0-9.]/g,''));
  if(isNaN(n))return'sc-neutral';
  if(higher)return n>=hi?'sc-ok':n>=lo?'sc-warn':'sc-bad';
  return n<=hi?'sc-ok':n<=lo?'sc-warn':'sc-bad';
}}
function scCard(val,label,tlArgs){{
  const cls=tlArgs?tlClass(val,...tlArgs):'sc-neutral';
  return`<div class="sc ${{cls}}"><div class="v">${{val}}</div><div class="l">${{label}}</div></div>`;
}}

// ── Overview ──────────────────────────────────────────────────────────────────
const OV_GROUPS=[
  {{label:'Read QC',color:'var(--amb)',metrics:[
    {{key:'total_raw_reads',label:'Raw reads (R1)',fmt:'int',tl:null}},
    {{key:'mean_qual',label:'Mean Phred (trimmed R1)',fmt:'str',tl:[30,25,true]}},
    {{key:'gc_pct',label:'GC content (trimmed R1)',fmt:'str',tl:null}},
    {{key:'mapping_rate',label:'Reads mapped to assembly',fmt:'pct',tl:[70,50,true]}},
  ]}},
  {{label:'Assembly',color:'var(--acc)',metrics:[
    {{key:'n_contigs',label:'Final contigs (deduplicated)',fmt:'int',tl:null}},
    {{key:'n50',label:'N50 (bp)',fmt:'int',tl:[10000,3000,true]}},
  ]}},
  {{label:'Viral Analysis',color:'var(--pur)',metrics:[
    {{key:'viral_consensus',label:'Viral consensus (≥{min_viral_tools} tools)',fmt:'int',tl:[10,1,true]}},
    {{key:'complete_viral',label:'Complete vOTUs (CheckV)',fmt:'int',tl:[1,1,true]}},
    {{key:'vmags',label:'vMAGs (vRhyme)',fmt:'int',tl:null}},
    {{key:'unbinned_viral',label:'Unbinned viral contigs',fmt:'int',tl:null}},
    {{key:'taxonomy_classified',label:'Viral seqs w/ taxonomy',fmt:'int',tl:[1,1,true]}},
    {{key:'host_pred_total',label:'Host predictions (PHIST)',fmt:'int',tl:null}},
  ]}},
  {{label:'Prokaryotic MAGs',color:'var(--grn)',metrics:[
    {{key:'total_bins',label:'Final bins (Binette)',fmt:'int',tl:[5,1,true]}},
    {{key:'hq_bins',label:'HQ MAGs (≥90% / ≤5%)',fmt:'int',tl:[5,1,true]}},
    {{key:'bacteria_bins',label:'Bacterial MAGs (GTDB)',fmt:'int',tl:null}},
    {{key:'archaea_bins',label:'Archaeal MAGs (GTDB)',fmt:'int',tl:null}},
  ]}},
];

function fmt(v){{if(v===undefined||v===null||v===""||v==="N/A")return"—";const s=String(v).replace(/,/g,"");const n=parseInt(s.replace(/[^0-9]/g,""));return isNaN(n)?s:n.toLocaleString();}}
function fmtv(v){{return(v===undefined||v===null||v===""||v==="N/A")?"—":v;}}

(function(){{
  try{{
  const sumEl=document.getElementById('run-summary');if(!sumEl)return;
  const ns=SAMPLES.length;
  const totViral=SAMPLES.reduce((a,s)=>a+(OVERVIEW[s]&&OVERVIEW[s].viral_consensus||0),0);
  const totHQ=SAMPLES.reduce((a,s)=>a+(OVERVIEW[s]&&OVERVIEW[s].hq_bins||0),0);
  const totBins=SAMPLES.reduce((a,s)=>a+(OVERVIEW[s]&&OVERVIEW[s].total_bins||0),0);
  const totComp=SAMPLES.reduce((a,s)=>a+(OVERVIEW[s]&&OVERVIEW[s].complete_viral||0),0);
  const totTax=SAMPLES.reduce((a,s)=>a+(OVERVIEW[s]&&OVERVIEW[s].taxonomy_classified||0),0);
  const hqPct=totBins>0?Math.round(100*totHQ/totBins):0;
  sumEl.innerHTML=`<span class="rs-label">Summary</span>`+
    `Analysed <b>${{ns}}</b> sample${{ns!==1?'s':''}}. `+
    `Viral analysis identified <b>${{totViral.toLocaleString()}}</b> consensus viral sequence${{totViral!==1?'s':''}} `+
    `(<b>${{totComp}}</b> complete; <b>${{totTax}}</b> classified). `+
    `Prokaryotic binning recovered <b>${{totBins}}</b> MAG${{totBins!==1?'s':''}} `+
    `(<b>${{totHQ}}</b> high-quality, <b>${{hqPct}}%</b> of total).`;
  }}catch(e){{console.warn('run summary error:',e);}}
}})();

(function(){{
  try{{
  const g=document.getElementById("ovg");if(!g)return;
  SAMPLES.forEach(s=>{{
    const o=OVERVIEW[s]||{{}};
    g.innerHTML+=`<div style="grid-column:1/-1;background:var(--bg2);border:1px solid var(--brd);border-left:3px solid var(--acc);border-radius:7px;padding:9px 16px;font-weight:600;font-size:.82rem;color:var(--acc);margin-top:8px">${{s}}</div>`;
    OV_GROUPS.forEach(grp=>{{
      g.innerHTML+=`<div class="ov-group-hdr" style="grid-column:1/-1;color:${{grp.color}}">${{grp.label}}</div>`;
      grp.metrics.forEach(m=>{{
        let val=m.fmt==='int'?fmt(o[m.key]):fmtv(o[m.key]);
        g.innerHTML+=scCard(val,m.label,m.tl);
      }});
    }});
  }});
  }}catch(e){{console.warn('overview grid error:',e);}}
}})();
</script></div></div></body></html>"""

    os.makedirs(os.path.dirname(out_html), exist_ok=True)
    with open(out_html,"w",encoding="utf-8") as f: f.write(html)
    print(f"[generate_report] done -> {out_html}")

except Exception:
    import datetime
    log_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '..', 'results', 'generate_report_error.log')
    with open(log_path, 'w') as _ef:
        _ef.write(f"[{datetime.datetime.now()}] FATAL ERROR\n")
        traceback.print_exc(file=_ef)
    print(f"\n[generate_report] FATAL ERROR - see {log_path}", file=sys.stderr)
    traceback.print_exc(file=sys.stderr)
    sys.exit(1)
