#!/usr/bin/env python3
"""
genome_map_universal.py — Circular genome maps for phages, viruses, and prokaryotes.

Three modes (--mode):
  phage      Phold/Pharokka GBK → PHROGS functional category colors
  virus      GeNomad genes TSV + FASTA → VOGDB category colors
  prok       Bakta GBK → COG-based colors (multi-contig aware)

Batch invocation (called by Snakemake rules):
  python3 genome_map_universal.py --mode phage \\
      --gbk annotation/phold/phold.gbk \\
      --fasta viral_nonredundant.fasta \\
      --checkv checkv/quality_summary.tsv \\
      --outdir annotation/genome_maps/phage \\
      --min-completeness 90 --top-n 5

  python3 genome_map_universal.py --mode virus \\
      --fasta viral_nonredundant.fasta \\
      --genomad-genes viral/genomad/sample_genes.tsv \\
      --checkv checkv/quality_summary.tsv \\
      --exclude-ids pharokka_hq_ids.txt \\
      --outdir annotation/genome_maps/virus \\
      --min-completeness 90 --top-n 5

  python3 genome_map_universal.py --mode prok \\
      --bakta-dir annotation/bakta \\
      --checkm2 bins/checkm2/quality_report.tsv \\
      --outdir annotation/genome_maps/prok \\
      --min-completeness 90 --max-contamination 5 --max-contigs 5 --top-n 5

Outputs per genome (in --outdir):
  {genome_id}_map.svg   inline-embeddable SVG for HTML report
  {genome_id}_map.pdf   vector PDF for publication
"""

import argparse
import csv
import os
import sys
from pathlib import Path

# ── Optional heavy deps — fail gracefully so pipeline can still touch done file ──
try:
    import matplotlib as mpl
    import matplotlib.patches as mpatches
    import matplotlib.pyplot as plt
    import numpy as np
    from Bio import SeqIO
    from pycirclize import Circos
    from pycirclize.parser import Genbank
    HAS_DEPS = True
except ImportError as exc:
    print(f"[genome_map] WARNING: Missing dependency ({exc}) — maps will be skipped",
          file=sys.stderr)
    HAS_DEPS = False


# ══════════════════════════════════════════════════════════════════════
#  Color schemes (colorblind-safe Okabe-Ito + Wong palette)
# ══════════════════════════════════════════════════════════════════════

# Phage mode — PHROGS functional categories
PHROGS_COLORS = {
    "head and packaging":                               "#0072B2",
    "tail":                                             "#009E73",
    "DNA, RNA and nucleotide metabolism":               "#E69F00",
    "lysis":                                            "#D55E00",
    "connector":                                        "#CC79A7",
    "transcription regulation":                         "#56B4E9",
    "moron, auxiliary metabolic gene and host takeover": "#F0E442",
    "other":                                            "#999999",
    "unknown function":                                 "#DDDDDD",
}

# Virus mode — VOGDB functional categories (two-letter codes)
VOGDB_COLORS = {
    "Xs": "#009E73",   # Structural/packaging — green
    "Xr": "#0072B2",   # Replication/repair — blue
    "Xp": "#56B4E9",   # Transcription — light blue
    "Xm": "#E69F00",   # Metabolism (AMG) — orange
    "Xh": "#D55E00",   # Host interaction — red-orange
    "Xl": "#CC79A7",   # Lysis/release — pink
    "Xi": "#F0E442",   # Integration/excision — yellow
    "Xu": "#DDDDDD",   # Unknown — gray
}

VOGDB_LABELS = {
    "Xs": "Xs — Structural/packaging",
    "Xr": "Xr — Replication/repair",
    "Xp": "Xp — Transcription",
    "Xm": "Xm — Metabolism (AMG)",
    "Xh": "Xh — Host interaction",
    "Xl": "Xl — Lysis/release",
    "Xi": "Xi — Integration/excision",
    "Xu": "Xu — Unknown",
}

# Keyword → VOGDB category inference from annotation_description + marker_activity.
# Applied when GeNomad only provides VOG accession IDs (no category code in TSV).
_VOGDB_KEYWORDS = [
    ("Xs", ["capsid", "portal", "terminase", "tail", "baseplate", "spike",
             "collar", "tube", "fiber", "sheath", "neck", "whisker",
             "packaging", "procapsid", "coat protein", "structural",
             "virion", "head protein", "major capsid", "minor capsid"]),
    ("Xr", ["polymerase", "helicase", "primase", "replication",
             "nuclease", "exonuclease", "endonuclease", "dna ligase",
             "methyltransferase", "dna repair", "dnab", "dnac"]),
    ("Xp", ["rna polymerase", "sigma factor", "transcription factor",
             "transcriptional regulator", "anti-sigma", "rna-binding"]),
    ("Xm", ["thymidylate", "ribonucleotide", "thymidine kinase",
             "dihydroorotate", "nicotinamide", "photosystem",
             "auxiliary metabolic", "coenzyme", "metabol"]),
    ("Xh", ["anti-crispr", "acr", "superinfection exclusion",
             "anti-restriction", "host takeover"]),
    ("Xl", ["endolysin", "lysin", "lysis", "holin", "spanin",
             "murein hydrolase", "amidase", "glucosaminidase"]),
    ("Xi", ["integrase", "excisionase", "site-specific recombinase",
             "excision", "excisase"]),
]


def _vogdb_category(description, marker_activity=""):
    """Infer VOGDB two-letter category from annotation description + marker_activity."""
    d = (description or "").lower()
    for cat, keywords in _VOGDB_KEYWORDS:
        if any(k in d for k in keywords):
            return cat
    # GeNomad marks viral hallmark genes (capsid, terminase, etc.) explicitly
    if "viral_hallmark" in (marker_activity or "").lower():
        return "Xs"
    return "Xu"


# Prokaryote mode — COG functional category letter codes
COG_COLORS = {
    "J": "#0072B2",   # Translation / ribosomal biogenesis — blue
    "K": "#56B4E9",   # Transcription — light blue
    "L": "#009E73",   # Replication, repair, recombination — green
    "D": "#E69F00",   # Cell cycle / cell division — orange
    "M": "#CC79A7",   # Cell wall / membrane biogenesis — pink
    "C": "#D55E00",   # Energy production and conversion — red-orange
    "E": "#F0E442",   # Amino acid transport and metabolism — yellow
    "G": "#8B4513",   # Carbohydrate metabolism — brown
    "P": "#FF8C00",   # Inorganic ion transport — dark orange
    "O": "#9400D3",   # Protein turnover / chaperones — purple
    "I": "#32CD32",   # Lipid metabolism — lime green
    "H": "#FF6B6B",   # Coenzyme metabolism — salmon
    "S": "#DDDDDD",   # Unknown / hypothetical — gray
    "R": "#BBBBBB",   # General function prediction — light gray
    "other": "#AAAAAA",
}

TRNA_COLOR   = "#9932CC"
CRISPR_COLOR = "#FF8C00"


# ══════════════════════════════════════════════════════════════════════
#  GC helpers
# ══════════════════════════════════════════════════════════════════════

def _gc_content_windows(seq, window=500, step=50):
    pos, vals = [], []
    n = len(seq)
    for s in range(0, n - window + 1, step):
        sub = seq[s:s + window].upper()
        gc = (sub.count("G") + sub.count("C")) / len(sub)
        pos.append(s + window // 2)
        vals.append(gc)
    pos.append(n)
    vals.append(vals[0] if vals else 0.0)
    return np.array(pos, dtype=float), np.array(vals, dtype=float)


def _gc_skew_windows(seq, window=500, step=50):
    pos, vals = [], []
    n = len(seq)
    for s in range(0, n - window + 1, step):
        sub = seq[s:s + window].upper()
        g, c = sub.count("G"), sub.count("C")
        skew = (g - c) / (g + c) if (g + c) > 0 else 0.0
        pos.append(s + window // 2)
        vals.append(skew)
    pos.append(n)
    vals.append(vals[0] if vals else 0.0)
    return np.array(pos, dtype=float), np.array(vals, dtype=float)


def _draw_gc_tracks(sector, seq, R, window=500):
    """Add GC content + GC skew tracks to a Circos sector."""
    gsize = len(seq)
    w = min(window, max(gsize // 20, 100))
    step = max(w // 10, 10)

    gc_pos, gc_val = _gc_content_windows(seq, w, step)
    sk_pos, sk_val = _gc_skew_windows(seq, w, step)
    gc_mean = float(np.mean(gc_val[:-1])) if len(gc_val) > 1 else 0.5

    # GC content — one track, draw rects with r_lim per window
    gc_track = sector.add_track((R["gc_in"], R["gc_out"]))
    gc_track.axis(fc="#F8F8F8", ec="#CCCCCC", lw=0.4)
    gc_dev = gc_val - gc_mean
    gc_abs = np.abs(gc_dev[:-1]).max() or 1e-6
    r_mid = (R["gc_in"] + R["gc_out"]) / 2
    r_rng = (R["gc_out"] - R["gc_in"]) / 2
    for i in range(len(gc_pos) - 1):
        v = gc_dev[i]
        if abs(v) < 1e-9:
            continue
        r_v = r_mid + (v / gc_abs) * r_rng
        col = "#D55E00" if v >= 0 else "#0072B2"
        lim = (r_mid, r_v) if v >= 0 else (r_v, r_mid)
        gc_track.rect(gc_pos[i], gc_pos[i + 1], r_lim=lim, fc=col, ec="none", alpha=0.7)
    gc_track.line(gc_pos, np.full_like(gc_pos, r_mid), lw=0.3, color="#888", ls="dashed")

    # GC skew — one track, same pattern
    sk_track = sector.add_track((R["sk_in"], R["sk_out"]))
    sk_track.axis(fc="#F8F8F8", ec="#CCCCCC", lw=0.4)
    sk_abs = np.abs(sk_val[:-1]).max() or 1e-6
    r_mid2 = (R["sk_in"] + R["sk_out"]) / 2
    r_rng2 = (R["sk_out"] - R["sk_in"]) / 2
    for i in range(len(sk_pos) - 1):
        v = sk_val[i]
        if abs(v) < 1e-9:
            continue
        r_v = r_mid2 + (v / sk_abs) * r_rng2
        col = "#009E73" if v >= 0 else "#CC79A7"
        lim = (r_mid2, r_v) if v >= 0 else (r_v, r_mid2)
        sk_track.rect(sk_pos[i], sk_pos[i + 1], r_lim=lim, fc=col, ec="none", alpha=0.7)
    sk_track.line(sk_pos, np.full_like(sk_pos, r_mid2), lw=0.3, color="#888", ls="dashed")

    return gc_mean


def _save_fasta(genome_id, seq, outdir):
    """Save individual FASTA file alongside the genome map SVG (used by HTML report 'Copy FASTA')."""
    try:
        fasta_path = os.path.join(outdir, f"{genome_id}.fasta")
        with open(fasta_path, "w") as fh:
            fh.write(f">{genome_id}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i+60] + "\n")
    except Exception as exc:
        print(f"[genome_map] WARNING: could not save FASTA for {genome_id}: {exc}", file=sys.stderr)


def _save_figure(fig, outdir, genome_id, title):
    """Save figure as SVG + PDF."""
    fig.suptitle(title, fontsize=9, fontweight="bold", y=0.97)
    svg_path = os.path.join(outdir, f"{genome_id}_map.svg")
    pdf_path = os.path.join(outdir, f"{genome_id}_map.pdf")
    try:
        fig.savefig(svg_path, format="svg", bbox_inches="tight")
        fig.savefig(pdf_path, format="pdf", bbox_inches="tight")
        print(f"[genome_map] Saved {svg_path}", file=sys.stderr)
    except Exception as exc:
        print(f"[genome_map] WARNING: could not save {genome_id}: {exc}", file=sys.stderr)
    finally:
        plt.close(fig)


def _mpl_setup():
    mpl.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "figure.dpi": 150,
        "savefig.dpi": 150,
        "savefig.bbox": "tight",
    })


# ══════════════════════════════════════════════════════════════════════
#  PHAGE MODE
# ══════════════════════════════════════════════════════════════════════

def draw_phage(gbk_record, seq, genome_id, name, outdir):
    """Draw circular genome map for a phage genome (PHROGS color scheme)."""
    gsize = len(seq)
    if gsize < 2000:
        print(f"[genome_map] Skipping {genome_id}: too short ({gsize} bp)", file=sys.stderr)
        return

    circos = Circos(sectors={genome_id: gsize})
    sector = circos.sectors[0]

    # Scale — placed just inside the GC-skew ring (sk_in=56) for visual proximity
    scale_track = sector.add_track((50, 54))
    scale_track.axis(fc="none", ec="none")
    scale_track.xticks(list(range(0, gsize, 1000)), labels=None,
                       tick_length=1.2, outer=False,
                       line_kws=dict(ec="#888888", lw=0.3))
    scale_track.xticks(list(range(0, gsize, 5000)),
                       labels=[f"{p//1000}kb" if p > 0 else "" for p in range(0, gsize, 5000)],
                       tick_length=2, outer=False, label_size=4,
                       label_orientation="vertical", label_margin=1,
                       line_kws=dict(ec="#333333", lw=0.5))

    R = {"fwd_out": 92, "fwd_in": 86,
         "rev_out": 85, "rev_in": 79,
         "feat_out": 78, "feat_in": 74,
         "gc_out": 73, "gc_in": 65,
         "sk_out": 64, "sk_in": 56}

    fwd_track = sector.add_track((R["fwd_in"], R["fwd_out"]))
    rev_track = sector.add_track((R["rev_in"], R["rev_out"]))
    fwd_track.axis(fc="none", ec="none")
    rev_track.axis(fc="none", ec="none")

    trna_pos, crispr_pos = [], []
    for feat in gbk_record.features:
        s, e = int(feat.location.start), int(feat.location.end)
        if feat.type == "tRNA":
            trna_pos.append((s, e))
        elif feat.type in ("repeat_region", "CRISPR"):
            note = " ".join(feat.qualifiers.get("note", []) + feat.qualifiers.get("rpt_type", []))
            if "CRISPR" in note.upper():
                crispr_pos.append((s, e))

    for feat in gbk_record.features:
        if feat.type != "CDS":
            continue
        s, e   = int(feat.location.start), int(feat.location.end)
        strand = feat.location.strand
        fn     = feat.qualifiers.get("function", ["unknown function"])[0]
        color  = PHROGS_COLORS.get(fn, PHROGS_COLORS["unknown function"])
        deg    = (e - s) / gsize * 360
        hdeg   = max(min(1.5, deg * 0.25), 0.3)
        if strand == 1:
            fwd_track.arrow(s, e, head_length=hdeg, shaft_ratio=0.5,
                            fc=color, ec="none", alpha=1.0)
        else:
            rev_track.arrow(e, s, head_length=hdeg, shaft_ratio=0.5,
                            fc=color, ec="none", alpha=0.85)

    feat_track = sector.add_track((R["feat_in"], R["feat_out"]))
    feat_track.axis(fc="none", ec="none")
    for s, e in trna_pos:
        feat_track.rect(s, e, fc=TRNA_COLOR, ec="none", alpha=0.9)
    for s, e in crispr_pos:
        feat_track.rect(s, e, fc=CRISPR_COLOR, ec="none", alpha=0.9)

    gc_mean = _draw_gc_tracks(sector, seq, R)

    fig = circos.plotfig(figsize=(7, 7))
    ax  = fig.axes[0]

    legend_items = [mpatches.Patch(fc=c, ec="#AAA", lw=0.3, label=k.capitalize())
                    for k, c in PHROGS_COLORS.items()]
    legend_items += [
        mpatches.Patch(fc=TRNA_COLOR,   ec="none", label="tRNA"),
        mpatches.Patch(fc=CRISPR_COLOR, ec="none", label="CRISPR"),
        mpatches.Patch(fc="#D55E00",    ec="none", label=f"GC > {gc_mean*100:.1f}%"),
        mpatches.Patch(fc="#0072B2",    ec="none", label=f"GC < {gc_mean*100:.1f}%"),
        mpatches.Patch(fc="#009E73",    ec="none", label="GC skew G>C"),
        mpatches.Patch(fc="#CC79A7",    ec="none", label="GC skew C>G"),
    ]
    fig.legend(handles=legend_items, loc="lower center", bbox_to_anchor=(0.5, -0.12),
               bbox_transform=fig.transFigure,
               frameon=True, framealpha=0.9, edgecolor="#CCC", fontsize=5.5,
               title="PHROGS category", title_fontsize=6, ncol=3,
               handlelength=1.0, handleheight=0.8)

    _save_figure(fig, outdir, genome_id, f"{name} — Phage genome map ({gsize//1000} kb)")


def batch_phage(args):
    """Batch: iterate phold GBK, select top-N HQ phages, generate maps."""
    gbk_path  = args.gbk
    fasta_path = args.fasta
    checkv_path = args.checkv
    outdir    = args.outdir
    min_comp  = args.min_completeness
    top_n     = args.top_n

    os.makedirs(outdir, exist_ok=True)

    if not gbk_path or not os.path.exists(gbk_path) or os.path.getsize(gbk_path) == 0:
        print("[genome_map] phage: phold GBK missing or empty — skipping", file=sys.stderr)
        return

    # Load CheckV completeness
    comp_map = {}
    if checkv_path and os.path.exists(checkv_path):
        with open(checkv_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                try:
                    comp = float(row.get("completeness", 0) or 0)
                except ValueError:
                    comp = 0.0
                comp_map[row.get("contig_id", "")] = comp

    # Load viral FASTA for sequences
    seqs = {}
    if fasta_path and os.path.exists(fasta_path):
        for rec in SeqIO.parse(fasta_path, "fasta"):
            seqs[rec.id] = str(rec.seq)

    # Parse GBK records and rank by CheckV completeness
    gbk_records = list(SeqIO.parse(gbk_path, "genbank"))
    if not gbk_records:
        print("[genome_map] phage: no records in GBK — skipping", file=sys.stderr)
        return

    ranked = []
    for rec in gbk_records:
        gid = rec.id
        comp = comp_map.get(gid, 0.0)
        if comp >= min_comp:
            ranked.append((gid, comp, rec))
    ranked.sort(key=lambda x: x[1], reverse=True)
    ranked = ranked[:top_n]

    print(f"[genome_map] phage: {len(ranked)} genomes qualify (≥{min_comp}% completeness)",
          file=sys.stderr)

    _mpl_setup()
    for gid, comp, rec in ranked:
        seq = seqs.get(gid, str(rec.seq))
        if not seq:
            seq = str(rec.seq)
        name = f"{gid} (completeness {comp:.1f}%)"
        try:
            draw_phage(rec, seq, gid, name, outdir)
            _save_fasta(gid, seq, outdir)
        except Exception as exc:
            print(f"[genome_map] phage: ERROR drawing {gid}: {exc}", file=sys.stderr)


# ══════════════════════════════════════════════════════════════════════
#  VIRUS MODE (non-phage, uses GeNomad genes TSV)
# ══════════════════════════════════════════════════════════════════════

def _parse_genomad_for_genome(genes_tsv, genome_id):
    """
    Return list of gene dicts for genome_id from GeNomad *_genes.tsv.
    GeNomad genes TSV columns: gene, start, end, length, strand, gc_content,
    genetic_code, rbs_motif, marker_activity, taxonomy, annotation_conjscan,
    annotation_pfam, annotation_vogdb, annotation_description, bitscore, evalue.
    Gene IDs have format: contig_id|N (pipe-delimited).
    """
    genes = []
    if not genes_tsv or not os.path.exists(genes_tsv) or os.path.getsize(genes_tsv) == 0:
        return genes
    with open(genes_tsv) as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            gene_id = row.get("gene", "")
            # gene ID: "contig_id|1" — extract contig
            contig = gene_id.rsplit("|", 1)[0] if "|" in gene_id else gene_id
            if contig != genome_id and not contig.startswith(genome_id + "_"):
                continue
            try:
                start  = int(row.get("start", 0))
                end    = int(row.get("end", 0))
                strand = 1 if row.get("strand", "+") in ("+", "1") else -1
            except (ValueError, TypeError):
                continue
            vogdb = row.get("annotation_vogdb", "") or ""
            marker_activity = row.get("marker_activity", "") or ""
            description = row.get("annotation_description", "") or vogdb
            category = _vogdb_category(description, marker_activity)
            genes.append({"start": start, "end": end, "strand": strand,
                          "category": category, "description": description})
    return genes


def draw_virus(genome_id, seq, genes, name, outdir):
    """Draw circular genome map for a non-phage virus (VOGDB color scheme)."""
    gsize = len(seq)
    if gsize < 2000:
        print(f"[genome_map] Skipping {genome_id}: too short ({gsize} bp)", file=sys.stderr)
        return

    circos = Circos(sectors={genome_id: gsize})
    scale_label = f"{gsize // 1000} kb" if gsize >= 1000 else f"{gsize} bp"
    sector = circos.sectors[0]

    # Scale — placed just inside the GC-skew ring (sk_in=60) for visual proximity
    tick_minor = max(gsize // 20, 500)
    tick_major = max(gsize // 5, 1000)
    scale_track = sector.add_track((53, 57))
    scale_track.axis(fc="none", ec="none")
    scale_track.xticks(list(range(0, gsize, tick_minor)), labels=None,
                       tick_length=1.2, outer=False,
                       line_kws=dict(ec="#888888", lw=0.3))
    major_labels = [f"{p//1000}kb" if p >= 1000 else str(p)
                    for p in range(0, gsize, tick_major)]
    scale_track.xticks(list(range(0, gsize, tick_major)), labels=major_labels,
                       tick_length=2, outer=False, label_size=4,
                       label_orientation="vertical", label_margin=1,
                       line_kws=dict(ec="#333", lw=0.5))

    R = {"fwd_out": 92, "fwd_in": 86,
         "rev_out": 85, "rev_in": 79,
         "gc_out": 77, "gc_in": 69,
         "sk_out": 68, "sk_in": 60}

    fwd_track = sector.add_track((R["fwd_in"], R["fwd_out"]))
    rev_track = sector.add_track((R["rev_in"], R["rev_out"]))
    fwd_track.axis(fc="none", ec="none")
    rev_track.axis(fc="none", ec="none")

    for g in genes:
        s, e  = g["start"], g["end"]
        if e <= s or e > gsize:
            continue
        strand  = g["strand"]
        cat     = g.get("category", "Xu")
        # map any unrecognized category to Xu
        color   = VOGDB_COLORS.get(cat, VOGDB_COLORS["Xu"])
        deg     = (e - s) / gsize * 360
        hdeg    = max(min(1.5, deg * 0.25), 0.3)
        if strand == 1:
            fwd_track.arrow(s, e, head_length=hdeg, shaft_ratio=0.5,
                            fc=color, ec="none", alpha=1.0)
        else:
            rev_track.arrow(e, s, head_length=hdeg, shaft_ratio=0.5,
                            fc=color, ec="none", alpha=0.85)

    gc_mean = _draw_gc_tracks(sector, seq, R)

    fig = circos.plotfig(figsize=(7, 7))

    legend_items = [mpatches.Patch(fc=c, ec="#AAA", lw=0.3, label=VOGDB_LABELS.get(k, k))
                    for k, c in VOGDB_COLORS.items()]
    legend_items += [
        mpatches.Patch(fc="#D55E00", ec="none", label=f"GC > {gc_mean*100:.1f}%"),
        mpatches.Patch(fc="#0072B2", ec="none", label=f"GC < {gc_mean*100:.1f}%"),
        mpatches.Patch(fc="#009E73", ec="none", label="GC skew G>C"),
        mpatches.Patch(fc="#CC79A7", ec="none", label="GC skew C>G"),
    ]
    fig.legend(handles=legend_items, loc="lower center", bbox_to_anchor=(0.5, -0.12),
               bbox_transform=fig.transFigure,
               frameon=True, framealpha=0.9, edgecolor="#CCC", fontsize=5.5,
               title="VOGDB category", title_fontsize=6, ncol=2,
               handlelength=1.0, handleheight=0.8)

    _save_figure(fig, outdir, genome_id, f"{name} — Viral genome map ({scale_label})")


# Taxonomic markers that identify bacteriophages in the viral taxonomy TSV.
# Checked against the 'lineage', 'final_order', and 'genomad_class' columns.
_PHAGE_LINEAGE_MARKERS = {
    "duplodnaviria",    # realm — tailed dsDNA phages
    "heunggongvirae",   # kingdom — Caudovirales clade
    "uroviricota",      # phylum
    "caudoviricetes",   # class (former Caudovirales order)
    "caudovirales",     # order (legacy but still used)
    "monodnaviria",     # realm — ssDNA phages (inoviruses, microviruses, etc.)
}


def _is_bacteriophage(row):
    """Return True if taxonomy columns indicate a bacteriophage."""
    lineage = (row.get("lineage") or "").lower()
    order   = (row.get("final_order") or "").lower()
    gclass  = (row.get("genomad_class") or "").lower()
    genus   = (row.get("final_genus") or "").lower()
    for marker in _PHAGE_LINEAGE_MARKERS:
        if marker in lineage or marker in order or marker in gclass:
            return True
    # genus names ending in 'phage' (e.g. "Lambdavirus" is NOT a phage genus name,
    # but literal "...phage" suffixes are — e.g. "T4-like phage" parsed as genus)
    if genus.endswith("phage"):
        return True
    return False


def _load_viral_taxonomy(tax_tsv):
    """Load viral_taxonomy_merged.tsv.

    Returns (tax_map, phage_ids):
      tax_map   — dict{seq_name: 'Order | Family | Genus'} for map titles
      phage_ids — set of seq_names identified as bacteriophages by taxonomy
    """
    tax_map   = {}
    phage_ids = set()
    if not tax_tsv or not os.path.exists(tax_tsv):
        return tax_map, phage_ids
    try:
        with open(tax_tsv) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                # TSV uses 'seq_name', not 'contig_id'
                gid = (row.get("seq_name") or row.get("contig_id") or "").strip()
                if not gid:
                    continue
                if _is_bacteriophage(row):
                    phage_ids.add(gid)
                # Build display label from deepest resolved rank
                _null = {"", "na", "unknown", "unclassified", "none", "nd"}
                parts = []
                for field in ("final_order", "final_family", "final_genus"):
                    v = (row.get(field) or "").strip()
                    if v and v.lower() not in _null:
                        parts.append(v)
                if parts:
                    tax_map[gid] = " | ".join(parts)
    except Exception as exc:
        print(f"[genome_map] WARNING: could not load viral taxonomy: {exc}", file=sys.stderr)
    return tax_map, phage_ids


def batch_virus(args):
    """Batch: iterate GeNomad genes TSV, select top-N non-phage viral sequences."""
    genomad_genes = args.genomad_genes
    fasta_path = args.fasta
    checkv_path = args.checkv
    outdir     = args.outdir
    min_comp   = args.min_completeness
    top_n      = args.top_n
    exclude_ids = set()

    tax_map, phage_ids = _load_viral_taxonomy(getattr(args, "viral_taxonomy", None))
    print(f"[genome_map] virus: taxonomy phage filter identified {len(phage_ids)} bacteriophages",
          file=sys.stderr)

    os.makedirs(outdir, exist_ok=True)

    # Load IDs to exclude: (1) already annotated by pharokka, (2) identified as
    # bacteriophages by taxonomy — virus mode is strictly for non-phage viruses.
    if args.exclude_ids and os.path.exists(args.exclude_ids):
        with open(args.exclude_ids) as f:
            exclude_ids = {line.strip() for line in f if line.strip()}
    exclude_ids |= phage_ids

    if not fasta_path or not os.path.exists(fasta_path) or os.path.getsize(fasta_path) == 0:
        print("[genome_map] virus: FASTA missing or empty — skipping", file=sys.stderr)
        return

    # Load CheckV completeness
    comp_map = {}
    if checkv_path and os.path.exists(checkv_path):
        with open(checkv_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                try:
                    comp = float(row.get("completeness", 0) or 0)
                except ValueError:
                    comp = 0.0
                comp_map[row.get("contig_id", "")] = comp

    # Rank by completeness, exclude phage IDs
    seqs = {r.id: str(r.seq) for r in SeqIO.parse(fasta_path, "fasta")}
    ranked = []
    for gid, seq in seqs.items():
        if gid in exclude_ids:
            continue
        comp = comp_map.get(gid, 0.0)
        if comp >= min_comp:
            ranked.append((gid, comp, seq))
    ranked.sort(key=lambda x: x[1], reverse=True)
    ranked = ranked[:top_n]

    print(f"[genome_map] virus: {len(ranked)} non-phage viruses qualify (≥{min_comp}%)",
          file=sys.stderr)

    _mpl_setup()
    for gid, comp, seq in ranked:
        genes = _parse_genomad_for_genome(genomad_genes, gid)
        tax_label = tax_map.get(gid, "")
        name = f"{gid} — {tax_label} ({comp:.1f}%)" if tax_label else f"{gid} — Unclassified ({comp:.1f}%)"
        try:
            draw_virus(gid, seq, genes, name, outdir)
            _save_fasta(gid, seq, outdir)
        except Exception as exc:
            print(f"[genome_map] virus: ERROR drawing {gid}: {exc}", file=sys.stderr)


# ══════════════════════════════════════════════════════════════════════
#  PROKARYOTE MODE (Bakta GBK + COG coloring)
# ══════════════════════════════════════════════════════════════════════

# Keyword-based fallback COG assignment from product description
_PRODUCT_KEYWORDS = [
    (["J"],  ["ribosom", "translation", "trna", "rrna", "16s ", "23s ", "5s "]),
    (["K"],  ["transcri", "rna polym", "sigma factor", "promoter"]),
    (["L"],  ["replicat", "repair", "recombina", "helicase", "gyrase",
              "topoisomer", "integrase", "transposase"]),
    (["D"],  ["cell division", "filamenting", "ftsz", "mincd", "parb"]),
    (["M"],  ["cell wall", "peptidoglycan", "lipopolysac", "membrane biosyn",
              "murein", "porin"]),
    (["C"],  ["atp synthase", "cytochrome", "nadh dehydro", "succinate dehyd",
              "oxidase", "reductase", "ferredoxin"]),
    (["E"],  ["amino acid", "transaminase", "aminotransfer", "synthetase",
              "aspartate", "glutamate", "arginine"]),
    (["G"],  ["carbohydrate", "glucose", "sugar", "glycoside", "phosphotransfer",
              "galactose", "xylose"]),
    (["P"],  ["metal transport", "iron", "zinc", "copper", "phosphate transport",
              "sulfate transport", "abc transporter"]),
    (["O"],  ["chaperone", "protease", "hsp70", "groel", "dnaj", "clpb",
              "proteasome", "ubiquitin"]),
    (["I"],  ["fatty acid", "lipid", "acyl-", "phospholipid"]),
    (["H"],  ["coenzyme", "folate", "thiamine", "biotin", "cobalamin",
              "riboflavin", "vitamin"]),
]


def _cog_from_product(product):
    p = product.lower()
    for cats, keywords in _PRODUCT_KEYWORDS:
        if any(k in p for k in keywords):
            return cats[0]
    return "S"


def _cog_from_dbxref(qualifiers):
    """Extract COG letter from Bakta db_xref qualifiers (e.g. 'COG:COG0001')."""
    for xref in qualifiers.get("db_xref", []):
        if xref.startswith("COG:"):
            # COG letter is embedded in COG description elsewhere;
            # for now return None and let product-based fallback handle it
            return None
    return None


def draw_prok(records, genome_id, name, outdir):
    """
    Draw circular genome map for a prokaryotic MAG (COG color scheme).
    records: list of BioPython SeqRecords (one per contig).
    Multi-contig handled as multiple Circos sectors.
    """
    if not records:
        return
    gsize_total = sum(len(r.seq) for r in records)
    if gsize_total < 1000:
        print(f"[genome_map] Skipping {genome_id}: too short ({gsize_total} bp)", file=sys.stderr)
        return

    sectors = {r.id: len(r.seq) for r in records}
    seq_map = {r.id: str(r.seq) for r in records}

    circos = Circos(sectors=sectors)
    scale_label = (f"{gsize_total // 1_000_000:.1f} Mb"
                   if gsize_total >= 1_000_000
                   else f"{gsize_total // 1000} kb")

    R = {"fwd_out": 92, "fwd_in": 86,
         "rev_out": 85, "rev_in": 79,
         "gc_out": 77, "gc_in": 69,
         "sk_out": 68, "sk_in": 60}

    # Build record lookup by ID; iterate over circos.sectors to draw each contig
    rec_by_id = {r.id: r for r in records}
    for sector in circos.sectors:
        rec    = rec_by_id.get(sector.name)
        if rec is None:
            continue
        gsize  = len(rec.seq)
        seq    = seq_map[rec.id]

        fwd_track = sector.add_track((R["fwd_in"], R["fwd_out"]))
        rev_track = sector.add_track((R["rev_in"], R["rev_out"]))
        fwd_track.axis(fc="none", ec="none")
        rev_track.axis(fc="none", ec="none")

        for feat in rec.features:
            if feat.type != "CDS":
                continue
            s, e   = int(feat.location.start), int(feat.location.end)
            strand = feat.location.strand
            product = feat.qualifiers.get("product", ["hypothetical protein"])[0]
            cog_cat = _cog_from_product(product)
            color   = COG_COLORS.get(cog_cat, COG_COLORS["other"])
            deg     = (e - s) / gsize * 360 if gsize > 0 else 0
            hdeg    = max(min(1.5, deg * 0.25), 0.3)
            if strand == 1:
                fwd_track.arrow(s, e, head_length=hdeg, shaft_ratio=0.5,
                                fc=color, ec="none", alpha=1.0)
            else:
                rev_track.arrow(e, s, head_length=hdeg, shaft_ratio=0.5,
                                fc=color, ec="none", alpha=0.85)

        _draw_gc_tracks(sector, seq, R)

    fig = circos.plotfig(figsize=(8, 8))

    legend_items = [mpatches.Patch(fc=c, ec="#AAA", lw=0.3, label=f"{k} - {_COG_LABELS.get(k, k)}")
                    for k, c in COG_COLORS.items() if k not in ("other",)]
    fig.legend(handles=legend_items, loc="lower center", bbox_to_anchor=(0.5, -0.12),
               bbox_transform=fig.transFigure,
               frameon=True, framealpha=0.9, edgecolor="#CCC", fontsize=5,
               title="COG category", title_fontsize=6, ncol=4,
               handlelength=1.0, handleheight=0.8)

    _save_figure(fig, outdir, genome_id, f"{name} — Prokaryote genome map ({scale_label})")


_COG_LABELS = {
    "J": "Translation", "K": "Transcription", "L": "DNA replication/repair",
    "D": "Cell cycle", "M": "Cell wall/membrane", "C": "Energy production",
    "E": "Amino acid metabolism", "G": "Carbohydrate metabolism",
    "P": "Inorganic ion transport", "O": "Protein turnover", "I": "Lipid metabolism",
    "H": "Coenzyme metabolism", "S": "Unknown", "R": "General function",
}


def _count_contigs_in_bakta_gbk(gbk_path):
    """Count GenBank records (= contigs) in a Bakta GBK file."""
    try:
        return sum(1 for _ in SeqIO.parse(gbk_path, "genbank"))
    except Exception:
        return 0


def batch_prok(args):
    """Batch: iterate CheckM2, select top-N HQ MAGs with ≤ max_contigs, generate maps."""
    bakta_dir  = args.bakta_dir
    checkm2    = args.checkm2
    outdir     = args.outdir
    min_comp   = args.min_completeness
    max_cont   = args.max_contamination
    max_ctg    = args.max_contigs
    top_n      = args.top_n

    os.makedirs(outdir, exist_ok=True)

    if not checkm2 or not os.path.exists(checkm2) or os.path.getsize(checkm2) == 0:
        print("[genome_map] prok: CheckM2 report missing — skipping", file=sys.stderr)
        return

    if not bakta_dir or not os.path.isdir(bakta_dir):
        print("[genome_map] prok: Bakta dir missing — skipping", file=sys.stderr)
        return

    # Parse CheckM2 report
    qualifying = []
    with open(checkm2) as f:
        for row in csv.DictReader(f, delimiter="\t"):
            name = row.get("Name", "").strip()
            try:
                comp = float(row.get("Completeness", 0) or 0)
                cont = float(row.get("Contamination", 100) or 100)
            except ValueError:
                continue
            if comp >= min_comp and cont <= max_cont:
                qualifying.append((name, comp, cont))

    qualifying.sort(key=lambda x: x[1], reverse=True)

    drawn = 0
    _mpl_setup()
    for mag_name, comp, cont in qualifying:
        if drawn >= top_n:
            break

        # Find Bakta GBK for this MAG
        gbk_path = os.path.join(bakta_dir, mag_name, f"{mag_name}.gbk")
        if not os.path.exists(gbk_path) or os.path.getsize(gbk_path) == 0:
            continue

        n_ctg = _count_contigs_in_bakta_gbk(gbk_path)
        if n_ctg == 0 or n_ctg > max_ctg:
            print(f"[genome_map] prok: {mag_name} has {n_ctg} contigs > {max_ctg} — skip",
                  file=sys.stderr)
            continue

        try:
            records = list(SeqIO.parse(gbk_path, "genbank"))
            name = f"{mag_name} ({comp:.1f}% comp, {cont:.1f}% cont, {n_ctg} ctg)"
            draw_prok(records, mag_name, name, outdir)
            drawn += 1
        except Exception as exc:
            print(f"[genome_map] prok: ERROR drawing {mag_name}: {exc}", file=sys.stderr)

    print(f"[genome_map] prok: {drawn} MAGs drawn", file=sys.stderr)


# ══════════════════════════════════════════════════════════════════════
#  CLI
# ══════════════════════════════════════════════════════════════════════

def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--mode", required=True, choices=["phage", "virus", "prok"],
                   help="Map mode: phage | virus | prok")
    p.add_argument("--outdir", required=True, help="Output directory for SVG/PDF files")

    # Phage-specific
    p.add_argument("--gbk",    help="GBK file (phold output for phage, Bakta for prok)")
    # Virus-specific
    p.add_argument("--genomad-genes",   help="GeNomad *_genes.tsv (virus mode)")
    p.add_argument("--exclude-ids",     help="File of genome IDs to exclude (virus mode)")
    p.add_argument("--viral-taxonomy",  help="viral_taxonomy_merged.tsv — adds class/family to title (optional, virus mode)")
    # Prok-specific
    p.add_argument("--bakta-dir",         help="Directory containing Bakta per-MAG subdirs")
    p.add_argument("--checkm2",           help="CheckM2 quality_report.tsv (prok mode)")
    p.add_argument("--max-contamination", type=float, default=5.0,
                   help="Max contamination %% for prok mode (default 5.0)")
    p.add_argument("--max-contigs",       type=int,   default=5,
                   help="Max contigs per MAG for prok mode (default 5)")
    # Shared
    p.add_argument("--fasta",            help="FASTA file for sequences")
    p.add_argument("--checkv",           help="CheckV quality_summary.tsv")
    p.add_argument("--min-completeness", type=float, default=90.0,
                   help="Minimum completeness %% (default 90.0)")
    p.add_argument("--top-n",            type=int,   default=5,
                   help="Maximum genomes to map per sample (default 5)")
    return p.parse_args()


def main():
    args = parse_args()

    if not HAS_DEPS:
        print("[genome_map] Dependencies missing — writing empty output dir", file=sys.stderr)
        os.makedirs(args.outdir, exist_ok=True)
        return

    if args.mode == "phage":
        batch_phage(args)
    elif args.mode == "virus":
        batch_virus(args)
    elif args.mode == "prok":
        batch_prok(args)


if __name__ == "__main__":
    main()
