"""data_loaders.py — All data-loading functions for the VAPOR HTML report.

Ported verbatim from generate_report.py with new loaders added for
alpha diversity, PCoA, KEGG, and eggnog aggregation.
"""
import os, re, glob, csv, json, random, zipfile
from collections import Counter, defaultdict


# ── Rule execution status ─────────────────────────────────────────────────────
#
# Rules that soft-fail write their real outcome into done.txt ('ok',
# 'skipped: <reason>' or 'failed: <reason>') instead of touching it empty.
# Without this the report cannot tell "the tool crashed" from "the tool found
# nothing", which is how a disk-full AMRFinderPlus run was once read as a
# biological zero across every sample. A rule whose status is 'failed' must be
# rendered as a gap, never as a count of 0.

STATUS_TRACKED_TOOLS = {
    "amrfinderplus": "bins/amrfinderplus/done.txt",
    "rgi":           "bins/rgi/done.txt",
    "galah_derep":   "bins/derep/done.txt",
    "gtdbtk":        "bins/gtdbtk/done.txt",
}

# Global (non-per-sample) rules tracked the same way, except their done.txt
# lives directly under {outdir}/ with no sample component -- e.g. the global
# vOTU catalog, built once over every sample and co-assembly group pooled
# together. Reported under GLOBAL_STATUS_LABEL, a pseudo-sample key that
# cannot collide with a real sample name (sample discovery never produces a
# name wrapped in parentheses).
STATUS_TRACKED_GLOBAL_TOOLS = {
    "votu_catalog_reps":       "votu_catalog/done.txt",
    "votu_catalog_matrices":   "votu_catalog/matrices_done.txt",
    "bacphlip_votu":           "votu_catalog/bacphlip/done.txt",
    "eggnog_viral":            "votu_catalog/eggnog_viral/done.txt",
    # Migradas para o catalogo global no item "(h)" (2026-08-18) e nao
    # rastreadas ate agora: uma falha em qualquer uma delas nao aparecia em
    # lugar nenhum do relatorio -- a aba simplesmente ficava vazia, que e a
    # mesma aparencia de "nao ha nada a mostrar".
    "votu_prodigal":           "votu_catalog/taxonomy/prodigal_done.txt",
    "votu_mmseqs_taxonomy":    "votu_catalog/taxonomy/mmseqs_inphared_done.txt",
    "votu_mmseqs_tax_custom":  "votu_catalog/taxonomy/mmseqs_custom_viral_done.txt",
    "votu_taxonomy":           "votu_catalog/taxonomy/taxonomy_done.txt",
    "votu_pharokka":           "votu_catalog/annotation/pharokka/done.txt",
    "votu_phold":              "votu_catalog/annotation/phold/done.txt",
    "votu_genome_map_phage":   "votu_catalog/annotation/genome_maps/phage_maps_done.txt",
    "votu_genome_map_virus":   "votu_catalog/annotation/genome_maps/virus_maps_done.txt",
    "votu_defensefinder_viral": "votu_catalog/defensefinder/done.txt",
    "votu_dbapis_viral":       "votu_catalog/dbapis/done.txt",
}

GLOBAL_STATUS_LABEL = "(global)"
# prefixo das chaves de grupo em load_tool_status; evita colisao com sample
GROUP_STATUS_PREFIX = "(coassembly) "


def _read_status_file(path):
    """Read a single done.txt into {'state', 'reason', 'raw'}. Shared by the
    per-sample and global branches of load_tool_status so both parse the
    'ok' / 'skipped: <reason>' / 'failed: <reason>' convention identically."""
    raw = ""
    if os.path.exists(path):
        try:
            with open(path) as fh:
                raw = fh.read().strip()
        except OSError:
            raw = ""
    if not raw:
        return {"state": "unknown", "reason": "no status recorded", "raw": raw}
    head, _, tail = raw.partition(":")
    head = head.strip().lower()
    state = head if head in ("ok", "skipped", "failed") else "unknown"
    reason = tail.strip() or ("" if state == "ok" else raw)
    return {"state": state, "reason": reason, "raw": raw}


def load_tool_status(outdir, samples, groups=()):
    """
    Read done.txt for every status-tracked tool, per sample, plus every
    status-tracked global (non-per-sample) rule under GLOBAL_STATUS_LABEL.

    Returns {sample: {tool: {'state': ..., 'reason': ..., 'raw': ...}}}, with
    an extra "(global)" pseudo-sample key holding the global rules' status.
    State is one of 'ok', 'skipped', 'failed', or 'unknown' — 'unknown' covers
    both a missing done.txt and a legacy empty one written before rules
    recorded status, and must not be presented as success.
    """
    if GLOBAL_STATUS_LABEL in samples:
        raise ValueError(
            f"sample name {GLOBAL_STATUS_LABEL!r} collides with the reserved "
            "pseudo-sample label used for global rule status"
        )

    status = {}
    for sample in samples:
        status[sample] = {}
        for tool, rel in STATUS_TRACKED_TOOLS.items():
            path = os.path.join(outdir, sample, rel)
            status[sample][tool] = _read_status_file(path)

    for group in (groups or ()):
        key = f"{GROUP_STATUS_PREFIX}{group}"
        if key in status:
            raise ValueError(f"co-assembly group key {key!r} collides with a sample name")
        status[key] = {}
        for tool, rel in STATUS_TRACKED_TOOLS.items():
            path = os.path.join(outdir, "coassembly", group, rel)
            status[key][tool] = _read_status_file(path)

    if STATUS_TRACKED_GLOBAL_TOOLS:
        status[GLOBAL_STATUS_LABEL] = {}
        for tool, rel in STATUS_TRACKED_GLOBAL_TOOLS.items():
            path = os.path.join(outdir, rel)
            status[GLOBAL_STATUS_LABEL][tool] = _read_status_file(path)

    return status


def tool_failed(status, sample, tool):
    """True when a tool's count must be shown as a gap rather than as zero."""
    return (status.get(sample, {}).get(tool, {}).get("state")
            in ("failed", "unknown"))


def summarize_tool_status(status):
    """Flat list of non-ok entries, for the report's execution-status section."""
    rows = []
    for sample in sorted(status):
        for tool in sorted(status[sample]):
            entry = status[sample][tool]
            if entry["state"] != "ok":
                rows.append({"sample": sample, "tool": tool,
                             "state": entry["state"], "reason": entry["reason"]})
    return rows


# ── Basic helpers ─────────────────────────────────────────────────────────────

def safe_float(v, d=0.0):
    try: return float(str(v).replace(",", "").replace(">", "").strip())
    except: return d


def safe_int(v, d=0):
    try: return int(str(v).replace(",", "").replace(">", "").strip())
    except: return d


def parse_tsv(path):
    rows = []
    if not path or not os.path.exists(path): return rows
    with open(path) as f:
        hdr = None
        for line in f:
            line = line.rstrip('\n').rstrip('\r')
            if not line: continue
            parts = line.split("\t")
            if hdr is None: hdr = parts; continue
            if not parts: continue
            if len(parts) < len(hdr):
                parts += [''] * (len(hdr) - len(parts))
            rows.append(dict(zip(hdr, parts[:len(hdr)])))
    return rows


def load_tsv(path, skip_empty=True):
    rows = []
    try:
        with open(path) as f:
            rdr = csv.DictReader(f, delimiter='\t')
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
            rdr = csv.DictReader(f)
            for row in rdr:
                rows.append(row)
    except Exception:
        pass
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
    counts = {1: 0, 2: 0, 3: 0, 4: 0}
    for row in parse_tsv(path):
        try:
            n = int(row.get("n_tools", 0))
            if n in counts: counts[n] += 1
        except: pass
    return counts


def parse_support_combos(path):
    """Which *combination* of detectors supports each contig, not just how many.

    The tool-support TSV written by rules/viral_detection.smk carries a `tools`
    column ("GeNomad,VirSorter2"), so the exact set intersections are already on
    disk — parse_support collapses them to a degree count, which cannot answer
    "which detector combinations agree, not just how many". Returns
    {"GeNomad,VirSorter2": n, …} with the tool names sorted inside each key."""
    combos = Counter()
    for row in parse_tsv(path):
        tools = (row.get("tools", "") or "").strip()
        if not tools:
            continue
        key = ",".join(sorted(t.strip() for t in tools.split(",") if t.strip()))
        if key:
            combos[key] += 1
    return dict(combos)


# ── fastp ─────────────────────────────────────────────────────────────────────

def parse_fastp_json(outdir, sample):
    """Parse fastp JSON → {"reads": [...], "trim": {...}}"""
    _empty = {"reads": [], "trim": {
        "reads_in": 0, "reads_written": 0,
        "adapter_r1_pct": 0.0, "adapter_r2_pct": 0.0, "bp_removed_pct": 0.0}}
    json_path = os.path.join(outdir, sample, "qc_raw", f"{sample}_fastp.json")
    if not os.path.exists(json_path): return _empty
    try:
        with open(json_path) as f: j = json.load(f)
    except Exception: return _empty

    reads = []
    for stage_key, stage_lbl in [("before_filtering", "raw"), ("after_filtering", "trimmed")]:
        # gc_content/q30_rate are only present in the aggregate summary block,
        # not in the per-read read{1,2}_{before,after}_filtering blocks.
        summ_stage = j.get("summary", {}).get(stage_key, {})
        gc_pct = safe_float(summ_stage.get("gc_content", 0)) * 100
        q30    = safe_float(summ_stage.get("q30_rate", 0) or 0)
        for read_key, read_lbl in [("read1", "R1"), ("read2", "R2")]:
            rd = j.get(f"{read_key}_{stage_key}") or summ_stage
            reads.append({
                "stage":           stage_lbl,
                "read":            read_lbl,
                "r_label":         f"{stage_lbl} {read_lbl}",
                "total_sequences": safe_int(rd.get("total_reads", 0)),
                "gc_percent":      gc_pct,
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


def parse_total_reads(outdir, sample):
    """Return total read count from flagstat (first line: 'N + 0 in total')."""
    flagstat_path = os.path.join(outdir, sample, "mapping", "flagstat.txt")
    if not os.path.exists(flagstat_path): return 0
    with open(flagstat_path) as f: content = f.read()
    m = re.search(r"^(\d+) \+ \d+ in total", content, re.MULTILINE)
    return int(m.group(1)) if m else 0


# ── Depth / contig data ───────────────────────────────────────────────────────

def collect_depth_data(depth_path):
    lengths, depths = [], []
    if not os.path.exists(depth_path): return lengths, depths
    with open(depth_path) as f:
        hdr = None
        for line in f:
            raw = line.rstrip("\n")
            p = raw.split("\t") if "\t" in raw else raw.split()
            if len(p) < 3: continue
            if hdr is None:
                hdr = p; continue
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
    return counts


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
        total_members += sum(safe_int(r.get("members", 0)) for r in rows)
    if n_bins == 0 and summary_rows:
        n_bins = len(summary_rows)
    return {"n_bins": n_bins, "total_members": total_members, "rows": summary_rows}


def collect_binner_counts(outdir, sample, checkm2_rows):
    b = os.path.join(outdir, sample, "bins")

    mb_count = len(glob.glob(os.path.join(b, "metabat2", "bin.*.fa")))

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

    das_dir  = os.path.join(b, "binette", "final_bins")
    das_bins = glob.glob(os.path.join(das_dir, "*.fa"))
    das = {"total": len(das_bins), "bacteria": 0, "archaea": 0, "unknown": len(das_bins)}

    return {
        "MetaBAT2":        {"total": mb_count,   "bacteria": 0, "archaea": 0, "unknown": 0},
        "VAMB":            {"total": vamb_count,  "bacteria": 0, "archaea": 0, "unknown": 0},
        "SemiBin2":        {"total": sb_count,    "bacteria": 0, "archaea": 0, "unknown": 0},
        "Binette (final)": das,
    }


def parse_checkm2_phyla(checkm2_rows):
    phyla = Counter()
    for row in checkm2_rows:
        tax = row.get("Taxonomic_lineage", row.get("taxonomic_lineage", "")).strip()
        if not tax: phyla["Unclassified"] += 1; continue
        phylum = "Unclassified"
        for p in tax.split(";"):
            if p.strip().startswith("p__"): phylum = p.strip()[3:] or "Unclassified"; break
        phyla[phylum] += 1
    return phyla


# ── Viral taxonomy ────────────────────────────────────────────────────────────

def load_viral_taxonomy(paths, samples):
    # valores nulos que qualquer fonte de taxonomia pode emitir
    _NULL_TAX = {"singleton", "unclassified", "nd", "none"}
    records = []
    for p, s in zip(paths, samples):
        for row in load_tsv(p):
            name = row.get('seq_name', '')
            if not name: continue
            final_family = row.get('final_family', '')
            final_genus  = row.get('final_genus', '')
            final_order  = row.get('final_order', '')
            lineage      = row.get('lineage', '')
            source       = row.get('source', '')
            if final_family.lower() in _NULL_TAX: final_family = ''
            if final_genus.lower()  in _NULL_TAX: final_genus  = ''
            if final_order.lower()  in _NULL_TAX: final_order  = ''
            if source in ('unclassified', ''):
                continue
            if not final_family or not final_genus:
                lf, lg, lo = _deepest_level(lineage)
                if not final_family: final_family = lf
                if not final_genus:  final_genus  = lg
                if not final_order:  final_order  = lo
            if source == 'genomad' and not final_family and not final_genus and not final_order:
                gc = row.get('genomad_class', '') or row.get('genomad_best', '')
                if gc and gc.lower() not in _NULL_TAX:
                    final_order = gc
            if not final_family and not final_genus and not final_order:
                continue
            best_tax = final_genus or final_family or final_order
            records.append({'sample': s,
                'Genome':        name,
                'final_family':  final_family,
                'final_genus':   final_genus,
                'final_order':   final_order,
                'Family':        final_family,
                'Genus':         final_genus,
                'Order':         final_order,
                'Best_taxonomy': best_tax,
                'GeNomad_best':  row.get('genomad_best', ''),
                'GeNomad_class': row.get('genomad_class', ''),
                'Source':        source,
                'Confidence':    row.get('confidence', ''),
                'Lineage':       lineage,
            })
    return records


def load_viral_source_distribution(paths, samples):
    """Per-sample distribution of classification sources for ALL viral contigs
    (including ones dropped from load_viral_taxonomy because they have no
    family/genus/order). Buckets: genomad, mmseqs_inphared, mmseqs_custom
    (IMG/VR or other user-supplied seqTaxDB), unknown (no hit)."""
    _NULL = {"singleton", "unclassified", "nd", "none", ""}
    dist = {}
    for p, s in zip(paths, samples):
        counts = Counter()
        for row in load_tsv(p):
            if not row.get('seq_name', ''):
                continue
            final_family = row.get('final_family', '').strip().lower()
            final_genus  = row.get('final_genus', '').strip().lower()
            final_order  = row.get('final_order', '').strip().lower()
            source       = row.get('source', '').strip()
            custom_lineage = row.get('custom_lineage', '').strip()
            genomad_class = (row.get('genomad_class', '') or row.get('genomad_best', '')).strip()
            if final_family not in _NULL or final_genus not in _NULL or final_order not in _NULL:
                bucket = source if source not in _NULL else 'genomad'
            elif genomad_class.lower() not in _NULL:
                bucket = 'genomad'
            elif custom_lineage.lower() not in _NULL:
                bucket = 'mmseqs_custom'
            else:
                bucket = 'unknown'
            counts[bucket] += 1
        dist[s] = dict(counts)
    return dist


def _deepest_level(lineage):
    if not lineage: return '', '', ''
    parts = [p.strip() for p in lineage.split(';')]
    family = ''; genus = ''; order = ''
    if len(parts) >= 8 and parts[7]: genus  = parts[7]
    if len(parts) >= 7 and parts[6]: family = parts[6]
    if len(parts) >= 6 and parts[5]: order  = parts[5]
    return family, genus, order


# ── GTDB-Tk ───────────────────────────────────────────────────────────────────

def load_gtdbtk(bac_paths, arc_paths, samples):
    records = []
    for bac_p, arc_p, s in zip(bac_paths, arc_paths, samples):
        for fpath in [bac_p, arc_p]:
            for row in load_tsv(fpath):
                classif = row.get('classification', '')
                if not classif or classif in ('N/A', 'NA', ''): continue
                def _g(prefix, c=classif):
                    if f';{prefix}__' in c:
                        v = c.split(f';{prefix}__')[-1].split(';')[0]
                        return v.strip() if v.strip() else ''
                    return ''
                domain = ('Archaea' if 'd__Archaea' in classif
                          else 'Bacteria' if 'd__Bacteria' in classif
                          else 'Unknown')
                bin_name = row.get('user_genome', '').strip()
                if not bin_name: continue
                records.append({'sample': s, 'Bin': bin_name,
                    'Domain': domain, 'Phylum': _g('p'), 'Class': _g('c'),
                    'Order': _g('o'), 'Family': _g('f'), 'Genus': _g('g'),
                    'Species': _g('s'), 'Full_classification': classif,
                    'RED_value': row.get('red_value', ''), 'Note': row.get('note', ''),
                })
    return records


# ── Custom prokaryote taxonomy (MMseqs2 LCA) ──────────────────────────────────

_PROK_RANK_LEVELS = ['Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species', 'Strain']


def _parse_mmseqs_lineage(lineage):
    """mmseqs `--tax-lineage 1` output, e.g.
    '-_Bacteria;p_Pseudomonadota;c_Betaproteobacteria;o_Burkholderiales' --
    rank-prefix codes are positional here (domain/strain have no single-
    letter NCBI rank code, hence '-'), so this is read by position against
    _PROK_RANK_LEVELS, not by the prefix letter. Only emits a path as deep
    as mmseqs itself resolved -- a shallow LCA naturally yields a short list."""
    if not lineage: return []
    names = []
    for tok in lineage.split(';'):
        tok = tok.strip()
        if not tok: continue
        _, _, name = tok.partition('_')
        names.append(name if name else tok)
    return names


def load_mmseqs_taxonomy_prok(paths_d, samples):
    """mmseqs_taxonomy_prok output (qseqid, taxid, rank, name, lineage) --
    mmseqs already computes a real per-PROTEIN lowest-common-ancestor, so
    this aggregates to genome level with a SECOND lowest-common-ancestor
    pass across each genome unit's own proteins, rather than a majority
    vote of best hits -- a vote would just be a softer version of the same
    "spurious specificity" (von Meijenfeldt et al. 2019, CAT/BAT) this path
    exists to avoid in the first place. Unclassified proteins (taxid 0) are
    dropped from the consensus instead of voting for 'no rank'."""
    records = []
    for s in samples:
        p = paths_d.get(s, '')
        if not p or not os.path.exists(p): continue
        lineages = defaultdict(list)
        for row in load_tsv(p):
            unit = _prok_genome_unit(row.get('qseqid', ''))
            if not unit: continue
            taxid = str(row.get('taxid', '0')).strip()
            if taxid in ('0', ''): continue
            lineage = _parse_mmseqs_lineage(row.get('lineage', ''))
            if lineage:
                lineages[unit].append(lineage)
        for unit, lins in lineages.items():
            consensus = list(lins[0])
            for lin in lins[1:]:
                agree = len(consensus)
                for i in range(len(consensus)):
                    if i >= len(lin) or lin[i] != consensus[i]:
                        agree = i
                        break
                consensus = consensus[:agree]
                if not consensus: break
            if not consensus: continue
            rec = {'sample': s, 'Bin': unit, 'Source': 'mmseqs_lca'}
            for i, level in enumerate(_PROK_RANK_LEVELS):
                rec[level] = consensus[i] if i < len(consensus) else ''
            rec['Organism'] = consensus[-1]
            rec['Sci_name'] = consensus[-1]
            records.append(rec)
    return records


# ── PHIST host prediction ─────────────────────────────────────────────────────

_FASTA_EXTS = ('.fasta', '.fna', '.ffn', '.fa')


def _strip_fasta_ext(name):
    """Remove UMA extensao de fasta do fim do nome, a mais longa primeiro.

    O codigo anterior era `.replace('.fa', '').replace('.fasta', '')` e errava
    de duas formas, ambas confirmadas contra dados reais em 2026-08-19:

    - `.fna` nao estava na lista. Os bins do VAMB (track co-assembly,
      rules/coassembly.smk coassembly_phist usa bin_ext=".fna") chamam-se
      "1.fna", entao o Host ficava "1.fna" e a juncao com o `user_genome` do
      GTDB-Tk ("1") falhava para o track inteiro -- taxonomia de hospedeiro
      vazia em todo grupo, sem erro nenhum. O track por amostra escapou por
      acaso: os bins do Binette sao ".fa".
    - `.replace` remove a substring em QUALQUER posicao e `.fa` vinha antes
      de `.fasta`, entao "x.fasta" virava "xsta".
    """
    for ext in _FASTA_EXTS:
        if name.endswith(ext):
            return name[:-len(ext)]
    return name


def load_phist(paths, samples):
    records = []
    for p, s in zip(paths, samples):
        for row in load_csv(p):
            virus_raw = row.get('phage', row.get('Phage', row.get('virus', row.get('Virus', ''))))
            host_raw  = row.get('host',  row.get('Host', ''))
            score     = row.get('#common-kmers', row.get('Score', row.get('score', '')))
            pval      = row.get('adj-pvalue', row.get('pvalue', ''))
            if not virus_raw or not host_raw: continue
            virus_clean = _strip_fasta_ext(os.path.basename(virus_raw))
            virus_clean = virus_clean.replace('contig_', '', 1)
            host_clean  = _strip_fasta_ext(os.path.basename(host_raw))
            records.append({'sample': s, 'Virus': virus_clean, 'Host': host_clean,
                            'Score': str(score), 'P_value': str(pval)})
    return records


def build_host_collapse(phist_data, votu_abund_by_sample, samples,
                        host_links=None):
    """Collapse viral RPKM by predicted host genus per sample.

    Returns {sample: [{genus, n_viruses, total_rpkm}]} sorted desc by total_rpkm.

    O genero vem do GTDB-Tk, via `host_links` (build_host_defense_links, que
    ja resolve Host -> Host_genus). Ate 2026-08-19 vinha de
    `host.split('_')[0]`, com o comentario "primeiro token do nome do arquivo
    do hospedeiro (ex. 'Bacteroides_fragilis.fa' -> 'Bacteroides')". Essa
    premissa vale para genomas de REFERENCIA nomeados por especie -- nao para
    os MAGs desta pipeline. Confirmado no report.html de ~/global/results:
    todo Host e "binette_binN", logo o "genero" de TODOS os hospedeiros era a
    string constante "binette", um unico feixe no grafico. No track
    co-assembly os bins sao numeros ("1.fna") e sairiam como generos "1",
    "136". O dado certo ja estava ao lado: o mesmo hospedeiro aparece em
    HOST_DEFENSE_LINKS com Host_genus "Neisseria".

    Sem `host_links` (chamada antiga) o genero fica 'Unknown' em vez de
    voltar a inventar um a partir do nome do arquivo.
    """
    genus_by_host = {}
    for link in (host_links or ()):
        g = (link.get('Host_genus') or '').strip()
        if g:
            genus_by_host[(link.get('sample', ''), link.get('Host', ''))] = g
    result = {}
    for s in samples:
        abund_by_rep = {}
        for row in votu_abund_by_sample.get(s, []):
            rep = row.get('representative', '')
            # votu_abundance.tsv last column is the coverm method name (e.g. 'rpkm')
            rpkm = 0.0
            for k, v in row.items():
                if k.lower() in ('rpkm', 'tpm', 'mean') and k != 'representative':
                    rpkm = safe_float(v); break
            if rep:
                abund_by_rep[rep] = rpkm

        genus_data = {}
        for row in phist_data:
            if row.get('sample') != s: continue
            virus = row.get('Virus', '')
            host  = row.get('Host', '')
            genus = genus_by_host.get((s, host), '') or 'Unknown'
            rpkm  = abund_by_rep.get(virus, 0.0)
            entry = genus_data.setdefault(genus, {'genus': genus, 'n_viruses': 0, 'total_rpkm': 0.0})
            entry['n_viruses']  += 1
            entry['total_rpkm'] += rpkm

        result[s] = sorted(genus_data.values(), key=lambda x: x['total_rpkm'], reverse=True)
    return result


# ── Defense / anti-defense systems (DefenseFinder) ───────────────────────────

def _split_proteins_in_syst(value):
    """DefenseFinder's 'protein_in_syst' column: comma-separated protein IDs
    of every gene assigned to a system. Used both to recover the originating
    contig (defense_amr.smk concatenates bins as '{genome}__{protein_id}',
    rules.prodigal_viral does not) and for defense-island detection."""
    return [x.strip() for x in (value or '').split(',') if x.strip()]


def _contig_from_protein_id(protein_id):
    """Recover the originating contig/replicon from a Prodigal protein ID.
    Prodigal always names genes '{seqid}_{gene_number}' -- a single trailing
    '_<digits>' -- so the rightmost split recovers seqid even when it itself
    contains underscores (e.g. 'k141_1234_3' -> 'k141_1234')."""
    if not protein_id:
        return ''
    return protein_id.rsplit('_', 1)[0]


def load_defensefinder(paths, samples):
    """DefenseFinder per-genome systems, merged with a 'genome' column by
    rules/defense_amr.smk (genome = bin name, or 'contigs_pseudogenome' in
    the low-depth fallback). 'Proteins' keeps the raw protein_in_syst list
    for defense-island detection (compute_defense_islands)."""
    records = []
    for p, s in zip(paths, samples):
        for row in load_tsv(p):
            genome = row.get('genome', '')
            sys_type = row.get('type', row.get('subtype', ''))
            if not genome or not sys_type: continue
            records.append({'sample': s, 'Bin': genome, 'System': sys_type,
                             'System_id': row.get('sys_id', sys_type),
                             'Genes': row.get('genes_count', ''),
                             'Proteins': _split_proteins_in_syst(row.get('protein_in_syst', ''))})
    return records


def load_antidefensefinder(paths, samples):
    """AntiDefenseFinder systems (DefenseFinder --antidefensefinder pass) —
    same table shape as load_defensefinder, kept in a separate file/loader
    so defense and anti-defense counts are never accidentally merged."""
    records = []
    for p, s in zip(paths, samples):
        for row in load_tsv(p):
            genome = row.get('genome', '')
            sys_type = row.get('type', row.get('subtype', ''))
            if not genome or not sys_type: continue
            records.append({'sample': s, 'Bin': genome, 'System': sys_type,
                             'System_id': row.get('sys_id', sys_type),
                             'Genes': row.get('genes_count', ''),
                             'Proteins': _split_proteins_in_syst(row.get('protein_in_syst', ''))})
    return records


# ── Viral-side anti-defense (Han et al. 2026 cold seep defensome paper) ──────
# DefenseFinder is the same tool/models/container as the host-side rule
# above, just pointed at viral proteins (rule defensefinder_viral). dbAPIS
# is a complementary sequence-similarity detector (rule dbapis_viral) —
# kept in its own loader/JS constant, never merged with the DefenseFinder
# calls (same "never merge tiers" rule as AMR curated/exploratory).

def load_antidefensefinder_viral(path, samples):
    """Anti-defense systems on viral proteins, from the global vOTU catalog's
    DefenseFinder --antidefensefinder pass (rule votu_defensefinder_viral,
    rules/votu_catalog.smk). Moved off per-sample paths on 2026-08-18
    (second half of "(h)", docs/ROADMAP_SIMPLIFICACAO.md): the per-sample
    rule fanned out N+G times over byte-identical input, and worse, its
    per-group output silently held the WHOLE catalog's systems (not just
    that group's), because both already read rules.votu_prodigal's global
    .faa. There is exactly ONE table to read, at
    {outdir}/votu_catalog/defensefinder/viral_antidefense_systems.tsv --
    the same records list is returned under every sample key, same
    duplication trade-off as load_phrogs/load_genome_maps (data_loaders.py),
    so the existing per-sample chart/table plumbing in hostdefense.js keeps
    working unmodified; the chart there is labeled "vOTU catalog --
    global, same across samples".

    The 'genome' column in the source TSV is just the catalog-wide protein
    set label (one DefenseFinder call across the whole catalog, not
    per-genome) -- the actual virus/contig (a NAMESPACED
    "{source_id}|{contig}" ID, since votu_prodigal runs on the namespaced
    catalog FASTA) is recovered from the first protein in protein_in_syst,
    which is reliable since every gene in one system sits on the same
    replicon."""
    records = []
    for row in load_tsv(path):
        sys_type = row.get('type', row.get('subtype', ''))
        proteins = _split_proteins_in_syst(row.get('protein_in_syst', ''))
        if not sys_type or not proteins: continue
        virus = _contig_from_protein_id(proteins[0])
        if not virus: continue
        records.append({'Virus': virus, 'System': sys_type,
                         'System_id': row.get('sys_id', sys_type),
                         'Genes': row.get('genes_count', ''), 'Source': 'DefenseFinder'})
    return [dict(r, sample=s) for s in samples for r in records]


def _load_apis_family_map(apis_db_dir):
    """seed_and_familyrep_all_infor.tsv (downloaded by rule dbapis_viral) --
    one row per APIS family, columns 'APIS families' (e.g. 'APIS001'),
    'APIS genes' (short characterized gene name, e.g. 'Apyc1'), 'Defense
    systems' (readable inhibited-defense-system label, e.g. 'pyrimidine
    cyclase system for antiphage resistance (Pycsar)') -- confirmed against
    a real download from pro.unl.edu/dbAPIS 2026-06-23. Indexed by BOTH the
    family ID and the gene name: dbAPIS protein headers (and therefore the
    'Family' field extracted from sseqid in load_dbapis_viral) use either
    convention depending on whether that family has a characterized name."""
    by_key = {}
    if not apis_db_dir:
        return by_key
    path = os.path.join(apis_db_dir, 'seed_and_familyrep_all_infor.tsv')
    if not os.path.exists(path):
        return by_key
    for row in load_tsv(path):
        fam    = row.get('APIS families', '').strip()
        gene   = row.get('APIS genes', '').strip()
        defsys = row.get('Defense systems', '').strip()
        if not defsys:
            continue
        if fam:  by_key[fam] = (gene, defsys)
        if gene: by_key[gene] = (gene, defsys)
    return by_key


def load_dbapis_viral(path, samples, apis_db_dir=''):
    """dbAPIS (Yan et al. 2023, NAR) DIAMOND blastp hits on viral proteins,
    from the global vOTU catalog's dbAPIS pass (rule votu_dbapis_viral,
    rules/votu_catalog.smk). Moved off per-sample paths on 2026-08-18, same
    move/rationale as load_antidefensefinder_viral above -- ONE table at
    {outdir}/votu_catalog/dbapis/dbapis_hits.tsv, the same records list
    duplicated under every sample key (see that function's docstring).
    Keeps only the best (lowest e-value) hit per query protein. sseqid is
    pipe-delimited: '{family_or_gene_id}|{IMGVR_UViG_id}|
    {genome_id}|{locus_with_coords}' (confirmed against a real run on
    litrp4, e.g. 'AcrIIA7|IMGVR_UViG_3300037418_004174|3300037418|
    Ga0395900_0000476_40112_40693') -- the first field alone is already a
    real, informative name (a dbAPIS family ID like 'APIS331', or a known
    gene name like 'AcrIIA7' for characterized Acr families). qseqid (and
    therefore 'Virus' via _contig_from_protein_id) is a NAMESPACED
    "{source_id}|{contig}" ID, since votu_prodigal runs on the namespaced
    catalog FASTA.

    'Gene'/'Defense_system_inhibited' add the readable translation via
    _load_apis_family_map (seed_and_familyrep_all_infor.tsv) -- e.g.
    'APIS331' -> gene name + 'restriction-modification system' instead of
    just the bare family ID. Falls back to the raw Family/empty string if
    apis_db_dir wasn't configured or the mapping file isn't there yet."""
    fam_map = _load_apis_family_map(apis_db_dir)
    best = {}
    for row in load_tsv(path):
        qseqid = row.get('qseqid', '')
        if not qseqid: continue
        try: evalue = float(row.get('evalue', '1') or '1')
        except ValueError: evalue = 1.0
        if qseqid not in best or evalue < best[qseqid][0]:
            best[qseqid] = (evalue, row)
    records = []
    for qseqid, (evalue, row) in best.items():
        virus = _contig_from_protein_id(qseqid)
        if not virus: continue
        sseqid = row.get('sseqid', '')
        family = sseqid.split('|', 1)[0] if sseqid else ''
        gene, defense_system = fam_map.get(family, ('', ''))
        records.append({'Virus': virus, 'Protein': qseqid,
                         'Family': family or sseqid, 'Gene': gene or family or sseqid,
                         'Defense_system_inhibited': defense_system,
                         'Hit': sseqid, 'Pident': row.get('pident', ''),
                         'Evalue': row.get('evalue', ''), 'Bitscore': row.get('bitscore', ''),
                         'Source': 'dbAPIS'})
    return [dict(r, sample=s) for s in samples for r in records]


# ── Defense islands (Han et al. 2026 / Beavogui et al. 2024 definition) ──────
# Arrays of defense genes separated by no more than 10 genes, containing
# >=5 genes from >=3 different defense systems.

def _read_protein_manifest(path):
    """Yield (name, mode, fna, faa, gff) tuples from a prok_bin_proteins
    manifest.txt — same format written by rules/defense_amr.smk's
    prok_bin_proteins rule, read independently here since the report side
    has no access to the Snakemake rules module's helper of the same name."""
    if not path or not os.path.exists(path):
        return
    with open(path) as f:
        for line in f:
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 5:
                continue
            yield parts[0], parts[1], parts[2], parts[3], parts[4]


def _ordered_proteins_by_contig(faa_path):
    """Protein IDs per contig, in genomic order — read straight from the
    .faa header order (Prodigal writes genes in scan order per sequence,
    numbering them '{seqid}_1', '{seqid}_2', ... — no GFF parsing needed)."""
    return {c: [g['Protein'] for g in genes]
            for c, genes in _genes_by_contig(faa_path).items()}


def _genes_by_contig(faa_path):
    """Genes per contig in genomic order, WITH coordinates and strand.

    Prodigal encodes the locus in the FASTA header itself:
        >{seqid}_{n} # {start} # {end} # {strand} # ID=...;partial=...
    so real bp coordinates and strand come for free from the .faa the
    pipeline already writes — no GFF parsing needed. Genes whose header
    lacks the coordinate fields (non-Prodigal input) still yield a record
    with Start/End/Strand set to None, and the report falls back to
    rendering by gene order."""
    by_contig = defaultdict(list)
    if not faa_path or not os.path.exists(faa_path):
        return by_contig
    try:
        with open(faa_path) as f:
            for line in f:
                if not line.startswith('>'):
                    continue
                header = line[1:].rstrip('\n')
                pid = header.split()[0].strip()
                contig = _contig_from_protein_id(pid)
                if not contig:
                    continue
                start = end = strand = None
                parts = [p.strip() for p in header.split('#')]
                # parts[0] is the id; 1..3 are start, end, strand (Prodigal)
                if len(parts) >= 4:
                    try:
                        start, end = int(parts[1]), int(parts[2])
                        strand = 1 if parts[3].lstrip('+') != '-1' else -1
                    except (ValueError, TypeError):
                        start = end = strand = None
                by_contig[contig].append(
                    {'Protein': pid, 'Start': start, 'End': end, 'Strand': strand})
    except Exception:
        pass
    return by_contig


def compute_defense_islands(manifest_paths, samples, defense_data, min_genes=5, min_systems=3, window=10):
    """One row per defense island: a run of defense genes on the same contig
    where consecutive defense-gene positions are never more than `window`
    genes apart, with >=min_genes genes from >=min_systems distinct systems.
    `defense_data` is load_defensefinder's output (needs the 'Proteins' field)."""
    islands = []
    for manifest_path, s in zip(manifest_paths, samples):
        if not manifest_path or not os.path.exists(manifest_path):
            continue
        # protein_id -> (System, System_id) for this sample's bins
        prot_to_sys = {}
        for rec in defense_data:
            if rec.get('sample') != s: continue
            for prot in rec.get('Proteins', []):
                prot_to_sys[prot] = (rec.get('Bin'), rec.get('System'), rec.get('System_id'))

        for name, mode, fna, faa, gff in _read_protein_manifest(manifest_path):
            by_contig = _genes_by_contig(faa)
            for contig, gene_recs in by_contig.items():
                ordered_prots = [g['Protein'] for g in gene_recs]
                hits = [(i, p, prot_to_sys[p]) for i, p in enumerate(ordered_prots) if p in prot_to_sys]
                if len(hits) < min_genes:
                    continue
                cluster = []
                def flush(cluster):
                    if len(cluster) < min_genes:
                        return
                    systems = {h[2][1] for h in cluster if h[2][1]}
                    if len(systems) < min_systems:
                        return
                    start_idx, end_idx = cluster[0][0], cluster[-1][0]
                    defense_idx = {h[0]: h[2][1] for h in cluster}
                    sysid_idx = {h[0]: h[2][2] for h in cluster}
                    window_genes = []
                    for i in range(start_idx, end_idx + 1):
                        rec = gene_recs[i]
                        window_genes.append({
                            'Protein': rec['Protein'], 'Index': i,
                            'System': defense_idx.get(i, ''),
                            'System_id': sysid_idx.get(i, ''),
                            'Start': rec.get('Start'), 'End': rec.get('End'),
                            'Strand': rec.get('Strand'),
                        })
                    coords = [g for g in window_genes
                              if g['Start'] is not None and g['End'] is not None]
                    islands.append({
                        'sample': s, 'Bin': name, 'Contig': contig,
                        'n_genes': len(cluster), 'n_systems': len(systems),
                        'Systems': sorted(systems),
                        'window_genes': window_genes,
                        'start_idx': start_idx, 'end_idx': end_idx,
                        # Genomic extent of the island (None when the input .faa
                        # carried no Prodigal coordinates — JS falls back to order).
                        'Start_bp': min(g['Start'] for g in coords) if coords else None,
                        'End_bp':   max(g['End']   for g in coords) if coords else None,
                    })
                for h in hits:
                    if cluster and h[0] - cluster[-1][0] > window:
                        flush(cluster)
                        cluster = []
                    cluster.append(h)
                flush(cluster)
    return islands


# ── AMR ──────────────────────────────────────────────────────────────────────

def _split_genome_prefix(value, sep="__"):
    """Split a '{genome}__{protein_id}' identifier produced by the
    defense_amr.smk protein-concatenation step into (genome, protein_id)."""
    if value and sep in value:
        genome, rest = value.split(sep, 1)
        return genome, rest
    return "", value or ""


_LOW_DEPTH_PSEUDO_GENOME = "contigs_pseudogenome"


def _prok_genome_unit(qseqid):
    """Resolve the taxonomic grouping unit for a prokaryote-side protein ID.

    Normal bins: the bin name (before '__'). low_depth_mode pools every
    contig's genes under one constant 'contigs_pseudogenome' prefix
    (rule prok_bin_proteins) -- grouping by that constant would blend every
    organism in the sample into a single majority-vote/LCA call. Regroup by
    the contig itself instead (recovered from the trailing
    '{contig}_{gene_n}' part via _contig_from_protein_id), so each contig
    gets its own taxonomy call, same granularity as the viral side."""
    prefix, rest = _split_genome_prefix(qseqid)
    if prefix == _LOW_DEPTH_PSEUDO_GENOME:
        return _contig_from_protein_id(rest)
    return prefix


def load_amr_consensus(paths, samples, *, low_depth_mode=False):
    """Consolidated AMR hits from consolidate_amr.py (AMRFinderPlus + RGI +
    DeepARG merged by CDS locus via ARO).  One row per locus; consensus_score
    = n_tools / 3.

    low_depth_mode: when True, contigs_pseudogenome is the analysis unit (no
    real bins), so those rows are kept (still filtered to n_tools >= 2).
    When False, contigs_pseudogenome is excluded entirely.
    """
    records = []
    for p, s in zip(paths, samples):
        for row in load_tsv(p):
            locus = row.get('locus', '')
            if not locus:
                continue
            genome, _ = _split_genome_prefix(locus)
            if genome == 'contigs_pseudogenome' and not low_depth_mode:
                continue
            if safe_int(row.get('n_tools', 0)) < 2:
                continue
            records.append({
                'sample':               s,
                'Bin':                  genome,
                'Gene':                 row.get('gene_name', ''),
                'ARO':                  row.get('aro_accession', ''),
                'Class':                row.get('drug_class', ''),
                'Resistance_mechanism': row.get('resistance_mechanism', ''),
                'n_tools':              safe_int(row.get('n_tools', 0)),
                'consensus_score':      safe_float(row.get('consensus_score', 0.0)),
                'tools_detected':       row.get('tools_detected', ''),
            })
    return records


# ── Host <-> Defense/AMR cross-link (Host & Defense report tab) ──────────────

def build_bin_annotation_summary(defense_data, antidefense_data, amr_consensus_data):
    """Defense/antidefense/AMR hits grouped once per (sample, bin) -- used as
    a small lookup table so the cross-link rows below don't have to embed
    (and duplicate) a host bin's full system/gene list in every row that
    references it."""
    def _add(target, key_field, records):
        for r in records:
            value = r.get(key_field)
            if not value:
                continue
            target.setdefault((r['sample'], r['Bin']), set()).add(value)

    defense_sets, antidefense_sets, amr_sets = {}, {}, {}
    _add(defense_sets,     'System', defense_data)
    _add(antidefense_sets, 'System', antidefense_data)
    _add(amr_sets,         'Gene',   amr_consensus_data)

    keys = set(defense_sets) | set(antidefense_sets) | set(amr_sets)
    summary = {}
    for sample, bin_name in keys:
        k = (sample, bin_name)
        summary[f"{sample}::{bin_name}"] = {
            'Defense_systems':    sorted(defense_sets.get(k, ())),
            'Antidefense_systems':sorted(antidefense_sets.get(k, ())),
            'AMR_genes':          sorted(amr_sets.get(k, ())),
        }
    return summary


def build_host_defense_links(phist_data, gtdb_data):
    """One row per predicted virus-host pair (PHIST), enriched with the
    host bin's GTDB-Tk taxonomy. Defense/antidefense/AMR detail is looked
    up client-side from BIN_ANNOTATIONS (build_bin_annotation_summary) by
    'sample::Host' instead of being embedded per row -- a host predicted
    for hundreds of viruses would otherwise duplicate its full system/gene
    list once per virus, which is what blew up report size 80MB -> 250MB."""
    gtdb_by_bin = {(r['sample'], r['Bin']): r for r in gtdb_data}

    links = []
    for row in phist_data:
        s, host = row['sample'], row.get('Host', '')
        if not host: continue
        tax = gtdb_by_bin.get((s, host), {})
        links.append({
            'sample': s, 'Virus': row.get('Virus', ''), 'Host': host,
            'Host_taxonomy': tax.get('Full_classification', ''),
            'Host_genus': tax.get('Genus', ''), 'Host_species': tax.get('Species', ''),
            'Score': row.get('Score', ''), 'P_value': row.get('P_value', ''),
        })
    return links


# ── Taxonomy enrichment / merge ───────────────────────────────────────────────

def enrich_taxonomy_with_checkv(tax_records, checkv_dict):
    cv_lookup = {}
    for s, rows in checkv_dict.items():
        cv_lookup[s] = {r.get('contig_id', r.get('contig', '')): r for r in rows}
    for rec in tax_records:
        cv_row = cv_lookup.get(rec.get('sample', ''), {}).get(rec.get('Genome', ''), {})
        rec['Completeness']   = cv_row.get('completeness', '')
        rec['Genome_length']  = cv_row.get('contig_length', '')
        rec['CheckV_quality'] = cv_row.get('checkv_quality', '')
    return tax_records


def collapse_taxonomy_to_votu(tax_records, outdir, samples):
    """One taxonomy row per vOTU representative, not per rep_seq contig.

    viral_nonredundant.fasta (MMseqs2, 95% identity, no aligned-fraction
    requirement) can keep multiple contigs that the stricter ICTV/Roux 2019
    vOTU definition (skani, 95% ANI + 85% AF — see rules/viral_binning.smk
    skani_cluster) considers the same viral population. Without collapsing,
    taxonomy tables/charts would count a single vOTU more than once.

    Groups records by (sample, representative) using vOTU_clusters.tsv, then
    picks the representative's own row if it was itself classified; otherwise
    falls back to the first classified member (the representative is chosen
    by CheckV completeness, not by confidence of taxonomic assignment, so it
    can occasionally be the one row in its cluster with no taxonomy hit).
    """
    membership = {}
    for s in samples:
        p = os.path.join(outdir, s, "viral", "votu", "vOTU_clusters.tsv")
        member_to_rep = {}
        for row in load_tsv(p):
            rep, member = row.get('representative', ''), row.get('member', '')
            if member:
                member_to_rep[member] = rep or member
        membership[s] = member_to_rep

    by_cluster = defaultdict(list)
    for rec in tax_records:
        s, genome = rec.get('sample', ''), rec.get('Genome', '')
        rep = membership.get(s, {}).get(genome, genome)
        by_cluster[(s, rep)].append(rec)

    collapsed = []
    for (s, rep), recs in by_cluster.items():
        chosen = next((r for r in recs if r.get('Genome') == rep), recs[0])
        chosen = dict(chosen)
        chosen['vOTU_members'] = len(recs)
        collapsed.append(chosen)
    return collapsed


def merge_prok_taxonomy(gtdb_records, mmseqs_prok_records, checkm2_dict, *, low_depth_mode=False):
    """Priority: GTDB-Tk (genome-level, most reliable) > MMseqs2 LCA (real
    per-genome lowest-common-ancestor against the custom IMG_NR DB -- see
    load_mmseqs_taxonomy_prok) > Unclassified. Replaces the old
    diamond_custom_prok best-hit/majority-vote source entirely (removed).

    Normally keyed off CheckM2 bins, since Completeness/Contamination only
    make sense for a real MAG. Under low_depth_mode there are no bins/CheckM2
    rows at all -- mmseqs records (already regrouped per-contig by
    _prok_genome_unit, not lumped into one sample-wide call) are surfaced
    directly instead of being silently dropped.
    """
    gtdb_bins   = {(r.get('sample', ''), r.get('Bin', '').replace('.fa', '')): r
                   for r in gtdb_records}
    mmseqs_bins = {(r.get('sample', ''), r.get('Bin', '')): r for r in mmseqs_prok_records}

    merged = []
    seen = set()
    for sample, rows in checkm2_dict.items():
        for cm_row in rows:
            bin_name = cm_row.get('Name', cm_row.get('name', '')).replace('.fa', '')
            if not bin_name: continue
            key = (sample, bin_name)
            seen.add(key)
            if key in gtdb_bins:
                base = dict(gtdb_bins[key])
                base['Source_tax'] = 'GTDB-Tk'
            elif key in mmseqs_bins:
                base = dict(mmseqs_bins[key])
                base['Bin'] = bin_name; base['sample'] = sample
                base['Source_tax'] = 'MMseqs2-LCA'
            else:
                base = {'sample': sample, 'Bin': bin_name,
                        'Domain': '', 'Phylum': '', 'Class': '', 'Order': '',
                        'Family': '', 'Genus': '', 'Species': '',
                        'Source_tax': 'Unclassified'}
            base['Completeness']  = cm_row.get('Completeness', '')
            base['Contamination'] = cm_row.get('Contamination', '')
            base['Genome_size']   = cm_row.get('Genome_Size', cm_row.get('genome_size', ''))
            merged.append(base)

    # low_depth_mode: no bins exist; aggregate per-contig MMseqs2-LCA records
    # by taxonomy string so the report stays manageable (one row per unique
    # lineage per sample instead of one row per contig).
    # Normal mode: skip this block to avoid flooding MERGED_PROK with raw contigs.
    if low_depth_mode and not merged:
        from collections import defaultdict
        groups = defaultdict(lambda: {'count': 0, 'rec': None})
        for key, rec in mmseqs_bins.items():
            if key in seen:
                continue
            sample = key[0]
            lineage = ';'.join(filter(None, [
                rec.get('Domain', ''), rec.get('Phylum', ''), rec.get('Class', ''),
                rec.get('Order', ''),  rec.get('Family', ''), rec.get('Genus', ''),
                rec.get('Species', ''),
            ])) or 'Unclassified'
            gk = (sample, lineage)
            groups[gk]['count'] += 1
            groups[gk]['rec'] = rec
        for (sample, lineage), val in groups.items():
            rec = val['rec']
            n   = val['count']
            merged.append({
                'sample':      sample,
                'Bin':         f'{n} contigs',
                'Domain':      rec.get('Domain', ''),
                'Phylum':      rec.get('Phylum', ''),
                'Class':       rec.get('Class', ''),
                'Order':       rec.get('Order', ''),
                'Family':      rec.get('Family', ''),
                'Genus':       rec.get('Genus', ''),
                'Species':     rec.get('Species', ''),
                'Full_classification': lineage,
                'Source_tax':  'MMseqs2-LCA',
                'Completeness': '', 'Contamination': '', 'Genome_size': '',
                'contig_count': n,
            })

    return merged


# ── Reads-only classification (Sylph + sylph-tax) ────────────────────────────

def load_reads_classify(abundance_path, host_path, samples):
    """Parse sylph-tax merged abundance and viral-by-host tables for the report.

    sylphmpa format uses | as rank separator and r__ (realm) for viral taxonomy.
    Merged TSV columns are full fastq paths — mapped to sample names by basename.

    Returns dict: viral, prok, archaea, host, has_data, samples
    """
    def _fastq_to_sample(col):
        """Strip directory + fastq extensions to recover sample name."""
        name = os.path.basename(col)
        for ext in ('.fastq.gz', '.fq.gz', '.fastq', '.fq'):
            if name.endswith(ext):
                name = name[:-len(ext)]
                break
        return name

    def _is_viral(clade):
        # ICTV realm taxonomy: r__Duplodnaviria, r__Monodnaviria, etc.
        first = clade.split('|')[0]
        if first.startswith('r__') and first not in ('r__', ):
            return True
        # UNKNOWN|...|t__IMGVR_* — IMG/VR hits with no realm annotation
        for part in clade.split('|'):
            if part.startswith('t__IMGVR_') or part.startswith('t__UHGV_'):
                return True
            if part.startswith('d__') and 'virus' in part.lower():
                return True
        return False

    def _is_archaea(clade):
        for part in clade.split('|'):
            if part == 'd__Archaea':
                return True
        return False

    # Standard rank prefixes for sylph-tax / sylphmpa format (ordered shallowest→deepest)
    _STD_RANKS = ('r__', 'k__', 'p__', 'c__', 'o__', 'f__', 'g__', 's__')

    def _last_std_rank_seg(clade):
        """Return the deepest standard-rank segment in the clade, or ''."""
        last = ''
        for seg in clade.split('|'):
            if any(seg.startswith(p) for p in _STD_RANKS):
                last = seg
        return last

    def _parent_has_same_eff(clade, cache, my_eff):
        """True if a shallower row shares the same deepest std-rank value.
        Identifies UNKNOWN-chain duplicates that repeat the parent's aggregate."""
        parts = clade.split('|')
        for i in range(1, len(parts)):
            parent = '|'.join(parts[:i])
            if cache.get(parent) == my_eff:
                return True
        return False

    viral, prok, archaea = [], [], []
    col_to_sample = {}  # merged TSV column header → sample name

    if abundance_path and os.path.exists(abundance_path):
        rows = list(load_tsv(abundance_path))
        if rows:
            taxon_col = list(rows[0].keys())[0]

            # Build col→sample map (merged headers are full fastq paths)
            for col in rows[0].keys():
                if col == taxon_col:
                    continue
                mapped = _fastq_to_sample(col)
                if mapped in samples:
                    col_to_sample[col] = mapped
                elif col in samples:
                    col_to_sample[col] = col

            all_clades = {row.get(taxon_col, '') for row in rows}
            eff_cache = {c: _last_std_rank_seg(c) for c in all_clades}

            for row in rows:
                clade = row.get(taxon_col, '').strip()
                if not clade:
                    continue

                # Keep "effective-rank rows": rows whose last std-rank segment is
                # their canonical rank label. Skip:
                # • all-UNKNOWN/t__ chains with no std rank → uninformative
                # • UNKNOWN-chain duplicates: a shallower row has the same eff rank
                #   (e.g., c__|UNKNOWN repeats c__ aggregate abundance)
                # This keeps one row per distinct recognized-rank node.
                # The JS then filters by _eff_rank to select the right level for
                # each chart without double-counting across ranks.
                my_eff = eff_cache.get(clade, '')
                if not my_eff:
                    continue
                if _parent_has_same_eff(clade, eff_cache, my_eff):
                    continue

                record = {'clade': clade, '_eff_rank': my_eff[:3]}  # e.g. 'c__'
                for col, sname in col_to_sample.items():
                    record[sname] = safe_float(row.get(col, 0.0))
                for s in samples:
                    if s not in record:
                        record[s] = 0.0

                if _is_viral(clade):
                    viral.append(record)
                elif _is_archaea(clade):
                    archaea.append(record)
                else:
                    prok.append(record)

    host = []
    if host_path and os.path.exists(host_path):
        for row in load_tsv(host_path):
            h = row.get('host_genus', '')
            if not h:
                continue
            record = {'host_genus': h, 'n_viral_taxa': safe_int(row.get('n_viral_taxa', 0))}
            for col, sname in col_to_sample.items():
                record[sname] = safe_float(row.get(col, 0.0))
            for s in samples:
                if s not in record:
                    record[s] = safe_float(row.get(s, 0.0))
            host.append(record)
        host.sort(key=lambda r: sum(r.get(s, 0) for s in samples), reverse=True)

    return {
        'viral':    viral,
        'prok':     prok,
        'archaea':  archaea,
        'host':     host,
        'samples':  samples,
        'has_data': bool(viral or prok or archaea),
    }


# ── New loaders: diversity + functional ──────────────────────────────────────

_ALPHA_DOMAIN_MAP = {'prokaryotic': 'prok', 'prok': 'prok', 'viral': 'viral', 'combined': 'combined'}
_ALPHA_INDEX_MAP  = {'richness': 'observed', 'shannon': 'shannon', 'simpson': 'simpson', 'chao1': 'chao1'}


def load_alpha_diversity(path):
    """Load alpha_diversity.tsv (wide: sample, domain, richness, shannon, simpson, chao1)
    and reshape into long format → list of {sample, domain, index, value}."""
    rows = []
    for row in load_tsv(path):
        sample = row.get('sample', '')
        domain = _ALPHA_DOMAIN_MAP.get(row.get('domain', ''), row.get('domain', 'combined'))
        if not sample:
            continue
        for col, index in _ALPHA_INDEX_MAP.items():
            if col in row:
                rows.append({'sample': sample, 'domain': domain, 'index': index,
                             'value': safe_float(row.get(col, 0))})
    return rows


def load_pcoord(path):
    """Load beta_pcoord_*.tsv → list of {sample, pc1, pc2, var_pc1, var_pc2}.
    PC1_var/PC2_var are written as percentages; stored here as 0-1 fractions
    (the JS multiplies by 100 when displaying)."""
    rows = []
    for row in load_tsv(path):
        sample  = row.get('sample', '')
        pc1     = safe_float(row.get('PC1', row.get('pc1', 0)))
        pc2     = safe_float(row.get('PC2', row.get('pc2', 0)))
        var_pc1 = safe_float(row.get('PC1_var', row.get('var_pc1', 0))) / 100.0
        var_pc2 = safe_float(row.get('PC2_var', row.get('var_pc2', 0))) / 100.0
        if sample:
            rows.append({'sample': sample, 'pc1': pc1, 'pc2': pc2,
                         'var_pc1': var_pc1, 'var_pc2': var_pc2})
    return rows


_COG_LABEL = {
    'J': 'Translation', 'K': 'Transcription', 'L': 'Replication/Repair',
    'D': 'Cell cycle', 'M': 'Cell membrane', 'C': 'Energy production',
    'E': 'Amino acid met.', 'G': 'Carbohydrate met.', 'P': 'Ion transport',
    'T': 'Signal transd.', 'V': 'Defense', 'O': 'Protein modif.',
    'U': 'Secretion', 'N': 'Cell motility', 'S': 'Unknown function',
}


def load_eggnog(outdir, samples):
    """Aggregate COG categories from eggnog_annotations.tsv per sample."""
    result = {}
    for s in samples:
        p = os.path.join(outdir, s, "annotation", "eggnog", "eggnog_annotations.tsv")
        cnt = Counter()
        for row in load_tsv(p):
            cog = str(row.get('COG_category', '') or row.get('cog_category', '') or 'S').strip()
            for c in cog:
                if c in _COG_LABEL:
                    cnt[_COG_LABEL[c]] += 1
        result[s] = dict(cnt)
    return result



def load_phrogs(outdir, samples):
    """PHROGS category counts from the global vOTU catalog's Pharokka run.

    Moved off per-sample paths on 2026-08-18 (second half of "(h)",
    docs/ROADMAP_SIMPLIFICACAO.md): pharokka now runs once over the whole
    vOTU catalog (rule votu_pharokka, rules/votu_catalog.smk), not once per
    sample. There is exactly ONE table to read, at
    {outdir}/votu_catalog/annotation/pharokka/pharokka_cds_final_merged_output.tsv
    -- the same dict is returned under every sample key so the existing
    per-sample chart plumbing in annotation.js keeps working unmodified;
    the chart title there was updated to say "vOTU catalog" instead of
    implying a per-sample count.
    """
    p = os.path.join(outdir, "votu_catalog", "annotation", "pharokka",
                     "pharokka_cds_final_merged_output.tsv")
    cnt = Counter()
    for row in load_tsv(p):
        cat = (row.get('phrog_category', '') or row.get('category', '') or
               'unknown function').lower().strip()
        if not cat: cat = 'unknown function'
        cnt[cat] += 1
    catalog_counts = dict(cnt)
    return {s: catalog_counts for s in samples}


def load_svg(svg_path):
    if not os.path.exists(svg_path): return ""
    try:
        with open(svg_path, encoding='utf-8') as f: return f.read()
    except Exception: return ""


def load_genome_maps(outdir, samples):
    """Load genome map SVGs for virus/prok modes, max 5 per category, from
    the global vOTU catalog.

    "virus" merges the phage/ and virus/ output subfolders into a single list
    -- the report shows one unified "Virus" view, with each genome tagged
    category="Phage" or category="Virus" (the backend split still exists on
    disk, decided by PHROGS hallmark-gene evidence in genome_map_universal.py,
    but the UI no longer forces the user to pick a mode to see all of them).

    Phage/virus maps moved to the global catalog on 2026-08-18 (second half
    of "(h)"): genome_map_phage/genome_map_virus now run once, over vOTU
    representatives (rules votu_genome_map_phage/votu_genome_map_virus,
    rules/votu_catalog.smk), at
    {outdir}/votu_catalog/annotation/genome_maps/{{phage,virus}}. The SAME
    catalog-wide "virus" list is returned under every sample key, same
    rationale as load_phrogs above. "prok" stays genuinely per-sample
    (genome_map_prok, rules/annotation.smk) -- prokaryotic MAGs are not
    deduplicated into a global catalog the way vOTUs are.
    """
    virus_base = os.path.join(outdir, "votu_catalog", "annotation", "genome_maps")
    virus_maps = []
    for mode, category in (("phage", "Phage"), ("virus", "Virus")):
        mdir = os.path.join(virus_base, mode)
        for svg_f in sorted(glob.glob(os.path.join(mdir, "*.svg")))[:5]:
            gid = os.path.basename(svg_f).replace("_map.svg", "")
            svg = load_svg(svg_f)
            if svg:
                seq = ""
                fasta_f = os.path.join(mdir, f"{gid}.fasta")
                if os.path.exists(fasta_f):
                    try:
                        with open(fasta_f) as ff:
                            seq = "".join(l.strip() for l in ff if not l.startswith(">"))
                    except Exception:
                        pass
                virus_maps.append({"id": gid, "svg": svg, "seq": seq, "category": category})

    result = {}
    for s in samples:
        result[s] = {"virus": virus_maps, "prok": []}
        base = os.path.join(outdir, s, "annotation", "genome_maps")
        mdir = os.path.join(base, "prok")
        for svg_f in sorted(glob.glob(os.path.join(mdir, "*.svg")))[:5]:
            gid = os.path.basename(svg_f).replace("_map.svg", "")
            svg = load_svg(svg_f)
            if svg:
                seq = ""
                fasta_f = os.path.join(mdir, f"{gid}.fasta")
                if os.path.exists(fasta_f):
                    try:
                        with open(fasta_f) as ff:
                            seq = "".join(l.strip() for l in ff if not l.startswith(">"))
                    except Exception:
                        pass
                result[s]["prok"].append({"id": gid, "svg": svg, "seq": seq})
    return result


# ── Path dict helper ──────────────────────────────────────────────────────────

def path_dict(paths, samples):
    """Map paths to sample names by matching sample name in path.

    Always returns a dict with ALL samples as keys; missing samples get empty string.
    """
    d = {s: '' for s in samples}
    for p in paths:
        for s in samples:
            if f"/{s}/" in p or p.endswith(f"/{s}"):
                d[s] = p; break
    # fallback: zip order if pattern matching failed to fill all samples
    if not all(d[s] for s in samples):
        for s, p in zip(samples, paths):
            if not d.get(s):
                d[s] = p
    return d


# ── Tool versions ─────────────────────────────────────────────────────────────

def collect_tool_versions():
    import subprocess as _sp
    import re as _re

    def _ver(cmd, line_idx=0):
        try:
            r = _sp.run(cmd, shell=True, capture_output=True, text=True, timeout=20)
            out = (r.stdout + r.stderr).strip()
            lines = [l for l in out.splitlines() if l.strip()]
            if not lines: return "n/a"
            target = lines[line_idx] if line_idx < len(lines) else lines[0]
            m = _re.search(r"[\d]+\.[\d]+[\d.]*", target)
            return m.group(0) if m else target[:40]
        except Exception: return "n/a"

    C = "conda run --no-capture-output -n"
    return {
        "fastp":         _ver(f"{C} env_qc fastp --version"),
        "MultiQC":       _ver(f"{C} env_qc multiqc --version"),
        "MEGAHIT":       _ver(f"{C} env_assembly megahit --version"),
        "MMseqs2":       _ver(f"{C} env_assembly mmseqs version"),
        "QUAST":         _ver(f"{C} env_qc quast.py --version"),
        "BWA-MEM2":      _ver(f"{C} env_mapping bwa-mem2 version"),
        "minimap2":      _ver(f"{C} env_mapping minimap2 --version"),
        "samtools":      _ver(f"{C} env_mapping samtools --version"),
        "CoverM":        _ver(f"{C} env_coverm coverm --version"),
        "VirSorter2":    _ver(f"{C} env_viral virsorter --version"),
        "GeNomad":       _ver(f"{C} env_genomad genomad --version"),
        "CheckV":        _ver(f"{C} env_viral checkv 2>&1 | head -1"),
        "vRhyme":        _ver(f"{C} env_vrhyme vRhyme --version"),
        "Diamond":       _ver(f"{C} env_viral diamond version"),
        "Prodigal":      _ver(f"{C} env_viral prodigal -v 2>&1 | head -2"),
        "MetaBAT2":      _ver(f"{C} env_binning metabat2 2>&1 | head -3"),
        "MaxBin2":       _ver(f"{C} env_binning run_MaxBin.pl 2>&1 | head -2"),
        "VAMB":          _ver(f"{C} env_binning vamb --version"),
        "SemiBin2":      _ver(f"{C} env_binning SemiBin2 --version"),
        "Binette":       _ver(f"{C} env_binette binette --version"),
        "CheckM2":       _ver(f"{C} env_checkm2 checkm2 --version"),
        "GTDB-Tk":       _ver(f"{C} env_gtdbtk gtdbtk --version"),
        "PHIST":         _ver(f"{C} env_phist phist 2>&1 | head -1"),
        "Pharokka":      _ver(f"{C} env_annotation pharokka --version"),
        "EggNOG-mapper": _ver(f"{C} env_annotation emapper.py --version"),
    }


# ── Co-assembly track (Plano 2) ────────────────────────────────────────────

def load_coassembly(outdir, groups):
    """Collect group-level CheckM2 + GTDB MAGs and CheckV/taxonomy vOTUs into
    report records (empty if none)."""
    out = []
    for g in groups:
        checkm2 = os.path.join(outdir, "coassembly", g, "checkm2", "quality_report.tsv")
        mags = []
        if os.path.exists(checkm2):
            for row in load_tsv(checkm2):
                mags.append({
                    "bin": row.get("Name", ""),
                    "completeness": safe_float(row.get("Completeness", 0)),
                    "contamination": safe_float(row.get("Contamination", 0)),
                    "classification": "",
                })
        # join GTDB classification if present
        for summ in glob.glob(os.path.join(outdir, "coassembly", g, "gtdbtk", "**", "*summary.tsv"), recursive=True):
            for row in load_tsv(summ):
                name = row.get("user_genome", "")
                for m in mags:
                    if m["bin"] == name:
                        m["classification"] = row.get("classification", "")

        # Group vOTUs (Plano 3): count vOTU representatives + quality tiers + taxonomy family counts
        n_votus = 0
        votu_reps = os.path.join(outdir, "coassembly", g, "viral", "votu", "votu_all_reps.fasta")
        if os.path.exists(votu_reps):
            with open(votu_reps) as fh:
                n_votus = sum(1 for line in fh if line.startswith(">"))

        quality_tiers = {}
        checkv_summary = os.path.join(outdir, "coassembly", g, "viral", "checkv", "quality_summary.tsv")
        if os.path.exists(checkv_summary):
            for row in load_tsv(checkv_summary):
                q = row.get("checkv_quality", "") or "Unknown"
                quality_tiers[q] = quality_tiers.get(q, 0) + 1

        family_counts = {}
        tax_merged = os.path.join(outdir, "coassembly", g, "viral", "taxonomy", "viral_taxonomy_merged.tsv")
        if os.path.exists(tax_merged):
            for row in load_tsv(tax_merged):
                fam = (row.get("final_family", "") or "").strip()
                if not fam:
                    fam = "Unclassified"
                family_counts[fam] = family_counts.get(fam, 0) + 1

        votu_families = sorted(
            [{"family": f, "count": c} for f, c in family_counts.items()],
            key=lambda x: -x["count"],
        )[:5]

        # Group vRhyme vMAGs (Plano 4): count CheckV-assessed vMAGs (short reads only)
        n_vmags = 0
        vrhyme_summary = os.path.join(outdir, "coassembly", g, "viral", "checkv_vrhyme", "quality_summary.tsv")
        if os.path.exists(vrhyme_summary):
            for row in load_tsv(vrhyme_summary):
                if (row.get("contig_id", "") or "").strip():
                    n_vmags += 1

        if mags or n_votus or n_vmags:
            out.append({
                "group": g,
                "mags": mags,
                "n_votus": n_votus,
                "n_vmags": n_vmags,
                "quality_tiers": quality_tiers,
                "votu_families": votu_families,
            })
    return {"groups": out, "has_data": bool(out)}


def load_votu_accumulation(outdir, groups, min_depth=1.0, n_perm=100, seed=0):
    """vOTU accumulation (collector) curve per co-assembly group.

    Computed per co-assembly group. Global accumulation over all samples is
    available from the vOTU catalog (presence_matrix.tsv), where vOTU identity
    is shared across every sample by construction.

    Presence: a vOTU counts as seen in a sample when any of its member contigs
    reaches `min_depth` mean coverage there, read from the group's VAMB
    abundance matrix (contig x sample totalAvgDepth) that co-binning already
    builds. Sample order is arbitrary, so the curve is averaged over `n_perm`
    random orders and reported with a 10-90 percentile band rather than one
    misleading ordering."""
    rng = random.Random(seed)
    out = {}
    for g in groups or []:
        matrix = os.path.join(outdir, "coassembly", g, "vamb", "abundance.tsv")
        clusters = os.path.join(outdir, "coassembly", g, "viral", "votu", "vOTU_clusters.tsv")
        if not (os.path.exists(matrix) and os.path.exists(clusters)):
            continue

        member_to_rep = {}
        for row in load_tsv(clusters):
            m, r = row.get("member", ""), row.get("representative", "")
            if m and r:
                member_to_rep[m] = r
        if not member_to_rep:
            continue

        # rep -> set of samples where any member reaches min_depth
        present = defaultdict(set)
        samples = []
        try:
            with open(matrix) as fh:
                header = fh.readline().rstrip("\n").split("\t")
                samples = header[1:]
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    rep = member_to_rep.get(parts[0])
                    if not rep:
                        continue
                    for i, s in enumerate(samples, start=1):
                        if i < len(parts) and safe_float(parts[i]) >= min_depth:
                            present[rep].add(s)
        except Exception:
            continue
        if not present or len(samples) < 2:
            continue

        by_sample = {s: {rep for rep, ss in present.items() if s in ss} for s in samples}
        n = len(samples)
        curves = []
        for _ in range(n_perm):
            order = samples[:]
            rng.shuffle(order)
            seen, row = set(), []
            for s in order:
                seen |= by_sample[s]
                row.append(len(seen))
            curves.append(row)

        def pct(vals, p):
            v = sorted(vals)
            return v[min(len(v) - 1, max(0, int(round((len(v) - 1) * p))))]

        cols = list(zip(*curves))
        out[g] = {
            "x":    list(range(1, n + 1)),
            "mean": [round(sum(c) / len(c), 1) for c in cols],
            "lo":   [pct(c, 0.10) for c in cols],
            "hi":   [pct(c, 0.90) for c in cols],
            "total": len(present),
            "n_samples": n,
            "min_depth": min_depth,
        }
    return out


# ── Global vOTU catalog ───────────────────────────────────────────────────────
#
# Richness is defined ONCE over the pooled viral sets of every sample and
# co-assembly group, so a per-sample count and a total are on the same scale.
# Summing per-sample counts is what inflated richness before this stage
# existed, and is never the right total.

def load_votu_catalog(outdir):
    """Global richness and how much the pool collapsed.

    Returns {'n_votus', 'n_pool', 'reduction_pct'}; zeros when absent.
    """
    d = os.path.join(outdir, "votu_catalog")
    clusters = os.path.join(d, "vOTU_clusters.tsv")
    provenance = os.path.join(d, "provenance.tsv")
    n_votus = n_pool = 0

    if os.path.exists(clusters):
        seen = set()
        with open(clusters) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                seen.add(row["votu_id"])
        n_votus = len(seen)

    if os.path.exists(provenance):
        with open(provenance) as fh:
            n_pool = max(sum(1 for _ in fh) - 1, 0)

    reduction = (100.0 * (1 - n_votus / n_pool)) if n_pool else 0.0
    return {"n_votus": n_votus, "n_pool": n_pool,
            "reduction_pct": round(reduction, 1)}


def load_votu_presence(outdir, samples):
    """Per-sample vOTU presence, keeping the two signals separate.

    'assembled' — a member contig came from that sample
    'recruited' — the sample's reads covered the representative past the cutoff
    'total'     — present by either signal; this is the sample's vOTU count
    """
    path = os.path.join(outdir, "votu_catalog", "presence_matrix.tsv")
    per_sample = {s: {"assembled": 0, "recruited": 0, "total": 0} for s in samples}
    votus = []

    if not os.path.exists(path):
        return {"votus": votus, "per_sample": per_sample}

    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            entry = {"votu_id": row["votu_id"],
                     "representative": row.get("representative", ""),
                     "samples": {}}
            for s in samples:
                state = (row.get(s) or "absent").strip()
                entry["samples"][s] = state
                if state in ("assembled", "both"):
                    per_sample[s]["assembled"] += 1
                if state in ("recruited", "both"):
                    per_sample[s]["recruited"] += 1
                if state != "absent":
                    per_sample[s]["total"] += 1
            votus.append(entry)

    return {"votus": votus, "per_sample": per_sample}


def load_coassembly_rich(outdir, groups):
    """Per-group rich data mirroring the per-sample viral + prokaryotic tabs,
    for the Co-assembly report tab. Reuses the same file-path-based loaders as
    the per-sample side (load_viral_taxonomy / load_gtdbtk / merge_prok_taxonomy
    / parse_tsv), passing GROUP names where those loaders expect sample names —
    every record it emits is keyed by group under the 'sample' field, so the JS
    render functions can filter by the selected group exactly like per-sample.

    Report-only: consumes what the co-assembly rules already produce. There is
    no group `make_votu_table`/VIBRANT aggregation, so vOTU lifestyle/abundance
    are intentionally absent (see the compact `load_coassembly` for counts)."""
    groups = list(groups or [])
    if not groups:
        return {"units": [], "checkv": {}, "checkv_vrh": {}, "checkm2": {},
                "tax": [], "source_dist": {}, "merged_prok": [], "mimag": {},
                "vlen": {}}

    def _g(g, *parts):
        return os.path.join(outdir, "coassembly", g, *parts)

    checkv     = {g: parse_tsv(_g(g, "viral", "checkv", "quality_summary.tsv")) for g in groups}
    checkv_vrh = {g: parse_tsv(_g(g, "viral", "checkv_vrhyme", "quality_summary.tsv")) for g in groups}
    checkm2    = {g: parse_tsv(_g(g, "checkm2", "quality_report.tsv")) for g in groups}

    tax_paths = [_g(g, "viral", "taxonomy", "viral_taxonomy_merged.tsv") for g in groups]
    tax = load_viral_taxonomy(tax_paths, groups)
    tax = enrich_taxonomy_with_checkv(tax, checkv)
    source_dist = load_viral_source_distribution(tax_paths, groups)

    bac_paths = [_g(g, "gtdbtk", "classify", "gtdbtk.bac120.summary.tsv") for g in groups]
    arc_paths = [_g(g, "gtdbtk", "classify", "gtdbtk.ar53.summary.tsv") for g in groups]
    gtdb = load_gtdbtk(bac_paths, arc_paths, groups)
    merged_prok = merge_prok_taxonomy(gtdb, [], checkm2, low_depth_mode=False)

    mimag = {}
    for g in groups:
        hq = mq = lq = 0
        for row in checkm2[g]:
            comp = safe_float(row.get('Completeness', 0))
            cont = safe_float(row.get('Contamination', 100))
            if comp >= 90 and cont <= 5:    hq += 1
            elif comp >= 50 and cont <= 10: mq += 1
            else:                           lq += 1
        mimag[g] = {'HQ': hq, 'MQ': mq, 'LQ': lq, 'total': hq + mq + lq}

    vlen = {}
    for g in groups:
        fa = _g(g, "viral", "votu", "votu_all_reps.fasta")
        vlen[g] = parse_fasta_lengths(fa) if os.path.exists(fa) else []

    tax_groups  = {r.get('sample', '') for r in tax}
    gtdb_groups = {r.get('sample', '') for r in gtdb}
    units = [g for g in groups
             if checkv[g] or checkm2[g] or vlen[g]
             or g in tax_groups or g in gtdb_groups]
    return {
        "units": units,
        "checkv": checkv, "checkv_vrh": checkv_vrh, "checkm2": checkm2,
        "tax": tax, "source_dist": source_dist,
        "merged_prok": merged_prok, "mimag": mimag, "vlen": vlen,
    }


def load_votu_lifestyle(outdir):
    """Global BACPHLIP lifestyle call per vOTU representative (bacphlip_votu).

    Reads votu_catalog/bacphlip/votu_lifestyle.tsv. Returns
    {'rows', 'counts', 'quality_tier_used'}; degrades silently to empty/zero
    when the file is absent or header-only -- a skipped/failed run must read
    as a gap in the report, not as a biological zero.
    """
    path = os.path.join(outdir, "votu_catalog", "bacphlip", "votu_lifestyle.tsv")
    rows = []
    counts = {"lytic": 0, "lysogenic": 0}
    tier_used = ""
    for row in load_tsv(path):
        lifestyle = (row.get("lifestyle", "") or "").strip()
        rows.append({
            "votu_id": row.get("votu_id", ""),
            "lifestyle": lifestyle,
            "virulent_score": row.get("virulent_score", ""),
            "checkv_quality": row.get("checkv_quality", ""),
        })
        if lifestyle in counts:
            counts[lifestyle] += 1
        if not tier_used:
            tier_used = (row.get("quality_tier_used", "") or "").strip()
    return {"rows": rows, "counts": counts, "quality_tier_used": tier_used}


def load_putative_amgs(outdir):
    """Global putative AMG candidates from eggNOG on vOTU representatives
    (eggnog_viral). Reads votu_catalog/eggnog_viral/putative_amgs.tsv.

    Named "putative" throughout -- AMG calls from annotation alone are
    recognized in the literature as prone to mis-annotation and require
    genomic-context inspection the pipeline does not perform. Returns
    {'rows', 'n_total', 'n_votus'}; degrades silently when absent.
    """
    path = os.path.join(outdir, "votu_catalog", "eggnog_viral", "putative_amgs.tsv")
    rows = []
    votus = set()
    for row in load_tsv(path):
        votu_id = row.get("votu_id", "")
        rows.append({
            "votu_id": votu_id,
            "protein_id": row.get("protein_id", ""),
            "KEGG_ko": row.get("KEGG_ko", ""),
            "KEGG_Pathway": row.get("KEGG_Pathway", ""),
            "COG_category": row.get("COG_category", ""),
            "Description": row.get("Description", ""),
        })
        if votu_id:
            votus.add(votu_id)
    return {"rows": rows, "n_total": len(rows), "n_votus": len(votus)}
