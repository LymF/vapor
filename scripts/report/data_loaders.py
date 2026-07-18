"""data_loaders.py — All data-loading functions for the VAPOR HTML report.

Ported verbatim from generate_report.py with new loaders added for
alpha diversity, PCoA, KEGG, and eggnog aggregation.
"""
import os, re, glob, csv, json, zipfile
from collections import Counter, defaultdict


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


# ── VIBRANT ───────────────────────────────────────────────────────────────────

def load_vibrant(outdir_base, samples):
    """Parse VIBRANT summary + AMG results for each sample."""
    scaffold_records = []
    amg_records = []
    for s in samples:
        vdir = os.path.join(outdir_base, s, 'viral', 'vibrant')
        for tsv in glob.glob(os.path.join(vdir, '**', 'VIBRANT_summary_results_*.tsv'), recursive=True):
            for row in load_tsv(tsv):
                scaffold = row.get('scaffold', '')
                if not scaffold: continue
                scaffold_records.append({'sample': s,
                    'Scaffold':    scaffold,
                    'Total_genes': row.get('total genes', ''),
                    'VOG_score':   row.get('VOG v-score', ''),
                    'Pfam_score':  row.get('Pfam v-score', ''),
                    'KEGG_score':  row.get('KEGG v-score', ''),
                })
        for tsv in glob.glob(os.path.join(vdir, '**', 'VIBRANT_AMG_pathways_*.tsv'), recursive=True):
            for row in load_tsv(tsv):
                pathway = row.get('Pathway', '') or row.get('pathway', '')
                if not pathway: continue
                n_amgs = row.get('Total AMGs', '') or row.get('total amgs', '1')
                kos    = row.get('Present AMG KOs', '') or row.get('present amg kos', '')
                amg_records.append({'sample':    s,
                    'Pathway':    pathway,
                    'Metabolism': row.get('Metabolism', '') or row.get('metabolism', ''),
                    'KEGG_map':   row.get('KEGG Entry', '') or row.get('kegg entry', ''),
                    'Total_AMGs': n_amgs,
                    'KOs':        kos,
                })
        for tsv in glob.glob(os.path.join(vdir, '**', 'VIBRANT_AMG_individuals_*.tsv'), recursive=True):
            for row in load_tsv(tsv):
                prot = row.get('protein', '') or row.get('gene', '')
                if not prot: continue
                amg_records.append({'sample':    s,
                    'Protein':    prot,
                    'Scaffold':   row.get('scaffold', ''),
                    'AMG_KO':     row.get('AMG KO', '') or row.get('KO', ''),
                    'AMG_gene':   row.get('AMG gene', '') or row.get('gene name', ''),
                    'Pathway':    row.get('metabolic pathway', '') or row.get('pathway', ''),
                    'Metabolism': row.get('AMG category', '') or row.get('category', ''),
                    'Total_AMGs': '1',
                    'KOs':        row.get('AMG KO', ''),
                })
    return scaffold_records, amg_records


# ── vConTACT3 ─────────────────────────────────────────────────────────────────

def load_vcontact3(paths_d, samples):
    """Load vConTACT3 final_assignments.csv (v3 output format)."""
    _VC3_NULL = {"singleton", "unclassified", "nd", "none", ""}
    records = []
    for s in samples:
        p = paths_d.get(s, '')
        rows = load_csv(p) if p.endswith('.csv') else load_tsv(p)
        for row in rows:
            genome = row.get('Genome', row.get('genome', ''))
            if not genome: continue
            ref_val = row.get('Reference', row.get('reference', ''))
            if str(ref_val).lower() in ('true', '1', 'yes'): continue
            fam_raw = row.get('family_prediction', row.get('Family', ''))
            gen_raw = row.get('genus_prediction',  row.get('Genus', ''))
            rlm     = row.get('realm_prediction',  row.get('realm_reference', ''))
            cls     = row.get('class_prediction',  '')
            ord_raw = row.get('order_prediction', '')
            fam  = fam_raw  if fam_raw  and fam_raw.lower()  not in _VC3_NULL else ''
            gen  = gen_raw  if gen_raw  and gen_raw.lower()  not in _VC3_NULL else ''
            ord_ = ord_raw  if ord_raw  and ord_raw.lower()  not in _VC3_NULL else ''
            is_novel = any(str(v).lower().startswith('novel_')
                           for v in [fam_raw, gen_raw, ord_raw] if v)
            novel_anchor = ''
            if is_novel:
                cls_clean = (cls if cls and not cls.lower().startswith('novel_')
                             and cls.lower() not in _VC3_NULL else '')
                if cls_clean:
                    novel_anchor = cls_clean
                else:
                    for v in [fam_raw, gen_raw, ord_raw]:
                        if v and 'of_' in v:
                            cand = v.rsplit('of_', 1)[-1].strip()
                            if not cand.lower().startswith('novel_') and cand.lower() not in _VC3_NULL:
                                novel_anchor = cand
                                break
                fam = ''
                gen = ''
                ord_ = ''
            is_singleton = not fam and not gen and not ord_
            best = novel_anchor if is_novel else (gen or fam or ord_ or cls or rlm or '')
            vc_status = ('novel' if is_novel
                         else 'singleton' if is_singleton
                         else 'classified')
            records.append({'sample':   s,
                'Genome':       genome,
                'VC':           row.get('VC', ''),
                'VC_status':    vc_status,
                'Family':       fam,
                'Genus':        gen,
                'Order':        ord_,
                'Realm':        rlm,
                'Best_taxonomy': best,
                'Novel_anchor': novel_anchor,
                'Is_novel':     str(is_novel),
            })
    return records


# ── Viral taxonomy ────────────────────────────────────────────────────────────

def load_viral_taxonomy(paths, samples):
    _VC3_NULL = {"singleton", "unclassified", "nd", "none"}
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
            if final_family.lower() in _VC3_NULL: final_family = ''
            if final_genus.lower()  in _VC3_NULL: final_genus  = ''
            if final_order.lower()  in _VC3_NULL: final_order  = ''
            if source == 'vcontact3' and not final_family and not final_genus:
                source = 'unclassified'
            if source in ('unclassified', ''):
                continue
            if not final_family or not final_genus:
                lf, lg, lo = _deepest_level(lineage)
                if not final_family: final_family = lf
                if not final_genus:  final_genus  = lg
                if not final_order:  final_order  = lo
            if source == 'genomad' and not final_family and not final_genus and not final_order:
                gc = row.get('genomad_class', '') or row.get('genomad_best', '')
                if gc and gc.lower() not in _VC3_NULL:
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
    (IMG/VR or other user-supplied seqTaxDB), vcontact3, unknown (no hit)."""
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

def load_phist(paths, samples):
    records = []
    for p, s in zip(paths, samples):
        for row in load_csv(p):
            virus_raw = row.get('phage', row.get('Phage', row.get('virus', row.get('Virus', ''))))
            host_raw  = row.get('host',  row.get('Host', ''))
            score     = row.get('#common-kmers', row.get('Score', row.get('score', '')))
            pval      = row.get('adj-pvalue', row.get('pvalue', ''))
            if not virus_raw or not host_raw: continue
            virus_clean = os.path.basename(virus_raw).replace('.fasta', '').replace('.fa', '')
            virus_clean = virus_clean.replace('contig_', '', 1)
            host_clean  = os.path.basename(host_raw).replace('.fa', '').replace('.fasta', '')
            records.append({'sample': s, 'Virus': virus_clean, 'Host': host_clean,
                            'Score': str(score), 'P_value': str(pval)})
    return records


def build_host_collapse(phist_data, votu_abund_by_sample, samples):
    """Collapse viral RPKM by predicted host genus per sample.

    Returns {sample: [{genus, n_viruses, total_rpkm}]} sorted desc by total_rpkm.
    Host genus is the first underscore-delimited token of the PHIST host filename
    (e.g. 'Bacteroides_fragilis.fa' → 'Bacteroides').
    """
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
            genus = host.split('_')[0] if host else 'Unknown'
            if not genus: genus = 'Unknown'
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

def load_antidefensefinder_viral(paths, samples):
    """Anti-defense systems on viral proteins. Unlike the bin-side loader,
    the 'genome' column here is just the sample (defensefinder_viral makes
    one call across the whole sample's viral protein set, not per-genome) —
    the actual virus/contig is recovered from the first protein in
    protein_in_syst, which is reliable since every gene in one system sits
    on the same replicon."""
    records = []
    for p, s in zip(paths, samples):
        for row in load_tsv(p):
            sys_type = row.get('type', row.get('subtype', ''))
            proteins = _split_proteins_in_syst(row.get('protein_in_syst', ''))
            if not sys_type or not proteins: continue
            virus = _contig_from_protein_id(proteins[0])
            if not virus: continue
            records.append({'sample': s, 'Virus': virus, 'System': sys_type,
                             'System_id': row.get('sys_id', sys_type),
                             'Genes': row.get('genes_count', ''), 'Source': 'DefenseFinder'})
    return records


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


def load_dbapis_viral(paths, samples, apis_db_dir=''):
    """dbAPIS (Yan et al. 2023, NAR) DIAMOND blastp hits on viral proteins
    (rule dbapis_viral). Keeps only the best (lowest e-value) hit per query
    protein. sseqid is pipe-delimited: '{family_or_gene_id}|{IMGVR_UViG_id}|
    {genome_id}|{locus_with_coords}' (confirmed against a real run on
    litrp4, e.g. 'AcrIIA7|IMGVR_UViG_3300037418_004174|3300037418|
    Ga0395900_0000476_40112_40693') -- the first field alone is already a
    real, informative name (a dbAPIS family ID like 'APIS331', or a known
    gene name like 'AcrIIA7' for characterized Acr families).

    'Gene'/'Defense_system_inhibited' add the readable translation via
    _load_apis_family_map (seed_and_familyrep_all_infor.tsv) -- e.g.
    'APIS331' -> gene name + 'restriction-modification system' instead of
    just the bare family ID. Falls back to the raw Family/empty string if
    apis_db_dir wasn't configured or the mapping file isn't there yet."""
    fam_map = _load_apis_family_map(apis_db_dir)
    records = []
    for p, s in zip(paths, samples):
        best = {}
        for row in load_tsv(p):
            qseqid = row.get('qseqid', '')
            if not qseqid: continue
            try: evalue = float(row.get('evalue', '1') or '1')
            except ValueError: evalue = 1.0
            if qseqid not in best or evalue < best[qseqid][0]:
                best[qseqid] = (evalue, row)
        for qseqid, (evalue, row) in best.items():
            virus = _contig_from_protein_id(qseqid)
            if not virus: continue
            sseqid = row.get('sseqid', '')
            family = sseqid.split('|', 1)[0] if sseqid else ''
            gene, defense_system = fam_map.get(family, ('', ''))
            records.append({'sample': s, 'Virus': virus, 'Protein': qseqid,
                             'Family': family or sseqid, 'Gene': gene or family or sseqid,
                             'Defense_system_inhibited': defense_system,
                             'Hit': sseqid, 'Pident': row.get('pident', ''),
                             'Evalue': row.get('evalue', ''), 'Bitscore': row.get('bitscore', ''),
                             'Source': 'dbAPIS'})
    return records


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
    by_contig = defaultdict(list)
    if not faa_path or not os.path.exists(faa_path):
        return by_contig
    try:
        with open(faa_path) as f:
            for line in f:
                if line.startswith('>'):
                    pid = line[1:].split()[0].strip()
                    contig = _contig_from_protein_id(pid)
                    if contig:
                        by_contig[contig].append(pid)
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
            by_contig = _ordered_proteins_by_contig(faa)
            for contig, ordered_prots in by_contig.items():
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
                    window_genes = [
                        {'Protein': ordered_prots[i], 'Index': i, 'System': defense_idx.get(i, '')}
                        for i in range(start_idx, end_idx + 1)
                    ]
                    islands.append({
                        'sample': s, 'Bin': name, 'Contig': contig,
                        'n_genes': len(cluster), 'n_systems': len(systems),
                        'Systems': sorted(systems),
                        'window_genes': window_genes,
                        'start_idx': start_idx, 'end_idx': end_idx,
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


def load_amr_consensus(paths, samples):
    """Consolidated AMR hits from consolidate_amr.py (AMRFinderPlus + RGI +
    DeepARG merged by CDS locus via ARO).  One row per locus; consensus_score
    = n_tools / 3."""
    records = []
    for p, s in zip(paths, samples):
        for row in load_tsv(p):
            locus = row.get('locus', '')
            if not locus:
                continue
            genome, _ = _split_genome_prefix(locus)
            # contigs_pseudogenome = all-contig fallback, not a real bin; skip entirely.
            # Only include consensus calls (≥2 tools) to keep the report manageable.
            if genome == 'contigs_pseudogenome':
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


def merge_prok_taxonomy(gtdb_records, mmseqs_prok_records, checkm2_dict):
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

    # low_depth_mode only: if no CheckM2 bins were found at all, fall back to
    # per-contig MMseqs2-LCA records. When real bins exist, skip this fallback
    # to avoid flooding MERGED_PROK with tens-of-thousands of individual contigs.
    if not merged:
        for key, rec in mmseqs_bins.items():
            if key in seen: continue
            seen.add(key)
            base = dict(rec)
            base['Source_tax'] = 'MMseqs2-LCA'
            base['Completeness'] = ''; base['Contamination'] = ''; base['Genome_size'] = ''
            merged.append(base)

    return merged


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
    """Load Pharokka CDS TSV and count PHROGS categories per sample."""
    result = {}
    for s in samples:
        p = os.path.join(outdir, s, "annotation", "pharokka",
                         "pharokka_cds_final_merged_output.tsv")
        cnt = Counter()
        for row in load_tsv(p):
            cat = (row.get('phrog_category', '') or row.get('category', '') or
                   'unknown function').lower().strip()
            if not cat: cat = 'unknown function'
            cnt[cat] += 1
        result[s] = dict(cnt)
    return result


def load_svg(svg_path):
    if not os.path.exists(svg_path): return ""
    try:
        with open(svg_path, encoding='utf-8') as f: return f.read()
    except Exception: return ""


def load_genome_maps(outdir, samples):
    """Load genome map SVGs for virus/prok modes, max 5 per category per sample.

    "virus" merges the phage/ and virus/ output subfolders into a single list
    -- the report shows one unified "Virus" view, with each genome tagged
    category="Phage" or category="Virus" (the backend split still exists on
    disk, decided by PHROGS hallmark-gene evidence in genome_map_universal.py,
    but the UI no longer forces the user to pick a mode to see all of them).
    """
    result = {}
    for s in samples:
        base = os.path.join(outdir, s, "annotation", "genome_maps")
        result[s] = {"virus": [], "prok": []}
        for mode, category in (("phage", "Phage"), ("virus", "Virus")):
            mdir = os.path.join(base, mode)
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
                    result[s]["virus"].append({"id": gid, "svg": svg, "seq": seq, "category": category})
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
        "metaSPAdes":    _ver(f"{C} env_assembly spades.py --version"),
        "MMseqs2":       _ver(f"{C} env_assembly mmseqs version"),
        "QUAST":         _ver(f"{C} env_qc quast.py --version"),
        "BWA-MEM2":      _ver(f"{C} env_mapping bwa-mem2 version"),
        "minimap2":      _ver(f"{C} env_mapping minimap2 --version"),
        "samtools":      _ver(f"{C} env_mapping samtools --version"),
        "CoverM":        _ver(f"{C} env_coverm coverm --version"),
        "VirSorter2":    _ver(f"{C} env_viral virsorter --version"),
        "GeNomad":       _ver(f"{C} env_genomad genomad --version"),
        "VIBRANT":       _ver(f"{C} env_viral VIBRANT_run.py 2>&1 | head -3"),
        "CheckV":        _ver(f"{C} env_viral checkv 2>&1 | head -1"),
        "vRhyme":        _ver(f"{C} env_vrhyme vRhyme --version"),
        "vConTACT3":     _ver(f"{C} env_vcontact3 vcontact3 --version"),
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
