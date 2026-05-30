# ══════════════════════════════════════════════════════════════════════
# rules/taxonomy.smk — BLOCK 9: Viral Taxonomy
#
# Hierarchical 3-tier priority:
#   Tier 1 — vConTACT3     : protein-sharing network (genus-level, most specific)
#   Tier 2 — Diamond/INPHARED : metagenome phage references (family/genus)
#   Tier 3 — Diamond/Custom + GeNomad fallback (class/order for unknowns)
#
# Diamond confidence thresholds (avg pident):
#   >= 70% → genus | >= 50% → family | >= 30% → order_putative
#
# Also contains diamond_custom_prok for prokaryote bin taxonomy.
# ══════════════════════════════════════════════════════════════════════


rule prodigal_viral:
    """Predict ORFs from the non-redundant viral genome set for Diamond taxonomy searches."""
    input:
        viral = rules.viral_nonredundant.output.fasta,
    output:
        faa  = f"{OUTDIR}/{{sample}}/viral/taxonomy/viral_proteins.faa",
        done = f"{OUTDIR}/{{sample}}/viral/taxonomy/prodigal_done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/prodigal_viral.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/prodigal_viral.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("prodigal")
    threads: 1
    shell:
        """
        mkdir -p $(dirname {output.faa})
        if [ ! -s {input.viral} ]; then
            touch {output.faa} {output.done}; exit 0
        fi
        prodigal -i {input.viral} -a {output.faa} -p meta -f gff > {log} 2>&1
        touch {output.done}
        """


rule diamond_inphared:
    """
    Diamond BLASTp vs INPHARED phage protein database.
    Builds the Diamond DB on first run from *_vConTACT2_proteins.faa.
    Tier 2 in viral_taxonomy (after vConTACT3).
    """
    input:
        faa  = rules.prodigal_viral.output.faa,
        done = rules.prodigal_viral.output.done,
    output:
        hits = f"{OUTDIR}/{{sample}}/viral/taxonomy/diamond_vs_inphared.tsv",
        done = f"{OUTDIR}/{{sample}}/viral/taxonomy/diamond_done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/diamond_inphared.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/diamond_inphared.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("diamond")
    threads: THREADS
    params:
        db = f"{INPHARED_DB}/inphared_proteins.dmnd"
    shell:
        """
        if [ ! -s {input.faa} ]; then
            printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n" > {output.hits}
            touch {output.done}; exit 0
        fi
        INPHARED_PROT=$(ls {INPHARED_DB}/*_vConTACT2_proteins.faa 2>/dev/null | sort | tail -1)
        if [ -z "$INPHARED_PROT" ]; then
            echo "[diamond] ERROR: INPHARED proteins not found" | tee {log}; exit 1
        fi
        if [ ! -f {params.db} ]; then
            diamond makedb --in "$INPHARED_PROT" --db {params.db} --threads {threads} >> {log} 2>&1
        fi
        printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n" \
            > {output.hits}
        diamond blastp \
            --query {input.faa} --db {params.db} --out /dev/stdout \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --max-target-seqs 1 --evalue 1e-5 --id 30 --query-cover 50 \
            --threads {threads} --sensitive --tmpdir /tmp \
            >> {output.hits} 2>> {log}
        touch {output.done}
        """


rule diamond_custom_viral:
    """
    Optional Diamond BLASTp vs user-supplied viral DB (e.g. IMG NR viral subset).
    Tier 3 in viral taxonomy: vConTACT3 > INPHARED > CUSTOM > GeNomad.
    Skipped gracefully when CUSTOM_VIRAL_DMND = "" or file missing.
    Prepare DB: python3 scripts/prepare_diamond_db.py ...
    """
    input:
        faa  = rules.prodigal_viral.output.faa,
        done = rules.prodigal_viral.output.done,
    output:
        hits = f"{OUTDIR}/{{sample}}/viral/taxonomy/diamond_vs_custom.tsv",
        done = f"{OUTDIR}/{{sample}}/viral/taxonomy/custom_viral_done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/diamond_custom_viral.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/diamond_custom_viral.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("diamond")
    threads: THREADS
    params:
        db = CUSTOM_VIRAL_DMND
    shell:
        """
        if [ -z "{params.db}" ] || [ ! -f "{params.db}" ]; then
            echo "[diamond_custom_viral] No custom DB configured — skipping" | tee {log}
            printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n" > {output.hits}
            touch {output.done}; exit 0
        fi
        if [ ! -s {input.faa} ]; then
            printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n" > {output.hits}
            touch {output.done}; exit 0
        fi
        printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n" \
            > {output.hits}
        diamond blastp \
            --query {input.faa} --db {params.db} --out /dev/stdout \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --max-target-seqs 1 --evalue 1e-5 --id 30 --query-cover 50 \
            --threads {threads} --sensitive --tmpdir /tmp \
            >> {output.hits} 2>> {log}
        touch {output.done}
        """


rule diamond_custom_prok:
    """
    Optional Diamond BLASTp vs user-supplied prokaryote DB (e.g. IMG NR).
    Classifies Binette MAGs not assigned by GTDB-Tk.
    Skipped gracefully when CUSTOM_PROK_DMND = "" or file missing.
    Reuses env_viral (prodigal + diamond are in that env).
    """
    input:
        faa      = rules.prodigal_viral.output.faa,  # reuse prodigal env; bins processed in shell
        bin_done = rules.binette.output.done,         # ensure bins exist before running
    output:
        hits = f"{OUTDIR}/{{sample}}/bins/diamond_custom_prok/diamond_vs_custom.tsv",
        done = f"{OUTDIR}/{{sample}}/bins/diamond_custom_prok/done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/diamond_custom_prok.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/diamond_custom_prok.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("diamond")
    threads: THREADS
    params:
        db       = CUSTOM_PROK_DMND,
        bins_dir = lambda wc: f"{OUTDIR}/{wc.sample}/bins/binette/final_bins",
        prok_faa = f"{OUTDIR}/{{sample}}/bins/diamond_custom_prok/all_bins.faa",
        outdir   = f"{OUTDIR}/{{sample}}/bins/diamond_custom_prok",
    shell:
        """
        if [ -z "{params.db}" ] || [ ! -f "{params.db}" ]; then
            echo "[diamond_custom_prok] No custom DB configured — skipping" | tee {log}
            printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n" > {output.hits}
            touch {output.done}; exit 0
        fi

        mkdir -p {params.outdir}

        # Predict proteins from all Binette bins (if not already done)
        if [ ! -s {params.prok_faa} ]; then
            for bin in {params.bins_dir}/*.fa; do
                [ -f "$bin" ] || continue
                bname=$(basename "$bin" .fa)
                prodigal -i "$bin" -a {params.outdir}/$bname.faa \
                         -p single -f gff -q 2>/dev/null || true
                [ -f {params.outdir}/$bname.faa ] && \
                    cat {params.outdir}/$bname.faa >> {params.prok_faa}
            done
        fi
        if [ ! -s {params.prok_faa} ]; then
            echo "[diamond_custom_prok] No bin proteins found" | tee -a {log}
            printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n" > {output.hits}
            touch {output.done}; exit 0
        fi
        printf "qseqid\tsseqid\tpident\tlength\tmismatch\tgapopen\tqstart\tqend\tsstart\tsend\tevalue\tbitscore\n" \
            > {output.hits}
        diamond blastp \
            --query {params.prok_faa} --db {params.db} --out /dev/stdout \
            --outfmt 6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore \
            --max-target-seqs 1 --evalue 1e-5 --id 50 --query-cover 70 \
            --threads {threads} --sensitive --tmpdir /tmp \
            >> {output.hits} 2>> {log}
        touch {output.done}
        """


rule vcontact3:
    """
    vConTACT3: protein-sharing network clustering, genus-level.
    Tier 1 in viral_taxonomy — most specific classification.
    Only runs on CheckV >= MQ (>= 50% complete) to reduce noise.
    Filtering: scripts/filter_checkv_hq.py (Tier1=MQ+, Tier2=length>=5kb).
    Options: SqRoot metric, 5 iterations, --reduce-memory, family+genus ranks.
    """
    input:
        viral  = rules.viral_nonredundant.output.fasta,
        checkv = rules.checkv.output.summary,
    output:
        done    = f"{OUTDIR}/{{sample}}/viral/vcontact3/done.txt",
        network = f"{OUTDIR}/{{sample}}/viral/vcontact3/genome_clusters.tsv",
    log:   f"{OUTDIR}/{{sample}}/logs/vcontact3.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/vcontact3.tsv"
    conda: "../envs/env_vcontact3.yaml"
    container:  CONTAINERS.get("vcontact3")
    threads: THREADS
    params:
        outdir      = f"{OUTDIR}/{{sample}}/viral/vcontact3",
        hq_fa       = f"{OUTDIR}/{{sample}}/viral/vcontact3/hq_viral.fasta",
        vc3_db      = VCONTACT3_DB,
        vc3_ver     = VCONTACT3_VER,
        scripts_dir = SCRIPTS_DIR,
    shell:
        """
        set -euo pipefail
        mkdir -p {params.outdir}

        python3 {params.scripts_dir}/filter_checkv_hq.py \
            {input.checkv} {input.viral} {params.hq_fa} \
            >> {log} 2>&1
        # Note: filter_checkv_hq.py keeps Complete/HQ/MQ + completeness>=50%
        # Unbinned short contigs (<50% complete) are still classified by
        # Diamond/INPHARED and GeNomad in rule viral_taxonomy

        if [ ! -s {params.hq_fa} ]; then
            echo "[vcontact3] No HQ/MQ genomes — skipping" | tee -a {log}
            mkdir -p {params.outdir}/vConTACT3_results
            printf "genome\tVC\tVC_status\tgenus\tfamily\torder\n" > {output.network}
            touch {output.done}; exit 0
        fi

        N=$(grep -c "^>" {params.hq_fa} || echo 0)
        echo "[vcontact3] Running on $N genomes" | tee -a {log}

        vcontact3 run \
            --nucleotide      {params.hq_fa} \
            --output          {params.outdir}/vConTACT3_results \
            --threads         {threads} \
            --db-path         {params.vc3_db} \
            --db-version      {params.vc3_ver} \
            --db-domain       prokaryotes \
            --distance-metric SqRoot \
            --max-iterations  5 \
            --reduce-memory \
            --target-rank family genus \
            --exports profiles completeness \
            --no-progress \
            --force-overwrite \
            >> {log} 2>&1 || true

        # vConTACT3 v3 writes final_assignments.csv in exports/
        VC3_OUT="{params.outdir}/vConTACT3_results/exports/final_assignments.csv"

        if [ -f "$VC3_OUT" ]; then
            cp "$VC3_OUT" {output.network}
            echo "[vcontact3] Output: $VC3_OUT ($(wc -l < $VC3_OUT) lines)" | tee -a {log}
        else
            echo "[vcontact3] final_assignments.csv not found — listing:" | tee -a {log}
            find {params.outdir}/vConTACT3_results -type f | tee -a {log}
            printf "Genome,family_prediction,genus_prediction,realm_prediction,order_prediction\n" > {output.network}
        fi
        touch {output.done}
        """


rule viral_taxonomy:
    """
    Merge taxonomy from all three tiers into one table per contig.
    Priority: vConTACT3 > Diamond/INPHARED > Diamond/Custom > GeNomad.
    GeNomad lineage parsed to deepest available level as final fallback.
    Output: viral_taxonomy_merged.tsv
    """
    input:
        genomad_done   = rules.genomad.output.done,
        diamond_hits   = rules.diamond_inphared.output.hits,
        diamond_done   = rules.diamond_inphared.output.done,
        custom_hits    = rules.diamond_custom_viral.output.hits,
        custom_done    = rules.diamond_custom_viral.output.done,
        vcontact3_net  = rules.vcontact3.output.network,
        vcontact3_done = rules.vcontact3.output.done,
        viral          = rules.viral_nonredundant.output.fasta,
    output:
        tsv  = f"{OUTDIR}/{{sample}}/viral/taxonomy/viral_taxonomy_merged.tsv",
        done = f"{OUTDIR}/{{sample}}/viral/taxonomy/taxonomy_done.txt",
    log:   f"{OUTDIR}/{{sample}}/logs/viral_taxonomy.log"
    benchmark: f"{OUTDIR}/{{sample}}/benchmarks/viral_taxonomy.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("diamond")
    threads: 1
    params:
        inphared_db = INPHARED_DB, custom_viral_meta = CUSTOM_VIRAL_META
    run:
        import csv, os, glob, collections
        from pathlib import Path

        lf = open(str(log[0]), "w")

        # All viral contigs
        contigs = []
        if os.path.exists(str(input.viral)):
            with open(str(input.viral)) as f:
                for line in f:
                    if line.startswith(">"): contigs.append(line[1:].split()[0])
        lf.write(f"Total viral contigs: {len(contigs)}\n")

        # ── Tier 1: vConTACT3 ─────────────────────────────────────────
        vc3_tax = {}
        vc3_path = str(input.vcontact3_net)
        if os.path.exists(vc3_path) and os.path.getsize(vc3_path) > 100:
            with open(vc3_path) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    g = row.get("genome", row.get("Genome","")).strip()
                    if not g: continue
                    vc3_tax[g] = {
                        "vc":     row.get("VC", row.get("cluster","")),
                        "status": row.get("VC_status", row.get("status","")),
                        "genus":  row.get("genus",  row.get("Genus","")),
                        "family": row.get("family", row.get("Family","")),
                        "order":  row.get("order",  row.get("Order","")),
                    }
        lf.write(f"vConTACT3: {len(vc3_tax)} genomes\n")

        # ── Tier 2: Diamond/INPHARED ──────────────────────────────────
        inphared_meta = {}
        data_files = sorted(glob.glob(os.path.join(str(params.inphared_db), "*_data.tsv")))
        if data_files:
            with open(data_files[-1]) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    acc = row.get("Accession","").strip()
                    if acc:
                        inphared_meta[acc] = {
                            "family": row.get("Family",""),
                            "genus":  row.get("Genus",""),
                            "order":  row.get("Order",""),
                        }
        lf.write(f"INPHARED entries: {len(inphared_meta)}\n")

        votes      = collections.defaultdict(collections.Counter)
        pident     = collections.defaultdict(list)
        ssci_map   = collections.defaultdict(list)
        skingd_map = collections.defaultdict(str)
        dpath  = str(input.diamond_hits)
        if os.path.exists(dpath) and os.path.getsize(dpath) > 0:
            with open(dpath) as f:
                for line in f:
                    if not line.strip() or line.startswith("qseqid"): continue
                    cols = line.strip().split("\t")
                    if len(cols) < 3: continue
                    contig = "_".join(cols[0].split("_")[:-1]) or cols[0]
                    acc    = "_".join(cols[1].split("_")[:-1]) or cols[1]
                    votes[contig][acc] += 1
                    pident[contig].append(float(cols[2]))
                    # cols[13]=sscinames, cols[14]=sskingdoms (new outfmt)
                    if len(cols) >= 14 and cols[13] not in ("N/A", "", "0"):
                        ssci_map[contig].append(cols[13])
                    if len(cols) >= 15 and cols[14] not in ("N/A", "", "0"):
                        skingd_map[contig] = cols[14]

        inphared_tax = {}
        for contig, v in votes.items():
            top_acc, top_votes = v.most_common(1)[0]
            meta = inphared_meta.get(top_acc, {})
            avg  = sum(pident[contig]) / len(pident[contig])
            if avg >= 70:   conf = "genus"
            elif avg >= 50: conf = "family"
            elif avg >= 30: conf = "order_putative"
            else:           conf = "unclassified"
            inphared_tax[contig] = {
                "family": meta.get("family",""), "genus": meta.get("genus",""),
                "order":  meta.get("order",""),  "avg_pident": round(avg, 2),
                "votes": top_votes, "top_acc": top_acc, "confidence": conf,
                "sscinames":  ";".join(dict.fromkeys(ssci_map.get(contig, []))),
                "sskingdoms": skingd_map.get(contig, ""),
            }
        lf.write(f"Diamond/INPHARED: {len(inphared_tax)} contigs\n")

        # ── Tier 3: Diamond/Custom DB ──────────────────────────────────
        custom_meta = {}
        if params.custom_viral_meta and os.path.exists(str(params.custom_viral_meta)):
            with open(str(params.custom_viral_meta)) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    acc = row.get("accession","").strip()
                    if acc:
                        custom_meta[acc] = {
                            "family": row.get("family",""),
                            "genus":  row.get("genus",""),
                            "order":  row.get("order",""),
                        }
        lf.write(f"Custom DB entries: {len(custom_meta)}\n")

        # Reuse same Diamond parse logic — majority vote per contig
        custom_votes  = collections.defaultdict(collections.Counter)
        custom_pident = collections.defaultdict(list)
        custom_ssci   = collections.defaultdict(list)
        cpath = str(input.custom_hits)
        if os.path.exists(cpath) and os.path.getsize(cpath) > 0:
            with open(cpath) as f:
                for line in f:
                    if not line.strip() or line.startswith("qseqid"): continue
                    cols = line.strip().split("\t")
                    if len(cols) < 3: continue
                    contig = "_".join(cols[0].split("_")[:-1]) or cols[0]
                    acc    = "_".join(cols[1].split("_")[:-1]) or cols[1]
                    custom_votes[contig][acc] += 1
                    custom_pident[contig].append(float(cols[2]))
                    if len(cols) >= 14 and cols[13] not in ("N/A", "", "0"):
                        custom_ssci[contig].append(cols[13])

        custom_tax = {}
        for contig, v in custom_votes.items():
            top_acc, top_votes = v.most_common(1)[0]
            meta = custom_meta.get(top_acc, {})
            avg  = sum(custom_pident[contig]) / len(custom_pident[contig])
            if avg >= 70:   conf = "genus"
            elif avg >= 50: conf = "family"
            elif avg >= 30: conf = "order_putative"
            else:           conf = "unclassified"
            custom_tax[contig] = {
                "family": meta.get("family",""), "genus": meta.get("genus",""),
                "order":  meta.get("order",""),  "avg_pident": round(avg,2),
                "top_acc": top_acc, "confidence": conf,
                "sscinames": ";".join(dict.fromkeys(custom_ssci.get(contig, []))),
            }

        # ── Tier 3: GeNomad fallback ──────────────────────────────────
        # Path derivado de input.genomad_done (idêntico ao Snakefile funcional)
        genomad_tax = {}
        gdir   = os.path.dirname(str(input.genomad_done))
        gfiles = glob.glob(os.path.join(gdir, "**", "*_virus_summary.tsv"), recursive=True)
        gpath  = gfiles[0] if gfiles else ""
        if gpath and os.path.getsize(gpath) > 0:
            with open(gpath) as f:
                for row in csv.DictReader(f, delimiter="\t"):
                    name  = row.get("seq_name","")
                    tax   = row.get("taxonomy","")
                    score = row.get("virus_score","0")
                    if not name or not tax: continue
                    parts = [p.strip() for p in tax.split(";")
                             if p.strip() and p.strip() != "Viruses"]
                    family=""; genus=""; order=""; cls=""; best=""
                    for p in parts:
                        if   p.endswith("viridae"):                          family = p
                        elif p.endswith("virales"):                          order  = p
                        elif p.endswith("viricetes") or p.endswith("viria"): cls   = p
                        elif p.endswith("virus") or p.endswith("phage"):     genus  = p
                    best = genus or family or order or cls or (parts[-1] if parts else "")
                    genomad_tax[name] = {
                        "family": family, "genus": genus, "order": order,
                        "class": cls, "best": best, "lineage": tax,
                        "score": float(score or 0),
                    }
        lf.write(f"GeNomad: {len(genomad_tax)} contigs\n")

        # ── Build final table ─────────────────────────────────────────
        rows = []; stats = collections.Counter()
        for contig in contigs:
            vc3 = vc3_tax.get(contig, {})
            iph = inphared_tax.get(contig, {})
            gmd = genomad_tax.get(contig, {})

            if vc3 and (vc3.get("family") or vc3.get("genus")):
                source = "vcontact3"
                ff, fg = vc3.get("family",""), vc3.get("genus","")
                fo     = vc3.get("order","")
                conf   = vc3.get("status","")
                lin    = ";".join(filter(None, [fo, ff, fg]))
                best   = fg or ff or fo

            elif (iph and iph.get("confidence") and
                  iph.get("confidence") not in ("unclassified", "order_putative")):
                source = "diamond_inphared"
                ff, fg = iph.get("family",""), iph.get("genus","")
                fo     = iph.get("order","")
                conf   = f"{iph['avg_pident']:.1f}% ({iph['confidence']})"
                lin    = ";".join(filter(None, [fo, ff, fg]))
                best   = fg or ff or fo

            elif (contig in custom_tax and
                  custom_tax[contig].get("confidence") not in ("unclassified", "order_putative", "")):
                _c     = custom_tax[contig]
                source = "diamond_custom"
                ff, fg = _c.get("family",""), _c.get("genus","")
                fo     = _c.get("order","")
                conf   = f"{float(_c.get('avg_pident',0) or 0):.1f}% ({_c.get('confidence','')})"
                lin    = ";".join(filter(None, [fo, ff, fg]))
                best   = fg or ff or fo

            elif gmd:
                source = "genomad"
                ff     = gmd.get("family","") or gmd.get("order","") or gmd.get("class","")
                fg     = gmd.get("genus","")
                fo     = gmd.get("order","")
                best   = gmd.get("best","")
                conf   = f"{gmd['score']:.3f}"
                lin    = gmd.get("lineage","")

            else:
                source = "unclassified"
                ff = fg = fo = conf = lin = best = ""

            stats[source] += 1
            rows.append({
                "seq_name":      contig,
                "final_family":  ff,
                "final_genus":   fg,
                "final_order":   fo,
                "best_taxonomy": best,
                "source":        source,
                "confidence":    conf,
                "lineage":       lin,
                "vc3_vc":        vc3.get("vc",""),
                "vc3_status":    vc3.get("status",""),
                "genomad_best":  gmd.get("best",""),
                "genomad_class": gmd.get("class",""),
                "genomad_score": gmd.get("score",""),
                "inphared_pident":    iph.get("avg_pident",""),
                "inphared_votes":     iph.get("votes",""),
                "inphared_sscinames": iph.get("sscinames",""),
                "custom_acc":         custom_tax.get(contig,{}).get("top_acc",""),
                "custom_pident":      custom_tax.get(contig,{}).get("avg_pident",""),
                "custom_sscinames":   custom_tax.get(contig,{}).get("sscinames",""),
            })

        lf.write("\nSummary:\n")
        for k, v in stats.most_common():
            lf.write(f"  {k}: {v}\n")
        total = len(rows)
        unclass = stats.get("unclassified", 0)
        lf.write(f"  Novel (unclassified): {unclass}/{total} = {100*unclass/total:.1f}%\n" if total else "")
        lf.close()

        fields = ["seq_name","final_family","final_genus","final_order","best_taxonomy",
                  "source","confidence","lineage",
                  "vc3_vc","vc3_status",
                  "genomad_best","genomad_class","genomad_score",
                  "inphared_pident","inphared_votes","inphared_sscinames",
                  "custom_acc","custom_pident","custom_sscinames"]
        with open(str(output.tsv), "w", newline="") as f:
            w = csv.DictWriter(f, fieldnames=fields, delimiter="\t")
            w.writeheader(); w.writerows(rows)
        Path(str(output.done)).write_text("ok\n")
