# ══════════════════════════════════════════════════════════════════════
# rules/pangenome.smk — Fase 1 do pangenoma dos clusters com ilha
#
# Ver docs/superpowers/specs/2026-08-19-pangenoma-clusters-defesa-design.md
#
# PADRAO OBRIGATORIO: selecao dependente de dados DENTRO de um job de
# numero fixo, como `mag_bakta` faz com qualifying_bins.txt. NUNCA um
# `checkpoint` -- o DAG dinamico quebraria a invariante de dry-run que o
# roadmap usa para verificar toda mudanca.
# ══════════════════════════════════════════════════════════════════════

import os as _os

PANGENOME_DIR = f"{MAG_CATALOG_DIR}/pangenome"


rule mag_pangenome_select:
    """Quais clusters merecem anotacao por membro, e por que.

    Portao: >= 3 membros E (ilha OU >= 3 sistemas OU ARG de consenso).
    O PlasmidFinder e registrado como sinal de mobilidade mas nao elege
    sozinho -- plasmidio sem defesa nem ARG nao motiva um pangenoma.
    """
    input:
        membership = rules.mag_catalog_membership.output.tsv,
        manifest   = rules.mag_catalog_proteins.output.manifest,
        df_systems = rules.mag_defensefinder.output.systems,
        consensus  = rules.mag_amr_consensus.output.consensus,
        plasmid    = rules.mag_abricate.output.plasmidfinder,
    output:
        candidates = f"{PANGENOME_DIR}/candidates.tsv",
        members    = f"{PANGENOME_DIR}/members.txt",
        done       = f"{PANGENOME_DIR}/select_done.txt",
    log:
        f"{OUTDIR}/logs/mag_pangenome_select.log"
    run:
        import csv as _csv
        import sys as _sys
        from collections import defaultdict

        _sys.path.insert(0, SCRIPTS_DIR)
        from pangenome_select import (load_membership, select_clusters)
        from defense_islands import find_islands, genes_by_contig
        from mag_catalog import resolve_prefixed_id

        _os.makedirs(f"{PANGENOME_DIR}", exist_ok=True)

        membership = load_membership(str(input.membership))
        known_reps = set(membership)

        # ── evidencia por representante ──────────────────────────────────
        # 1. sistemas de defesa e as proteinas de cada um (para a ilha)
        n_systems = defaultdict(int)
        prot_to_sys = defaultdict(dict)
        with open(str(input.df_systems), newline="") as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                genome = (row.get("genome") or "").strip()
                stype  = (row.get("type") or row.get("subtype") or "").strip()
                if not genome or not stype:
                    continue
                n_systems[genome] += 1
                for prot in (row.get("protein_in_syst") or "").split(","):
                    prot = prot.strip()
                    if prot:
                        prot_to_sys[genome][prot] = (
                            genome, stype, row.get("sys_id", stype))

        # 2. ilhas: precisa do .faa de cada representante, via manifesto
        n_islands = defaultdict(int)
        with open(str(input.manifest)) as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 5:
                    continue
                name, _mode, _fna, faa, _gff = parts[:5]
                if name not in prot_to_sys:
                    continue
                n_islands[name] = len(find_islands(genes_by_contig(faa),
                                                   prot_to_sys[name]))

        # 3. ARG de consenso (n_tools >= 2) e hits de plasmidio
        n_args = defaultdict(int)
        with open(str(input.consensus), newline="") as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                try:
                    if int(row.get("n_tools") or 0) < 2:
                        continue
                except ValueError:
                    continue
                locus = (row.get("locus") or "").strip()
                # NUNCA cortar em separador: o locus e {genome}__{protein},
                # e o genome ja contem "__" ({source}__{bin}). Casar contra
                # os nomes conhecidos, como as vistas fazem.
                genome, _rest = resolve_prefixed_id(locus, known_reps)
                if genome:
                    n_args[genome] += 1

        n_plasmid = defaultdict(int)
        with open(str(input.plasmid), newline="") as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                genome = (row.get("genome") or "").strip()
                if genome:
                    n_plasmid[genome] += 1

        evidence = {rep: {"n_islands": n_islands.get(rep, 0),
                          "n_systems": n_systems.get(rep, 0),
                          "n_args":    n_args.get(rep, 0),
                          "n_plasmid": n_plasmid.get(rep, 0)}
                    for rep in membership}

        rows = select_clusters(membership, evidence)

        cols = ["representative_id", "n_members", "n_islands", "n_systems",
                "n_args", "n_plasmid", "criterio", "eligible"]
        with open(str(output.candidates), "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for r in rows:
                w.writerow(r)

        eligible = [r for r in rows if r["eligible"]]
        with open(str(output.members), "w") as f:
            for r in eligible:
                for m in membership[r["representative_id"]]:
                    f.write(m["member_id"] + "\n")

        n_mem = sum(len(membership[r["representative_id"]]) for r in eligible)
        with open(str(log[0]), "w") as lf:
            lf.write(f"[pangenome_select] {len(rows)} clusters avaliados, "
                     f"{len(eligible)} elegiveis, {n_mem} membros a anotar\n")
            if not eligible:
                lf.write("[pangenome_select] nenhum cluster elegivel -- as "
                         "regras seguintes vao pular. Nao e erro: com poucas "
                         "especies compartilhadas entre amostras o catalogo "
                         "pode nao ter cluster com 3+ membros.\n")

        write_status(str(output.done),
                     "ok" if eligible else "skipped: no eligible clusters")


rule mag_pangenome_proteins:
    """Prodigal nos MEMBROS dos clusters elegiveis.

    Gemea de `mag_catalog_proteins` com outro diretorio de entrada: os
    genomas do POOL (`mag_catalog/genomes/`), filtrados por members.txt.
    Mesmo formato de manifesto, de proposito -- e o que deixa
    `mag_pangenome_defensefinder` e `mag_pangenome_amr` serem herdadas com
    `use rule` sem alterar uma linha do corpo delas.
    """
    input:
        members = rules.mag_pangenome_select.output.members,
        select  = rules.mag_pangenome_select.output.done,
    output:
        manifest = f"{PANGENOME_DIR}/proteins/manifest.txt",
        done     = f"{PANGENOME_DIR}/proteins/done.txt",
    log:
        f"{OUTDIR}/logs/mag_pangenome_proteins.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_pangenome_proteins.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("prodigal")
    threads: 1
    params:
        genomes_dir = f"{MAG_CATALOG_DIR}/genomes",
        outdir      = f"{PANGENOME_DIR}/proteins",
    run:
        _os.makedirs(params.outdir, exist_ok=True)
        rows = []
        with open(str(log[0]), "w") as lf:
            names = [ln.strip() for ln in open(str(input.members))
                     if ln.strip()]
            lf.write(f"[pangenome_proteins] {len(names)} membros a anotar\n")
            for name in names:
                fna = _os.path.join(params.genomes_dir, f"{name}.fa")
                if not _os.path.exists(fna):
                    lf.write(f"[pangenome_proteins] AUSENTE: {fna}\n")
                    continue
                faa = _os.path.join(params.outdir, f"{name}.faa")
                gff = _os.path.join(params.outdir, f"{name}.gff")
                shell("prodigal -i {fna} -a {faa} -f gff -o {gff} "
                      "-p single -q >> {log} 2>&1 || true")
                if _os.path.exists(faa) and _os.path.getsize(faa) > 0:
                    rows.append((name, "bin", fna, faa, gff))
                else:
                    lf.write(f"[pangenome_proteins] FALHOU: {name} "
                             f"(prodigal nao gerou {faa})\n")

            with open(str(output.manifest), "w") as mf:
                for r in rows:
                    mf.write("\t".join(r) + "\n")
            lf.write(f"[pangenome_proteins] manifesto com {len(rows)} genomas\n")

        if not names:
            write_status(str(output.done), "skipped: no eligible members")
        elif not rows:
            write_status(str(output.done),
                         "failed: %d membros, 0 proteomas" % len(names))
        else:
            write_status(str(output.done), "ok")


# Defesa e AMR por membro, herdadas das regras globais do catalogo.
# use rule ... as ... with: substitui apenas input/output/log/benchmark; o
# corpo (run:) vem inteiro da regra original e referencia input.<nome>
# literalmente. Os corpos ja consomem manifesto e nao conhecem wildcard,
# por isso apontam sem alteracao para o manifesto dos MEMBROS em vez do
# manifesto das representantes.
#
# mag_pangenome_argnorm e mag_pangenome_amr_consensus encadeiam com as
# regras mag_pangenome_* anteriores (amrfinderplus/rgi/deeparg deste
# arquivo), NUNCA com as globais mag_* do catalogo -- reapontar para as
# globais misturaria resultado de membro com o de representante, em
# silencio.

use rule mag_defensefinder as mag_pangenome_defensefinder with:
    input:
        manifest = rules.mag_pangenome_proteins.output.manifest,
        done     = rules.mag_pangenome_proteins.output.done,
    output:
        done        = f"{PANGENOME_DIR}/defensefinder/done.txt",
        systems     = f"{PANGENOME_DIR}/defensefinder/defensefinder_systems.tsv",
        antisystems = f"{PANGENOME_DIR}/defensefinder/antidefensefinder_systems.tsv",
    log:
        f"{OUTDIR}/logs/mag_pangenome_defensefinder.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_pangenome_defensefinder.tsv"


use rule mag_amrfinderplus as mag_pangenome_amrfinderplus with:
    input:
        manifest = rules.mag_pangenome_proteins.output.manifest,
        done     = rules.mag_pangenome_proteins.output.done,
    output:
        done    = f"{PANGENOME_DIR}/amrfinderplus/done.txt",
        results = f"{PANGENOME_DIR}/amrfinderplus/amrfinder_results.tsv",
    log:
        f"{OUTDIR}/logs/mag_pangenome_amrfinderplus.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_pangenome_amrfinderplus.tsv"


use rule mag_rgi_card as mag_pangenome_rgi_card with:
    input:
        manifest = rules.mag_pangenome_proteins.output.manifest,
        done     = rules.mag_pangenome_proteins.output.done,
    output:
        done    = f"{PANGENOME_DIR}/rgi/done.txt",
        results = f"{PANGENOME_DIR}/rgi/rgi_results.txt",
    log:
        f"{OUTDIR}/logs/mag_pangenome_rgi.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_pangenome_rgi.tsv"


use rule mag_deeparg as mag_pangenome_deeparg with:
    input:
        manifest = rules.mag_pangenome_proteins.output.manifest,
        done     = rules.mag_pangenome_proteins.output.done,
    output:
        done    = f"{PANGENOME_DIR}/deeparg/done.txt",
        results = f"{PANGENOME_DIR}/deeparg/deeparg_results.mapping.ARG",
    log:
        f"{OUTDIR}/logs/mag_pangenome_deeparg.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_pangenome_deeparg.tsv"


use rule mag_argnorm_normalize as mag_pangenome_argnorm with:
    input:
        amrfinder      = rules.mag_pangenome_amrfinderplus.output.results,
        amrfinder_done = rules.mag_pangenome_amrfinderplus.output.done,
        deeparg        = rules.mag_pangenome_deeparg.output.results,
        deeparg_done   = rules.mag_pangenome_deeparg.output.done,
    output:
        done             = f"{PANGENOME_DIR}/argnorm/done.txt",
        amrfinder_normed = f"{PANGENOME_DIR}/argnorm/amrfinderplus_normed.tsv",
        deeparg_normed   = f"{PANGENOME_DIR}/argnorm/deeparg_normed.tsv",
    log:
        f"{OUTDIR}/logs/mag_pangenome_argnorm.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_pangenome_argnorm.tsv"


use rule mag_amr_consensus as mag_pangenome_amr_consensus with:
    input:
        argnorm_done     = rules.mag_pangenome_argnorm.output.done,
        rgi_done         = rules.mag_pangenome_rgi_card.output.done,
        amrfinder_normed = rules.mag_pangenome_argnorm.output.amrfinder_normed,
        deeparg_normed   = rules.mag_pangenome_argnorm.output.deeparg_normed,
        rgi_results      = rules.mag_pangenome_rgi_card.output.results,
    output:
        done      = f"{PANGENOME_DIR}/amr_consensus/done.txt",
        consensus = f"{PANGENOME_DIR}/amr_consensus/amr_consensus.tsv",
    log:
        f"{OUTDIR}/logs/mag_pangenome_amr_consensus.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_pangenome_amr_consensus.tsv"


rule mag_pangenome_matrix:
    """Matriz gene x membro e sumario por cluster.

    O sumario e o que decide a fase 2: ela so se justifica se houver
    cluster com >= 5 membros AVALIAVEIS e variacao real
    (n_genes_variaveis > 0). Com 3-6 membros o PPanGGOLiN roda, mas a
    separacao shell/cloud nao tem sustentacao (recomendado: >= 15).
    """
    input:
        candidates = rules.mag_pangenome_select.output.candidates,
        membership = rules.mag_catalog_membership.output.tsv,
        quality    = rules.mag_catalog_quality.output.tsv,
        df_systems = rules.mag_pangenome_defensefinder.output.systems,
        consensus  = rules.mag_pangenome_amr_consensus.output.consensus,
    output:
        matrix  = f"{PANGENOME_DIR}/gene_by_member.tsv",
        summary = f"{PANGENOME_DIR}/cluster_summary.tsv",
        done    = f"{PANGENOME_DIR}/done.txt",
    log:
        f"{OUTDIR}/logs/mag_pangenome_matrix.log"
    run:
        import csv as _csv
        import sys as _sys
        from collections import defaultdict

        _sys.path.insert(0, SCRIPTS_DIR)
        from pangenome_select import load_membership, load_completeness
        from pangenome_matrix import build_matrix, summarize_clusters
        from mag_catalog import resolve_prefixed_id

        membership   = load_membership(str(input.membership))
        completeness = load_completeness(str(input.quality))

        clusters = []
        with open(str(input.candidates), newline="") as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                if (row.get("eligible") or "").strip() in ("True", "true", "1"):
                    clusters.append(row["representative_id"])

        members_by_rep = {rep: [m["member_id"] for m in membership.get(rep, [])]
                          for rep in clusters}
        known_members = {m for ms in members_by_rep.values() for m in ms}

        gene_hits = defaultdict(set)
        with open(str(input.df_systems), newline="") as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                genome = (row.get("genome") or "").strip()
                stype  = (row.get("type") or row.get("subtype") or "").strip()
                if genome and stype:
                    gene_hits[genome].add(stype)

        with open(str(input.consensus), newline="") as f:
            for row in _csv.DictReader(f, delimiter="\t"):
                try:
                    if int(row.get("n_tools") or 0) < 2:
                        continue
                except ValueError:
                    continue
                locus = (row.get("locus") or "").strip()
                gene  = (row.get("gene_name") or "").strip()
                # Idem: casar contra os MEMBROS conhecidos, nunca cortar.
                genome, _rest = resolve_prefixed_id(locus, known_members)
                if genome and gene:
                    gene_hits[genome].add(gene)

        rows = build_matrix(clusters, members_by_rep, gene_hits, completeness)
        summary = summarize_clusters(rows, members_by_rep, completeness)

        # A completude viaja no cabecalho: "4/6" so e interpretavel com o
        # denominador a vista.
        all_members = sorted({m for ms in members_by_rep.values() for m in ms})
        with open(str(output.matrix), "w") as f:
            f.write("# completude: " + ", ".join(
                f"{m}={completeness.get(m, 0.0):.1f}" for m in all_members) + "\n")
            f.write("cluster\tgene\tfreq\tn_present\tn_evaluable\t"
                    + "\t".join(all_members) + "\n")
            for r in rows:
                states = [r["states"].get(m, "") for m in all_members]
                f.write("\t".join([r["representative_id"], r["gene"], r["freq"],
                                   str(r["n_present"]), str(r["n_evaluable"])]
                                  + states) + "\n")

        cols = ["representative_id", "n_members", "n_members_avaliaveis",
                "n_genes_core", "n_genes_variaveis", "n_genes_singleton"]
        with open(str(output.summary), "w", newline="") as f:
            w = _csv.DictWriter(f, fieldnames=cols, delimiter="\t")
            w.writeheader()
            for s in summary:
                w.writerow(s)

        fase2 = [s for s in summary
                 if s["n_members_avaliaveis"] >= 5 and s["n_genes_variaveis"] > 0]
        with open(str(log[0]), "w") as lf:
            lf.write(f"[pangenome_matrix] {len(clusters)} clusters, "
                     f"{len(rows)} linhas de gene\n")
            lf.write(f"[pangenome_matrix] clusters que sustentariam a fase 2 "
                     f"(>= 5 membros avaliaveis e variacao): {len(fase2)}\n")
            if not fase2:
                lf.write("[pangenome_matrix] nenhum. A matriz acima e a "
                         "resposta honesta; o PPanGGOLiN sobre 3-4 membros "
                         "produziria particao sem sustentacao.\n")

        write_status(str(output.done),
                     "ok" if clusters else "skipped: no eligible clusters")
