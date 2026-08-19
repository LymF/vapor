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
        quality    = rules.mag_catalog_quality.output.tsv,
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
