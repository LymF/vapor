# ══════════════════════════════════════════════════════════════════════
# rules/mag_catalog.smk — Catálogo GLOBAL de MAGs procarióticos
#
# O análogo procariótico do catálogo de vOTU. Mesmo princípio do item "(h)"
# em docs/ROADMAP_SIMPLIFICACAO.md: **computar no representante, herdar no
# membro**. 95% ANI é nível de espécie para procarioto também, então um MAG
# da mesma espécie recuperado em três amostras não precisa ser anotado três
# vezes — e, mais importante, não DEVE: nada garantia que as três recebessem
# a mesma anotação, e o relatório apresentaria as divergências como se
# fossem biologia.
#
# O que MUDA em relação ao desenho anterior:
#   antes  `galah_derep` era POR AMOSTRA e alimentava só o GTDB-Tk. Bins de
#          amostras diferentes nunca se encontravam.
#   agora  um galah único sobre o pool namespaced de todas as amostras e
#          grupos, com provenance.tsv, e as análises rodando apenas sobre as
#          representantes.
#
# O que NÃO muda: binning, CheckM2 e GUNC seguem por amostra. São
# inerentemente por amostra (a cobertura que separa os bins é daquela
# amostra), e o CheckM2 precisa existir ANTES da desreplicação, porque é o
# critério com que o galah escolhe o representante de cada cluster.
# ══════════════════════════════════════════════════════════════════════

import os as _os
import sys as _sys

_sys.path.insert(0, SCRIPTS_DIR)
from mag_catalog import (
    build_pool as _mag_build_pool,
    member_map as _mag_member_map,
    resolve_prefixed_id as _mag_resolve_prefixed,
    merge_checkm2 as _mag_merge_checkm2,
    parse_galah_clusters as _mag_parse_clusters,
    read_provenance as _mag_read_provenance,
    representative_view as _mag_representative_view,
)

MAG_CATALOG_DIR = f"{OUTDIR}/mag_catalog"


def _mag_sources():
    """(source_type, source_id, bins_dir, bin_ext) para cada fonte de MAG.

    Binette (por amostra) produz .fa; VAMB (co-assembly) produz .fna. As duas
    extensões são declaradas aqui e o pool normaliza tudo para .fa, para que
    nenhum consumidor a jusante precise conhecer a origem do bin.
    """
    srcs = [("sample", s, f"{OUTDIR}/{s}/bins/binette/final_bins", ".fa")
            for s in SAMPLES]
    if COASSEMBLY_ENABLED and COASSEMBLY_BINNING and not LONG_READS:
        srcs += [("group", g, f"{OUTDIR}/coassembly/{g}/vamb/run/bins", ".fna")
                 for g in GROUPS]
    return srcs


def _mag_checkm2_pairs():
    """(source_id, quality_report) na MESMA enumeração de _mag_sources().

    Pares, não duas listas paralelas: re-chavear a qualidade sobre os IDs
    namespaced exige saber de qual fonte cada relatório veio, e casar duas
    listas por posição é como se produz um catálogo silenciosamente trocado.
    """
    pairs = [(s, f"{OUTDIR}/{s}/bins/checkm2/quality_report.tsv") for s in SAMPLES]
    if COASSEMBLY_ENABLED and COASSEMBLY_BINNING and not LONG_READS:
        pairs += [(g, f"{OUTDIR}/coassembly/{g}/checkm2/quality_report.tsv")
                  for g in GROUPS]
    return pairs


def _mag_binning_done(wildcards):
    """Sentinelas de binning — o pool depende de que todo bin já exista."""
    deps = [f"{OUTDIR}/{s}/bins/binette/done.txt" for s in SAMPLES]
    if COASSEMBLY_ENABLED and COASSEMBLY_BINNING and not LONG_READS:
        deps += [f"{OUTDIR}/coassembly/{g}/vamb/done.txt" for g in GROUPS]
    return deps


def _mag_checkm2_reports(wildcards):
    return [p for _, p in _mag_checkm2_pairs()]


rule mag_catalog_pool:
    """Reúne todo bin num diretório só, com nome namespaced, + provenance.

    O prefixo NÃO é cosmético: o Binette emite `binette_bin1` em toda amostra
    e o VAMB emite números nus (`1`, `136`). Juntar sem prefixar sobrescreveria
    bins de organismos diferentes em silêncio.
    """
    input:
        binning = _mag_binning_done,
    output:
        provenance = f"{MAG_CATALOG_DIR}/provenance.tsv",
        done       = f"{MAG_CATALOG_DIR}/pool_done.txt",
    log:
        f"{OUTDIR}/logs/mag_catalog_pool.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_catalog_pool.tsv"
    params:
        genomes = f"{MAG_CATALOG_DIR}/genomes",
    run:
        stats = _mag_build_pool(_mag_sources(), str(params.genomes),
                                str(output.provenance))
        with open(str(log[0]), "w") as lf:
            lf.write(f"[mag_catalog_pool] genomas: {stats['n_genomes']}\n")
            lf.write(f"[mag_catalog_pool] fontes com bins: {stats['n_sources']}"
                     f" de {len(_mag_sources())}\n")
        if stats["n_genomes"] == 0:
            write_status(str(output.done), "skipped: nenhum MAG em nenhuma fonte")
        else:
            write_status(str(output.done), "ok")


rule mag_catalog_quality:
    """Um quality_report do CheckM2 para o catálogo inteiro, re-chaveado.

    O galah casa `--checkm2-quality-report` pelo NOME do genoma. Como o pool
    renomeia todo bin para o ID namespaced, a coluna `Name` tem de ser
    reescrita para o mesmo ID. Concatenar os relatórios crus deixaria o galah
    sem qualidade para nenhum genoma, e ele escolheria o representante por
    outro critério sem avisar — trocando "o melhor MAG da espécie" por "algum
    MAG da espécie".
    """
    input:
        pool    = rules.mag_catalog_pool.output.done,
        reports = _mag_checkm2_reports,
    output:
        tsv = f"{MAG_CATALOG_DIR}/checkm2_quality_report.tsv",
    log:
        f"{OUTDIR}/logs/mag_catalog_quality.log"
    run:
        n_rows, n_src = _mag_merge_checkm2(_mag_checkm2_pairs(), str(output.tsv))
        with open(str(log[0]), "w") as lf:
            lf.write(f"[mag_catalog_quality] {n_rows} genomas com qualidade, "
                     f"de {n_src} relatorios\n")
            if n_rows == 0:
                lf.write("[mag_catalog_quality] AVISO: nenhuma linha. O galah "
                         "vai escolher representantes sem criterio de "
                         "qualidade.\n")


rule mag_catalog_derep:
    """galah — UM cluster global a MAG_DEREP_ANI sobre todas as fontes.

    Substitui o `galah_derep` por amostra. Aquele desreplicava os bins da
    própria amostra e alimentava só o GTDB-Tk, então dois MAGs da mesma
    espécie vindos de amostras diferentes nunca se encontravam: eram anotados
    duas vezes, e podiam receber anotações divergentes.
    """
    input:
        pool    = rules.mag_catalog_pool.output.done,
        quality = rules.mag_catalog_quality.output.tsv,
    output:
        cluster = f"{MAG_CATALOG_DIR}/galah_clusters.tsv",
        done    = f"{MAG_CATALOG_DIR}/derep_done.txt",
    log:
        f"{OUTDIR}/logs/mag_catalog_derep.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_catalog_derep.tsv"
    conda: "../envs/env_derep.yaml"
    container: CONTAINERS.get("galah")
    threads: THREADS
    params:
        genomes = f"{MAG_CATALOG_DIR}/genomes",
        repdir  = f"{MAG_CATALOG_DIR}/representatives",
        ani     = MAG_DEREP_ANI,
        enabled = MAG_DEREP_ENABLED,
    shell:
        """
        # O galah aborta se o diretorio de representantes existir e nao estiver
        # vazio, entao um resume do Snakemake sempre falharia aqui.
        rm -rf {params.repdir}
        mkdir -p {params.repdir}
        shopt -s nullglob
        BINS=({params.genomes}/*.fa)

        if [ ${{#BINS[@]}} -eq 0 ]; then
            echo "[mag_catalog_derep] nenhum MAG no pool" | tee {log}
            printf "representative\tmember\n" > {output.cluster}
            printf "skipped: nenhum MAG no pool\n" > {output.done}; exit 0
        fi

        if [ "{params.enabled}" != "True" ]; then
            echo "[mag_catalog_derep] desativado via config — todo MAG e seu proprio representante" | tee {log}
            for fa in "${{BINS[@]}}"; do ln -sf "$(readlink -f "$fa")" {params.repdir}/; done
            printf "representative\tmember\n" > {output.cluster}
            for fa in "${{BINS[@]}}"; do printf "%s\t%s\n" "$fa" "$fa" >> {output.cluster}; done
            printf "skipped: desativado via config\n" > {output.done}; exit 0
        fi

        GALAH_EXIT=0
        galah cluster \
            --genome-fasta-files "${{BINS[@]}}" \
            --ani {params.ani} \
            --checkm2-quality-report {input.quality} \
            --output-cluster-definition {output.cluster} \
            --output-representative-fasta-directory {params.repdir} \
            --threads {threads} \
            > {log} 2>&1 || GALAH_EXIT=$?

        # Fallback: sem representantes utilizaveis, espelha o pool inteiro.
        # Isso NAO e desreplicacao -- o status registra a diferenca, senao um
        # galah quebrado passaria por "esta especie tem N representantes".
        if [ -z "$(ls {params.repdir}/*.fa 2>/dev/null)" ]; then
            echo "[mag_catalog_derep] AVISO: galah nao produziu representantes — usando o pool" | tee -a {log}
            for fa in "${{BINS[@]}}"; do ln -sf "$(readlink -f "$fa")" {params.repdir}/; done
            [ -s {output.cluster} ] || printf "representative\tmember\n" > {output.cluster}
        fi

        if [ $GALAH_EXIT -ne 0 ]; then
            printf "failed: galah cluster exit %s (usando os bins originais)\n" "$GALAH_EXIT" > {output.done}
        else
            printf "ok\n" > {output.done}
        fi
        """


rule mag_catalog_membership:
    """Vista: qual amostra contribuiu qual MAG, e para qual representante.

    É isto que cumpre "as análises só nas representantes, mas sabendo de qual
    amostra cada MAG veio". A anotação vive no representante; a leitura por
    amostra é uma junção — exatamente como `viral_taxonomy` virou uma vista
    sobre o catálogo de vOTU em vez de uma computação por amostra.
    """
    input:
        derep      = rules.mag_catalog_derep.output.done,
        cluster    = rules.mag_catalog_derep.output.cluster,
        provenance = rules.mag_catalog_pool.output.provenance,
    output:
        tsv = f"{MAG_CATALOG_DIR}/mag_membership.tsv",
    log:
        f"{OUTDIR}/logs/mag_catalog_membership.log"
    run:
        clusters = _mag_parse_clusters(str(input.cluster))
        prov     = _mag_read_provenance(str(input.provenance))
        view     = _mag_representative_view(clusters, prov)
        with open(str(output.tsv), "w") as out:
            out.write("source_id\toriginal_bin_id\tmember_id\trepresentative_id\n")
            for row in view:
                out.write("\t".join(row) + "\n")
        reps = {r[3] for r in view}
        with open(str(log[0]), "w") as lf:
            lf.write(f"[mag_catalog_membership] {len(view)} MAGs -> "
                     f"{len(reps)} representantes\n")
            if view and len(reps) == len(view):
                lf.write("[mag_catalog_membership] nenhum cluster com mais de "
                         "um membro: ou nao ha especies compartilhadas entre "
                         "amostras, ou a desreplicacao nao rodou (ver "
                         "derep_done.txt).\n")


rule mag_catalog_gtdbtk:
    """GTDB-Tk classify_wf sobre as REPRESENTANTES do catálogo, uma vez.

    Antes rodava por amostra (e um gêmeo por grupo), sobre bins desreplicados
    apenas dentro da própria amostra. Duas amostras com a mesma espécie
    pagavam o classify_wf duas vezes — e podiam receber classificações
    diferentes, porque o MAG representante de cada amostra era um genoma
    diferente. O GTDB-Tk é a etapa mais cara do lado procariótico, então este
    é também o maior ganho de compute do catálogo.
    """
    input:
        derep = rules.mag_catalog_derep.output.done,
    output:
        done    = f"{MAG_CATALOG_DIR}/gtdbtk/done.txt",
        bac_tsv = f"{MAG_CATALOG_DIR}/gtdbtk/classify/gtdbtk.bac120.summary.tsv",
        ar_tsv  = f"{MAG_CATALOG_DIR}/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
    log:
        f"{OUTDIR}/logs/mag_catalog_gtdbtk.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_catalog_gtdbtk.tsv"
    conda: "../envs/env_gtdbtk.yaml"
    container:  CONTAINERS.get("gtdbtk")
    threads: THREADS
    params:
        bins_dir = f"{MAG_CATALOG_DIR}/representatives",
        outdir   = f"{MAG_CATALOG_DIR}/gtdbtk",
    shell:
        """
        mkdir -p {params.outdir}/classify
        N_BINS=$(find {params.bins_dir} -maxdepth 1 -name "*.fa" 2>/dev/null | wc -l)
        if [ "$N_BINS" -eq 0 ]; then
            echo "[mag_catalog_gtdbtk] nenhuma representante — skipping" | tee {log}
            printf "user_genome\tclassification\n" > {output.bac_tsv}
            printf "user_genome\tclassification\n" > {output.ar_tsv}
            printf "skipped: nenhuma representante\n" > {output.done}; exit 0
        fi

        export GTDBTK_DATA_PATH={GTDBTK_DB}
        GTDBTK_EXIT=0
        gtdbtk classify_wf \
            --genome_dir {params.bins_dir} \
            --out_dir    {params.outdir} \
            --cpus       {threads} \
            --extension  fa \
            >> {log} 2>&1 || GTDBTK_EXIT=$?

        mkdir -p {params.outdir}/classify
        [ -f {output.bac_tsv} ] || printf "user_genome\tclassification\n" > {output.bac_tsv}
        [ -f {output.ar_tsv}  ] || printf "user_genome\tclassification\n" > {output.ar_tsv}
        if [ $GTDBTK_EXIT -ne 0 ]; then
            printf "failed: classify_wf exit %s\n" "$GTDBTK_EXIT" > {output.done}
        else
            printf "ok\n" > {output.done}
        fi
        """


def _mag_gtdbtk_view(source_id, membership_path, bac_in, ar_in, bac_out, ar_out,
                     log_path):
    """Escreve a vista por fonte das tabelas GLOBAIS do GTDB-Tk.

    A classificação vive no representante; cada MAG da fonte herda a do seu.
    A coluna `user_genome` sai com o nome ORIGINAL do bin (`binette_bin1`,
    `136`) e não com o ID namespaced, porque é assim que o relatório e o
    `final/` sempre a leram — o catálogo muda quem computa, não o contrato de
    quem lê. Mesmo desenho de `viral_taxonomy` sobre o catálogo de vOTU.
    """
    import csv as _csv

    members = []
    with open(membership_path, newline="") as fh:
        for row in _csv.DictReader(fh, delimiter="\t"):
            if (row.get("source_id") or "").strip() == source_id:
                members.append(((row.get("original_bin_id") or "").strip(),
                                (row.get("representative_id") or "").strip()))

    n_written = 0
    for src, dst in ((bac_in, bac_out), (ar_in, ar_out)):
        header, by_rep = None, {}
        if _os.path.exists(src):
            with open(src, newline="") as fh:
                r = _csv.reader(fh, delimiter="\t")
                header = next(r, None)
                if header and "user_genome" in header:
                    gi = header.index("user_genome")
                    for row in r:
                        if gi < len(row):
                            by_rep[row[gi].strip()] = row
        header = header or ["user_genome", "classification"]
        gi = header.index("user_genome") if "user_genome" in header else 0
        with open(dst, "w", newline="") as out:
            w = _csv.writer(out, delimiter="\t")
            w.writerow(header)
            for bin_name, rep in members:
                row = by_rep.get(rep)
                if row is None:
                    continue          # representante sem classificação: omite
                row = list(row)
                row[gi] = bin_name
                w.writerow(row)
                n_written += 1

    with open(log_path, "w") as lf:
        lf.write(f"[mag_gtdbtk_view] fonte={source_id}: {len(members)} MAGs, "
                 f"{n_written} linhas herdadas do representante\n")
        if members and n_written == 0:
            lf.write("[mag_gtdbtk_view] AVISO: nenhum representante desta "
                     "fonte tem classificacao. Se o GTDB-Tk global rodou, "
                     "isso e ausencia real; se falhou, ver "
                     "mag_catalog/gtdbtk/done.txt.\n")


rule gtdbtk:
    """Vista por amostra do GTDB-Tk global (o nome é mantido de propósito).

    Era a regra que rodava o `classify_wf` por amostra. Virou vista em
    2026-08-19, quando a classificação passou para o catálogo global. O nome
    e os caminhos de saída foram preservados porque `rule phist`,
    `finalize.smk` e o relatório os referenciam — o catálogo muda quem
    computa, não o contrato de quem lê.
    """
    input:
        gtdbtk_done = rules.mag_catalog_gtdbtk.output.done,
        bac         = rules.mag_catalog_gtdbtk.output.bac_tsv,
        ar          = rules.mag_catalog_gtdbtk.output.ar_tsv,
        membership  = rules.mag_catalog_membership.output.tsv,
    output:
        done    = f"{OUTDIR}/{{sample}}/bins/gtdbtk/done.txt",
        bac_tsv = f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.bac120.summary.tsv",
        ar_tsv  = f"{OUTDIR}/{{sample}}/bins/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/gtdbtk.log"
    params:
        source_id = lambda wc: wc.sample,
    run:
        _os.makedirs(_os.path.dirname(str(output.bac_tsv)), exist_ok=True)
        _mag_gtdbtk_view(params.source_id, str(input.membership),
                         str(input.bac), str(input.ar),
                         str(output.bac_tsv), str(output.ar_tsv), str(log[0]))
        write_status(str(output.done), "ok")


if COASSEMBLY_ENABLED and COASSEMBLY_BINNING and not LONG_READS:

    rule gtdbtk_group:
        """Vista por grupo do GTDB-Tk global. Gêmea de `rule gtdbtk`.

        Nome preservado porque `coassembly_phist` e o `ruleorder` contra
        `multisplit_gtdbtk` o referenciam.

        NOTA: a trilha `multisplit` (VAMB sobre todos os grupos de uma vez)
        NÃO entra no catálogo — seus bins não estão no pool. Segue computando
        o próprio GTDB-Tk, o que é coerente: é um experimento de co-binning
        alternativo, não uma fonte de MAGs finais.
        """
        input:
            gtdbtk_done = rules.mag_catalog_gtdbtk.output.done,
            bac         = rules.mag_catalog_gtdbtk.output.bac_tsv,
            ar          = rules.mag_catalog_gtdbtk.output.ar_tsv,
            membership  = rules.mag_catalog_membership.output.tsv,
        output:
            done    = f"{OUTDIR}/coassembly/{{group}}/gtdbtk/done.txt",
            bac_tsv = f"{OUTDIR}/coassembly/{{group}}/gtdbtk/classify/gtdbtk.bac120.summary.tsv",
            ar_tsv  = f"{OUTDIR}/coassembly/{{group}}/gtdbtk/classify/gtdbtk.ar53.summary.tsv",
        log:
            f"{OUTDIR}/coassembly/{{group}}/logs/gtdbtk.log"
        params:
            source_id = lambda wc: wc.group,
        run:
            _os.makedirs(_os.path.dirname(str(output.bac_tsv)), exist_ok=True)
            _mag_gtdbtk_view(params.source_id, str(input.membership),
                             str(input.bac), str(input.ar),
                             str(output.bac_tsv), str(output.ar_tsv), str(log[0]))
            write_status(str(output.done), "ok")


# ══════════════════════════════════════════════════════════════════════
# Análises nas REPRESENTANTES
#
# O catálogo é o caminho único: toda análise que depende só da sequência roda
# uma vez na representante e o bin membro herda o resultado por join no
# `mag_membership.tsv` (mesmo princípio do item (h) do lado viral). Só o
# `defense_amr_enabled` desliga esse bloco, porque sem ele não há consumidor.
# ══════════════════════════════════════════════════════════════════════

def _mag_source_reps(membership_path, source_id):
    """Representantes que ESTA fonte contribuiu -> nome original do bin.

    {representative_id: [bin_name, ...]}. Um representante pode responder por
    mais de um bin da mesma amostra (duas linhagens da mesma espécie no mesmo
    metagenoma), por isso lista e não escalar.
    """
    import csv as _csv
    if not _os.path.exists(membership_path):
        return {}
    with open(membership_path, newline="") as fh:
        return _mag_member_map(list(_csv.DictReader(fh, delimiter="\t")), source_id)


def _mag_view_by_genome(global_tsv, membership_path, source_id, out_tsv,
                        log_path, genome_col="genome"):
    """Vista por fonte de uma tabela global com UMA LINHA POR GENOMA.

    Implementa HERANÇA, não identidade: todo MAG desta fonte recebe a linha do
    SEU representante, e não só os MAGs que por acaso são representantes. Foi
    exatamente essa distinção que o `viral_taxonomy` errou em 2026-08-18 --
    filtrar a tabela global pelo prefixo da fonte recupera apenas o membro que
    É representante, e com 32 amostras o representante da maioria dos clusters
    pertence a outra amostra, deixando a coluna vazia em silêncio.

    A coluna de genoma sai com o nome ORIGINAL do bin, porque é assim que o
    relatório e o `final/` sempre a leram.
    """
    import csv as _csv

    reps = _mag_source_reps(membership_path, source_id)
    header, by_rep = None, {}
    if _os.path.exists(global_tsv):
        with open(global_tsv, newline="") as fh:
            r = _csv.reader(fh, delimiter="\t")
            header = next(r, None)
            if header and genome_col in header:
                gi = header.index(genome_col)
                for row in r:
                    if gi < len(row):
                        by_rep.setdefault(row[gi].strip(), []).append(row)

    header = header or [genome_col]
    gi = header.index(genome_col) if genome_col in header else 0
    n = 0
    _os.makedirs(_os.path.dirname(out_tsv) or ".", exist_ok=True)
    with open(out_tsv, "w", newline="") as out:
        w = _csv.writer(out, delimiter="\t")
        w.writerow(header)
        for rep, bin_names in reps.items():
            for row in by_rep.get(rep, []):
                for bin_name in bin_names:
                    row = list(row)
                    row[gi] = bin_name
                    w.writerow(row)
                    n += 1
    with open(log_path, "a") as lf:
        lf.write(f"[mag_view] {_os.path.basename(out_tsv)} fonte={source_id}: "
                 f"{len(reps)} representantes, {n} linhas herdadas\n")
        if reps and n == 0:
            lf.write("[mag_view] AVISO: nenhum representante desta fonte tem "
                     "linha na tabela global. Ausencia real ou falha da regra "
                     "global -- ver o done.txt dela, nao assuma zero.\n")
    return n


rule mag_catalog_proteins:
    """Prodigal por representante — o hub de proteínas do catálogo.

    Gêmeo global do `prok_bin_proteins`, com o MESMO formato de manifesto
    (`name / mode / fna / faa / gff`), para que as regras de defesa e AMR
    possam ser herdadas com `use rule` sem alterar uma linha do corpo
    delas. `-p single`, igual ao caminho de bins do original: são genomas,
    não metagenomas.
    """
    input:
        derep = rules.mag_catalog_derep.output.done,
    output:
        manifest = f"{MAG_CATALOG_DIR}/proteins/manifest.txt",
        done     = f"{MAG_CATALOG_DIR}/proteins/done.txt",
    log:
        f"{OUTDIR}/logs/mag_catalog_proteins.log"
    benchmark:
        f"{OUTDIR}/benchmarks/mag_catalog_proteins.tsv"
    conda: "../envs/env_viral.yaml"
    container:  CONTAINERS.get("prodigal")
    threads: 1
    params:
        reps_dir = f"{MAG_CATALOG_DIR}/representatives",
        outdir   = f"{MAG_CATALOG_DIR}/proteins",
        enabled  = DEFENSE_AMR_ENABLED,
    run:
        import glob

        _os.makedirs(params.outdir, exist_ok=True)
        rows = []
        with open(str(log[0]), "w") as lf:
            reps = (sorted(glob.glob(_os.path.join(params.reps_dir, "*.fa")))
                    if params.enabled else [])
            if not params.enabled:
                lf.write("[mag_catalog_proteins] defense_amr_enabled=False -- pulando\n")
            lf.write(f"[mag_catalog_proteins] {len(reps)} representantes\n")
            for rep_fa in reps:
                name = _os.path.splitext(_os.path.basename(rep_fa))[0]
                faa = _os.path.join(params.outdir, f"{name}.faa")
                gff = _os.path.join(params.outdir, f"{name}.gff")
                shell("prodigal -i {rep_fa} -a {faa} -f gff -o {gff} "
                      "-p single -q >> {log} 2>&1 || true")
                if _os.path.exists(faa) and _os.path.getsize(faa) > 0:
                    rows.append((name, "bins", rep_fa, faa, gff))
            with open(str(output.manifest), "w") as mf:
                for row in rows:
                    mf.write("\t".join(row) + "\n")
            lf.write(f"[mag_catalog_proteins] {len(rows)} genomas no manifesto\n")
        if reps and not rows:
            write_status(str(output.done),
                         "failed: %d representantes, 0 proteomas" % len(reps))
        else:
            write_status(str(output.done), "ok")


# ══════════════════════════════════════════════════════════════════════
# VISTAS por fonte — o único lugar que escreve em {sample}/bins/... e
# {coassembly}/{group}/bins/...
#
# As ferramentas discordam sobre ONDE fica o genoma na saída, e a vista
# precisa das duas variantes (verificado nas saídas reais de ~/global/results):
#
#   coluna `genome`   DefenseFinder, ABRicate
#   prefixo do ID     AMRFinderPlus (`Protein identifier`), RGI (`ORF_ID`),
#                     DeepARG (`#ARG`), MMseqs2 (`qseqid`), argNorm (herda a
#                     coluna da ferramenta de origem), consenso (`locus`)
#
# O prefixo NÃO pode ser cortado no primeiro `__`: o ID do catálogo já é
# `{source}__{bin}`, então uma proteína sai `S1__binette_bin1__k141_1_5` e o
# corte devolveria `S1` — todo hit de AMR atribuído à AMOSTRA em vez do MAG.
# `resolve_prefixed_id` casa contra os representantes conhecidos. Depois da
# reescrita para o nome original do bin o corte no primeiro `__` volta a
# estar certo (nomes de bin não contêm o separador), que é o que o relatório
# faz em `_split_genome_prefix`.
# ══════════════════════════════════════════════════════════════════════


def _mag_view_by_prefix(global_tsv, membership_path, source_id, out_tsv,
                        log_path, id_col):
    """Vista por fonte de uma tabela global cujo genoma vive no PREFIXO do ID."""
    import csv as _csv

    reps = _mag_source_reps(membership_path, source_id)
    header, by_rep = None, {}
    if _os.path.exists(global_tsv):
        # Descarta o preambulo de comentario ("## ..." do eggNOG-mapper,
        # "# argNorm version: ..." do argNorm) sem tocar em cabecalhos que
        # comecam com "#" sem espaco, como o "#query"/"#ARG".
        with open(global_tsv, newline="") as fh:
            lines = [ln for ln in fh
                     if not (ln.startswith("##") or ln.startswith("# "))]
        r = _csv.reader(lines, delimiter="\t")
        header = next(r, None)
        idx = None
        if header:
            for cand in (id_col, id_col.lstrip("#")):
                if cand in header:
                    idx = header.index(cand)
                    break
        if idx is not None:
            for row in r:
                if idx >= len(row):
                    continue
                rep, rest = _mag_resolve_prefixed(row[idx], reps)
                if rep is not None:
                    by_rep.setdefault(rep, []).append((row, idx, rest))

    header = header or [id_col]
    n = 0
    _os.makedirs(_os.path.dirname(out_tsv) or ".", exist_ok=True)
    with open(out_tsv, "w", newline="") as out:
        w = _csv.writer(out, delimiter="\t")
        w.writerow(header)
        for rep, bin_names in reps.items():
            for row, idx, rest in by_rep.get(rep, []):
                for bin_name in bin_names:
                    row_out = list(row)
                    row_out[idx] = f"{bin_name}{'__'}{rest}"
                    w.writerow(row_out)
                    n += 1
    with open(log_path, "a") as lf:
        lf.write(f"[mag_view] {_os.path.basename(out_tsv)} fonte={source_id}: "
                 f"{len(reps)} representantes, {n} linhas herdadas\n")
        if reps and n == 0:
            lf.write("[mag_view] AVISO: nenhuma linha casou com um "
                     "representante desta fonte. Ausencia real ou falha da "
                     "regra global -- ver o done.txt dela, nao assuma zero.\n")
    return n


def _mag_manifest_view(global_manifest, membership_path, source_id,
                       out_manifest, log_path):
    """Manifesto de proteínas por fonte, apontando para o FASTA do representante.

    O nome do genoma vira o nome ORIGINAL do bin, mas `fna`/`faa`/`gff`
    continuam sendo os do representante — é dele que vêm os IDs de proteína
    nas tabelas de defesa, e é contra esses IDs que
    `compute_defense_islands` (relatório) casa os genes. Apontar para um
    proteoma próprio do membro daria zero ilhas em silêncio.
    """
    reps = _mag_source_reps(membership_path, source_id)
    rows_by_rep = {}
    if _os.path.exists(global_manifest):
        with open(global_manifest) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 5:
                    rows_by_rep[parts[0]] = parts

    n = 0
    _os.makedirs(_os.path.dirname(out_manifest) or ".", exist_ok=True)
    with open(out_manifest, "w") as out:
        for rep, bin_names in reps.items():
            parts = rows_by_rep.get(rep)
            if parts is None:
                continue
            for bin_name in bin_names:
                out.write("\t".join([bin_name] + parts[1:]) + "\n")
                n += 1
    with open(log_path, "a") as lf:
        lf.write(f"[mag_view] manifest fonte={source_id}: {len(reps)} "
                 f"representantes, {n} genomas\n")
    return n


def _mag_status_view(global_done, out_done):
    """Repassa o status da regra GLOBAL para o done.txt da vista.

    Uma vista vazia porque a ferramenta global falhou não pode aparecer como
    zero biológico no relatório — `load_tool_status` já distingue
    `ok`/`skipped:`/`failed:`, e essa distinção só sobrevive se a vista
    propagar o status em vez de escrever "ok" por ter conseguido escrever um
    arquivo vazio.
    """
    status = "ok"
    if _os.path.exists(global_done):
        txt = open(global_done).read().strip()
        if txt:
            status = txt.splitlines()[0]
    write_status(str(out_done), status)


# (tabela global, coluna de ID) por saída da vista. `None` = coluna `genome`.
_MAG_VIEW_PREFIX_COLS = {
    "amrfinder":        "Protein identifier",
    "rgi":              "ORF_ID",
    # `read_id`, NAO `#ARG`: no DeepARG a primeira coluna e o nome do gene
    # ("MDFA"), nao o identificador da proteina. Usar `#ARG` aqui daria uma
    # vista vazia em silencio, porque nenhum nome de gene casa com um
    # representante do catalogo.
    "deeparg":          "read_id",
    "amrfinder_normed": "Protein identifier",
    "deeparg_normed":   "read_id",
    "consensus":        "locus",
    "mmseqs":           "qseqid",
    "eggnog":           "#query",
}


def _mag_write_views(source_id, membership, inp, outp, log_path):
    """Escreve TODAS as vistas de uma fonte. `inp`/`outp` são dicts nomeados."""
    open(log_path, "w").close()

    _mag_manifest_view(str(inp["manifest"]), membership, source_id,
                       str(outp["manifest"]), log_path)
    _mag_status_view(str(inp["proteins_done"]), str(outp["proteins_done"]))

    for key, genome_col in (("df_systems", "genome"), ("df_anti", "genome"),
                            ("vfdb", "genome"), ("plasmidfinder", "genome"),
                            ("bakta", "bin"), ("kegg", "mag")):
        if key in outp:
            _mag_view_by_genome(str(inp[key]), membership, source_id,
                                str(outp[key]), log_path, genome_col=genome_col)

    for key, id_col in _MAG_VIEW_PREFIX_COLS.items():
        if key in outp:
            _mag_view_by_prefix(str(inp[key]), membership, source_id,
                                str(outp[key]), log_path, id_col)

    for key in ("df_done", "amr_done", "rgi_done", "deeparg_done",
                "abricate_done", "argnorm_done", "consensus_done",
                "mmseqs_done", "bakta_done", "eggnog_done", "kegg_done"):
        if key in outp:
            _mag_status_view(str(inp[key]), str(outp[key]))


def _mag_view_io(base):
    """Saídas de uma vista sob `base` ({sample} ou coassembly/{group})."""
    return {
        "manifest":       f"{base}/bins/proteins/manifest.txt",
        "proteins_done":  f"{base}/bins/proteins/done.txt",
        "df_systems":     f"{base}/bins/defensefinder/defensefinder_systems.tsv",
        "df_anti":        f"{base}/bins/defensefinder/antidefensefinder_systems.tsv",
        "df_done":        f"{base}/bins/defensefinder/done.txt",
        "amrfinder":      f"{base}/bins/amrfinderplus/amrfinder_results.tsv",
        "amr_done":       f"{base}/bins/amrfinderplus/done.txt",
        "rgi":            f"{base}/bins/rgi/rgi_results.txt",
        "rgi_done":       f"{base}/bins/rgi/done.txt",
        "deeparg":        f"{base}/bins/deeparg/deeparg_results.mapping.ARG",
        "deeparg_done":   f"{base}/bins/deeparg/done.txt",
        "vfdb":           f"{base}/bins/abricate/vfdb_results.tsv",
        "plasmidfinder":  f"{base}/bins/abricate/plasmidfinder_results.tsv",
        "abricate_done":  f"{base}/bins/abricate/done.txt",
        "amrfinder_normed": f"{base}/bins/argnorm/amrfinderplus_normed.tsv",
        "deeparg_normed":   f"{base}/bins/argnorm/deeparg_normed.tsv",
        "argnorm_done":     f"{base}/bins/argnorm/done.txt",
        "consensus":        f"{base}/bins/amr_consensus/amr_consensus.tsv",
        "consensus_done":   f"{base}/bins/amr_consensus/done.txt",
        "bakta":            f"{base}/annotation/bakta/bakta_summary.tsv",
        "bakta_done":       f"{base}/annotation/bakta/done.txt",
        "eggnog":           f"{base}/annotation/eggnog/eggnog_annotations.tsv",
        "eggnog_done":      f"{base}/annotation/eggnog/done.txt",
        "kegg":             f"{base}/annotation/kegg_decoder/ko_per_mag.tsv",
        "kegg_done":        f"{base}/annotation/kegg_decoder/done.txt",
    }


_MAG_VIEW_GLOBAL = lambda: {
    "manifest":         rules.mag_catalog_proteins.output.manifest,
    "proteins_done":    rules.mag_catalog_proteins.output.done,
    "df_systems":       rules.mag_defensefinder.output.systems,
    "df_anti":          rules.mag_defensefinder.output.antisystems,
    "df_done":          rules.mag_defensefinder.output.done,
    "amrfinder":        rules.mag_amrfinderplus.output.results,
    "amr_done":         rules.mag_amrfinderplus.output.done,
    "rgi":              rules.mag_rgi_card.output.results,
    "rgi_done":         rules.mag_rgi_card.output.done,
    "deeparg":          rules.mag_deeparg.output.results,
    "deeparg_done":     rules.mag_deeparg.output.done,
    "vfdb":             rules.mag_abricate.output.vfdb,
    "plasmidfinder":    rules.mag_abricate.output.plasmidfinder,
    "abricate_done":    rules.mag_abricate.output.done,
    "amrfinder_normed": rules.mag_argnorm_normalize.output.amrfinder_normed,
    "deeparg_normed":   rules.mag_argnorm_normalize.output.deeparg_normed,
    "argnorm_done":     rules.mag_argnorm_normalize.output.done,
    "consensus":        rules.mag_amr_consensus.output.consensus,
    "consensus_done":   rules.mag_amr_consensus.output.done,
    "mmseqs":           rules.mag_mmseqs_taxonomy_prok.output.hits,
    "mmseqs_done":      rules.mag_mmseqs_taxonomy_prok.output.done,
    "bakta":            rules.mag_bakta.output.summary,
    "bakta_done":       rules.mag_bakta.output.done,
    "eggnog":           rules.mag_eggnog_prok.output.annot_tsv,
    "eggnog_done":      rules.mag_eggnog_prok.output.done,
    "kegg":             rules.mag_extract_kegg_kos.output.ko_table,
    "kegg_done":        rules.mag_extract_kegg_kos.output.done,
}
