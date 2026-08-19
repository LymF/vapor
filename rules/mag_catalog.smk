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
