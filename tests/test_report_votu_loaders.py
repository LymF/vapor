import os
import pytest
from report.data_loaders import (
    STATUS_TRACKED_GLOBAL_TOOLS, STATUS_TRACKED_TOOLS,
    load_votu_catalog, load_votu_presence,
    load_tool_status, summarize_tool_status, GLOBAL_STATUS_LABEL,
)


def _make_catalog(tmp_path, presence_rows, n_pool=10):
    d = tmp_path / "votu_catalog"
    d.mkdir()
    with open(d / "vOTU_clusters.tsv", "w") as fh:
        fh.write("votu_id\trepresentative\tmember\n")
        for i in range(1, len(presence_rows) + 1):
            fh.write(f"vOTU_{i:05d}\trep{i}\trep{i}\n")
    with open(d / "provenance.tsv", "w") as fh:
        fh.write("member_id\tsource_type\tsource_id\toriginal_contig_id\n")
        for i in range(n_pool):
            fh.write(f"S1|c{i}\tsample\tS1\tc{i}\n")
    with open(d / "presence_matrix.tsv", "w") as fh:
        fh.write("votu_id\trepresentative\tS1\tS2\n")
        for i, (a, b) in enumerate(presence_rows, start=1):
            fh.write(f"vOTU_{i:05d}\trep{i}\t{a}\t{b}\n")
    return str(tmp_path)


def test_load_votu_catalog_reports_global_richness(tmp_path):
    outdir = _make_catalog(tmp_path, [("both", "absent"), ("absent", "recruited")],
                           n_pool=10)
    cat = load_votu_catalog(outdir)
    assert cat["n_votus"] == 2
    assert cat["n_pool"] == 10
    assert cat["reduction_pct"] == pytest.approx(80.0)


def test_load_votu_presence_counts_both_signals(tmp_path):
    outdir = _make_catalog(tmp_path, [
        ("both", "absent"),
        ("assembled", "recruited"),
        ("absent", "absent"),
    ])
    pres = load_votu_presence(outdir, ["S1", "S2"])
    assert pres["per_sample"]["S1"]["assembled"] == 2   # 'both' + 'assembled'
    assert pres["per_sample"]["S1"]["recruited"] == 1   # 'both'
    assert pres["per_sample"]["S1"]["total"] == 2       # presente por qualquer sinal
    assert pres["per_sample"]["S2"]["recruited"] == 1
    assert pres["per_sample"]["S2"]["total"] == 1


def test_load_votu_catalog_missing_returns_empty(tmp_path):
    cat = load_votu_catalog(str(tmp_path))
    assert cat["n_votus"] == 0
    assert cat["reduction_pct"] == 0.0


def test_load_votu_presence_missing_returns_zeros(tmp_path):
    pres = load_votu_presence(str(tmp_path), ["S1"])
    assert pres["per_sample"]["S1"]["total"] == 0
    assert pres["votus"] == []


def test_load_tool_status_reports_failed_global_rule(tmp_path):
    outdir = tmp_path
    catalog_dir = outdir / "votu_catalog"
    catalog_dir.mkdir()
    (catalog_dir / "done.txt").write_text("failed: empty pool\n")
    # matrices_done.txt intentionally absent -> reported as 'unknown'.

    status = load_tool_status(str(outdir), ["S1", "S2"])

    assert GLOBAL_STATUS_LABEL in status
    assert status[GLOBAL_STATUS_LABEL]["votu_catalog_reps"]["state"] == "failed"
    assert status[GLOBAL_STATUS_LABEL]["votu_catalog_reps"]["reason"] == "empty pool"
    assert status[GLOBAL_STATUS_LABEL]["votu_catalog_matrices"]["state"] == "unknown"

    rows = summarize_tool_status(status)
    global_rows = [r for r in rows if r["sample"] == GLOBAL_STATUS_LABEL]
    # Uma linha por ferramenta global rastreada. Nao fixar o numero: o conjunto
    # cresce a cada regra que migra para o catalogo (bacphlip_votu e
    # eggnog_viral entraram em 1c66c6c e deixaram este teste quebrado ate
    # 18/08). Fixar so as duas que este teste de fato exercita.
    assert len(global_rows) == len(STATUS_TRACKED_GLOBAL_TOOLS)
    tools_reported = {r["tool"] for r in global_rows}
    assert {"votu_catalog_reps", "votu_catalog_matrices"} <= tools_reported


def test_load_tool_status_global_rule_does_not_corrupt_per_sample_counts(tmp_path):
    outdir = tmp_path
    catalog_dir = outdir / "votu_catalog"
    catalog_dir.mkdir()
    (catalog_dir / "done.txt").write_text("failed: empty pool\n")
    (catalog_dir / "matrices_done.txt").write_text("ok\n")

    status = load_tool_status(str(outdir), ["S1", "S2"])

    # Per-sample entries are untouched by the global rules: only the
    # sample-scoped tools appear under each real sample key.
    # Nao fixar o conjunto: ele encolhe quando uma regra migra para o
    # catalogo global (o galah_derep saiu em 2026-08-19, ao virar
    # mag_catalog_derep). Fixar so a propriedade que o teste exercita --
    # as chaves por amostra sao exatamente as de STATUS_TRACKED_TOOLS.
    assert set(status["S1"].keys()) == set(STATUS_TRACKED_TOOLS)
    assert "galah_derep" not in status["S1"], (
        "a desreplicacao virou global; um done.txt por amostra em bins/derep "
        "nao existe mais e apareceria como lacuna permanente no relatorio")
    assert set(status["S2"].keys()) == set(status["S1"].keys())
    assert set(status.keys()) == {"S1", "S2", GLOBAL_STATUS_LABEL}

    rows = summarize_tool_status(status)
    per_sample_rows = [r for r in rows if r["sample"] != GLOBAL_STATUS_LABEL]
    # Nenhum done.txt por amostra existe em disco -> toda ferramenta por
    # amostra fica 'unknown', |STATUS_TRACKED_TOOLS| por amostra, sem ser
    # afetada pela falha separada da regra global. Contagem derivada, nao
    # fixa: o conjunto encolhe a cada regra que migra para um catalogo
    # global.
    assert len(per_sample_rows) == 2 * len(STATUS_TRACKED_TOOLS)
    assert {(r["sample"], r["state"]) for r in per_sample_rows} == {
        ("S1", "unknown"), ("S2", "unknown"),
    }


def test_load_tool_status_rejects_colliding_sample_name(tmp_path):
    with pytest.raises(ValueError):
        load_tool_status(str(tmp_path), ["S1", GLOBAL_STATUS_LABEL])


def test_load_tool_status_reads_coassembly_groups(tmp_path):
    """Group-scoped done.txt files must be read and keyed separately.

    Before this, load_tool_status only iterated over samples, so a tool that
    failed on the co-assembly track appeared neither as a gap nor as an error.
    """
    from report.data_loaders import (
        load_tool_status, tool_failed, GROUP_STATUS_PREFIX,
    )

    outdir = tmp_path / "results"
    grp = outdir / "coassembly" / "G1" / "bins"
    (grp / "amrfinderplus").mkdir(parents=True)
    (grp / "amrfinderplus" / "done.txt").write_text("failed: amrfinder exit 1\n")
    (grp / "rgi").mkdir(parents=True)
    (grp / "rgi" / "done.txt").write_text("ok\n")

    status = load_tool_status(str(outdir), ["S1"], ["G1"])

    key = f"{GROUP_STATUS_PREFIX}G1"
    assert key in status
    assert status[key]["amrfinderplus"]["state"] == "failed"
    assert status[key]["amrfinderplus"]["reason"] == "amrfinder exit 1"
    assert status[key]["rgi"]["state"] == "ok"
    # a real failure must render as a gap, never as a zero
    assert tool_failed(status, key, "amrfinderplus")
    assert not tool_failed(status, key, "rgi")
    # an absent group done.txt stays 'unknown', which also counts as a gap
    assert status[key]["gtdbtk"]["state"] == "unknown"
    assert tool_failed(status, key, "gtdbtk")
    # samples are unaffected
    assert set(status["S1"]) == set(status[key])


def test_load_tool_status_group_key_cannot_collide_with_sample(tmp_path):
    from report.data_loaders import load_tool_status, GROUP_STATUS_PREFIX
    outdir = tmp_path / "results"
    outdir.mkdir()
    colliding = f"{GROUP_STATUS_PREFIX}G1"
    with pytest.raises(ValueError):
        load_tool_status(str(outdir), [colliding], ["G1"])


def test_every_global_done_file_is_status_tracked():
    """
    Toda regra global do catalogo que escreve um done.txt tem de estar em
    STATUS_TRACKED_GLOBAL_TOOLS.

    Este teste existe porque o item "(h)" (18/08/2026) migrou dez regras para
    o catalogo global e nenhuma foi registrada: pharokka, phold, os dois
    genome_map, defensefinder e dbapis viral podiam falhar sem que aparecesse
    uma linha sequer no relatorio -- a aba so ficava vazia, que e exatamente
    a aparencia de "nao havia nada a mostrar". Sem uma checagem automatica, a
    proxima regra migrada repete o esquecimento.
    """
    import re

    here = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    smk = open(os.path.join(here, "rules", "votu_catalog.smk")).read()

    # Apenas as f-strings de saida ancoradas em CATALOG_DIR; referencias de
    # input (`done = rules.votu_prodigal.output.done`) nao sao f-strings e
    # nao casam.
    declared = {
        "votu_catalog/" + m
        for m in re.findall(r'f"\{CATALOG_DIR\}/([^"]*done[^"]*\.txt)"', smk)
    }
    assert declared, "nenhum done.txt encontrado -- regex do teste quebrou"

    untracked = declared - set(STATUS_TRACKED_GLOBAL_TOOLS.values())
    assert not untracked, (
        "regras globais escrevem done.txt sem rastreio no relatorio: "
        + ", ".join(sorted(untracked))
    )


# ── eggNOG-mapper: o preambulo "##" nao pode virar cabecalho ──────────────

def test_emapper_tsv_ignora_o_preambulo_de_comentarios(tmp_path):
    """Com o preambulo lido como cabecalho, COG_category some e o relatorio
    contava toda proteina como 'Function unknown'."""
    from report.data_loaders import load_emapper_tsv
    p = tmp_path / "eggnog_annotations.tsv"
    p.write_text(
        "## Thu Jul  9 10:18:39 2026\n"
        "## emapper-2.1.13\n"
        "#query\tseed_ortholog\tCOG_category\n"
        "S1__binette_bin1__LLOGBO_00001\t1158612.I580_00126\tC\n"
        "## 1 queries scanned\n"
    )
    rows = load_emapper_tsv(str(p))
    assert len(rows) == 1
    assert rows[0]["COG_category"] == "C"
    assert rows[0]["#query"] == "S1__binette_bin1__LLOGBO_00001"


# ── Colapso da taxonomia por vOTU ────────────────────────────────────────

def _votu_table(tmp_path, sample, rows):
    d = tmp_path / sample / "viral" / "votu"
    d.mkdir(parents=True)
    p = d / f"{sample}_vOTU_table.tsv"
    p.write_text("votu_id\trepresentative\tmember\n" +
                 "".join("\t".join(r) + "\n" for r in rows))
    return str(tmp_path)


def test_colapso_agrupa_membros_do_mesmo_votu(tmp_path):
    """A fonte e o {sample}_vOTU_table.tsv. Ate 2026-08-19 lia-se um
    vOTU_clusters.tsv por amostra que ninguem escreve mais: nada colapsava e
    a mesma populacao viral era contada varias vezes."""
    from report.data_loaders import collapse_taxonomy_to_votu
    outdir = _votu_table(tmp_path, "S1", [
        ["vOTU_1", "S1|k141_10", "k141_10"],
        ["vOTU_1", "S1|k141_10", "k141_55"],
    ])
    recs = [
        {"sample": "S1", "Genome": "k141_10", "Family": "Siphoviridae"},
        {"sample": "S1", "Genome": "k141_55", "Family": "Myoviridae"},
    ]
    out = collapse_taxonomy_to_votu(recs, outdir, ["S1"])
    assert len(out) == 1
    assert out[0]["vOTU_members"] == 2


def test_colapso_prefere_a_linha_do_representante(tmp_path):
    """O representante vem NAMESPACED na tabela e o Genome vem nu — comparar
    os dois cruamente nunca casa e a escolha cai sempre no primeiro membro."""
    from report.data_loaders import collapse_taxonomy_to_votu
    outdir = _votu_table(tmp_path, "S1", [
        ["vOTU_1", "S1|k141_10", "k141_55"],
        ["vOTU_1", "S1|k141_10", "k141_10"],
    ])
    recs = [
        {"sample": "S1", "Genome": "k141_55", "Family": "Myoviridae"},
        {"sample": "S1", "Genome": "k141_10", "Family": "Siphoviridae"},
    ]
    out = collapse_taxonomy_to_votu(recs, outdir, ["S1"])
    assert out[0]["Genome"] == "k141_10"


def test_colapso_com_representante_de_outra_amostra(tmp_path):
    """O representante pode ser de outra amostra: agrupa-se por ele mesmo
    assim, e a linha escolhida e a do primeiro membro local."""
    from report.data_loaders import collapse_taxonomy_to_votu
    outdir = _votu_table(tmp_path, "S1", [
        ["vOTU_1", "S2|k141_99", "k141_10"],
        ["vOTU_1", "S2|k141_99", "k141_55"],
    ])
    recs = [
        {"sample": "S1", "Genome": "k141_10", "Family": "Siphoviridae"},
        {"sample": "S1", "Genome": "k141_55", "Family": "Myoviridae"},
    ]
    out = collapse_taxonomy_to_votu(recs, outdir, ["S1"])
    assert len(out) == 1 and out[0]["Genome"] == "k141_10"


# ── Trilha de reads: coluna do sylph -> amostra ──────────────────────────

def _merged_table(tmp_path, cols, rows):
    p = tmp_path / "merged_relative_abundance.tsv"
    p.write_text("clade_name\t" + "\t".join(cols) + "\n" +
                 "".join("\t".join(r) + "\n" for r in rows))
    return str(p)


def test_reads_classify_casa_colunas_paired_end(tmp_path):
    """Em PE o sylph registra o arquivo do -1 ("{sample}_R1.fastq.gz"). Tirar
    so a extensao devolve "{sample}_R1", que nao e amostra nenhuma: a coluna
    era descartada e TODA a trilha de reads saia zerada, com has_data True."""
    from report.data_loaders import load_reads_classify
    p = _merged_table(tmp_path, ["/data/S1_R1.fastq.gz", "/data/S2_1.fastq.gz"],
                      [["r__Duplodnaviria", "12.5", "3.5"]])
    out = load_reads_classify(p, "", ["S1", "S2"])
    assert out["has_data"]
    v = out["viral"][0]
    assert v["S1"] == 12.5 and v["S2"] == 3.5


def test_reads_classify_casa_colunas_single_end(tmp_path):
    from report.data_loaders import load_reads_classify
    p = _merged_table(tmp_path, ["/data/S1.fastq.gz"],
                      [["r__Duplodnaviria", "9.0"]])
    out = load_reads_classify(p, "", ["S1"])
    assert out["viral"][0]["S1"] == 9.0


def test_reads_classify_nao_confunde_amostras_com_prefixo_comum(tmp_path):
    """'S1' prefixa 'S1_extra': a coluna de S1_extra tem de ir para ela, nao
    para S1."""
    from report.data_loaders import load_reads_classify
    p = _merged_table(tmp_path,
                      ["/data/S1_extra_R1.fastq.gz", "/data/S1_R1.fastq.gz"],
                      [["r__Duplodnaviria", "1.0", "2.0"]])
    out = load_reads_classify(p, "", ["S1", "S1_extra"])
    v = out["viral"][0]
    assert v["S1_extra"] == 1.0 and v["S1"] == 2.0


# ── Curva de acumulacao de vOTU do grupo ─────────────────────────────────

def test_acumulacao_conta_votu_de_provirus(tmp_path):
    """A matriz de abundancia e chaveada pelo contig da montagem
    ("k141_10"); o membro do vOTU vem do conjunto aparado pelo CheckV
    ("k141_10_1"). Sem desfazer o sufixo, todo vOTU de provirus some da
    curva -- em silencio."""
    from report.data_loaders import load_votu_accumulation
    g = tmp_path / "coassembly" / "G1"
    (g / "vamb").mkdir(parents=True)
    (g / "viral" / "votu").mkdir(parents=True)
    (g / "vamb" / "abundance.tsv").write_text(
        "contigname\tS1\tS2\n"
        "k141_10\t5.0\t0.0\n"
        "k141_77\t0.0\t4.0\n")
    (g / "viral" / "votu" / "vOTU_clusters.tsv").write_text(
        "votu_id\trepresentative\tmember\n"
        "vOTU_1\tk141_10_1\tk141_10_1\n"
        "vOTU_2\tk141_77\tk141_77\n")
    out = load_votu_accumulation(str(tmp_path), ["G1"], min_depth=1.0, n_perm=5)
    assert out["G1"]["total"] == 2, "o vOTU do provirus tem de aparecer"
    assert out["G1"]["mean"][-1] == 2.0
