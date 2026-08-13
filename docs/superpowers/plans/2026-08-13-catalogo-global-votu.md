# Catálogo global de vOTU — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Substituir o clustering de vOTU por amostra — que nunca funcionou — por um catálogo global único, com presença/abundância por amostra derivada de recrutamento de leituras e anotação viral rodando uma vez sobre os representantes do catálogo.

**Architecture:** Um novo estágio global `{OUTDIR}/votu_catalog/` concatena os conjuntos virais não-redundantes de todas as amostras e grupos com IDs prefixados, roda `skani triangle --sparse` uma vez sobre o pool, agrupa por single-linkage com o critério ICTV (ANI ≥ 95 % e AF ≥ 85 %), e extrai três camadas de representantes. Cada amostra tem suas leituras remapeadas contra o catálogo para produzir matrizes de presença e abundância. Taxonomia, PHIST e anotação perdem o wildcard `{sample}` e passam a consumir os representantes globais.

**Tech Stack:** Snakemake 9.21, Python 3 (stdlib apenas na lógica de clustering), skani 0.3.2, bwa-mem2 / minimap2, CoverM, pytest.

**Spec:** `docs/superpowers/specs/2026-08-13-catalogo-global-votu-design.md`

## Global Constraints

- **Critério de vOTU:** ANI ≥ `votu_ani` (padrão 95.0) **e** `max(af_q, af_r) ≥ votu_af` (padrão 85.0). ICTV / Roux et al. 2019.
- **Colunas do `skani triangle --sparse`** (verificado com skani 0.3.2, ordem fixa): `Ref_file`(0), `Query_file`(1), `ANI`(2), `Align_fraction_ref`(3), `Align_fraction_query`(4), `Ref_name`(5), `Query_name`(6). Nomes de genoma são as colunas **5 e 6**.
- **Prefixo de ID no pool:** `{source_id}|{contig_id}`. Obrigatório — 3.668 IDs colidem entre amostras nos dados de referência.
- **Corte de presença por recrutamento:** ≥ 75 % do comprimento do representante coberto (`votu_presence_min_coverage`, padrão 75.0). Roux et al. 2017.
- **Filtros de mapeamento herdados de `coverm_viral`:** `--min-read-aligned-length 45`, `--contig-end-exclusion 75`.
- **Identidade mínima de leitura no recrutamento** (`votu_recruit_min_identity`): 95 para leituras curtas e PacBio HiFi; **85 para ONT**. O 95 fixo do `coverm_viral` atual rejeitaria quase todas as leituras ONT, cuja taxa de erro por leitura fica acima disso — usar o mesmo valor para todas as tecnologias zeraria a matriz de presença em corridas ONT sem nenhum aviso.
- **Suporte a leituras longas:** toda regra de mapeamento do catálogo é definida nos dois ramos (`if LONG_READS:` / `else:`) produzindo o **mesmo caminho de saída**, como `rules/mapping.smk` já faz com `bwa_mem`/`minimap2_lr`.
- **Status de regra:** toda regra nova grava o desfecho real no `done.txt` (`ok` / `skipped: <motivo>` / `failed: <motivo>`) via `write_status()` (`Snakefile:92`) ou `printf` equivalente em `shell:`. Nunca `touch` vazio.
- **Testes:** `pytest`, arquivos em `tests/`, import do repo raiz já resolvido por `tests/conftest.py`. Rodar com o env conda `snakemake` (`pytest` precisa estar instalado nele — Task 1 cuida disso).
- **Idioma:** docstrings e comentários de código em inglês (padrão do repo); mensagens de commit em português.

### Correção à spec

A spec §2.1 diz que os conjuntos não-redundantes da co-montagem estão em `coassembly/{group}/viral/consensus/*_viral_nonredundant.fasta`. O caminho real é **`coassembly/{group}/viral/checkv/{group}_viral_trimmed.fasta`** (`rules/coassembly.smk:1516`). Este plano usa o caminho real.

---

## Estrutura de arquivos

| arquivo | responsabilidade |
|---|---|
| `scripts/votu_catalog.py` **(criar)** | Lógica pura de pool, parsing de skani e clustering. Sem dependência de Snakemake, testável direto. |
| `rules/votu_catalog.smk` **(criar)** | Regras Snakemake do estágio global: pool, skani, cluster, reps, map, abundance. |
| `tests/test_votu_catalog.py` **(criar)** | Testes unitários de `scripts/votu_catalog.py`. |
| `rules/viral_binning.smk` (modificar) | Remove `skani_votu`, `skani_cluster`, `viral_votu_reps`; repontar `make_votu_table`. |
| `rules/coassembly.smk` (modificar) | Corrigir o parser espelhado nas regras de grupo. |
| `rules/taxonomy.smk`, `rules/host_prediction.smk`, `rules/annotation.smk` (modificar) | Repontar os 8 consumidores para o catálogo. |
| `rules/abundance.smk` (modificar) | Repontar `coverm_viral` / `votu_abundance`. |
| `Snakefile` (modificar) | Novas variáveis de config, `include:` da nova regra, alvos em `_t_viral()`. |
| `config.yaml` (modificar) | Chaves `votu_catalog_enabled`, `votu_presence_min_coverage`. |
| `scripts/report/data_loaders.py`, `components/*.js` (modificar) | Carregadores e visões do catálogo. |

A lógica pesada vai para `scripts/votu_catalog.py` de propósito: o `run:` block de uma regra Snakemake não é testável isoladamente, e o defeito que este trabalho corrige passou despercebido justamente por estar enterrado num `run:` sem teste.

---

## Task 1: Módulo de pool com prefixação de IDs

**Files:**
- Create: `scripts/votu_catalog.py`
- Create: `tests/test_votu_catalog.py`

**Interfaces:**
- Consumes: nada.
- Produces:
  - `build_pool(sources, pool_path, provenance_path) -> dict` onde `sources` é `list[tuple[str, str, str]]` de `(source_type, source_id, fasta_path)` com `source_type` em `{"sample", "group"}`. Retorna `{"n_sequences": int, "n_sources": int, "n_skipped": int}`.
  - `prefixed_id(source_id, contig_id) -> str`.

- [ ] **Step 1: Instalar pytest no env snakemake**

```bash
source /home/alumnos/lmelo/miniforge3/etc/profile.d/conda.sh
conda activate snakemake
mamba install -y -c conda-forge pytest
pytest --version
```

- [ ] **Step 2: Escrever o teste que falha**

Criar `tests/test_votu_catalog.py`:

```python
import os
import pytest
from votu_catalog import build_pool, prefixed_id


def _write_fasta(path, records):
    with open(path, "w") as fh:
        for name, seq in records:
            fh.write(f">{name}\n{seq}\n")


def test_prefixed_id_joins_with_pipe():
    assert prefixed_id("P01_RNG", "MEGAHIT_k141_10006") == "P01_RNG|MEGAHIT_k141_10006"


def test_build_pool_prefixes_and_records_provenance(tmp_path):
    a = tmp_path / "a.fasta"
    b = tmp_path / "b.fasta"
    # Same contig ID in both samples -- the collision the prefix exists to prevent.
    _write_fasta(a, [("MEGAHIT_k141_10006", "ACGT"), ("MEGAHIT_k141_2", "TTTT")])
    _write_fasta(b, [("MEGAHIT_k141_10006", "GGGG")])

    pool = tmp_path / "pool.fasta"
    prov = tmp_path / "provenance.tsv"
    stats = build_pool(
        [("sample", "S1", str(a)), ("sample", "S2", str(b))],
        str(pool), str(prov),
    )

    assert stats["n_sequences"] == 3
    assert stats["n_sources"] == 2

    names = [l[1:].strip() for l in open(pool) if l.startswith(">")]
    assert names == ["S1|MEGAHIT_k141_10006", "S1|MEGAHIT_k141_2", "S2|MEGAHIT_k141_10006"]
    assert len(set(names)) == 3          # no collision survived

    rows = [l.rstrip("\n").split("\t") for l in open(prov)]
    assert rows[0] == ["member_id", "source_type", "source_id", "original_contig_id"]
    assert rows[1] == ["S1|MEGAHIT_k141_10006", "sample", "S1", "MEGAHIT_k141_10006"]
    assert rows[3] == ["S2|MEGAHIT_k141_10006", "sample", "S2", "MEGAHIT_k141_10006"]


def test_build_pool_keeps_sequence_content(tmp_path):
    a = tmp_path / "a.fasta"
    _write_fasta(a, [("c1", "ACGTACGT")])
    pool = tmp_path / "pool.fasta"
    build_pool([("sample", "S1", str(a))], str(pool), str(tmp_path / "p.tsv"))
    assert "ACGTACGT" in open(pool).read()


def test_build_pool_skips_missing_and_empty_sources(tmp_path):
    a = tmp_path / "a.fasta"
    _write_fasta(a, [("c1", "ACGT")])
    empty = tmp_path / "empty.fasta"
    empty.write_text("")

    stats = build_pool(
        [("sample", "S1", str(a)),
         ("sample", "S2", str(empty)),
         ("group", "G1", str(tmp_path / "missing.fasta"))],
        str(tmp_path / "pool.fasta"), str(tmp_path / "p.tsv"),
    )
    assert stats["n_sequences"] == 1
    assert stats["n_sources"] == 1
    assert stats["n_skipped"] == 2


def test_build_pool_uses_first_whitespace_token_as_id(tmp_path):
    a = tmp_path / "a.fasta"
    with open(a, "w") as fh:
        fh.write(">c1 length=500 cov=3.2\nACGT\n")
    pool = tmp_path / "pool.fasta"
    build_pool([("sample", "S1", str(a))], str(pool), str(tmp_path / "p.tsv"))
    names = [l[1:].strip() for l in open(pool) if l.startswith(">")]
    assert names == ["S1|c1"]
```

Adicionar `scripts/` ao path em `tests/conftest.py`:

```python
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
```

- [ ] **Step 3: Rodar o teste e confirmar que falha**

Run: `conda activate snakemake && python -m pytest tests/test_votu_catalog.py -v`
Expected: FAIL com `ModuleNotFoundError: No module named 'votu_catalog'`

- [ ] **Step 4: Implementar `scripts/votu_catalog.py`**

```python
"""votu_catalog.py -- pure logic for the global vOTU catalog.

Kept out of Snakemake `run:` blocks on purpose: the per-sample clustering
this module replaces was broken for its entire lifetime precisely because
it lived in an untested `run:` block. Everything here is importable and
covered by tests/test_votu_catalog.py.
"""
import os


def prefixed_id(source_id, contig_id):
    """Namespace a contig ID by its source.

    Contig IDs are only unique within an assembly -- 'MEGAHIT_k141_10006'
    exists in many samples. Pooling without a prefix silently merges
    unrelated contigs into the same catalog entry.
    """
    return f"{source_id}|{contig_id}"


def build_pool(sources, pool_path, provenance_path):
    """Concatenate viral FASTAs into one pool with namespaced IDs.

    sources: list of (source_type, source_id, fasta_path); source_type is
             'sample' or 'group'.
    Returns {'n_sequences', 'n_sources', 'n_skipped'}.
    Missing or empty sources are skipped and counted, not fatal -- a sample
    with no viral contigs is a real outcome.
    """
    os.makedirs(os.path.dirname(pool_path) or ".", exist_ok=True)
    n_sequences = n_sources = n_skipped = 0

    with open(pool_path, "w") as out, open(provenance_path, "w") as prov:
        prov.write("member_id\tsource_type\tsource_id\toriginal_contig_id\n")
        for source_type, source_id, fasta_path in sources:
            if not fasta_path or not os.path.exists(fasta_path) \
                    or os.path.getsize(fasta_path) == 0:
                n_skipped += 1
                continue
            wrote_any = False
            with open(fasta_path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        contig_id = line[1:].strip().split()[0]
                        member_id = prefixed_id(source_id, contig_id)
                        out.write(f">{member_id}\n")
                        prov.write(
                            f"{member_id}\t{source_type}\t{source_id}\t{contig_id}\n")
                        n_sequences += 1
                        wrote_any = True
                    else:
                        out.write(line)
            if wrote_any:
                n_sources += 1
            else:
                n_skipped += 1

    return {"n_sequences": n_sequences, "n_sources": n_sources,
            "n_skipped": n_skipped}
```

- [ ] **Step 5: Rodar os testes e confirmar que passam**

Run: `conda activate snakemake && python -m pytest tests/test_votu_catalog.py -v`
Expected: PASS, 5 testes.

- [ ] **Step 6: Commit**

```bash
git add scripts/votu_catalog.py tests/test_votu_catalog.py tests/conftest.py
git commit -m "feat(votu): pool global com IDs prefixados por origem"
```

---

## Task 2: Parser do skani esparso

**Files:**
- Modify: `scripts/votu_catalog.py`
- Modify: `tests/test_votu_catalog.py`

**Interfaces:**
- Consumes: `prefixed_id` (Task 1).
- Produces: `parse_skani_sparse(path, ani_min, af_min, valid_ids) -> list[tuple[str, str]]` — arestas `(name_a, name_b)` que passam o critério. `valid_ids` é um `set[str]`; pares com nome fora do conjunto são descartados.

- [ ] **Step 1: Escrever o teste que falha**

Acrescentar a `tests/test_votu_catalog.py`:

```python
from votu_catalog import parse_skani_sparse

SKANI_HEADER = ("Ref_file\tQuery_file\tANI\tAlign_fraction_ref\t"
                "Align_fraction_query\tRef_name\tQuery_name\n")


def _skani_row(ani, af_ref, af_query, ref_name, query_name):
    return (f"/path/pool.fasta\t/path/pool.fasta\t{ani}\t{af_ref}\t{af_query}\t"
            f"{ref_name}\t{query_name}\n")


def test_parse_uses_columns_five_and_six_for_names(tmp_path):
    """Names live in Ref_name/Query_name, NOT in the first two file-path columns."""
    p = tmp_path / "ani.tsv"
    p.write_text(SKANI_HEADER + _skani_row(99.0, 95.0, 95.0, "S1|a", "S1|b"))
    edges = parse_skani_sparse(str(p), 95.0, 85.0, {"S1|a", "S1|b"})
    assert edges == [("S1|a", "S1|b")]


def test_parse_applies_ani_and_af_thresholds(tmp_path):
    p = tmp_path / "ani.tsv"
    p.write_text(
        SKANI_HEADER
        + _skani_row(99.0, 95.0, 95.0, "a", "b")   # passa
        + _skani_row(94.9, 95.0, 95.0, "a", "c")   # ANI baixo
        + _skani_row(99.0, 80.0, 84.9, "a", "d")   # AF baixo nos dois lados
        + _skani_row(99.0, 84.0, 90.0, "a", "e")   # max(AF) passa
    )
    ids = {"a", "b", "c", "d", "e"}
    edges = parse_skani_sparse(str(p), 95.0, 85.0, ids)
    assert sorted(edges) == [("a", "b"), ("a", "e")]


def test_parse_drops_self_comparisons(tmp_path):
    p = tmp_path / "ani.tsv"
    p.write_text(SKANI_HEADER + _skani_row(100.0, 100.0, 100.0, "a", "a"))
    assert parse_skani_sparse(str(p), 95.0, 85.0, {"a"}) == []


def test_parse_drops_unknown_names(tmp_path):
    p = tmp_path / "ani.tsv"
    p.write_text(SKANI_HEADER + _skani_row(99.0, 95.0, 95.0, "a", "ghost"))
    assert parse_skani_sparse(str(p), 95.0, 85.0, {"a"}) == []


def test_parse_rejects_dense_matrix_format(tmp_path):
    """The dense PHYLIP matrix must yield zero edges, not garbage ones.

    This is the exact format the old parser was silently fed.
    """
    p = tmp_path / "dense.tsv"
    p.write_text("3\nS1|a\nS1|b\t0.00\nS1|c\t0.00\t0.00\n")
    edges = parse_skani_sparse(str(p), 95.0, 85.0, {"S1|a", "S1|b", "S1|c"})
    assert edges == []


def test_parse_tolerates_missing_file(tmp_path):
    assert parse_skani_sparse(str(tmp_path / "nope.tsv"), 95.0, 85.0, {"a"}) == []


def test_parse_skips_malformed_rows(tmp_path):
    p = tmp_path / "ani.tsv"
    p.write_text(SKANI_HEADER
                 + "só\tdois\n"
                 + _skani_row("NA", 95.0, 95.0, "a", "b")
                 + _skani_row(99.0, 95.0, 95.0, "a", "c"))
    assert parse_skani_sparse(str(p), 95.0, 85.0, {"a", "b", "c"}) == [("a", "c")]
```

- [ ] **Step 2: Rodar e confirmar falha**

Run: `python -m pytest tests/test_votu_catalog.py -v -k parse`
Expected: FAIL com `ImportError: cannot import name 'parse_skani_sparse'`

- [ ] **Step 3: Implementar**

Acrescentar a `scripts/votu_catalog.py`:

```python
# skani triangle --sparse column layout (verified against skani 0.3.2):
#   0 Ref_file  1 Query_file  2 ANI  3 Align_fraction_ref
#   4 Align_fraction_query  5 Ref_name  6 Query_name
# The genome names are columns 5 and 6 -- columns 0 and 1 are FILE PATHS.
# Reading names from 0/1 (as the removed per-sample rule did) produces a
# parser that silently never matches anything.
_COL_ANI = 2
_COL_AF_REF = 3
_COL_AF_QUERY = 4
_COL_REF_NAME = 5
_COL_QUERY_NAME = 6
_N_SPARSE_COLS = 7


def parse_skani_sparse(path, ani_min, af_min, valid_ids):
    """Read `skani triangle --sparse` output into vOTU edges.

    An edge is kept when ANI >= ani_min AND max(af_ref, af_query) >= af_min
    (ICTV / Roux et al. 2019). Self-comparisons and names absent from
    valid_ids are dropped.

    Returns a list of (name_a, name_b) tuples.
    """
    edges = []
    if not path or not os.path.exists(path):
        return edges

    with open(path) as fh:
        first = fh.readline()
        # A real sparse file starts with the tab-separated header; anything
        # else (e.g. the dense PHYLIP matrix, which opens with a bare count)
        # simply produces no parseable rows below.
        if first.startswith("Ref_file"):
            pass
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < _N_SPARSE_COLS:
                continue
            try:
                ani = float(parts[_COL_ANI])
                af_ref = float(parts[_COL_AF_REF])
                af_query = float(parts[_COL_AF_QUERY])
            except ValueError:
                continue
            a = parts[_COL_REF_NAME]
            b = parts[_COL_QUERY_NAME]
            if a == b:
                continue
            if a not in valid_ids or b not in valid_ids:
                continue
            if ani >= ani_min and max(af_ref, af_query) >= af_min:
                edges.append((a, b))
    return edges
```

- [ ] **Step 4: Rodar e confirmar que passam**

Run: `python -m pytest tests/test_votu_catalog.py -v`
Expected: PASS, 12 testes.

- [ ] **Step 5: Commit**

```bash
git add scripts/votu_catalog.py tests/test_votu_catalog.py
git commit -m "feat(votu): parser do skani esparso com os indices de coluna corretos"
```

---

## Task 3: Clustering single-linkage com IDs estáveis

**Files:**
- Modify: `scripts/votu_catalog.py`
- Modify: `tests/test_votu_catalog.py`

**Interfaces:**
- Consumes: `parse_skani_sparse` (Task 2).
- Produces:
  - `cluster_votus(ids, edges, completeness) -> list[list[str]]` — componentes conexos, ordenados por tamanho decrescente e depois pelo ID do representante. `completeness` é `dict[str, float]`.
  - `pick_representative(members, completeness) -> str`.
  - `assign_votu_ids(clusters) -> list[tuple[str, str, str]]` — linhas `(votu_id, representative, member)` prontas para escrita.
  - `ClusteringCollapseError` — exceção levantada quando nenhum cluster se forma.

- [ ] **Step 1: Escrever o teste que falha**

Acrescentar a `tests/test_votu_catalog.py`:

```python
from votu_catalog import (
    cluster_votus, pick_representative, assign_votu_ids,
    ClusteringCollapseError, write_clusters,
)


def test_cluster_merges_transitively():
    ids = ["a", "b", "c", "d"]
    edges = [("a", "b"), ("b", "c")]
    clusters = cluster_votus(ids, edges, {})
    assert sorted(len(c) for c in clusters) == [1, 3]
    big = [c for c in clusters if len(c) == 3][0]
    assert set(big) == {"a", "b", "c"}


def test_cluster_keeps_singletons():
    clusters = cluster_votus(["a", "b"], [], {})
    assert sorted(clusters) == [["a"], ["b"]]


def test_pick_representative_prefers_highest_completeness():
    assert pick_representative(["a", "b", "c"], {"a": 40.0, "b": 95.0, "c": 70.0}) == "b"


def test_pick_representative_breaks_ties_by_member_order():
    assert pick_representative(["a", "b"], {"a": 50.0, "b": 50.0}) == "a"


def test_pick_representative_handles_missing_completeness():
    assert pick_representative(["a", "b"], {"b": 10.0}) == "b"


def test_clusters_sorted_by_size_then_representative():
    ids = ["z", "y", "m", "n", "o"]
    edges = [("m", "n"), ("n", "o")]          # cluster de 3
    clusters = cluster_votus(ids, edges, {})
    assert len(clusters[0]) == 3               # maior primeiro
    reps = [pick_representative(c, {}) for c in clusters[1:]]
    assert reps == sorted(reps)                # empates em ordem estavel


def test_assign_votu_ids_is_deterministic():
    ids = ["a", "b", "c", "d"]
    edges = [("a", "b")]
    rows1 = assign_votu_ids(cluster_votus(ids, edges, {}))
    rows2 = assign_votu_ids(cluster_votus(ids, edges, {}))
    assert rows1 == rows2
    votu_ids = sorted({r[0] for r in rows1})
    assert votu_ids == ["vOTU_00001", "vOTU_00002", "vOTU_00003"]


def test_assign_votu_ids_emits_one_row_per_member():
    rows = assign_votu_ids([["a", "b"], ["c"]])
    assert len(rows) == 3
    assert all(len(r) == 3 for r in rows)


def test_write_clusters_raises_when_nothing_collapsed(tmp_path):
    """N sequences -> N clusters is the signature of the historical bug."""
    clusters = [["a"], ["b"], ["c"]]
    with pytest.raises(ClusteringCollapseError, match="no clusters formed"):
        write_clusters(clusters, n_input=3, path=str(tmp_path / "c.tsv"))


def test_write_clusters_allows_single_input_sequence(tmp_path):
    """One genome cannot collapse -- that is not the bug."""
    out = tmp_path / "c.tsv"
    write_clusters([["a"]], n_input=1, path=str(out))
    assert out.read_text().startswith("votu_id\trepresentative\tmember\n")


def test_write_clusters_writes_header_and_rows(tmp_path):
    out = tmp_path / "c.tsv"
    write_clusters([["a", "b"], ["c"]], n_input=3, path=str(out),
                   completeness={"a": 10.0, "b": 90.0})
    rows = [l.rstrip("\n").split("\t") for l in open(out)]
    assert rows[0] == ["votu_id", "representative", "member"]
    assert rows[1] == ["vOTU_00001", "b", "a"]
    assert rows[2] == ["vOTU_00001", "b", "b"]
    assert rows[3] == ["vOTU_00002", "c", "c"]
```

- [ ] **Step 2: Rodar e confirmar falha**

Run: `python -m pytest tests/test_votu_catalog.py -v -k "cluster or representative or votu_ids"`
Expected: FAIL com `ImportError: cannot import name 'cluster_votus'`

- [ ] **Step 3: Implementar**

Acrescentar a `scripts/votu_catalog.py`:

```python
class ClusteringCollapseError(RuntimeError):
    """Raised when clustering produced exactly one cluster per input sequence.

    That outcome is the signature of the defect this module replaces: the
    old rule fed skani's dense matrix to an edge-list parser, so no edge was
    ever built and every contig stayed a singleton for the pipeline's entire
    history. Failing loudly here is the point -- a silent N-in-N-out catalog
    reports assembly redundancy as viral richness.
    """


def pick_representative(members, completeness):
    """Cluster representative: highest CheckV completeness, ties by member order."""
    return max(
        members,
        key=lambda m: (completeness.get(m, 0.0), -members.index(m)),
    )


def cluster_votus(ids, edges, completeness):
    """Single-linkage connected components over the kept edges.

    Returns clusters sorted by size (descending) then by representative ID,
    so that a re-run over the same pool yields the same vOTU labels.
    """
    neigh = {i: set() for i in ids}
    for a, b in edges:
        if a in neigh and b in neigh:
            neigh[a].add(b)
            neigh[b].add(a)

    seen = set()
    clusters = []
    for start in ids:
        if start in seen:
            continue
        comp = []
        stack = [start]
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            comp.append(node)
            stack.extend(n for n in neigh[node] if n not in seen)
        clusters.append(comp)

    clusters.sort(key=lambda c: (-len(c), pick_representative(c, completeness)))
    return clusters


def assign_votu_ids(clusters, completeness=None):
    """Turn clusters into (votu_id, representative, member) rows."""
    completeness = completeness or {}
    rows = []
    for idx, members in enumerate(clusters, start=1):
        votu_id = f"vOTU_{idx:05d}"
        rep = pick_representative(members, completeness)
        for member in members:
            rows.append((votu_id, rep, member))
    return rows


def write_clusters(clusters, n_input, path, completeness=None):
    """Write vOTU_clusters.tsv, refusing to emit a catalog that collapsed nothing."""
    if n_input > 1 and len(clusters) == n_input:
        raise ClusteringCollapseError(
            f"no clusters formed: {n_input} input sequences produced "
            f"{len(clusters)} clusters. This is the signature of a skani "
            f"output-format or parser mismatch -- check that "
            f"`skani triangle --sparse` was used and that genome names are "
            f"read from columns 5 and 6."
        )
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write("votu_id\trepresentative\tmember\n")
        for votu_id, rep, member in assign_votu_ids(clusters, completeness):
            fh.write(f"{votu_id}\t{rep}\t{member}\n")
    return len(clusters)
```

- [ ] **Step 4: Rodar e confirmar que passam**

Run: `python -m pytest tests/test_votu_catalog.py -v`
Expected: PASS, 23 testes.

- [ ] **Step 5: Commit**

```bash
git add scripts/votu_catalog.py tests/test_votu_catalog.py
git commit -m "feat(votu): clustering single-linkage com IDs estaveis e guarda anti-colapso"
```

---

## Task 4: Teste de regressão contra os dados reais

**Files:**
- Create: `tests/test_votu_catalog_regression.py`
- Create: `tests/data/votu_regression/README.md`

**Interfaces:**
- Consumes: `build_pool`, `parse_skani_sparse`, `cluster_votus` (Tasks 1-3).
- Produces: nada consumido por tarefas posteriores.

Este teste ancora o comportamento no baseline medido durante o design: pool de 6 amostras, 9.653 contigs → 5.524 vOTUs (42,8 % de redução). Roda apenas quando os dados de referência existem; caso contrário, é pulado.

- [ ] **Step 1: Documentar como reproduzir o insumo**

Criar `tests/data/votu_regression/README.md`:

```markdown
# Baseline de regressão do catálogo de vOTU

O teste `tests/test_votu_catalog_regression.py` roda sobre um recorte real da
corrida de junho/2026 e é pulado quando os dados não estão presentes.

Para gerar o insumo (requer skani >= 0.3.2 e os resultados em
`~/amazon/results`):

```bash
R=~/amazon/results
for s in P01_RNG_08_947 P01_RNG_08_948 P01_RNG_3_947 P01_RNG_3_948 \
         P06_TAP_3_957 P06_TAP_3_958; do
  awk -v S="$s" '/^>/{print ">"S"|"substr($1,2); next}{print}' \
      "$R/$s/viral/consensus/${s}_viral_nonredundant.fasta"
done > pool.fasta

skani triangle -i pool.fasta -o pool_sparse.tsv -t 8 --slow --sparse
```

Valores medidos em 2026-08-13 com skani 0.3.2:

| métrica | valor |
|---|---|
| sequências no pool | 9653 |
| vOTUs | 5524 |
| redução | 42,8 % |
| vOTUs em >1 amostra | 577 (10,4 %) |
```

- [ ] **Step 2: Escrever o teste**

Criar `tests/test_votu_catalog_regression.py`:

```python
"""Regression test against a real 6-sample slice of the June/2026 run.

Skipped when the reference data is absent. See
tests/data/votu_regression/README.md for how to regenerate it.
"""
import os
import pytest
from votu_catalog import parse_skani_sparse, cluster_votus

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "votu_regression")
POOL = os.path.join(DATA_DIR, "pool.fasta")
ANI = os.path.join(DATA_DIR, "pool_sparse.tsv")

pytestmark = pytest.mark.skipif(
    not (os.path.exists(POOL) and os.path.exists(ANI)),
    reason="reference data absent -- see tests/data/votu_regression/README.md",
)


@pytest.fixture(scope="module")
def clusters():
    ids = [l[1:].split()[0] for l in open(POOL) if l.startswith(">")]
    edges = parse_skani_sparse(ANI, 95.0, 85.0, set(ids))
    return ids, cluster_votus(ids, edges, {})


def test_pool_size_matches_baseline(clusters):
    ids, _ = clusters
    assert len(ids) == 9653


def test_clustering_collapses_to_baseline(clusters):
    _, cl = clusters
    # Exact count is pinned: any drift means the criterion or parser changed.
    assert len(cl) == 5524


def test_reduction_is_substantial(clusters):
    ids, cl = clusters
    reduction = 1 - len(cl) / len(ids)
    assert reduction > 0.40


def test_shared_votus_are_detected(clusters):
    """~10% of vOTUs span more than one sample -- the whole point of pooling."""
    _, cl = clusters
    multi = sum(1 for c in cl if len({m.split("|")[0] for m in c}) > 1)
    assert multi == 577
```

- [ ] **Step 3: Gerar o insumo e rodar**

```bash
source /home/alumnos/lmelo/miniforge3/etc/profile.d/conda.sh
conda activate skani_test
mkdir -p tests/data/votu_regression
cd tests/data/votu_regression
R=~/amazon/results
for s in P01_RNG_08_947 P01_RNG_08_948 P01_RNG_3_947 P01_RNG_3_948 \
         P06_TAP_3_957 P06_TAP_3_958; do
  awk -v S="$s" '/^>/{print ">"S"|"substr($1,2); next}{print}' \
      "$R/$s/viral/consensus/${s}_viral_nonredundant.fasta"
done > pool.fasta
skani triangle -i pool.fasta -o pool_sparse.tsv -t 8 --slow --sparse
cd -
conda activate snakemake
python -m pytest tests/test_votu_catalog_regression.py -v
```

Expected: PASS, 4 testes.

- [ ] **Step 4: Excluir os dados grandes do versionamento**

Acrescentar a `.gitignore`:

```
tests/data/votu_regression/*.fasta
tests/data/votu_regression/*.tsv
```

- [ ] **Step 5: Commit**

```bash
git add tests/test_votu_catalog_regression.py tests/data/votu_regression/README.md .gitignore
git commit -m "test(votu): regressao ancorada no baseline real de 6 amostras"
```

---

## Task 5: Regras do catálogo — pool, skani, cluster, reps

**Files:**
- Create: `rules/votu_catalog.smk`
- Modify: `Snakefile` (variáveis de config, `include:`, alvos)
- Modify: `config.yaml`

**Interfaces:**
- Consumes: `build_pool`, `parse_skani_sparse`, `cluster_votus`, `write_clusters` (Tasks 1-3).
- Produces (caminhos que Tasks 6-9 consomem):
  - `{OUTDIR}/votu_catalog/pool.fasta`
  - `{OUTDIR}/votu_catalog/provenance.tsv`
  - `{OUTDIR}/votu_catalog/skani_ani.tsv`
  - `{OUTDIR}/votu_catalog/vOTU_clusters.tsv`
  - `{OUTDIR}/votu_catalog/catalog_all_reps.fasta`
  - `{OUTDIR}/votu_catalog/catalog_mq_reps.fasta`
  - `{OUTDIR}/votu_catalog/catalog_hq_10kb_reps.fasta`
  - `{OUTDIR}/votu_catalog/done.txt`

- [ ] **Step 1: Acrescentar as chaves de config**

Em `config.yaml`, logo após `votu_af`:

```yaml
votu_catalog_enabled:       true    # catálogo global de vOTU (substitui o clustering por amostra)
votu_presence_min_coverage: 75.0    # % do representante coberto p/ contar presença (Roux 2017)
votu_recruit_min_identity:  null    # % identidade mín. da leitura; null = 95 (SR/HiFi) ou 85 (ONT)
```

Em `Snakefile`, junto de `VOTU_AF` (linha ~231):

```python
VOTU_CATALOG_ENABLED   = config.get("votu_catalog_enabled", True)
VOTU_PRESENCE_MIN_COV  = config.get("votu_presence_min_coverage", 75.0)

# ONT reads carry a per-read error rate well above 5%, so the 95% identity
# filter used for Illumina/HiFi would reject nearly every ONT alignment and
# silently empty the presence matrix. Default per technology; overridable.
_recruit_id = config.get("votu_recruit_min_identity")
if _recruit_id is None:
    _recruit_id = 85 if (LONG_READS and LR_TECH == "ont") else 95
VOTU_RECRUIT_MIN_ID = _recruit_id
```

Esta atribuição precisa vir **depois** de `LONG_READS` e `LR_TECH` (linhas 141-143).

- [ ] **Step 2: Criar `rules/votu_catalog.smk`**

```python
# ══════════════════════════════════════════════════════════════════════
# rules/votu_catalog.smk — BLOCK 7.5: Global vOTU catalog
#
# Replaces the former per-sample skani chain (skani_votu / skani_cluster /
# viral_votu_reps). A vOTU is defined ONCE over the pooled viral sets of
# every sample and co-assembly group, so richness is comparable across
# samples and against the literature. Per-sample presence comes from
# read recruitment against the catalog (see votu_catalog_abundance).
#
#   votu_catalog_pool    — concatenate all viral sets, namespacing IDs
#   votu_catalog_skani   — skani triangle --sparse over the pool
#   votu_catalog_cluster — single-linkage vOTUs (ANI >= 95, AF >= 85)
#   votu_catalog_reps    — all / MQ+ / HQ+>=10kb representative FASTAs
# ══════════════════════════════════════════════════════════════════════

import sys as _sys
_sys.path.insert(0, SCRIPTS_DIR)
from votu_catalog import (
    build_pool, parse_skani_sparse, cluster_votus, write_clusters,
)

CATALOG_DIR = f"{OUTDIR}/votu_catalog"


def _catalog_sources():
    """(source_type, source_id, fasta_path) for every viral set entering the pool."""
    sources = [
        ("sample", s,
         f"{OUTDIR}/{s}/viral/consensus/{s}_viral_nonredundant.fasta")
        for s in SAMPLES
    ]
    if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:
        sources += [
            ("group", g,
             f"{OUTDIR}/coassembly/{g}/viral/checkv/{g}_viral_trimmed.fasta")
            for g in COASSEMBLY_GROUPS
        ]
    return sources


def _catalog_input_fastas(wildcards):
    return [path for _, _, path in _catalog_sources()]


def _catalog_checkv_summaries(wildcards):
    paths = [f"{OUTDIR}/{s}/viral/checkv/quality_summary.tsv" for s in SAMPLES]
    if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:
        paths += [f"{OUTDIR}/coassembly/{g}/viral/checkv/quality_summary.tsv"
                  for g in COASSEMBLY_GROUPS]
    return paths


rule votu_catalog_pool:
    """Concatenate every viral non-redundant set into one namespaced pool.

    Contig IDs are only unique within an assembly, so they are prefixed with
    their source. Without this the pool silently merges unrelated contigs.
    """
    input:
        fastas = _catalog_input_fastas,
    output:
        pool       = f"{CATALOG_DIR}/pool.fasta",
        provenance = f"{CATALOG_DIR}/provenance.tsv",
    log:
        f"{OUTDIR}/logs/votu_catalog_pool.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_pool.tsv"
    run:
        stats = build_pool(_catalog_sources(),
                           str(output.pool), str(output.provenance))
        with open(str(log[0]), "w") as lf:
            lf.write(f"[votu_catalog_pool] sequences: {stats['n_sequences']}\n")
            lf.write(f"[votu_catalog_pool] sources used: {stats['n_sources']}\n")
            lf.write(f"[votu_catalog_pool] sources skipped (missing/empty): "
                     f"{stats['n_skipped']}\n")


rule votu_catalog_skani:
    """skani triangle over the pooled catalog.

    --sparse is REQUIRED, not an optimisation: the default dense matrix
    emits ANI only, with no aligned fraction, so the ICTV AF >= 85 criterion
    cannot be evaluated from it at all.
    """
    input:
        pool = rules.votu_catalog_pool.output.pool,
    output:
        ani = f"{CATALOG_DIR}/skani_ani.tsv",
    log:
        f"{OUTDIR}/logs/votu_catalog_skani.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_skani.tsv"
    conda: "../envs/env_derep.yaml"
    container: CONTAINERS.get("skani")
    threads: THREADS
    shell:
        """
        mkdir -p $(dirname {output.ani})
        N_SEQ=$(grep -c '^>' {input.pool} 2>/dev/null || echo 0)
        if [ "$N_SEQ" -eq 0 ]; then
            echo "[votu_catalog_skani] Empty pool" | tee {log}
            printf "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n" > {output.ani}
            exit 0
        fi
        echo "[votu_catalog_skani] $N_SEQ genomes in pool" | tee {log}
        skani triangle \
            -i {input.pool} \
            -o {output.ani} \
            -t {threads} \
            --slow \
            --sparse \
            >> {log} 2>&1
        """


rule votu_catalog_cluster:
    """Single-linkage vOTU clustering over the global pool.

    Fails loudly when nothing collapses -- N sequences producing N clusters
    is the exact signature of the format/parser mismatch this stage replaces.
    """
    input:
        ani    = rules.votu_catalog_skani.output.ani,
        pool   = rules.votu_catalog_pool.output.pool,
        checkv = _catalog_checkv_summaries,
    output:
        clusters = f"{CATALOG_DIR}/vOTU_clusters.tsv",
    log:
        f"{OUTDIR}/logs/votu_catalog_cluster.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_cluster.tsv"
    params:
        ani_min = VOTU_ANI,
        af_min  = VOTU_AF,
    run:
        import csv

        ids = []
        with open(str(input.pool)) as fh:
            for line in fh:
                if line.startswith(">"):
                    ids.append(line[1:].strip().split()[0])

        # CheckV completeness, re-keyed with the pool's namespaced IDs.
        completeness = {}
        for src, path in zip(_catalog_sources(), list(input.checkv)):
            source_id = src[1]
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    contig = (row.get("contig_id") or "").strip()
                    if not contig:
                        continue
                    try:
                        completeness[f"{source_id}|{contig}"] = float(
                            row.get("completeness") or 0)
                    except ValueError:
                        continue

        edges = parse_skani_sparse(str(input.ani), params.ani_min,
                                   params.af_min, set(ids))
        clusters = cluster_votus(ids, edges, completeness)
        n_clusters = write_clusters(clusters, len(ids), str(output.clusters),
                                    completeness)

        with open(str(log[0]), "w") as lf:
            lf.write(f"[votu_catalog_cluster] genomes={len(ids)} "
                     f"edges={len(edges)} clusters={n_clusters} "
                     f"ani>={params.ani_min} af>={params.af_min}\n")
            if len(ids):
                lf.write(f"[votu_catalog_cluster] reduction="
                         f"{100 * (1 - n_clusters / len(ids)):.1f}%\n")


rule votu_catalog_reps:
    """Extract the three representative tiers used downstream.

    Same quality gates as the removed per-sample viral_votu_reps, applied
    once over the global catalog:
      all      — one per vOTU; recruitment reference and report base
      mq       — MQ+ (Complete/HQ/MQ or completeness >= 50%); taxonomy,
                 PHIST, annotation
      hq_10kb  — HQ+/Complete and >= 10 kb; vConTACT3
    """
    input:
        pool     = rules.votu_catalog_pool.output.pool,
        clusters = rules.votu_catalog_cluster.output.clusters,
        checkv   = _catalog_checkv_summaries,
    output:
        all_fasta     = f"{CATALOG_DIR}/catalog_all_reps.fasta",
        mq_fasta      = f"{CATALOG_DIR}/catalog_mq_reps.fasta",
        hq_10kb_fasta = f"{CATALOG_DIR}/catalog_hq_10kb_reps.fasta",
        done          = f"{CATALOG_DIR}/done.txt",
    log:
        f"{OUTDIR}/logs/votu_catalog_reps.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_reps.tsv"
    params:
        keep_tiers = VIRAL_KEEP_TIERS,
    run:
        import csv

        reps = set()
        with open(str(input.clusters)) as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            for row in rdr:
                reps.add(row["representative"])

        quality = {}
        completeness = {}
        for src, path in zip(_catalog_sources(), list(input.checkv)):
            source_id = src[1]
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    contig = (row.get("contig_id") or "").strip()
                    if not contig:
                        continue
                    key = f"{source_id}|{contig}"
                    quality[key] = (row.get("checkv_quality") or "").strip()
                    try:
                        completeness[key] = float(row.get("completeness") or 0)
                    except ValueError:
                        completeness[key] = 0.0

        def is_mq(rid):
            return (quality.get(rid, "") in params.keep_tiers
                    or completeness.get(rid, 0.0) >= 50.0)

        def is_hq(rid):
            return quality.get(rid, "") in ("Complete", "High-quality")

        seqs = {}
        cur = None
        with open(str(input.pool)) as fh:
            for line in fh:
                if line.startswith(">"):
                    cur = line[1:].strip().split()[0]
                    if cur in reps:
                        seqs[cur] = []
                    else:
                        cur = None
                elif cur is not None:
                    seqs[cur].append(line.strip())

        n_all = n_mq = n_hq = 0
        with open(str(output.all_fasta), "w") as fa, \
             open(str(output.mq_fasta), "w") as fm, \
             open(str(output.hq_10kb_fasta), "w") as fh10:
            for rid, chunks in seqs.items():
                seq = "".join(chunks)
                record = f">{rid}\n{seq}\n"
                fa.write(record)
                n_all += 1
                if is_mq(rid):
                    fm.write(record)
                    n_mq += 1
                if is_hq(rid) and len(seq) >= 10000:
                    fh10.write(record)
                    n_hq += 1

        with open(str(log[0]), "w") as lf:
            lf.write(f"[votu_catalog_reps] all: {n_all}\n")
            lf.write(f"[votu_catalog_reps] MQ+ (taxonomy/PHIST/annotation): {n_mq}\n")
            lf.write(f"[votu_catalog_reps] HQ+/>=10kb (vConTACT3): {n_hq}\n")
        write_status(str(output.done), "ok")
```

- [ ] **Step 3: Registrar o módulo e os alvos no `Snakefile`**

Acrescentar o `include:` junto dos demais (após o de `viral_binning.smk`):

```python
include: "rules/votu_catalog.smk"
```

Em `_t_viral()`, substituir o bloco `if VOTU_CLUSTERING_ENABLED:` (linhas 438-440) por:

```python
    if VOTU_CATALOG_ENABLED:
        t.append(f"{OUTDIR}/votu_catalog/vOTU_clusters.tsv")
        t.append(f"{OUTDIR}/votu_catalog/catalog_all_reps.fasta")
        t.append(f"{OUTDIR}/votu_catalog/done.txt")
```

Confirmar que `VIRAL_KEEP_TIERS`, `COASSEMBLY_ENABLED`, `COASSEMBLY_VIRAL` e `COASSEMBLY_GROUPS` existem como globais antes do `include:`; se algum tiver outro nome no `Snakefile`, usar o nome real em `rules/votu_catalog.smk`.

- [ ] **Step 4: Verificar que o DAG monta**

```bash
source /home/alumnos/lmelo/miniforge3/etc/profile.d/conda.sh
conda activate snakemake
python - <<'PY'
import yaml
c = yaml.safe_load(open('config.yaml'))
c['tracks'] = {'reads': False, 'viral': True, 'prok': True}
yaml.safe_dump(c, open('/tmp/cfg_dag.yaml','w'), sort_keys=False, allow_unicode=True)
PY
snakemake -n --cores 1 --configfile /tmp/cfg_dag.yaml 2>&1 | grep -E "votu_catalog|Error|Exception|total"
```

Expected: as quatro regras `votu_catalog_*` aparecem com count 1 cada; nenhum erro.

- [ ] **Step 5: Commit**

```bash
git add rules/votu_catalog.smk Snakefile config.yaml
git commit -m "feat(votu): estagio global do catalogo (pool, skani, cluster, reps)"
```

---

## Task 6: Recrutamento de leituras e matrizes de presença/abundância

**Files:**
- Modify: `rules/votu_catalog.smk`
- Modify: `Snakefile` (alvos)

**Interfaces:**
- Consumes: `catalog_all_reps.fasta`, `vOTU_clusters.tsv`, `provenance.tsv` (Task 5).
- Produces:
  - `{OUTDIR}/votu_catalog/mapping/{sample}.catalog.sorted.bam`
  - `{OUTDIR}/votu_catalog/coverm/{sample}.tsv`
  - `{OUTDIR}/votu_catalog/presence_matrix.tsv` — colunas `votu_id` + uma por amostra, valores `absent` / `assembled` / `recruited` / `both`
  - `{OUTDIR}/votu_catalog/votu_abundance_matrix.tsv` — colunas `votu_id` + uma por amostra, valor da métrica `COVERM_METHOD`

- [ ] **Step 1: Acrescentar as regras de mapeamento e abundância**

Ao final de `rules/votu_catalog.smk`:

```python
> **Padrão de dois ramos.** As duas variantes abaixo produzem exatamente o mesmo
> caminho de saída (`{CATALOG_DIR}/mapping/{sample}.catalog.sorted.bam`), então
> `votu_catalog_coverm` e `votu_catalog_matrices` funcionam sem saber qual
> tecnologia foi usada — mesma estratégia de `rules/mapping.smk` com
> `bwa_mem` / `minimap2_lr`.

```python
if LONG_READS:

    rule votu_catalog_map:
        """Long-read recruitment of one sample against the whole catalog.

        minimap2 indexes the reference on the fly, so there is no separate
        index rule on this branch. ONT: -ax map-ont; HiFi: -ax map-hifi.
        """
        input:
            reads = _clean_lr,
            fasta = rules.votu_catalog_reps.output.all_fasta,
        output:
            bam = f"{CATALOG_DIR}/mapping/{{sample}}.catalog.sorted.bam",
            bai = f"{CATALOG_DIR}/mapping/{{sample}}.catalog.sorted.bam.bai",
        log:
            f"{OUTDIR}/{{sample}}/logs/votu_catalog_map.log"
        benchmark:
            f"{OUTDIR}/{{sample}}/benchmarks/votu_catalog_map.tsv"
        conda: "../envs/env_mapping.yaml"
        container: CONTAINERS.get("minimap2")
        threads: THREADS
        shell:
            """
            mkdir -p $(dirname {output.bam})
            if [ ! -s {input.fasta} ]; then
                echo "[votu_catalog_map] Empty catalog -- writing empty BAM" | tee {log}
                printf "@HD\tVN:1.6\tSO:coordinate\n" \
                    | samtools view -bS -o {output.bam} - 2>> {log}
                samtools index {output.bam} 2>> {log}
                exit 0
            fi
            if [ "{LR_TECH}" = "hifi" ]; then
                PRESET="map-hifi"
            else
                PRESET="map-ont"
            fi
            minimap2 -ax $PRESET \
                -t {threads} \
                {input.fasta} {input.reads} 2> {log} \
                | samtools sort -@ {threads} -o {output.bam} - 2>> {log}
            samtools index {output.bam} 2>> {log}
            """

else:

    rule votu_catalog_index:
        """Index the catalog once; every sample maps against the same reference."""
        input:
            fasta = rules.votu_catalog_reps.output.all_fasta,
        output:
            idx = f"{CATALOG_DIR}/mapping/catalog_index.bwt.2bit.64",
        log:
            f"{OUTDIR}/logs/votu_catalog_index.log"
        conda: "../envs/env_mapping.yaml"
        container: CONTAINERS.get("bwa_mem2")
        params:
            prefix = f"{CATALOG_DIR}/mapping/catalog_index",
        shell:
            """
            mkdir -p $(dirname {output.idx})
            bwa-mem2 index -p {params.prefix} {input.fasta} > {log} 2>&1
            """

    rule votu_catalog_map:
        """Competitive short-read recruitment of one sample against the catalog.

        This is what makes per-sample presence assembly-independent: a virus
        present but too low-coverage to assemble in a sample is still detected
        here, which assembly-derived presence cannot do.
        """
        input:
            tr1   = _clean_r1,
            tr2   = _clean_r2,
            idx   = rules.votu_catalog_index.output.idx,
            fasta = rules.votu_catalog_reps.output.all_fasta,
        output:
            bam = f"{CATALOG_DIR}/mapping/{{sample}}.catalog.sorted.bam",
            bai = f"{CATALOG_DIR}/mapping/{{sample}}.catalog.sorted.bam.bai",
        log:
            f"{OUTDIR}/{{sample}}/logs/votu_catalog_map.log"
        benchmark:
            f"{OUTDIR}/{{sample}}/benchmarks/votu_catalog_map.tsv"
        conda: "../envs/env_mapping.yaml"
        container: CONTAINERS.get("bwa_mem2")
        threads: THREADS
        params:
            prefix     = f"{CATALOG_DIR}/mapping/catalog_index",
            single_end = SINGLE_END,
        shell:
            """
            mkdir -p $(dirname {output.bam})
            if [ ! -s {input.fasta} ]; then
                echo "[votu_catalog_map] Empty catalog -- writing empty BAM" | tee {log}
                printf "@HD\tVN:1.6\tSO:coordinate\n" \
                    | samtools view -bS -o {output.bam} - 2>> {log}
                samtools index {output.bam} 2>> {log}
                exit 0
            fi
            if [ "{params.single_end}" = "True" ]; then
                bwa-mem2 mem -t {threads} {params.prefix} {input.tr1} 2> {log} \
                    | samtools sort -@ {threads} -o {output.bam} - 2>> {log}
            else
                bwa-mem2 mem -t {threads} {params.prefix} {input.tr1} {input.tr2} 2> {log} \
                    | samtools sort -@ {threads} -o {output.bam} - 2>> {log}
            fi
            samtools index {output.bam} 2>> {log}
            """
```

As demais regras desta tarefa ficam **fora** do `if/else`, no nível do módulo:


rule votu_catalog_coverm:
    """CoverM over the catalog BAM, with the same filters as coverm_viral."""
    input:
        bam   = rules.votu_catalog_map.output.bam,
        fasta = rules.votu_catalog_reps.output.all_fasta,
    output:
        tsv = f"{CATALOG_DIR}/coverm/{{sample}}.tsv",
    log:
        f"{OUTDIR}/{{sample}}/logs/votu_catalog_coverm.log"
    benchmark:
        f"{OUTDIR}/{{sample}}/benchmarks/votu_catalog_coverm.tsv"
    conda: "../envs/env_coverm.yaml"
    container: CONTAINERS.get("coverm")
    threads: THREADS
    params:
        method   = COVERM_METHOD,
        min_id   = VOTU_RECRUIT_MIN_ID,
        long_reads = LONG_READS,
    shell:
        """
        mkdir -p $(dirname {output.tsv})
        if [ ! -s {input.fasta} ]; then
            echo "[votu_catalog_coverm] Empty catalog" | tee {log}
            printf "Contig\t{params.method}\tcovered_fraction\n" > {output.tsv}
            exit 0
        fi
        # --min-read-aligned-length is a short-read filter: a 45 bp floor is
        # meaningless for reads averaging kilobases, and on ONT it would only
        # discard genuinely short alignments that carry real signal.
        if [ "{params.long_reads}" = "True" ]; then
            LEN_FILTER=""
        else
            LEN_FILTER="--min-read-aligned-length 45"
        fi
        echo "[votu_catalog_coverm] min identity: {params.min_id}%" | tee {log}
        coverm contig \
            --bam-files {input.bam} \
            --min-read-percent-identity {params.min_id} \
            $LEN_FILTER \
            --contig-end-exclusion 75 \
            --methods {params.method} covered_fraction \
            --threads {threads} \
            --output-file {output.tsv} \
            >> {log} 2>&1
        """


rule votu_catalog_matrices:
    """Build the vOTU x sample presence and abundance matrices.

    Two independent presence signals, reported side by side and never merged:
      assembled — the vOTU has a member contig from that sample (provenance)
      recruited — >= votu_presence_min_coverage of the representative is
                  covered by that sample's reads (Roux et al. 2017)
    """
    input:
        clusters   = rules.votu_catalog_cluster.output.clusters,
        provenance = rules.votu_catalog_pool.output.provenance,
        coverm     = expand(f"{CATALOG_DIR}/coverm/{{sample}}.tsv", sample=SAMPLES),
    output:
        presence  = f"{CATALOG_DIR}/presence_matrix.tsv",
        abundance = f"{CATALOG_DIR}/votu_abundance_matrix.tsv",
        done      = f"{CATALOG_DIR}/matrices_done.txt",
    log:
        f"{OUTDIR}/logs/votu_catalog_matrices.log"
    benchmark:
        f"{OUTDIR}/benchmarks/votu_catalog_matrices.tsv"
    params:
        min_cov = VOTU_PRESENCE_MIN_COV,
        method  = COVERM_METHOD,
    run:
        import csv
        from collections import defaultdict

        member_to_votu = {}
        votu_rep = {}
        with open(str(input.clusters)) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                member_to_votu[row["member"]] = row["votu_id"]
                votu_rep[row["votu_id"]] = row["representative"]

        # assembled presence: straight from provenance x clusters, free.
        assembled = defaultdict(set)
        with open(str(input.provenance)) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                votu = member_to_votu.get(row["member_id"])
                if votu and row["source_type"] == "sample":
                    assembled[votu].add(row["source_id"])

        rep_to_votu = {rep: votu for votu, rep in votu_rep.items()}
        recruited = defaultdict(set)
        abundance = defaultdict(dict)
        for sample, path in zip(SAMPLES, list(input.coverm)):
            if not os.path.exists(path):
                continue
            with open(path) as fh:
                rdr = csv.reader(fh, delimiter="\t")
                header = next(rdr, None)
                if not header:
                    continue
                # CoverM names columns "<bam> <method>"; positions are stable:
                # 0 = Contig, 1 = configured method, 2 = covered_fraction.
                for parts in rdr:
                    if len(parts) < 3:
                        continue
                    votu = rep_to_votu.get(parts[0])
                    if not votu:
                        continue
                    try:
                        value = float(parts[1])
                        covered = float(parts[2])
                    except ValueError:
                        continue
                    # CoverM reports covered_fraction in [0,1]; the config
                    # threshold is a percentage.
                    if covered * 100.0 >= params.min_cov:
                        recruited[votu].add(sample)
                        abundance[votu][sample] = value

        votus = sorted(votu_rep)
        with open(str(output.presence), "w") as fh:
            fh.write("votu_id\trepresentative\t" + "\t".join(SAMPLES) + "\n")
            for votu in votus:
                cells = []
                for s in SAMPLES:
                    a = s in assembled[votu]
                    r = s in recruited[votu]
                    cells.append("both" if (a and r) else
                                 "assembled" if a else
                                 "recruited" if r else "absent")
                fh.write(f"{votu}\t{votu_rep[votu]}\t" + "\t".join(cells) + "\n")

        with open(str(output.abundance), "w") as fh:
            fh.write("votu_id\t" + "\t".join(SAMPLES) + "\n")
            for votu in votus:
                vals = [f"{abundance[votu].get(s, 0.0):.6g}" for s in SAMPLES]
                fh.write(f"{votu}\t" + "\t".join(vals) + "\n")

        with open(str(log[0]), "w") as lf:
            lf.write(f"[votu_catalog_matrices] vOTUs: {len(votus)}\n")
            lf.write(f"[votu_catalog_matrices] metric: {params.method}\n")
            lf.write(f"[votu_catalog_matrices] presence cutoff: "
                     f"{params.min_cov}% of representative covered\n")
            for s in SAMPLES:
                n_a = sum(1 for v in votus if s in assembled[v])
                n_r = sum(1 for v in votus if s in recruited[v])
                lf.write(f"  {s}: assembled={n_a} recruited={n_r}\n")
        write_status(str(output.done), "ok")
```

- [ ] **Step 2: Acrescentar os alvos**

Em `_t_viral()`, dentro do bloco `if VOTU_CATALOG_ENABLED:`:

```python
        t.append(f"{OUTDIR}/votu_catalog/presence_matrix.tsv")
        t.append(f"{OUTDIR}/votu_catalog/votu_abundance_matrix.tsv")
        t.append(f"{OUTDIR}/votu_catalog/matrices_done.txt")
```

- [ ] **Step 3: Verificar o DAG no modo leituras curtas**

```bash
conda activate snakemake
snakemake -n --cores 1 --configfile /tmp/cfg_dag.yaml 2>&1 | grep -E "votu_catalog_(map|coverm|matrices|index)|Error|Exception"
```

Expected: `votu_catalog_map` e `votu_catalog_coverm` com count = nº de amostras; `votu_catalog_index` e `votu_catalog_matrices` com count 1. Nenhum erro.

- [ ] **Step 4: Verificar o DAG no modo leituras longas**

O ramo ONT/HiFi é código diferente e precisa da sua própria verificação — um `if LONG_READS:` que nunca é exercitado é exatamente como um defeito silencioso entra.

```bash
conda activate snakemake
python - <<'PY'
import yaml
c = yaml.safe_load(open('config.yaml'))
c['tracks'] = {'reads': False, 'viral': True, 'prok': True}
c['long_reads'] = True
c['lr_tech'] = 'ont'
yaml.safe_dump(c, open('/tmp/cfg_lr.yaml','w'), sort_keys=False, allow_unicode=True)
PY
snakemake -n --cores 1 --configfile /tmp/cfg_lr.yaml 2>&1 | grep -E "votu_catalog|minimap2|Error|Exception"
```

Expected: `votu_catalog_map` presente com count = nº de amostras; **`votu_catalog_index` ausente** (minimap2 não usa índice separado); nenhum erro.

Se o `find_samples()` não achar amostras long-read no `FASTQ_DIR` atual, o DAG fica vazio — nesse caso, validar apenas que o workflow **parseia** sem erro:

```bash
snakemake --configfile /tmp/cfg_lr.yaml --list-target-rules 2>&1 | grep -E "votu_catalog|Error"
```

Expected: `votu_catalog_map`, `votu_catalog_coverm`, `votu_catalog_matrices` listadas; `votu_catalog_index` não.

- [ ] **Step 5: Confirmar o limiar de identidade por tecnologia**

```bash
conda activate snakemake
python - <<'PY'
import yaml
for tech, expected in [("ont", 85), ("hifi", 95)]:
    c = yaml.safe_load(open('config.yaml'))
    lr = c.get('long_reads', False)
    print(f"{tech}: config long_reads={lr}")
PY
grep -n "VOTU_RECRUIT_MIN_ID\|_recruit_id" Snakefile
```

Expected: a lógica em `Snakefile` resolve para 85 quando `long_reads: true` e `lr_tech: ont`, e 95 nos demais casos.

- [ ] **Step 6: Commit**

```bash
git add rules/votu_catalog.smk Snakefile config.yaml
git commit -m "feat(votu): recrutamento de leituras (SR + ONT/HiFi) e matrizes de presenca"
```

---

## Task 7: Remover a cadeia por amostra e repontar os consumidores

**Files:**
- Modify: `rules/viral_binning.smk` (remove linhas 248-290, 291-392, 393-501; ajusta `make_votu_table` em 502+)
- Modify: `rules/taxonomy.smk:80,346,438`
- Modify: `rules/host_prediction.smk:23`
- Modify: `rules/annotation.smk:35,461,507`
- Modify: `rules/abundance.smk:11,63`
- Modify: `Snakefile` (remove `VOTU_CLUSTERING_ENABLED` se ficar órfã)

**Interfaces:**
- Consumes: as saídas de `votu_catalog_reps` (Task 5).
- Produces: nenhuma interface nova; é uma migração de entradas.

- [ ] **Step 1: Remover as três regras**

Em `rules/viral_binning.smk`, apagar os blocos completos `rule skani_votu:` (linha 248), `rule skani_cluster:` (291) e `rule viral_votu_reps:` (393), preservando `rule make_votu_table:` (502).

- [ ] **Step 2: Repontar os oito consumidores**

Substituições literais, uma por linha:

| arquivo:linha | de | para |
|---|---|---|
| `taxonomy.smk:80` | `viral = rules.viral_votu_reps.output.mq_fasta,` | `viral = rules.votu_catalog_reps.output.mq_fasta,` |
| `taxonomy.smk:346` | `viral = rules.viral_votu_reps.output.hq_10kb_fasta,` | `viral = rules.votu_catalog_reps.output.hq_10kb_fasta,` |
| `taxonomy.smk:438` | `viral           = rules.viral_votu_reps.output.mq_fasta,` | `viral           = rules.votu_catalog_reps.output.mq_fasta,` |
| `host_prediction.smk:23` | `viral  = rules.viral_votu_reps.output.mq_fasta,` | `viral  = rules.votu_catalog_reps.output.mq_fasta,` |
| `annotation.smk:35` | `viral_nr = rules.viral_votu_reps.output.mq_fasta,` | `viral_nr = rules.votu_catalog_reps.output.mq_fasta,` |
| `annotation.smk:461` | `viral_nr     = rules.viral_votu_reps.output.mq_fasta,` | `viral_nr     = rules.votu_catalog_reps.output.mq_fasta,` |
| `annotation.smk:507` | `viral_nr   = rules.viral_votu_reps.output.mq_fasta,` | `viral_nr   = rules.votu_catalog_reps.output.mq_fasta,` |
| `viral_binning.smk:522` | `votu_reps     = rules.viral_votu_reps.output.all_fasta,` | `votu_reps     = rules.votu_catalog_reps.output.all_fasta,` |

- [ ] **Step 3: Repontar `abundance.smk`**

`rule votu_abundance` (linha 63) tem `clusters = rules.skani_cluster.output.clusters,` — trocar por:

```python
        clusters  = rules.votu_catalog_cluster.output.clusters,
```

Os IDs em `vOTU_clusters.tsv` agora são prefixados (`S1|contig`) enquanto o `coverm_viral` da amostra usa IDs crus. Acrescentar a tradução no `run:` de `votu_abundance`, logo após a leitura do arquivo de clusters:

```python
        # Catalog IDs are namespaced ("<sample>|<contig>"); this rule's CoverM
        # table uses the sample's own bare contig IDs. Strip the prefix of THIS
        # sample's members and ignore members from other samples.
        prefix = f"{wildcards.sample}|"
        member_to_rep = {
            m[len(prefix):]: r[len(prefix):] if r.startswith(prefix) else r
            for m, r in member_to_rep.items()
            if m.startswith(prefix)
        }
```

- [ ] **Step 4: Limpar globais órfãs**

```bash
grep -rn "VOTU_CLUSTERING_ENABLED\|viral_votu_reps\|skani_cluster\|skani_votu" Snakefile rules/*.smk | grep -v coassembly
```

Remover `VOTU_CLUSTERING_ENABLED` do `Snakefile:229` e a chave `use_votu` de `config.yaml` se nenhuma ocorrência sobrar fora de `coassembly.smk`.

- [ ] **Step 5: Verificar o DAG**

```bash
conda activate snakemake
snakemake -n --cores 1 --configfile /tmp/cfg_dag.yaml 2>&1 | tail -30
```

Expected: DAG monta sem erro; `grep -c "viral_votu_reps"` na saída = 0.

- [ ] **Step 6: Commit**

```bash
git add rules/ Snakefile config.yaml
git commit -m "refactor(votu): remove a cadeia por amostra e reponta os 8 consumidores"
```

---

## Task 8: Corrigir o parser espelhado da co-montagem

**Files:**
- Modify: `rules/coassembly.smk:1745-1750` (comando skani), `rules/coassembly.smk:1800-1830` (parser)

**Interfaces:**
- Consumes: `parse_skani_sparse`, `cluster_votus`, `write_clusters` (Tasks 2-3).
- Produces: nenhuma interface nova.

A aba de co-montagem continua exibindo vOTUs por grupo, mas hoje reporta `N contigs = N vOTUs` pelo mesmo defeito. Sem esta tarefa a aba fica errada mesmo com o catálogo correto.

- [ ] **Step 1: Acrescentar `--sparse` ao comando**

Em `rules/coassembly.smk`, no `shell:` de `coassembly_skani_votu` (~linha 1745), acrescentar `--sparse \` após `--slow \`, e trocar o cabeçalho dos caminhos de saída vazia para o cabeçalho esparso:

```
            printf "Ref_file\tQuery_file\tANI\tAlign_fraction_ref\tAlign_fraction_query\tRef_name\tQuery_name\n" > {output.ani}
```

(nas duas guardas: "Disabled via config" e "Empty viral set")

- [ ] **Step 2: Substituir o parser pelo módulo compartilhado**

No `run:` de `coassembly_skani_cluster`, substituir todo o bloco de leitura da matriz e de componentes conexos pelo módulo:

```python
                from votu_catalog import (
                    parse_skani_sparse, cluster_votus, write_clusters,
                )

                edges = parse_skani_sparse(str(input.ani), params.ani_min,
                                           params.af_min, set(ids))
                clusters = cluster_votus(ids, edges, completeness)
                n_clusters = write_clusters(clusters, len(ids),
                                            str(output.clusters), completeness)
                msg = (f"[skani_cluster] genomes={len(ids)} edges={len(edges)} "
                       f"clusters={n_clusters} ani>={params.ani_min} "
                       f"af>={params.af_min}\n")
                _lf.write(msg)
                print(msg, end="")
```

O formato de `vOTU_clusters.tsv` da co-montagem passa a ter três colunas (`votu_id`, `representative`, `member`) em vez de duas. Ajustar `data_loaders.py:1700` (`load_votu_accumulation`), que lê esse arquivo, para usar `csv.DictReader` com as chaves `representative` e `member` — que continuam existindo.

- [ ] **Step 3: Verificar o DAG com co-montagem ligada**

```bash
conda activate snakemake
python - <<'PY'
import yaml
c = yaml.safe_load(open('config.yaml'))
c['tracks'] = {'reads': False, 'viral': True, 'prok': True}
c['coassembly']['enabled'] = True
c['coassembly']['viral'] = True
yaml.safe_dump(c, open('/tmp/cfg_coas.yaml','w'), sort_keys=False, allow_unicode=True)
PY
snakemake -n --cores 1 --configfile /tmp/cfg_coas.yaml 2>&1 | grep -E "coassembly_skani|Error|Exception"
```

Expected: as duas regras aparecem, sem erro.

- [ ] **Step 4: Commit**

```bash
git add rules/coassembly.smk scripts/report/data_loaders.py
git commit -m "fix(coassembly): parser do skani esparso na track de grupo"
```

---

## Task 9: Carregadores do relatório

**Files:**
- Modify: `scripts/report/data_loaders.py`
- Modify: `scripts/report/renderer.py`
- Create: `tests/test_report_votu_loaders.py`

**Interfaces:**
- Consumes: `presence_matrix.tsv`, `votu_abundance_matrix.tsv`, `vOTU_clusters.tsv` (Tasks 5-6).
- Produces:
  - `load_votu_catalog(outdir) -> dict` com `{'n_votus': int, 'n_pool': int, 'reduction_pct': float}`
  - `load_votu_presence(outdir, samples) -> dict` com `{'votus': [...], 'per_sample': {sample: {'assembled': int, 'recruited': int, 'total': int}}}`
  - Constantes JS `VOTU_CATALOG`, `VOTU_PRESENCE`

- [ ] **Step 1: Escrever o teste que falha**

Criar `tests/test_report_votu_loaders.py`:

```python
import os
import pytest
from report.data_loaders import load_votu_catalog, load_votu_presence


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
```

- [ ] **Step 2: Rodar e confirmar falha**

Run: `python -m pytest tests/test_report_votu_loaders.py -v`
Expected: FAIL com `ImportError: cannot import name 'load_votu_catalog'`

- [ ] **Step 3: Implementar os carregadores**

Acrescentar a `scripts/report/data_loaders.py`:

```python
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
```

- [ ] **Step 4: Rodar e confirmar que passam**

Run: `python -m pytest tests/test_report_votu_loaders.py -v`
Expected: PASS, 4 testes.

- [ ] **Step 5: Ligar ao renderer**

Em `scripts/report/renderer.py`, acrescentar aos imports de `.data_loaders`:

```python
    load_votu_catalog, load_votu_presence,
```

Após `tool_status_issues = summarize_tool_status(tool_status)`:

```python
    votu_catalog  = load_votu_catalog(outdir)
    votu_presence = load_votu_presence(outdir, samples)
```

No dicionário de dados, junto de `"TOOL_STATUS_ISSUES"`:

```python
        "VOTU_CATALOG":  votu_catalog,
        "VOTU_PRESENCE": votu_presence,
```

- [ ] **Step 6: Commit**

```bash
git add scripts/report/data_loaders.py scripts/report/renderer.py tests/test_report_votu_loaders.py
git commit -m "feat(report): carregadores do catalogo global de vOTU"
```

---

## Task 10: Visões do relatório — por amostra e total

**Files:**
- Modify: `scripts/report/components/overview.js:46,121,160`
- Modify: `scripts/report/components/viral.js`
- Modify: `scripts/report/components/shell.html`
- Modify: `scripts/report/data_loaders.py:1684-1688` (comentário obsoleto)

**Interfaces:**
- Consumes: `VOTU_CATALOG`, `VOTU_PRESENCE` (Task 9).
- Produces: nenhuma interface nova.

- [ ] **Step 1: Trocar o KPI inflado pela riqueza global**

Em `scripts/report/components/overview.js`, a linha 46 hoje é:

```javascript
      { val: fmt(totalViral),   label: 'Viral vOTUs',       sub: 'consensus' },
```

Substituir por:

```javascript
      { val: fmt(typeof VOTU_CATALOG !== 'undefined' ? VOTU_CATALOG.n_votus : 0),
        label: 'Viral vOTUs', sub: 'catálogo global' },
```

`totalViral` (soma de `viral_consensus` por amostra) deixa de ser riqueza. Onde ela ainda for usada como "contigs virais", renomear o rótulo para `Viral contigs` — nunca `vOTUs`.

Na linha 160 (cartão por amostra), substituir:

```javascript
      { val: fmt(ov.viral_consensus),  label: 'Viral vOTUs' },
```

por:

```javascript
      { val: fmt((VOTU_PRESENCE.per_sample[sample] || {}).total || 0),
        label: 'vOTUs presentes' },
      { val: fmt((VOTU_PRESENCE.per_sample[sample] || {}).assembled || 0),
        label: 'vOTUs montados' },
      { val: fmt(ov.viral_consensus), label: 'Viral contigs' },
```

- [ ] **Step 2: Acrescentar o painel de presença**

Em `scripts/report/components/shell.html`, dentro da seção viral, antes do primeiro `chart-card`:

```html
        <div id="votu-catalog-summary"></div>
        <div class="chart-card"><div id="votu-presence-chart" class="echarts-box"></div></div>
        <div id="votu-presence-table"></div>
```

Em `scripts/report/components/viral.js`, acrescentar antes do `window.renderViral`:

```javascript
  // ── Global vOTU catalog ──────────────────────────────────────────────────
  // Per-sample counts and the total come from the SAME catalog, so they are
  // on one scale. Summing per-sample counts is not the total and never was.
  function buildVotuCatalog() {
    const cat = typeof VOTU_CATALOG !== 'undefined' ? VOTU_CATALOG : null;
    const box = document.getElementById('votu-catalog-summary');
    if (!box || !cat || !cat.n_votus) { if (box) box.innerHTML = ''; return; }
    box.innerHTML =
      `<div class="chart-card"><h3 style="margin-top:0">Catálogo global de vOTU</h3>` +
      `<p><strong>${cat.n_votus.toLocaleString()}</strong> vOTUs a partir de ` +
      `<strong>${cat.n_pool.toLocaleString()}</strong> contigs virais ` +
      `(${cat.reduction_pct}% de redundância removida). Clusterização ICTV: ` +
      `95% ANI + 85% AF, num único passo sobre o conjunto completo.</p></div>`;
  }

  function buildVotuPresence() {
    const pres = typeof VOTU_PRESENCE !== 'undefined' ? VOTU_PRESENCE : null;
    if (!pres) return;
    const samples = typeof SAMPLES !== 'undefined' ? SAMPLES : [];
    const assembled = samples.map(s => (pres.per_sample[s] || {}).assembled || 0);
    const recruited = samples.map(s => (pres.per_sample[s] || {}).recruited || 0);

    mkChart('votu-presence-chart', {
      title: { text: 'vOTUs por amostra — montagem vs recrutamento' },
      tooltip: { trigger: 'axis' },
      legend: { data: ['Montados', 'Recrutados'], top: 28 },
      xAxis: { type: 'category', data: samples, axisLabel: { rotate: 45 } },
      yAxis: { type: 'value', name: 'vOTUs' },
      grid: { bottom: 110, top: 70 },
      series: [
        { name: 'Montados',  type: 'bar', data: assembled },
        { name: 'Recrutados', type: 'bar', data: recruited },
      ],
    });

    // "Em qual amostra está cada vírus" -- capped for page weight.
    const box = document.getElementById('votu-presence-table');
    if (!box) return;
    const shown = pres.votus.slice(0, 500);
    const head = '<th>vOTU</th>' + samples.map(s => `<th>${s}</th>`).join('');
    const rows = shown.map(v => {
      const cells = samples.map(s => {
        const st = v.samples[s];
        const mark = st === 'both' ? '●' : st === 'assembled' ? '◐'
                   : st === 'recruited' ? '○' : '';
        return `<td title="${st}">${mark}</td>`;
      }).join('');
      return `<tr><td>${v.votu_id}</td>${cells}</tr>`;
    }).join('');
    box.innerHTML =
      `<div class="chart-card"><h3 style="margin-top:0">Presença por amostra</h3>` +
      `<p>● montado e recrutado · ◐ só montado · ○ só recrutado. ` +
      `Mostrando ${shown.length} de ${pres.votus.length} vOTUs.</p>` +
      `<div class="table-wrap"><table class="vapor-table"><thead><tr>${head}</tr></thead>` +
      `<tbody>${rows}</tbody></table></div></div>`;
  }
```

Chamar as duas em `window.renderViral`:

```javascript
    buildVotuCatalog();
    buildVotuPresence();
```

- [ ] **Step 3: Corrigir o comentário obsoleto**

Em `scripts/report/data_loaders.py:1684-1688`, o texto diz que a curva de acumulação só é calculável na co-montagem porque os vOTUs por amostra são independentes. Substituir por:

```python
    Computed per co-assembly group. Global accumulation over all samples is
    available from the vOTU catalog (presence_matrix.tsv), where vOTU identity
    is shared across every sample by construction.
```

- [ ] **Step 4: Verificar a sintaxe JS e rodar os testes**

```bash
node --check scripts/report/components/overview.js
node --check scripts/report/components/viral.js
conda activate snakemake && python -m pytest tests/ -v
```

Expected: JS OK; todos os testes passam.

- [ ] **Step 5: Commit**

```bash
git add scripts/report/components/ scripts/report/data_loaders.py
git commit -m "feat(report): riqueza global, vOTUs por amostra e matriz de presenca"
```

---

## Task 11: Documentação e verificação final

**Files:**
- Modify: `CLAUDE.md`
- Modify: `docs/PIPELINE_METHODS.md`
- Modify: `README.md`

**Interfaces:**
- Consumes: tudo.
- Produces: nada.

- [ ] **Step 1: Atualizar a descrição dos módulos**

Em `CLAUDE.md`, na lista de `rules/`, acrescentar após a linha de `viral_binning.smk`:

```markdown
  - `votu_catalog.smk`: Catálogo global de vOTU — pool de todos os conjuntos virais com IDs prefixados por origem, `skani triangle --sparse` num passo único, clustering ICTV (95% ANI + 85% AF), representantes em três camadas, e matrizes de presença/abundância vOTU×amostra por recrutamento de leituras. Substitui o clustering por amostra.
```

Em `docs/PIPELINE_METHODS.md`, substituir qualquer descrição de vOTU por amostra pela definição global, citando Roux et al. 2019 (critério) e Roux et al. 2017 (corte de 75 % de cobertura para presença).

- [ ] **Step 2: Verificação completa do DAG**

Os três modos têm que montar: leituras curtas, co-montagem e leituras longas.

```bash
conda activate snakemake
snakemake -n --cores 1 --configfile /tmp/cfg_dag.yaml  2>&1 | tail -20   # short reads
snakemake -n --cores 1 --configfile /tmp/cfg_coas.yaml 2>&1 | tail -20   # co-assembly
snakemake --configfile /tmp/cfg_lr.yaml --list-target-rules 2>&1 | tail -20  # long reads
grep -rn "viral_votu_reps\|rules.skani_cluster\|rules.skani_votu" rules/ Snakefile | grep -v coassembly_skani
```

Expected: os três montam/parseiam sem erro; o `grep` não retorna nada.

- [ ] **Step 3: Rodar a suíte completa**

```bash
conda activate snakemake && python -m pytest tests/ -v
```

Expected: todos passam (os de regressão são pulados se os dados não estiverem presentes).

- [ ] **Step 4: Commit**

```bash
git add CLAUDE.md docs/PIPELINE_METHODS.md README.md
git commit -m "docs(votu): documenta o catalogo global e a definicao de vOTU"
```

---

## Auto-revisão

**Cobertura da spec:**

| seção da spec | tarefa |
|---|---|
| §1.1 três defeitos do parser | Tasks 2 (índices), 5 (`--sparse`), 3 (guarda anti-colapso) |
| §2.1 pool prefixado + provenance | Task 1, Task 5 |
| §2.2 skani `--sparse` | Task 5 |
| §2.3 clustering, IDs estáveis, validação | Task 3, Task 5 |
| §2.4 três camadas de representantes | Task 5 |
| §3.1 presença por montagem | Task 6 |
| §3.2 presença por recrutamento | Task 6 |
| §3.3 matrizes de saída | Task 6 |
| §4 migração dos 8 consumidores | Task 7 |
| §4.1 regras removidas | Task 7 |
| §4.2 co-montagem preservada e corrigida | Task 8 |
| §5 relatório | Tasks 9, 10 |
| §6 configuração | Task 5 |
| §7 verificação | Tasks 4, 5, 6, 7, 8, 11 |

**Riscos conhecidos, deliberadamente aceitos:**

1. **Ordem das colunas do CoverM** (Task 6) — o código assume `Contig`, método, `covered_fraction` nas posições 0/1/2, como o `coverm_viral` atual já assume. Se a versão do CoverM mudar a ordem, a matriz sai zerada. Quem implementar deve conferir o cabeçalho real do primeiro `coverm/{sample}.tsv` gerado antes de dar a Task 6 por concluída.
2. **`_clean_r1` / `_clean_r2` / `_clean_lr`** (Task 6) são funções do `Snakefile` (linhas 163-180). Se `rules/votu_catalog.smk` for incluído antes delas existirem, o `include:` falha — inserir o `include:` depois de `viral_binning.smk`, como a Task 5 especifica.
3. **Ramo long-read não exercitado por dados reais.** A Task 6 define o ramo minimap2 e a Task 11 valida que ele parseia e entra no DAG, mas o conjunto de referência disponível é Illumina PE — nenhuma execução real ONT/HiFi do catálogo foi feita. O limiar de identidade de 85 % para ONT é uma escolha fundamentada na taxa de erro por leitura da tecnologia, não um valor calibrado contra estes dados. Primeira corrida ONT: conferir no log do `votu_catalog_coverm` que a matriz de presença não saiu vazia antes de confiar nos números.
4. **Co-montagem sem recrutamento próprio.** Os contigs dos grupos entram no catálogo (Task 5), mas o recrutamento é sempre por amostra — grupos não têm leituras próprias, e as leituras de um grupo são as das amostras que o compõem. Um vOTU recuperado apenas na co-montagem aparece no catálogo e ganha presença por recrutamento nas amostras onde suas leituras mapeiam; sua presença "por montagem" fica registrada como do grupo, não das amostras. Isso é intencional e está refletido no `presence_matrix.tsv`, que só tem colunas de amostra.
