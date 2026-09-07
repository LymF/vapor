# Calibração LOO de cutoffs para os rollups de taxonomia mmseqs2 — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Calibrar, via leave-one-out sobre um self-alignment all-vs-all de
cada seqTaxDB, uma fração-de-concordância mínima por rank/taxon, e usá-la
nos 3 rollups de consenso mmseqs2 (INPHARED, custom viral, custom prok) no
lugar da unanimidade estrita atual — ganhando profundidade de classificação
sem reintroduzir "spurious specificity".

**Architecture:** Um módulo Python novo e puro
(`scripts/mmseqs_taxonomy_calibration.py`) calcula os cutoffs a partir de
`mapping.tsv` + `taxdump/` (já produzidos por `prepare_mmseqs_taxdb.py`)
mais um novo `self_align.tsv` (all-vs-all, gerado uma vez via um flag
`--calibrate` nesse mesmo script) e persiste `rank_cutoffs.json` ao lado da
seqTaxDB. Os dois pontos de rollup existentes
(`scripts/votu_catalog.py::mmseqs_lca_rollup`, movido de dentro de
`rules/votu_catalog.smk`, e
`scripts/report/data_loaders.py::load_mmseqs_taxonomy_prok`) passam a
aceitar um `cutoffs: dict | None` opcional e trocam a checagem de
unanimidade por fração-vs-cutoff, caindo em unanimidade quando não há
cutoff calibrado para aquele rank/taxon.

**Tech Stack:** Python 3 puro (sem pandas nesta parte), pytest, mmseqs2
(`search`, `convertalis`) via `subprocess` — mesmo padrão já usado em
`scripts/prepare_mmseqs_taxdb.py`.

**Spec:** `docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md`

## Global Constraints

- Zero regressão quando `rank_cutoffs.json` não existe: comportamento
  idêntico ao atual (unanimidade, `fraction == 1.0`).
- Nenhuma das 3 rules de taxonomia (`votu_mmseqs_taxonomy`,
  `votu_mmseqs_taxonomy_custom`, `mag_mmseqs_taxonomy_prok`) muda de
  comportamento ou de output — só o consumo em `votu_taxonomy` e no report
  muda.
- `rank_cutoffs.json` vive ao lado da própria seqTaxDB (mesmo diretório do
  path passado a `--out` em `prepare_mmseqs_taxdb.py`), nunca em
  `CATALOG_DIR`.
- Funções puras (sem `subprocess`, sem argumento implícito de working
  directory) vivem em `scripts/mmseqs_taxonomy_calibration.py` e
  `scripts/votu_catalog.py`; orquestração (chamadas a `mmseqs`, leitura de
  `snakemake.*`) fica nas rules/scripts que já orquestravam antes.

---

## Task 1: Módulo de calibração — parsing puro (taxdump, mapping, self-align)

**Files:**
- Create: `scripts/mmseqs_taxonomy_calibration.py`
- Test: `tests/test_mmseqs_taxonomy_calibration.py`

**Interfaces:**
- Produces: `load_taxdump(taxdump_dir) -> (parent_of: dict[str,str], rank_of: dict[str,str], name_of: dict[str,str])`
- Produces: `lineage_by_rank(taxid: str, parent_of, rank_of, name_of) -> dict[str,str]`
- Produces: `load_protein_taxid_map(mapping_path: str) -> dict[str,str]`
- Produces: `iter_self_align_rows(self_align_path: str) -> Iterator[tuple[str,str,float]]`

- [ ] **Step 1: Write the failing tests for `load_taxdump` and `lineage_by_rank`**

```python
# tests/test_mmseqs_taxonomy_calibration.py
import os
import pytest
from mmseqs_taxonomy_calibration import load_taxdump, lineage_by_rank


def _write_taxdump(tmp_path, nodes_lines, names_lines):
    d = tmp_path / "taxdump"
    d.mkdir()
    (d / "nodes.dmp").write_text("\n".join(nodes_lines) + "\n")
    (d / "names.dmp").write_text("\n".join(names_lines) + "\n")
    return str(d)


def test_load_taxdump_le_parent_rank_name(tmp_path):
    taxdump_dir = _write_taxdump(
        tmp_path,
        nodes_lines=[
            "1\t|\t1\t|\tno rank\t|\t-\t|",
            "2\t|\t1\t|\trealm\t|\t-\t|",
            "3\t|\t2\t|\tfamily\t|\t-\t|",
        ],
        names_lines=[
            "1\t|\troot\t|\t-\t|\tscientific name\t|",
            "2\t|\tDuplodnaviria\t|\t-\t|\tscientific name\t|",
            "3\t|\tMariniviridae\t|\t-\t|\tscientific name\t|",
        ],
    )
    parent_of, rank_of, name_of = load_taxdump(taxdump_dir)
    assert parent_of["3"] == "2"
    assert rank_of["3"] == "family"
    assert name_of["3"] == "Mariniviridae"


def test_lineage_by_rank_anda_ate_a_raiz(tmp_path):
    taxdump_dir = _write_taxdump(
        tmp_path,
        nodes_lines=[
            "1\t|\t1\t|\tno rank\t|\t-\t|",
            "2\t|\t1\t|\trealm\t|\t-\t|",
            "3\t|\t2\t|\tfamily\t|\t-\t|",
        ],
        names_lines=[
            "1\t|\troot\t|\t-\t|\tscientific name\t|",
            "2\t|\tDuplodnaviria\t|\t-\t|\tscientific name\t|",
            "3\t|\tMariniviridae\t|\t-\t|\tscientific name\t|",
        ],
    )
    parent_of, rank_of, name_of = load_taxdump(taxdump_dir)
    lineage = lineage_by_rank("3", parent_of, rank_of, name_of)
    assert lineage == {"realm": "Duplodnaviria", "family": "Mariniviridae"}


def test_lineage_by_rank_taxid_desconhecido_devolve_vazio(tmp_path):
    taxdump_dir = _write_taxdump(
        tmp_path,
        nodes_lines=["1\t|\t1\t|\tno rank\t|\t-\t|"],
        names_lines=["1\t|\troot\t|\t-\t|\tscientific name\t|"],
    )
    parent_of, rank_of, name_of = load_taxdump(taxdump_dir)
    assert lineage_by_rank("999", parent_of, rank_of, name_of) == {}
```

- [ ] **Step 2: Run tests to verify they fail**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_mmseqs_taxonomy_calibration.py -v`
Expected: FAIL with `ModuleNotFoundError: No module named 'mmseqs_taxonomy_calibration'`

- [ ] **Step 3: Implement `load_taxdump` + `lineage_by_rank`**

```python
# scripts/mmseqs_taxonomy_calibration.py
"""mmseqs_taxonomy_calibration.py -- leave-one-out calibration of a
minimum agreement-fraction per rank/taxon for the mmseqs2 taxonomy
rollups (scripts/votu_catalog.py::mmseqs_lca_rollup and
scripts/report/data_loaders.py::load_mmseqs_taxonomy_prok).

See docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md
for the full design and why this calibrates a FRACTION, not a bitscore
cutoff like VITAP does: the mmseqs taxonomy output the rollups consume
(qseqid, taxid, rank, name, lineage) carries no per-protein score, so a
bitscore cutoff has nowhere to plug in without changing the 3 taxonomy
rules themselves (out of scope here).
"""
import json
import os


def load_taxdump(taxdump_dir):
    """Parse nodes.dmp/names.dmp (format written by
    prepare_mmseqs_taxdb.py::build_taxdump) into
    (parent_of, rank_of, name_of), each {taxid: value}, all as strings."""
    parent_of, rank_of, name_of = {}, {}, {}
    nodes_path = os.path.join(taxdump_dir, 'nodes.dmp')
    names_path = os.path.join(taxdump_dir, 'names.dmp')
    with open(nodes_path) as f:
        for line in f:
            parts = [p.strip() for p in line.split('\t|\t')]
            if len(parts) < 3:
                continue
            taxid, parent, rank = parts[0], parts[1], parts[2]
            parent_of[taxid] = parent
            rank_of[taxid] = rank
    with open(names_path) as f:
        for line in f:
            parts = [p.strip() for p in line.split('\t|\t')]
            if len(parts) < 2:
                continue
            taxid, name = parts[0], parts[1]
            name_of[taxid] = name
    return parent_of, rank_of, name_of


def lineage_by_rank(taxid, parent_of, rank_of, name_of):
    """Walk taxid up to the root (id '1'), returning {rank: name} for
    every ancestor including taxid itself. Unknown taxid -> {}."""
    if taxid not in parent_of:
        return {}
    lineage = {}
    current = taxid
    seen = set()
    while current != '1' and current not in seen:
        seen.add(current)
        rank = rank_of.get(current)
        name = name_of.get(current)
        if rank and rank != 'no rank' and name:
            lineage[rank] = name
        current = parent_of.get(current, '1')
    return lineage
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_mmseqs_taxonomy_calibration.py -v`
Expected: PASS (3 tests)

- [ ] **Step 5: Write failing tests for `load_protein_taxid_map` and `iter_self_align_rows`**

```python
# append to tests/test_mmseqs_taxonomy_calibration.py
from mmseqs_taxonomy_calibration import load_protein_taxid_map, iter_self_align_rows


def test_load_protein_taxid_map(tmp_path):
    mapping = tmp_path / "mapping.tsv"
    mapping.write_text("prot_1\t3\nprot_2\t5\n")
    result = load_protein_taxid_map(str(mapping))
    assert result == {"prot_1": "3", "prot_2": "5"}


def test_iter_self_align_rows_parseia_query_target_score(tmp_path):
    align = tmp_path / "self_align.tsv"
    align.write_text("prot_1\tprot_2\t150.3\nprot_2\tprot_1\t150.3\n")
    rows = list(iter_self_align_rows(str(align)))
    assert rows == [("prot_1", "prot_2", 150.3), ("prot_2", "prot_1", 150.3)]


def test_iter_self_align_rows_ignora_linhas_vazias(tmp_path):
    align = tmp_path / "self_align.tsv"
    align.write_text("prot_1\tprot_2\t150.3\n\n")
    rows = list(iter_self_align_rows(str(align)))
    assert rows == [("prot_1", "prot_2", 150.3)]
```

- [ ] **Step 6: Run to verify failure**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_mmseqs_taxonomy_calibration.py -v -k "protein_taxid or self_align"`
Expected: FAIL with `ImportError`

- [ ] **Step 7: Implement `load_protein_taxid_map` and `iter_self_align_rows`**

```python
# append to scripts/mmseqs_taxonomy_calibration.py

def load_protein_taxid_map(mapping_path):
    """Load a mapping.tsv ('{protein_id}\\t{taxid}' per line, no header --
    same file write_mapping() in prepare_mmseqs_taxdb.py already writes)."""
    protein_taxid = {}
    with open(mapping_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            protein_id, taxid = line.split('\t')
            protein_taxid[protein_id] = taxid
    return protein_taxid


def iter_self_align_rows(self_align_path):
    """Yield (query, target, score) from an mmseqs convertalis
    '--format-output query,target,bits' TSV (no header)."""
    with open(self_align_path) as f:
        for line in f:
            line = line.rstrip('\n')
            if not line:
                continue
            query, target, score = line.split('\t')
            yield query, target, float(score)
```

- [ ] **Step 8: Run tests to verify they pass**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_mmseqs_taxonomy_calibration.py -v`
Expected: PASS (6 tests)

- [ ] **Step 9: Commit**

```bash
git add scripts/mmseqs_taxonomy_calibration.py tests/test_mmseqs_taxonomy_calibration.py
git commit -m "feat(mmseqs_taxonomy_calibration): taxdump/mapping/self-align parsing

Pure parsing helpers for the LOO calibration module -- reads the same
taxdump/mapping.tsv prepare_mmseqs_taxdb.py already writes, plus a new
self_align.tsv (mmseqs convertalis --format-output query,target,bits).
See docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md.

No behavior change yet -- compute_rank_cutoffs (next commit) is what
actually turns this into cutoffs."
```

---

## Task 2: Módulo de calibração — `compute_rank_cutoffs` + JSON I/O

**Files:**
- Modify: `scripts/mmseqs_taxonomy_calibration.py`
- Test: `tests/test_mmseqs_taxonomy_calibration.py`

**Interfaces:**
- Consumes: nothing from Task 1 directly (takes already-loaded dicts as args, decoupled from I/O).
- Produces: `compute_rank_cutoffs(self_align_rows: Iterable[tuple[str,str,float]], protein_taxid: dict[str,str], taxid_lineage: dict[str,dict[str,str]], ranks: list[str]) -> dict[str, dict[str, float]]`
- Produces: `write_rank_cutoffs(cutoffs: dict, out_path: str) -> None`
- Produces: `load_rank_cutoffs(cutoffs_path: str) -> dict | None`

- [ ] **Step 1: Write failing tests for `compute_rank_cutoffs`**

```python
# append to tests/test_mmseqs_taxonomy_calibration.py
from mmseqs_taxonomy_calibration import compute_rank_cutoffs


def test_compute_rank_cutoffs_taxon_bem_separado_da_cutoff_alta(tmp_path=None):
    # Two units of family "Mariniviridae" (taxids 10, 11), two of
    # "Kyanoviridae" (20, 21). Every protein's best cross-unit hit lands
    # correctly on its own family -- f_same == 1.0 for everyone, no
    # cross-family confusion recorded, so cutoff == 1.0 * 0.75 (no cross
    # data) per the "no cross-taxon data" fallback.
    taxid_lineage = {
        "10": {"family": "Mariniviridae"},
        "11": {"family": "Mariniviridae"},
        "20": {"family": "Kyanoviridae"},
        "21": {"family": "Kyanoviridae"},
    }
    protein_taxid = {
        "p10a": "10", "p10b": "10",
        "p11a": "11",
        "p20a": "20",
        "p21a": "21",
    }
    rows = [
        ("p10a", "p11a", 200.0),  # best cross-unit hit for p10a: unit 11, same family
        ("p10b", "p11a", 190.0),  # same
        ("p11a", "p10a", 200.0),
        ("p20a", "p21a", 180.0),  # same family (Kyanoviridae)
        ("p21a", "p20a", 180.0),
    ]
    cutoffs = compute_rank_cutoffs(rows, protein_taxid, taxid_lineage, ["family"])
    assert cutoffs["family"]["Mariniviridae"] == pytest.approx(0.75)
    assert cutoffs["family"]["Kyanoviridae"] == pytest.approx(0.75)


def test_compute_rank_cutoffs_com_confusao_cruzada_usa_ponto_medio():
    # p10a's best cross-unit hit is unit 20 (WRONG family) -- this is a
    # false attribution to Kyanoviridae. p10b still hits correctly (unit
    # 11). So for Mariniviridae: f_same values are {p10a: 0/1=0.0 (its
    # only cross-unit hit disagreed), p10b: 1.0} -> same_fracs = [0.0, 1.0]
    # (one value per protein here since each has exactly one best hit;
    # the aggregation is per QUERY PROTEIN, not per unit, when a unit has
    # only one protein contributing evidence in this toy example).
    taxid_lineage = {
        "10": {"family": "Mariniviridae"},
        "11": {"family": "Mariniviridae"},
        "20": {"family": "Kyanoviridae"},
    }
    protein_taxid = {"p10a": "10", "p10b": "10", "p11a": "11", "p20a": "20"}
    rows = [
        ("p10a", "p20a", 300.0),  # best hit for p10a is WRONG family
        ("p10b", "p11a", 150.0),  # best hit for p10b is correct family
        ("p11a", "p10b", 150.0),
        ("p20a", "p10a", 300.0),  # best hit for p20a is WRONG family too
    ]
    cutoffs = compute_rank_cutoffs(rows, protein_taxid, taxid_lineage, ["family"])
    # Mariniviridae: same_fracs come from unit 10 (1 of 1 proteins correct
    # -> wait, p10a disagreed and p10b agreed, both belong to unit 10 ->
    # f_same(unit=10) = 1/2 = 0.5; unit 11 has only p11a, whose best hit
    # (p10b, unit 10, Mariniviridae) agrees -> f_same(unit=11) = 1.0.
    # same_fracs(Mariniviridae) = [0.5, 1.0] -> min = 0.5.
    # cross_fracs(Mariniviridae): unit 20 (true family Kyanoviridae) has
    # 1 protein (p20a), whose best hit falsely lands on Mariniviridae ->
    # f_cross(unit=20, target=Mariniviridae) = 1/1 = 1.0.
    # cutoff = (0.5 + 1.0) / 2 = 0.75
    assert cutoffs["family"]["Mariniviridae"] == pytest.approx(0.75)


def test_compute_rank_cutoffs_ignora_self_hit_literal():
    # A protein's hit against itself must not count as evidence -- if it
    # did, every unit would trivially get f_same == 1.0 regardless of
    # real cross-unit signal.
    taxid_lineage = {"10": {"family": "Mariniviridae"}}
    protein_taxid = {"p10a": "10"}
    rows = [("p10a", "p10a", 999.0)]  # only row is a literal self-hit
    cutoffs = compute_rank_cutoffs(rows, protein_taxid, taxid_lineage, ["family"])
    assert cutoffs == {}


def test_compute_rank_cutoffs_ignora_hit_mesma_unidade_taxonomica():
    # Two different proteins of the SAME unit (same species-level taxid)
    # hitting each other must not count as cross-unit evidence either --
    # they are trivially "the same taxon" by construction.
    taxid_lineage = {"10": {"family": "Mariniviridae"}}
    protein_taxid = {"p10a": "10", "p10b": "10"}
    rows = [("p10a", "p10b", 500.0), ("p10b", "p10a", 500.0)]
    cutoffs = compute_rank_cutoffs(rows, protein_taxid, taxid_lineage, ["family"])
    assert cutoffs == {}
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_mmseqs_taxonomy_calibration.py -v -k compute_rank_cutoffs`
Expected: FAIL with `ImportError: cannot import name 'compute_rank_cutoffs'`

- [ ] **Step 3: Implement `compute_rank_cutoffs`**

```python
# append to scripts/mmseqs_taxonomy_calibration.py
import collections


def compute_rank_cutoffs(self_align_rows, protein_taxid, taxid_lineage, ranks):
    """Leave-one-out calibration of a minimum agreement fraction per
    rank/taxon, from an all-vs-all self-alignment of a seqTaxDB against
    itself.

    For each protein Q, find its single best-scoring hit among rows whose
    target belongs to a DIFFERENT taxonomic unit (taxid) than Q's own --
    this excludes both the trivial literal self-hit (query == target) and
    any hit against another protein of the same species-level unit, with
    no extra genome-key bookkeeping needed (two units sharing a taxid
    already share an identical lineage by construction, see
    prepare_mmseqs_taxdb.py::build_taxdump).

    For every rank, group by the TRUE taxon of each protein's own unit at
    that rank: `f_same` is the fraction of a unit's proteins whose best
    cross-unit hit landed on that unit's own taxon (correct); `f_cross` is
    the fraction of a *different* unit's proteins whose best cross-unit
    hit falsely landed on this taxon instead. The cutoff is the VITAP
    leave-one-out formula (same_taxon.min() vs cross_taxon.max(),
    midpoint; same.min() * 0.75 with no cross-taxon evidence) applied to
    these fractions instead of raw bitscores -- see the design doc for
    why a bitscore cutoff can't plug into the rollups this feeds.

    Returns {rank: {taxon_name: cutoff}} -- only rank/taxon pairs with at
    least one contributing unit are present; anything absent means "not
    enough data to calibrate," which callers treat as "fall back to
    requiring unanimity."
    """
    best_hit = {}  # query_protein -> (score, target_protein)
    for query, target, score in self_align_rows:
        if query == target:
            continue
        q_unit = protein_taxid.get(query)
        t_unit = protein_taxid.get(target)
        if q_unit is None or t_unit is None or q_unit == t_unit:
            continue
        current = best_hit.get(query)
        if current is None or score > current[0]:
            best_hit[query] = (score, target)

    # unit -> rank -> [bool agree, ...] one per protein with a valid best hit
    agreement = collections.defaultdict(lambda: collections.defaultdict(list))
    # (rank, false_taxon) -> [bool, ...] one per protein of a unit whose
    # true taxon at that rank is NOT false_taxon, recording whether its
    # best hit falsely landed there
    false_attr = collections.defaultdict(list)

    for query, (score, target) in best_hit.items():
        q_unit = protein_taxid[query]
        t_unit = protein_taxid[target]
        q_lineage = taxid_lineage.get(q_unit, {})
        t_lineage = taxid_lineage.get(t_unit, {})
        for rank in ranks:
            q_taxon = q_lineage.get(rank)
            t_taxon = t_lineage.get(rank)
            if q_taxon is None or t_taxon is None:
                continue
            agree = q_taxon == t_taxon
            agreement[q_unit][rank].append(agree)
            if not agree:
                false_attr[(rank, t_taxon)].append(q_unit)

    # unit -> rank -> true taxon at that rank (for aggregating same_fracs by taxon)
    unit_taxon_at_rank = {
        unit: {rank: taxid_lineage.get(unit, {}).get(rank) for rank in ranks}
        for unit in agreement
    }

    same_fracs = collections.defaultdict(lambda: collections.defaultdict(list))
    for unit, by_rank in agreement.items():
        for rank, flags in by_rank.items():
            taxon = unit_taxon_at_rank[unit][rank]
            if taxon is None or not flags:
                continue
            same_fracs[rank][taxon].append(sum(flags) / len(flags))

    cross_fracs = collections.defaultdict(lambda: collections.defaultdict(list))
    for (rank, false_taxon), offending_units in false_attr.items():
        counts = collections.Counter(offending_units)
        for unit, n_false in counts.items():
            n_total = len(agreement[unit][rank])
            if n_total == 0:
                continue
            cross_fracs[rank][false_taxon].append(n_false / n_total)

    result = {}
    for rank, by_taxon in same_fracs.items():
        for taxon, fracs in by_taxon.items():
            same_min = min(fracs)
            cross_list = cross_fracs.get(rank, {}).get(taxon, [])
            cutoff = (same_min + max(cross_list)) / 2 if cross_list else same_min * 0.75
            result.setdefault(rank, {})[taxon] = cutoff
    return result


def write_rank_cutoffs(cutoffs, out_path):
    with open(out_path, 'w') as f:
        json.dump(cutoffs, f, indent=2, sort_keys=True)


def load_rank_cutoffs(cutoffs_path):
    """Return the cutoffs dict from cutoffs_path, or None if it doesn't
    exist -- signals 'this DB has not been calibrated yet' to callers,
    which then fall back to requiring unanimity."""
    if not os.path.exists(cutoffs_path):
        return None
    with open(cutoffs_path) as f:
        return json.load(f)
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_mmseqs_taxonomy_calibration.py -v`
Expected: PASS (10 tests)

- [ ] **Step 5: Write and run failing tests for `write_rank_cutoffs`/`load_rank_cutoffs`**

```python
# append to tests/test_mmseqs_taxonomy_calibration.py
from mmseqs_taxonomy_calibration import write_rank_cutoffs, load_rank_cutoffs


def test_load_rank_cutoffs_arquivo_ausente_devolve_none(tmp_path):
    assert load_rank_cutoffs(str(tmp_path / "nope.json")) is None


def test_write_then_load_rank_cutoffs_roundtrip(tmp_path):
    cutoffs = {"family": {"Mariniviridae": 0.75}}
    out_path = str(tmp_path / "rank_cutoffs.json")
    write_rank_cutoffs(cutoffs, out_path)
    assert load_rank_cutoffs(out_path) == cutoffs
```

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_mmseqs_taxonomy_calibration.py -v`
Expected: PASS (12 tests) -- these two are already implemented by Step 3,
so no further implementation step is needed here, just confirmation.

- [ ] **Step 6: Commit**

```bash
git add scripts/mmseqs_taxonomy_calibration.py tests/test_mmseqs_taxonomy_calibration.py
git commit -m "feat(mmseqs_taxonomy_calibration): compute_rank_cutoffs (LOO fraction calibration)

Implements the fraction-based leave-one-out calibration described in
docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md:
per rank/taxon, min(same-taxon agreement fraction) vs max(cross-taxon
false-attribution fraction), VITAP's midpoint formula applied to
agreement fractions instead of bitscores (the rollups this feeds have no
per-protein score available -- see design doc for why).

Not wired into anything yet -- Task 3 calls this from prepare_mmseqs_taxdb.py."
```

---

## Task 3: `--calibrate` em `scripts/prepare_mmseqs_taxdb.py`

**Files:**
- Modify: `scripts/prepare_mmseqs_taxdb.py`

**Interfaces:**
- Consumes: `compute_rank_cutoffs`, `write_rank_cutoffs`, `load_protein_taxid_map`, `load_taxdump`, `lineage_by_rank`, `iter_self_align_rows` from `scripts/mmseqs_taxonomy_calibration.py` (Tasks 1-2).
- Consumes: `FORMAT_RANKS` (module-level dict already in this file, `scripts/prepare_mmseqs_taxdb.py:91-96`).

No automated test for this task: `prepare_mmseqs_taxdb.py` has no existing
test file (it is a manual setup script that calls `mmseqs`, `wget`-free,
run by the user per INSTALL.md) and this task only wires already-tested
pure functions into it plus two `subprocess.run` calls -- consistent with
the file's existing untested `build()`/`main()` orchestration.

- [ ] **Step 1: Add the `--calibrate` flag**

In `scripts/prepare_mmseqs_taxdb.py`, inside `parse_args()`
(`scripts/prepare_mmseqs_taxdb.py:106`), add after the existing
`p.add_argument('--threads', ...)` line:

```python
    p.add_argument('--calibrate', action='store_true',
                    help='After building/verifying the seqTaxDB, run an all-vs-all '
                         'self-search and calibrate rank_cutoffs.json via leave-one-out '
                         '(see docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md). '
                         'Safe to re-run alone against an already-built seqTaxDB.')
```

- [ ] **Step 2: Add the `calibrate()` orchestration function**

Add this function right after `build()` (`scripts/prepare_mmseqs_taxdb.py`,
after line 356, before `def main():`):

```python
def calibrate(seqdb_path, out_dir, ranks, threads):
    """Run an all-vs-all self-search of seqdb_path against itself and
    calibrate rank_cutoffs.json from it. Idempotent: re-running just
    redoes the search and overwrites the cutoffs file."""
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from mmseqs_taxonomy_calibration import (
        load_taxdump, lineage_by_rank, load_protein_taxid_map,
        iter_self_align_rows, compute_rank_cutoffs, write_rank_cutoffs,
    )

    tmp_dir        = os.path.join(out_dir, 'tmp_calibrate')
    result_path    = os.path.join(out_dir, 'self_align_result')
    self_align_path = os.path.join(out_dir, 'self_align.tsv')
    cutoffs_path   = os.path.join(out_dir, 'rank_cutoffs.json')
    mapping_path   = os.path.join(out_dir, 'mapping.tsv')
    taxdump_dir    = os.path.join(out_dir, 'taxdump')

    subprocess.run(['rm', '-rf', tmp_dir] + glob.glob(result_path + '*'))
    os.makedirs(tmp_dir, exist_ok=True)

    print("[prepare_mmseqs_taxdb] Running all-vs-all self-search for calibration "
          "(mmseqs search) -- this can take a while...")
    subprocess.run(['mmseqs', 'search', seqdb_path, seqdb_path, result_path, tmp_dir,
                     '--threads', str(threads)], check=True)

    print("[prepare_mmseqs_taxdb] Converting self-alignment to TSV (mmseqs convertalis)...")
    subprocess.run(['mmseqs', 'convertalis', seqdb_path, seqdb_path, result_path,
                     self_align_path, '--format-output', 'query,target,bits'], check=True)

    print("[prepare_mmseqs_taxdb] Computing leave-one-out rank cutoffs...")
    parent_of, rank_of, name_of = load_taxdump(taxdump_dir)
    protein_taxid = load_protein_taxid_map(mapping_path)
    taxid_lineage = {
        taxid: lineage_by_rank(taxid, parent_of, rank_of, name_of)
        for taxid in set(protein_taxid.values())
    }
    cutoffs = compute_rank_cutoffs(
        iter_self_align_rows(self_align_path), protein_taxid, taxid_lineage, ranks)
    write_rank_cutoffs(cutoffs, cutoffs_path)
    n_taxa = sum(len(v) for v in cutoffs.values())
    print(f"[prepare_mmseqs_taxdb] Wrote {cutoffs_path} "
          f"({len(cutoffs)} ranks, {n_taxa} taxa calibrated)")
    return cutoffs_path
```

- [ ] **Step 3: Call `calibrate()` from `main()` when the flag is set**

In `scripts/prepare_mmseqs_taxdb.py`, `main()` (currently ending around
line 385), right after the existing `seqdb_path = build(...)` line and
before the `print(f"\n[prepare_mmseqs_taxdb] Done!")` line, add:

```python
    if args.calibrate:
        calibrate(seqdb_path, args.out, FORMAT_RANKS[args.format], args.threads)
```

- [ ] **Step 4: Verify the flag parses and the new function is importable**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -c "
import sys; sys.path.insert(0, 'scripts')
import prepare_mmseqs_taxdb as p
print(hasattr(p, 'calibrate'))
sys.argv = ['prepare_mmseqs_taxdb.py', '--faa', 'x', '--format', 'inphared',
            '--inphared-tax', 'y', '--out', 'z', '--calibrate']
args = p.parse_args()
print(args.calibrate)
"`
Expected: prints `True` then `True`

- [ ] **Step 5: Commit**

```bash
git add scripts/prepare_mmseqs_taxdb.py
git commit -m "feat(prepare_mmseqs_taxdb): add --calibrate (LOO rank_cutoffs.json)

Wires scripts/mmseqs_taxonomy_calibration.py into the DB-prep script: an
all-vs-all mmseqs search + convertalis, then compute_rank_cutoffs writes
rank_cutoffs.json beside the seqTaxDB. Manual, opt-in step (same pattern
as the rest of this script) -- re-running an existing DB with --calibrate
just redoes the calibration, build()'s own idempotency is untouched.

Not consumed by anything yet -- Tasks 4-5 wire the rollups to read it."
```

---

## Task 4: `mmseqs_lca_rollup` — mover para `scripts/votu_catalog.py` e calibrar

**Files:**
- Modify: `scripts/votu_catalog.py`
- Modify: `rules/votu_catalog.smk`
- Test: `tests/test_votu_catalog.py`

**Interfaces:**
- Consumes: `load_rank_cutoffs` from `scripts/mmseqs_taxonomy_calibration.py` (Task 2).
- Produces: `mmseqs_lca_rollup(hits_path: str, ranks: list[str], cutoffs: dict | None = None) -> dict` (replaces the module-level `_mmseqs_lca_rollup` currently inline in `rules/votu_catalog.smk:224-272`; same return shape plus a new `"confidence"` float key per contig).

- [ ] **Step 1: Write failing tests for the new fraction/cutoff behavior**

```python
# append to tests/test_votu_catalog.py
from votu_catalog import mmseqs_lca_rollup


def _write_hits(path, rows):
    """rows: list of (qseqid, taxid, rank, name, lineage)."""
    with open(path, "w") as f:
        f.write("qseqid\ttaxid\trank\tname\tlineage\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def test_rollup_sem_cutoffs_mantem_unanimidade(tmp_path):
    hits = tmp_path / "hits.tsv"
    _write_hits(str(hits), [
        ("k1_1", "10", "family", "Mariniviridae",
         "r_Duplodnaviria;k_X;p_X;c_Caudoviricetes;o_X;f_Mariniviridae"),
        ("k1_2", "10", "family", "Mariniviridae",
         "r_Duplodnaviria;k_X;p_X;c_Caudoviricetes;o_X;f_Mariniviridae"),
        ("k1_3", "11", "family", "Kyanoviridae",
         "r_Duplodnaviria;k_X;p_X;c_Caudoviricetes;o_X;f_Kyanoviridae"),
    ])
    ranks = ["realm", "kingdom", "phylum", "class", "order", "family", "subfamily", "genus"]
    result = mmseqs_lca_rollup(str(hits), ranks)
    # class-level (index 3) is unanimous across all 3 proteins (Caudoviricetes);
    # family-level (index 5) splits 2-vs-1 -- with no cutoffs, unanimity is
    # required, so the rollup stops at class, not family.
    assert result["k1"]["family"] == ""
    assert result["k1"]["rank"] == "class"


def test_rollup_com_cutoff_calibrado_aceita_maioria(tmp_path):
    hits = tmp_path / "hits.tsv"
    _write_hits(str(hits), [
        ("k1_1", "10", "family", "Mariniviridae",
         "r_Duplodnaviria;k_X;p_X;c_Caudoviricetes;o_X;f_Mariniviridae"),
        ("k1_2", "10", "family", "Mariniviridae",
         "r_Duplodnaviria;k_X;p_X;c_Caudoviricetes;o_X;f_Mariniviridae"),
        ("k1_3", "11", "family", "Kyanoviridae",
         "r_Duplodnaviria;k_X;p_X;c_Caudoviricetes;o_X;f_Kyanoviridae"),
    ])
    ranks = ["realm", "kingdom", "phylum", "class", "order", "family", "subfamily", "genus"]
    cutoffs = {"family": {"Mariniviridae": 0.6}}
    result = mmseqs_lca_rollup(str(hits), ranks, cutoffs=cutoffs)
    # 2/3 = 0.667 >= 0.6 cutoff -- Mariniviridae wins at family level now.
    assert result["k1"]["family"] == "Mariniviridae"
    assert result["k1"]["rank"] == "family"
    assert result["k1"]["confidence"] == pytest.approx(2 / 3)


def test_rollup_fracao_abaixo_do_cutoff_ainda_trunca(tmp_path):
    hits = tmp_path / "hits.tsv"
    _write_hits(str(hits), [
        ("k1_1", "10", "family", "Mariniviridae",
         "r_Duplodnaviria;k_X;p_X;c_Caudoviricetes;o_X;f_Mariniviridae"),
        ("k1_2", "11", "family", "Kyanoviridae",
         "r_Duplodnaviria;k_X;p_X;c_Caudoviricetes;o_X;f_Kyanoviridae"),
        ("k1_3", "12", "family", "Straboviridae",
         "r_Duplodnaviria;k_X;p_X;c_Caudoviricetes;o_X;f_Straboviridae"),
    ])
    ranks = ["realm", "kingdom", "phylum", "class", "order", "family", "subfamily", "genus"]
    cutoffs = {"family": {"Mariniviridae": 0.6}}
    result = mmseqs_lca_rollup(str(hits), ranks, cutoffs=cutoffs)
    # 1/3 = 0.333 < 0.6 cutoff -- stays at class.
    assert result["k1"]["family"] == ""
    assert result["k1"]["rank"] == "class"
```

- [ ] **Step 2: Run to verify failure**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_votu_catalog.py -v -k rollup`
Expected: FAIL with `ImportError: cannot import name 'mmseqs_lca_rollup'`

- [ ] **Step 3: Move and rewrite `_mmseqs_lca_rollup` into `scripts/votu_catalog.py`**

First, in `rules/votu_catalog.smk`, delete the entire `_mmseqs_lca_rollup`
function (lines 224-272, from `def _mmseqs_lca_rollup(hits_path, ranks):`
through the closing `return result`).

Then add this to `scripts/votu_catalog.py` (append at the end of the file):

```python
def mmseqs_lca_rollup(hits_path, ranks, cutoffs=None):
    """Roll up per-PROTEIN mmseqs2 taxonomy LCA calls (qseqid, taxid,
    rank, name, lineage -- from `mmseqs taxonomy` + `createtsv`) to
    per-CONTIG, by walking rank levels and, at each level, requiring the
    majority taxon's support fraction among proteins reaching that level
    to meet a calibrated cutoff.

    Shared by votu_mmseqs_taxonomy (INPHARED) and
    votu_mmseqs_taxonomy_custom (e.g. IMG/VR) results -- both produce the
    same shape over the same 8-level ICTV rank scheme (realm..genus),
    just against different seqTaxDBs, hence different `cutoffs`.

    cutoffs: {rank: {taxon: min_fraction}} from
    mmseqs_taxonomy_calibration.load_rank_cutoffs(), or None. Missing
    cutoffs (None, or no entry for a given rank/taxon) means "not
    calibrated" -- the fraction must then be exactly 1.0 (unanimity,
    the pre-calibration behavior) to keep that rank. See
    docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md.

    Returns {contig: {order, family, subfamily, genus, rank, lineage,
    confidence, n_proteins}}.
    """
    import csv, os, re, collections
    _RANK_PREFIX = re.compile(r'^[a-z]+_')
    contig_protein_lineages = collections.defaultdict(list)
    if os.path.exists(hits_path) and os.path.getsize(hits_path) > 0:
        with open(hits_path) as f:
            for row in csv.DictReader(f, delimiter="\t"):
                protein_id = row.get("qseqid", "")
                rank       = (row.get("rank") or "").strip()
                lineage    = (row.get("lineage") or "").strip()
                if not protein_id or not lineage or rank in ("", "no rank", "root"):
                    continue
                contig = "_".join(protein_id.split("_")[:-1]) or protein_id
                names = [_RANK_PREFIX.sub("", p).strip() for p in lineage.split(";") if p.strip()]
                if names:
                    contig_protein_lineages[contig].append(names)

    result = {}
    for contig, lineages in contig_protein_lineages.items():
        common = []
        confidence = 1.0
        depth = 0
        while True:
            level_values = [lin[depth] for lin in lineages if len(lin) > depth]
            if not level_values:
                break
            counts = collections.Counter(level_values)
            taxon, n_votes = counts.most_common(1)[0]
            fraction = n_votes / len(level_values)
            rank_name = ranks[depth] if depth < len(ranks) else ranks[-1]
            threshold = 1.0
            if cutoffs is not None:
                threshold = cutoffs.get(rank_name, {}).get(taxon, 1.0)
            if fraction < threshold:
                break
            common.append(taxon)
            confidence = fraction
            depth += 1
        if not common:
            continue
        d = len(common)
        result[contig] = {
            "order":      common[4] if d >= 5 else "",
            "family":     common[5] if d >= 6 else "",
            "subfamily":  common[6] if d >= 7 else "",
            "genus":      common[7] if d >= 8 else "",
            "rank":       ranks[d - 1] if d <= len(ranks) else ranks[-1],
            "lineage":    ";".join(common),
            "confidence": confidence,
            "n_proteins": len(lineages),
        }
    return result
```

Note: the previous implementation used `zip(*lineages)`, which silently
stops at the length of the SHORTEST protein lineage in a contig -- a
contig with one protein resolved to genus (8 names) and another only to
order (4 names) never even looked past index 4, regardless of agreement.
The `while True` / `level_values = [... if len(lin) > depth]` loop above
fixes this as a side effect: a level only requires votes from proteins
that actually reach it, matching the "majority fraction among proteins
reaching that level" semantics the cutoff needs to mean anything. This is
an intentional, in-scope fix, not a silent behavior change to call out
separately -- note it in the commit message.

- [ ] **Step 4: Wire the import and cutoffs-loading into `rules/votu_catalog.smk`**

At the top of `rules/votu_catalog.smk` (line 17-19), change:

```python
from votu_catalog import (
    build_pool, parse_skani_sparse, cluster_votus, write_clusters,
)
```

to:

```python
from votu_catalog import (
    build_pool, parse_skani_sparse, cluster_votus, write_clusters,
    mmseqs_lca_rollup,
)
from mmseqs_taxonomy_calibration import load_rank_cutoffs
```

In the `votu_taxonomy` rule (`rules/votu_catalog.smk`, starts at line 682),
add a `params:` block (there isn't one currently) right before `output:`:

```python
    params:
        inphared_seqtaxdb = f"{INPHARED_DB}/inphared_mmseqs_taxdb/seqTaxDB",
        custom_seqtaxdb   = CUSTOM_VIRAL_MMSEQS_DB,
```

Then inside the `run:` block, replace the two rollup calls (currently at
lines 758 and 761):

```python
        _RANKS = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus']
        mmseqs_tax = _mmseqs_lca_rollup(str(input.mmseqs_hits), _RANKS)
        lf.write(f"MMseqs2/INPHARED: {len(mmseqs_tax)} contigs (namespaced keys)\n")

        custom_tax = _mmseqs_lca_rollup(str(input.custom_hits), _RANKS)
        lf.write(f"MMseqs2/Custom: {len(custom_tax)} contigs (namespaced keys)\n")
```

with:

```python
        _RANKS = ['realm', 'kingdom', 'phylum', 'class', 'order', 'family', 'subfamily', 'genus']

        inphared_cutoffs_path = os.path.join(os.path.dirname(str(params.inphared_seqtaxdb)), "rank_cutoffs.json")
        inphared_cutoffs = load_rank_cutoffs(inphared_cutoffs_path)
        lf.write(f"INPHARED rank cutoffs: {'calibrated (' + inphared_cutoffs_path + ')' if inphared_cutoffs else 'not calibrated -- unanimity'}\n")
        mmseqs_tax = mmseqs_lca_rollup(str(input.mmseqs_hits), _RANKS, cutoffs=inphared_cutoffs)
        lf.write(f"MMseqs2/INPHARED: {len(mmseqs_tax)} contigs (namespaced keys)\n")

        custom_cutoffs = None
        if params.custom_seqtaxdb:
            custom_cutoffs_path = os.path.join(os.path.dirname(str(params.custom_seqtaxdb)), "rank_cutoffs.json")
            custom_cutoffs = load_rank_cutoffs(custom_cutoffs_path)
            lf.write(f"Custom viral rank cutoffs: {'calibrated (' + custom_cutoffs_path + ')' if custom_cutoffs else 'not calibrated -- unanimity'}\n")
        custom_tax = mmseqs_lca_rollup(str(input.custom_hits), _RANKS, cutoffs=custom_cutoffs)
        lf.write(f"MMseqs2/Custom: {len(custom_tax)} contigs (namespaced keys)\n")
```

Finally, update the two confidence-text format strings later in the same
`run:` block (currently, around line 843-851):

```python
            if mms and (mms.get("family") or mms.get("genus") or mms.get("order")):
                ff, fg, fo = mms.get("family",""), mms.get("genus",""), mms.get("order","")
                candidates.append(("mmseqs_inphared", ff, fg, fo, mms.get("lineage",""),
                                    f"{mms.get('rank','')} ({mms.get('n_proteins',0)} proteins)",
                                    fg or ff or fo))

            if cms and (cms.get("family") or cms.get("genus") or cms.get("order")):
                ff, fg, fo = cms.get("family",""), cms.get("genus",""), cms.get("order","")
                candidates.append(("mmseqs_custom", ff, fg, fo, cms.get("lineage",""),
                                    f"{cms.get('rank','')} ({cms.get('n_proteins',0)} proteins)",
                                    fg or ff or fo))
```

to include the confidence fraction:

```python
            if mms and (mms.get("family") or mms.get("genus") or mms.get("order")):
                ff, fg, fo = mms.get("family",""), mms.get("genus",""), mms.get("order","")
                candidates.append(("mmseqs_inphared", ff, fg, fo, mms.get("lineage",""),
                                    f"{mms.get('rank','')} ({mms.get('confidence',1.0):.2f} frac, {mms.get('n_proteins',0)} proteins)",
                                    fg or ff or fo))

            if cms and (cms.get("family") or cms.get("genus") or cms.get("order")):
                ff, fg, fo = cms.get("family",""), cms.get("genus",""), cms.get("order","")
                candidates.append(("mmseqs_custom", ff, fg, fo, cms.get("lineage",""),
                                    f"{cms.get('rank','')} ({cms.get('confidence',1.0):.2f} frac, {cms.get('n_proteins',0)} proteins)",
                                    fg or ff or fo))
```

- [ ] **Step 5: Run tests to verify they pass**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_votu_catalog.py -v`
Expected: PASS, including the 3 new rollup tests and every pre-existing
test in this file (no regressions).

- [ ] **Step 6: Sanity-check the Snakefile still parses**

Run: `cd /home/alumnos/lmelo/vapor-pc && snakemake -n --use-conda 2>&1 | tail -30`
Expected: no `SyntaxError`/`NameError` from `rules/votu_catalog.smk`; the
dry-run plan should list `votu_taxonomy` same as before this change (a
full green dry-run needs real sample data and databases configured, which
this environment may not have -- clean output with no traceback pointing
at `votu_catalog.smk` is sufficient confirmation here).

- [ ] **Step 7: Commit**

```bash
git add scripts/votu_catalog.py rules/votu_catalog.smk tests/test_votu_catalog.py
git commit -m "feat(votu_catalog): calibrated fraction cutoff in mmseqs_lca_rollup

Moves _mmseqs_lca_rollup out of rules/votu_catalog.smk into
scripts/votu_catalog.py as mmseqs_lca_rollup (testable, matches the
project's pure-logic-in-scripts/-glue-in-rules/ convention already used
for genomad_base_contig/resolve_genomad_key in this same file).

Behavior change: instead of requiring 100% protein agreement to keep a
rank, computes the majority taxon's support fraction among proteins that
reach that rank and compares it to a calibrated cutoff from
scripts/mmseqs_taxonomy_calibration.py (rank_cutoffs.json beside the
seqTaxDB, produced by 'prepare_mmseqs_taxdb.py --calibrate'). No
rank_cutoffs.json -> falls back to the old unanimity requirement, zero
regression for uncalibrated DBs.

Side effect fix: the old zip(*lineages)-based walk silently stopped at
the length of a contig's SHORTEST resolved protein lineage, so a single
shallow protein call could mask real agreement at deeper ranks among the
contig's other proteins. The rewritten walk only requires votes from
proteins that actually reach each rank -- a real, in-scope correctness
fix needed for the fraction/cutoff logic to mean what it says, not a
silent side change.

votu_taxonomy now loads rank_cutoffs.json (if present) for both the
INPHARED and custom viral seqTaxDBs and reports it in its own log; the
report's confidence text now shows the support fraction, not just rank
and protein count.

See docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md."
```

---

## Task 5: `load_mmseqs_taxonomy_prok` — mesma calibração no report

**Files:**
- Modify: `scripts/report/data_loaders.py`
- Modify: `scripts/report/renderer.py`
- Modify: `rules/report.smk`
- Test: `tests/test_report_votu_loaders.py`

**Interfaces:**
- Consumes: `load_rank_cutoffs` from `scripts/mmseqs_taxonomy_calibration.py` (Task 2).
- Produces: `load_mmseqs_taxonomy_prok(paths_d: dict, samples: list, cutoffs: dict | None = None) -> list[dict]` (adds one new optional parameter to the existing signature; return shape unchanged except each row gains a `'Confidence'` float key).

- [ ] **Step 1: Write failing tests**

```python
# append to tests/test_report_votu_loaders.py
from report.data_loaders import load_mmseqs_taxonomy_prok


def _write_prok_hits(path, rows):
    with open(path, "w") as f:
        f.write("qseqid\ttaxid\trank\tname\tlineage\n")
        for r in rows:
            f.write("\t".join(r) + "\n")


def test_prok_rollup_sem_cutoffs_mantem_unanimidade(tmp_path):
    hits = tmp_path / "taxonomy.tsv"
    _write_prok_hits(str(hits), [
        ("bin1__k1_1", "10", "family", "Rhodobacteraceae", "d_Bacteria;p_X;c_X;o_X;f_Rhodobacteraceae"),
        ("bin1__k1_2", "10", "family", "Rhodobacteraceae", "d_Bacteria;p_X;c_X;o_X;f_Rhodobacteraceae"),
        ("bin1__k1_3", "11", "family", "Other", "d_Bacteria;p_X;c_X;o_X;f_Other"),
    ])
    records = load_mmseqs_taxonomy_prok({"S1": str(hits)}, ["S1"])
    assert len(records) == 1
    assert records[0]["Family"] == ""  # 2-vs-1 split, no cutoffs -> unanimity required -> stops at Order


def test_prok_rollup_com_cutoff_aceita_maioria(tmp_path):
    hits = tmp_path / "taxonomy.tsv"
    _write_prok_hits(str(hits), [
        ("bin1__k1_1", "10", "family", "Rhodobacteraceae", "d_Bacteria;p_X;c_X;o_X;f_Rhodobacteraceae"),
        ("bin1__k1_2", "10", "family", "Rhodobacteraceae", "d_Bacteria;p_X;c_X;o_X;f_Rhodobacteraceae"),
        ("bin1__k1_3", "11", "family", "Other", "d_Bacteria;p_X;c_X;o_X;f_Other"),
    ])
    cutoffs = {"family": {"Rhodobacteraceae": 0.6}}
    records = load_mmseqs_taxonomy_prok({"S1": str(hits)}, ["S1"], cutoffs=cutoffs)
    assert records[0]["Family"] == "Rhodobacteraceae"
    assert records[0]["Confidence"] == pytest.approx(2 / 3)
```

Note: `_prok_genome_unit` (`scripts/report/data_loaders.py:1085`) delegates
to `_split_genome_prefix` (`scripts/report/data_loaders.py:1076`), which
splits on the FIRST `__` only (`value.split(sep, 1)`). The fixture IDs
above use a plain `{bin}__{contig}_{n}` shape (one `__`) so the unit comes
out as `"bin1"`, matching what this task actually needs to verify (the
aggregation logic, not the grouping key).

**Separate finding, out of scope for this task**: a real
`mag_mmseqs_taxonomy_prok` qseqid is 3-segment
(`{sample}__{bin}__{contig}_{n}`, per this project's CLAUDE.md — e.g.
`S1__binette_bin1__k141_1_5`), and `_split_genome_prefix`'s first-`__` cut
on a 3-segment ID returns just `{sample}` (`"S1"`), not `{sample}__{bin}`
(`"S1__binette_bin1"`) — the exact bug pattern CLAUDE.md documents
elsewhere as `resolve_prefixed_id`'s reason for existing ("Never cut the
ID at the first `__`"). If confirmed, `load_mmseqs_taxonomy_prok` would
be aggregating real hits.tsv files by SAMPLE instead of by MAG today,
independent of anything in this plan. Flag this to the user/reviewer
rather than fixing it here — it's a pre-existing correctness question
about the grouping key, not the aggregation rule this task changes, and
deserves its own investigation and spec if confirmed.

- [ ] **Step 2: Run to verify failure**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_report_votu_loaders.py -v -k prok_rollup`
Expected: FAIL (assertion or TypeError on the unexpected `cutoffs` kwarg)

- [ ] **Step 3: Rewrite `load_mmseqs_taxonomy_prok`**

Replace the function body in `scripts/report/data_loaders.py:708-747`
entirely:

```python
def load_mmseqs_taxonomy_prok(paths_d, samples, cutoffs=None):
    """mmseqs_taxonomy_prok output (qseqid, taxid, rank, name, lineage) --
    mmseqs already computes a real per-PROTEIN lowest-common-ancestor, so
    this aggregates to genome level with a SECOND pass across each genome
    unit's own proteins: at each rank, the majority taxon's support
    fraction among proteins reaching that rank is compared against a
    calibrated cutoff (from
    scripts/mmseqs_taxonomy_calibration.load_rank_cutoffs(), rank_cutoffs.json
    beside CUSTOM_PROK_MMSEQS_DB). No cutoffs (None, or no entry for a
    given rank/taxon) falls back to requiring unanimity -- the original
    behavior, chosen specifically to avoid a majority vote reintroducing
    the "spurious specificity" problem (von Meijenfeldt et al. 2019,
    CAT/BAT) with an arbitrary, uncalibrated threshold. A calibrated
    fraction is not that same risk: it's validated against this DB's own
    leave-one-out self-alignment, not a number picked by hand. See
    docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md.
    Unclassified proteins (taxid 0) are dropped instead of voting for
    'no rank'."""
    import collections as _collections
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
            common = []
            confidence = 1.0
            depth = 0
            while True:
                level_values = [lin[depth] for lin in lins if len(lin) > depth]
                if not level_values:
                    break
                counts = _collections.Counter(level_values)
                taxon, n_votes = counts.most_common(1)[0]
                fraction = n_votes / len(level_values)
                rank_name = _PROK_RANK_LEVELS[depth].lower() if depth < len(_PROK_RANK_LEVELS) else _PROK_RANK_LEVELS[-1].lower()
                threshold = 1.0
                if cutoffs is not None:
                    threshold = cutoffs.get(rank_name, {}).get(taxon, 1.0)
                if fraction < threshold:
                    break
                common.append(taxon)
                confidence = fraction
                depth += 1
            if not common: continue
            rec = {'sample': s, 'Bin': unit, 'Source': 'mmseqs_lca', 'Confidence': confidence}
            for i, level in enumerate(_PROK_RANK_LEVELS):
                rec[level] = common[i] if i < len(common) else ''
            rec['Organism'] = common[-1]
            rec['Sci_name'] = common[-1]
            records.append(rec)
    return records
```

- [ ] **Step 4: Run tests to verify they pass**

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_report_votu_loaders.py -v -k prok_rollup`
Expected: PASS

Then run the full report test file to check for regressions:

Run: `cd /home/alumnos/lmelo/vapor-pc && python -m pytest tests/test_report_votu_loaders.py -v`
Expected: PASS, all tests including pre-existing ones.

- [ ] **Step 5: Wire `cutoffs` through to the call site**

In `rules/report.smk`, in the `generate_report` rule's `params:` block
(starts at `rules/report.smk:97`), add after the existing
`gtdbtk_db = GTDBTK_DB,` line:

```python
        custom_prok_mmseqs_db = CUSTOM_PROK_MMSEQS_DB,
```

In `scripts/report/renderer.py`, add the import at the top (wherever the
existing `from report.data_loaders import (..., load_mmseqs_taxonomy_prok, load_phist, ...)`
block is, `scripts/report/renderer.py:19`) — no new import needed, just add
one line right before the existing call at line 144:

```python
    mmseqs_prok_cutoffs_path = os.path.join(
        os.path.dirname(getattr(snakemake.params, 'custom_prok_mmseqs_db', '') or '.'),
        "rank_cutoffs.json")
    mmseqs_prok_cutoffs = load_rank_cutoffs(mmseqs_prok_cutoffs_path) \
        if getattr(snakemake.params, 'custom_prok_mmseqs_db', '') else None
    mmseqs_prok_data = load_mmseqs_taxonomy_prok(mmseqs_prok_paths, samples, cutoffs=mmseqs_prok_cutoffs)
```

replacing the current single-line call:

```python
    mmseqs_prok_data = load_mmseqs_taxonomy_prok(mmseqs_prok_paths, samples)
```

Add the needed import near the top of `scripts/report/renderer.py` (same
place `_sys.path`/`SCRIPTS_DIR`-style imports for this package already
live — check the top of the file for the existing `sys.path` setup used
to import `report.data_loaders`, and add alongside it):

```python
from mmseqs_taxonomy_calibration import load_rank_cutoffs
```

- [ ] **Step 6: Commit**

```bash
git add scripts/report/data_loaders.py scripts/report/renderer.py rules/report.smk tests/test_report_votu_loaders.py
git commit -m "feat(report): calibrated fraction cutoff in load_mmseqs_taxonomy_prok

Same change as the votu_catalog rollup (previous commit): the genome-
level consensus pass now compares the majority taxon's support fraction
against a calibrated cutoff (rank_cutoffs.json beside
CUSTOM_PROK_MMSEQS_DB, from 'prepare_mmseqs_taxdb.py --calibrate')
instead of requiring 100% protein agreement. No rank_cutoffs.json ->
falls back to unanimity, zero regression for an uncalibrated custom prok DB.

Wired through rules/report.smk's params -> renderer.py -> the loader's
new optional `cutoffs` kwarg.

See docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md."
```

---

## Task 6: Documentação (INSTALL.md, CLAUDE.md)

**Files:**
- Modify: `INSTALL.md`
- Modify: `CLAUDE.md`

- [ ] **Step 1: Add `--calibrate` to INSTALL.md**

Find the section of `INSTALL.md` documenting `prepare_mmseqs_taxdb.py`
(search for `prepare_mmseqs_taxdb.py` in the file) and add, right after
the existing usage example for building a seqTaxDB, a new paragraph:

```markdown
After building a seqTaxDB, optionally calibrate it with `--calibrate`:

    python scripts/prepare_mmseqs_taxdb.py --format inphared \
        --inphared-tax <path> --faa <path> --out <seqtaxdb-dir> --calibrate

This runs an all-vs-all self-search once and writes `rank_cutoffs.json`
beside the seqTaxDB. The taxonomy rollups (`votu_taxonomy`, the report's
prokaryote taxonomy panel) automatically use it when present, falling
back to requiring unanimous protein agreement when it's absent -- so this
step is optional and safe to add later to an already-built DB. See
`docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md`
for why this exists.
```

- [ ] **Step 2: Update CLAUDE.md's `votu_catalog.smk` and `taxonomy.smk` bullets**

In `CLAUDE.md`, find the `votu_catalog.smk` bullet (mentions
"prodigal-gv gene calling", "MMseqs2 taxonomy", etc.) and the
`taxonomy.smk` bullet, and append one sentence to each noting the
calibration exists, e.g. for `votu_catalog.smk`:

```
`votu_taxonomy`'s mmseqs rollup (`scripts/votu_catalog.py::mmseqs_lca_rollup`)
accepts a calibrated per-rank/taxon agreement-fraction cutoff from
`scripts/mmseqs_taxonomy_calibration.py` (`rank_cutoffs.json` beside a
seqTaxDB, produced by `prepare_mmseqs_taxdb.py --calibrate`) instead of
requiring unanimous protein agreement -- falls back to unanimity when a
DB hasn't been calibrated.
```

and analogously for the prokaryote custom-taxonomy consumer in
`taxonomy.smk`'s bullet or `annotation.smk`/wherever `mag_mmseqs_taxonomy_prok`
is described.

- [ ] **Step 3: Commit**

```bash
git add INSTALL.md CLAUDE.md
git commit -m "docs: document --calibrate and the calibrated taxonomy rollup cutoff

See docs/superpowers/specs/2026-09-06-mmseqs-taxonomy-loo-calibration-design.md."
```
