# Pangenoma dos clusters com ilha — Fase 1 (matriz gene × membro)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Construir a matriz gene × membro (defesa + AMR) para os clusters de MAG que passam um portão de interesse, medindo de quebra se a fase 2 (PPanGGOLiN) tem sustentação nos dados.

**Architecture:** Nove regras globais novas em `rules/pangenome.smk` (select, proteins, defensefinder, e a cadeia de cinco do AMR, mais a matriz), todas seguindo o padrão `mag_bakta` (seleção dependente de dados DENTRO de um job de número fixo, DAG estático — nunca `checkpoint`). A anotação por membro reaproveita máquina existente: `mag_pangenome_proteins` é `mag_catalog_proteins` apontado para outro diretório, e as regras de defesa/AMR são herdadas com `use rule ... as ... with:` sobre o manifesto novo. A lógica pura sai para três scripts testáveis.

**Tech Stack:** Snakemake, Python 3.11, pytest, prodigal (env_viral), DefenseFinder (env_defense), AMRFinderPlus/RGI/DeepARG/argNorm (env_defense/env_argnorm).

**Spec:** `docs/superpowers/specs/2026-08-19-pangenoma-clusters-defesa-design.md`

## Global Constraints

- **Mínimo de membros por cluster: 3.** Valor fixo, não configurável nesta fase.
- **Piso de completude para avaliabilidade: CheckM2 >= 70.0%** — o mesmo de `mag_bakta` (`BAKTA_MIN_COMPLETENESS`, `config.yaml`).
- **A matriz tem TRÊS estados: `x` presente, `.` ausente, `?` não avaliável.** Membro `?` nunca entra no denominador de frequência.
- **Sem flag de config.** A regra roda sempre; quando nenhum cluster passa, escreve `skipped: no eligible clusters` no `done.txt`. Nunca `touch` vazio.
- **Proibido `checkpoint` do Snakemake.** Quebraria a invariante de dry-run do roadmap. Padrão obrigatório: `mag_bakta` (`rules/annotation.smk:25`).
- **Proibido cortar ID no primeiro `__`.** Usar `resolve_prefixed_id` de `scripts/mag_catalog.py` contra o conjunto de nomes conhecidos.
- **AMR entra só por consenso `n_tools >= 2`**, com nome normalizado pelo argNorm.
- **Commits sem rodapé `Co-Authored-By`** — convenção deste repositório.
- **Testes rodam com:** `conda run -n snakemake python -m pytest tests/ -q` (não há pytest no python base). Piso: os 171 testes atuais continuam passando.
- **Delta de DAG esperado ao final: +9 jobs globais fixos, ZERO por amostra.** Qualquer job por amostra é bug de desenho. (A spec diz "+5" porque agrupou defesa+AMR em duas linhas de tabela; o AMR são cinco regras encadeadas. A Task 7 corrige a spec.)

---

## Estrutura de arquivos

| arquivo | responsabilidade |
|---|---|
| `scripts/defense_islands.py` | **criar** — núcleo puro da detecção de ilha, extraído do relatório |
| `scripts/report/data_loaders.py` | **modificar** — passa a importar de `defense_islands` em vez de definir |
| `scripts/pangenome_select.py` | **criar** — portão puro: quem é elegível e por quê |
| `scripts/pangenome_matrix.py` | **criar** — matriz de três estados e sumário por cluster |
| `rules/pangenome.smk` | **criar** — as nove regras |
| `Snakefile` | **modificar** — `include:` do módulo novo |
| `tests/test_defense_islands.py` | **criar** |
| `tests/test_pangenome_select.py` | **criar** |
| `tests/test_pangenome_matrix.py` | **criar** |

---

### Task 1: Extrair a detecção de ilha para um módulo compartilhado

Hoje `compute_defense_islands` vive em `scripts/report/data_loaders.py:959` — camada de relatório, que a pipeline não alcança. Sem esta extração passariam a existir duas definições de ilha divergindo em silêncio.

**Files:**
- Create: `scripts/defense_islands.py`
- Modify: `scripts/report/data_loaders.py:897-1022`
- Test: `tests/test_defense_islands.py`

**Interfaces:**
- Consumes: nada (primeira tarefa)
- Produces:
  - `find_islands(genes_by_contig, prot_to_sys, min_genes=5, min_systems=3, window=10) -> list[dict]` — núcleo puro. `genes_by_contig` é `{contig: [{'Protein','Start','End','Strand'}, ...]}` em ordem genômica; `prot_to_sys` é `{protein_id: (bin, system, system_id)}`. Cada dict devolvido tem as chaves `Contig, n_genes, n_systems, Systems, window_genes, start_idx, end_idx, Start_bp, End_bp` (sem `sample` nem `Bin` — quem chama acrescenta).
  - `genes_by_contig(faa_path) -> dict` — mesma leitura de header do Prodigal já usada pelo relatório.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_defense_islands.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from defense_islands import find_islands


def _genes(n, contig="k141_1"):
    """n genes consecutivos com coordenadas de 1000 bp cada."""
    return {contig: [{"Protein": f"{contig}_{i+1}", "Start": i * 1000 + 1,
                      "End": (i + 1) * 1000, "Strand": 1} for i in range(n)]}


class TestFindIslands:
    def test_run_of_defense_genes_is_an_island(self):
        genes = _genes(10)
        prot_to_sys = {f"k141_1_{i}": ("bin1", f"Sys{i}", f"id{i}")
                       for i in range(1, 6)}
        islands = find_islands(genes, prot_to_sys, min_genes=5, min_systems=3)
        assert len(islands) == 1
        assert islands[0]["n_genes"] == 5
        assert islands[0]["n_systems"] == 5
        assert islands[0]["Contig"] == "k141_1"

    def test_too_few_systems_is_not_an_island(self):
        # 5 genes, mas todos do MESMO sistema: nao e ilha.
        genes = _genes(10)
        prot_to_sys = {f"k141_1_{i}": ("bin1", "RM", "id1") for i in range(1, 6)}
        assert find_islands(genes, prot_to_sys, min_genes=5, min_systems=3) == []

    def test_genes_further_apart_than_window_split_into_two_clusters(self):
        genes = _genes(40)
        # 3 genes no inicio, 3 no fim: nenhum grupo alcanca min_genes=5
        hits = [1, 2, 3, 30, 31, 32]
        prot_to_sys = {f"k141_1_{i}": ("bin1", f"Sys{i}", f"id{i}") for i in hits}
        assert find_islands(genes, prot_to_sys, min_genes=5, min_systems=3) == []

    def test_island_carries_genomic_extent(self):
        genes = _genes(10)
        prot_to_sys = {f"k141_1_{i}": ("bin1", f"Sys{i}", f"id{i}")
                       for i in range(1, 6)}
        island = find_islands(genes, prot_to_sys, min_genes=5, min_systems=3)[0]
        assert island["Start_bp"] == 1
        assert island["End_bp"] == 5000

    def test_no_defense_genes_yields_nothing(self):
        assert find_islands(_genes(10), {}, min_genes=5, min_systems=3) == []
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n snakemake python -m pytest tests/test_defense_islands.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'defense_islands'`

- [ ] **Step 3: Write minimal implementation**

Criar `scripts/defense_islands.py` movendo o corpo de `compute_defense_islands` (`data_loaders.py:959-1022`) e de `_genes_by_contig`/`_contig_from_protein_id` (`data_loaders.py:920-957`). O loop externo por amostra/manifesto NÃO vem junto — fica em quem chama.

```python
#!/usr/bin/env python3
"""Deteccao de ilha de defesa -- nucleo puro, compartilhado.

Ate 2026-08-19 isto vivia so em scripts/report/data_loaders.py, ou seja, na
camada de RELATORIO: a pipeline nao alcancava. O portao do pangenoma
(rules/pangenome.smk) precisa da mesma definicao, e duas copias divergiriam
em silencio -- a familia de defeito que docs/ROADMAP_SIMPLIFICACAO.md
persegue. Uma definicao, dois consumidores.
"""

import os
from collections import defaultdict


def _contig_from_protein_id(pid):
    """'k141_219139_5' -> 'k141_219139'. O Prodigal numera '{seqid}_{n}'."""
    if not pid or "_" not in pid:
        return ""
    return pid.rsplit("_", 1)[0]


def genes_by_contig(faa_path):
    """Genes por contig em ordem genomica, com coordenadas e fita.

    O Prodigal codifica o locus no proprio header do FASTA:
        >{seqid}_{n} # {start} # {end} # {strand} # ID=...;partial=...
    entao coordenadas reais saem do .faa que a pipeline ja escreve.
    """
    by_contig = defaultdict(list)
    if not faa_path or not os.path.exists(faa_path):
        return by_contig
    with open(faa_path) as f:
        for line in f:
            if not line.startswith(">"):
                continue
            header = line[1:].rstrip("\n")
            pid = header.split()[0].strip()
            contig = _contig_from_protein_id(pid)
            if not contig:
                continue
            start = end = strand = None
            parts = [p.strip() for p in header.split("#")]
            if len(parts) >= 4:
                try:
                    start, end, strand = int(parts[1]), int(parts[2]), int(parts[3])
                except (ValueError, IndexError):
                    start = end = strand = None
            by_contig[contig].append({"Protein": pid, "Start": start,
                                      "End": end, "Strand": strand})
    return by_contig


def find_islands(genes_by_contig_map, prot_to_sys,
                 min_genes=5, min_systems=3, window=10):
    """Uma entrada por ilha: corrida de genes de defesa no mesmo contig onde
    posicoes consecutivas nunca distam mais que `window` genes, com
    >= min_genes genes de >= min_systems sistemas distintos.

    `prot_to_sys`: {protein_id: (bin, system, system_id)}.
    """
    islands = []
    for contig, gene_recs in genes_by_contig_map.items():
        ordered = [g["Protein"] for g in gene_recs]
        hits = [(i, p, prot_to_sys[p]) for i, p in enumerate(ordered)
                if p in prot_to_sys]
        if len(hits) < min_genes:
            continue

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
                    "Protein": rec["Protein"], "Index": i,
                    "System": defense_idx.get(i, ""),
                    "System_id": sysid_idx.get(i, ""),
                    "Start": rec.get("Start"), "End": rec.get("End"),
                    "Strand": rec.get("Strand"),
                })
            coords = [g for g in window_genes
                      if g["Start"] is not None and g["End"] is not None]
            islands.append({
                "Contig": contig,
                "n_genes": len(cluster), "n_systems": len(systems),
                "Systems": sorted(systems), "window_genes": window_genes,
                "start_idx": start_idx, "end_idx": end_idx,
                "Start_bp": min(g["Start"] for g in coords) if coords else None,
                "End_bp": max(g["End"] for g in coords) if coords else None,
            })

        cluster = []
        for h in hits:
            if cluster and h[0] - cluster[-1][0] > window:
                flush(cluster)
                cluster = []
            cluster.append(h)
        flush(cluster)
    return islands
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n snakemake python -m pytest tests/test_defense_islands.py -q`
Expected: PASS (5 testes)

- [ ] **Step 5: Reescrever `compute_defense_islands` como casca sobre o módulo**

Em `scripts/report/data_loaders.py`, substituir o corpo de `compute_defense_islands` (linhas 959-1022) e remover `_genes_by_contig`/`_contig_from_protein_id` de lá, importando do módulo novo. Manter a assinatura pública EXATA — o `renderer.py` chama com os mesmos argumentos.

```python
def compute_defense_islands(manifest_paths, samples, defense_data,
                            min_genes=5, min_systems=3, window=10):
    """One row per defense island. O nucleo vive em scripts/defense_islands.py
    desde 2026-08-19, compartilhado com o portao do pangenoma -- aqui fica
    apenas o laco por amostra/manifesto, que e especifico do relatorio."""
    islands = []
    for manifest_path, s in zip(manifest_paths, samples):
        if not manifest_path or not os.path.exists(manifest_path):
            continue
        prot_to_sys = {}
        for rec in defense_data:
            if rec.get('sample') != s:
                continue
            for prot in rec.get('Proteins', []):
                prot_to_sys[prot] = (rec.get('Bin'), rec.get('System'),
                                     rec.get('System_id'))

        for name, mode, fna, faa, gff in _read_protein_manifest(manifest_path):
            for isl in _find_islands(_genes_by_contig(faa), prot_to_sys,
                                     min_genes=min_genes,
                                     min_systems=min_systems, window=window):
                isl['sample'] = s
                isl['Bin'] = name
                islands.append(isl)
    return islands
```

E no topo de `data_loaders.py`, junto dos outros imports:

```python
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from defense_islands import find_islands as _find_islands, genes_by_contig as _genes_by_contig
```

`_ordered_proteins_by_contig` (`data_loaders.py:912`) continua funcionando: ela chama `_genes_by_contig`, que agora é o nome importado.

- [ ] **Step 6: Run the full suite — nada do relatório pode ter quebrado**

Run: `conda run -n snakemake python -m pytest tests/ -q`
Expected: PASS, 176 testes (171 atuais + 5 novos)

- [ ] **Step 7: Commit**

```bash
git add scripts/defense_islands.py tests/test_defense_islands.py scripts/report/data_loaders.py
git commit -m "refactor(defesa): extrai a deteccao de ilha para um modulo compartilhado

A logica vivia so em scripts/report/data_loaders.py, na camada de
relatorio, onde a pipeline nao alcanca. O portao do pangenoma precisa da
mesma definicao; duas copias divergiriam em silencio."
```

---

### Task 2: O portão — quem é elegível e por quê

**Files:**
- Create: `scripts/pangenome_select.py`
- Test: `tests/test_pangenome_select.py`

**Interfaces:**
- Consumes: `defense_islands.find_islands` (Task 1)
- Produces:
  - `load_membership(path) -> dict[str, list[dict]]` — `{representative_id: [{'source_id','original_bin_id','member_id'}, ...]}`, lido de `mag_membership.tsv` (colunas `source_id / original_bin_id / member_id / representative_id`).
  - `load_completeness(path) -> dict[str, float]` — `{genome: completeness}` do `checkm2_quality_report.tsv` do catálogo (colunas `Name`, `Completeness`).
  - `select_clusters(membership, evidence, min_members=3) -> list[dict]` — cada dict com `representative_id, n_members, n_islands, n_systems, n_args, n_plasmid, criterio, eligible`. `evidence` é `{representative_id: {'n_islands': int, 'n_systems': int, 'n_args': int, 'n_plasmid': int}}`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_pangenome_select.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from pangenome_select import select_clusters

MEMB = {
    "S1__binette_bin1": [{"member_id": f"S{i}__binette_bin1"} for i in range(1, 5)],
    "S9__binette_bin2": [{"member_id": "S9__binette_bin2"}, {"member_id": "S8__binette_bin7"}],
}
NO_EVIDENCE = {"n_islands": 0, "n_systems": 0, "n_args": 0, "n_plasmid": 0}


def _sel(evidence, membership=MEMB, **kw):
    return {c["representative_id"]: c
            for c in select_clusters(membership, evidence, **kw)}


class TestSelectClusters:
    def test_island_with_enough_members_is_eligible(self):
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_islands=1)}
        out = _sel(ev)
        assert out["S1__binette_bin1"]["eligible"] is True
        assert out["S1__binette_bin1"]["criterio"] == "ilha"

    def test_two_members_is_never_eligible_however_strong_the_evidence(self):
        # O piso de 3 e absoluto: com 2 nao ha frequencia que signifique nada.
        ev = {"S9__binette_bin2": dict(NO_EVIDENCE, n_islands=5, n_args=9)}
        assert _sel(ev)["S9__binette_bin2"]["eligible"] is False

    def test_three_defense_systems_without_island_is_eligible(self):
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_systems=3)}
        out = _sel(ev)["S1__binette_bin1"]
        assert out["eligible"] is True
        assert out["criterio"] == "sistemas"

    def test_two_defense_systems_without_island_is_not(self):
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_systems=2)}
        assert _sel(ev)["S1__binette_bin1"]["eligible"] is False

    def test_consensus_arg_alone_is_eligible(self):
        # O AMR entra por direito proprio: um cluster sem ilha e com ARG e
        # exatamente um caso que um portao so-de-defesa perderia.
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_args=1)}
        out = _sel(ev)["S1__binette_bin1"]
        assert out["eligible"] is True
        assert out["criterio"] == "amr"

    def test_plasmid_alone_is_NOT_eligible(self):
        # Plasmidio sem defesa nem ARG e sinal de mobilidade, nao motivo.
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_plasmid=4)}
        assert _sel(ev)["S1__binette_bin1"]["eligible"] is False

    def test_plasmid_reinforces_but_island_names_the_criterion(self):
        ev = {"S1__binette_bin1": dict(NO_EVIDENCE, n_islands=1, n_plasmid=2)}
        out = _sel(ev)["S1__binette_bin1"]
        assert out["eligible"] is True
        assert out["criterio"] == "ilha"

    def test_cluster_without_evidence_row_is_reported_not_dropped(self):
        # Toda linha aparece no candidates.tsv: a selecao tem de ser
        # auditavel, nao um numero que apareceu.
        out = _sel({})
        assert set(out) == {"S1__binette_bin1", "S9__binette_bin2"}
        assert all(c["eligible"] is False for c in out.values())

    def test_n_members_counts_the_membership_rows(self):
        assert _sel({})["S1__binette_bin1"]["n_members"] == 4
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n snakemake python -m pytest tests/test_pangenome_select.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'pangenome_select'`

- [ ] **Step 3: Write minimal implementation**

```python
#!/usr/bin/env python3
"""Portao do pangenoma: quais clusters de MAG merecem anotacao por membro.

O principio (h) (docs/ROADMAP_SIMPLIFICACAO.md) computa na representante e
herda no membro. Comparar conteudo genico DENTRO de um cluster exige o
contrario, e por isso este portao existe: reintroduz anotacao por membro
apenas onde a pergunta de core/acessorio tem sentido -- dezenas de
genomas, nao todos.
"""

import csv

MIN_MEMBERS = 3
MIN_SYSTEMS = 3

# Ordem = precedencia na hora de nomear o criterio disparado. A ilha e a
# evidencia mais forte, o AMR entra por direito proprio, e o plasmidio NAO
# aparece aqui: e sinal de mobilidade que reforca, nunca criterio isolado.
_CRITERIA = (
    ("ilha",      lambda e: e.get("n_islands", 0) >= 1),
    ("sistemas",  lambda e: e.get("n_systems", 0) >= MIN_SYSTEMS),
    ("amr",       lambda e: e.get("n_args", 0) >= 1),
)


def load_membership(path):
    """{representative_id: [{'source_id','original_bin_id','member_id'}, ...]}"""
    by_rep = {}
    with open(path, newline="") as f:
        for row in csv.DictReader(f, delimiter="\t"):
            rep = (row.get("representative_id") or "").strip()
            member = (row.get("member_id") or "").strip()
            if not rep or not member:
                continue
            by_rep.setdefault(rep, []).append({
                "source_id": (row.get("source_id") or "").strip(),
                "original_bin_id": (row.get("original_bin_id") or "").strip(),
                "member_id": member,
            })
    return by_rep


def load_completeness(path):
    """{genome: completeness} do checkm2_quality_report.tsv do catalogo."""
    out = {}
    try:
        with open(path, newline="") as f:
            for row in csv.DictReader(f, delimiter="\t"):
                name = (row.get("Name") or "").strip()
                if not name:
                    continue
                try:
                    out[name] = float(row.get("Completeness") or 0)
                except ValueError:
                    continue
    except OSError:
        pass
    return out


def select_clusters(membership, evidence, min_members=MIN_MEMBERS):
    """Uma linha por cluster, elegivel ou nao -- a selecao tem de ser
    auditavel. `criterio` nomeia o primeiro criterio satisfeito, ou o
    motivo da recusa."""
    rows = []
    for rep, members in sorted(membership.items()):
        ev = evidence.get(rep, {})
        n_members = len(members)
        criterio, eligible = "", False
        if n_members < min_members:
            criterio = f"poucos membros ({n_members} < {min_members})"
        else:
            for name, test in _CRITERIA:
                if test(ev):
                    criterio, eligible = name, True
                    break
            if not eligible:
                criterio = "sem evidencia de defesa/amr"
        rows.append({
            "representative_id": rep,
            "n_members":  n_members,
            "n_islands":  ev.get("n_islands", 0),
            "n_systems":  ev.get("n_systems", 0),
            "n_args":     ev.get("n_args", 0),
            "n_plasmid":  ev.get("n_plasmid", 0),
            "criterio":   criterio,
            "eligible":   eligible,
        })
    return rows
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n snakemake python -m pytest tests/test_pangenome_select.py -q`
Expected: PASS (9 testes)

- [ ] **Step 5: Commit**

```bash
git add scripts/pangenome_select.py tests/test_pangenome_select.py
git commit -m "feat(pangenoma): portao de selecao de clusters por ilha/defesa/AMR

Plasmidio reforca mas nao elege; 3 membros e piso absoluto; todo cluster
aparece no relatorio da selecao, elegivel ou nao."
```

---

### Task 3: A regra `mag_pangenome_select`

**Files:**
- Create: `rules/pangenome.smk`
- Modify: `Snakefile:395` (adicionar `include:` logo após `rules/annotation.smk`)

**Interfaces:**
- Consumes: `pangenome_select.load_membership/load_completeness/select_clusters` (Task 2), `defense_islands.find_islands` (Task 1)
- Produces:
  - `rules.mag_pangenome_select.output.candidates` → `{MAG_CATALOG_DIR}/pangenome/candidates.tsv`
  - `rules.mag_pangenome_select.output.members` → `{MAG_CATALOG_DIR}/pangenome/members.txt`, um `member_id` por linha (sem extensão) apenas dos clusters elegíveis
  - `rules.mag_pangenome_select.output.done`

- [ ] **Step 1: Criar o módulo com a regra**

`rules/pangenome.smk`. Note que a regra **não tem wildcard nenhum** — é global, como toda `mag_*`.

```python
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
```

- [ ] **Step 2: Registrar o include**

Em `Snakefile`, logo após a linha `include: "rules/annotation.smk"`:

```python
include: "rules/pangenome.smk"
```

Precisa vir **depois** de `annotation.smk` e `defense_amr.smk`? Não: ela referencia `rules.mag_defensefinder`, `rules.mag_amr_consensus` e `rules.mag_abricate`, todas definidas em `defense_amr.smk` (include 396). Portanto **o include tem de vir depois de `rules/defense_amr.smk`** — inserir entre `defense_amr.smk` e `finalize.smk`.

- [ ] **Step 3: Verificar que o Snakefile ainda parseia**

Run: `conda run -n snakemake snakemake -n --configfile /tmp/config_dagtest.yaml --quiet 2>&1 | tail -5`
(usar o `config_dagtest.yaml` gerado como no roadmap: `sed 's|^outdir:.*|outdir: "/tmp/dagout"|' config_amazon_18-08-26.yaml`)
Expected: sem `NameError` / `AttributeError`; a regra aparece no dry-run

- [ ] **Step 4: Confirmar o delta de DAG**

Run: comparar contagens antes/depois com o método do roadmap:
```bash
conda run -n snakemake snakemake -n --forceall --configfile /tmp/config_dagtest.yaml \
  | grep -E "^[a-z_]+ +[0-9]+$" > /tmp/dag_task3.txt
```
Expected: `mag_pangenome_select  1` presente; total = anterior **+1**. Zero jobs por amostra.

- [ ] **Step 5: Commit**

```bash
git add rules/pangenome.smk Snakefile
git commit -m "feat(pangenoma): regra de selecao de clusters, sem checkpoint

Segue o padrao mag_bakta: selecao dependente de dados dentro de um job de
numero fixo, para nao quebrar a invariante de dry-run do roadmap."
```

---

### Task 4: Proteínas dos membros selecionados

**Files:**
- Modify: `rules/pangenome.smk` (acrescentar a regra ao fim)

**Interfaces:**
- Consumes: `rules.mag_pangenome_select.output.members` (Task 3)
- Produces: `rules.mag_pangenome_proteins.output.manifest` → `{PANGENOME_DIR}/proteins/manifest.txt`, no formato `name \t mode \t fna \t faa \t gff` — **o mesmo de `mag_catalog_proteins`**, que é o que permite herdar as regras de defesa/AMR sem tocar no corpo delas.

- [ ] **Step 1: Escrever a regra**

Cópia paramétrica de `mag_catalog_proteins` (`rules/mag_catalog.smk:504`), lendo `members.txt` em vez de varrer `representatives/`.

```python
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
                      "-p single -q >> {log} 2>&1")
                rows.append((name, "bin", fna, faa, gff))

            with open(str(output.manifest), "w") as mf:
                for r in rows:
                    mf.write("\t".join(r) + "\n")
            lf.write(f"[pangenome_proteins] manifesto com {len(rows)} genomas\n")

        write_status(str(output.done),
                     "ok" if rows else "skipped: no eligible members")
```

- [ ] **Step 2: Verificar o formato do manifesto contra o do catálogo**

Run:
```bash
grep -n "mf.write" rules/mag_catalog.smk rules/pangenome.smk
```
Expected: as duas escrevem `"\t".join(...)` de uma tupla de 5 campos na mesma ordem. Se divergirem, a Task 5 quebra em silêncio (o leitor pula linhas com menos de 5 campos).

- [ ] **Step 3: Confirmar o delta de DAG**

Run: dry-run como na Task 3.
Expected: `mag_pangenome_proteins  1`; total anterior **+1**.

- [ ] **Step 4: Commit**

```bash
git add rules/pangenome.smk
git commit -m "feat(pangenoma): prodigal nos membros dos clusters elegiveis"
```

---

### Task 5: Defesa e AMR por membro, herdando as regras globais

**Files:**
- Modify: `rules/pangenome.smk` (acrescentar ao fim)

**Interfaces:**
- Consumes: `rules.mag_pangenome_proteins.output.manifest` (Task 4)
- Produces:
  - `rules.mag_pangenome_defensefinder.output.systems` → `{PANGENOME_DIR}/defensefinder/defensefinder_systems.tsv` (coluna `genome` = `member_id`)
  - `rules.mag_pangenome_amr_consensus.output.consensus` → `{PANGENOME_DIR}/amr_consensus/amr_consensus.tsv` (coluna `locus` prefixada por `{member_id}__`)

- [ ] **Step 1: Herdar as regras**

`use rule ... as ... with:` só substitui `input`/`output`/`log`/`benchmark`; o corpo (`run:`) vem inteiro. Os corpos consomem manifesto e não conhecem `{sample}` — foi para isso que o formato foi preservado na migração de 2026-08-19.

```python
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
```

- [ ] **Step 2: Confirmar os nomes de `input` herdados**

Os corpos herdados referenciam `input.<nome>` literalmente: um nome trocado levanta `AttributeError` **em tempo de execução**, não no dry-run. Os nomes acima já foram conferidos contra `rules/defense_amr.smk` em 2026-08-19 e são estes:

| regra original | linha | nomes de `input` |
|---|---|---|
| `mag_defensefinder` | `defense_amr.smk:63` | `manifest`, `done` |
| `mag_amrfinderplus` | `defense_amr.smk:205` | `manifest`, `done` |
| `mag_rgi_card` | `defense_amr.smk:278` | `manifest`, `done` |
| `mag_deeparg` | `defense_amr.smk:392` | `manifest`, `done` |
| `mag_argnorm_normalize` | `defense_amr.smk:573` | `amrfinder`, `amrfinder_done`, `deeparg`, `deeparg_done` |
| `mag_amr_consensus` | `defense_amr.smk:642` | `argnorm_done`, `rgi_done`, `amrfinder_normed`, `deeparg_normed`, `rgi_results` |

Run, para confirmar que nada mudou desde então:
```bash
for r in mag_defensefinder mag_amrfinderplus mag_rgi_card mag_deeparg \
         mag_argnorm_normalize mag_amr_consensus; do
  echo "=== $r"
  grep -n "rule $r:" -A 22 rules/defense_amr.smk | sed -n '/input:/,/output:/p'
done
```
Expected: os nomes da tabela acima, sem divergência.

- [ ] **Step 3: Confirmar o delta de DAG**

Run: dry-run como antes.
Expected: mais 6 jobs globais (`defensefinder`, `amrfinderplus`, `rgi_card`, `deeparg`, `argnorm`, `amr_consensus`). Total anterior **+6**.

> **Nota sobre o total:** a spec previu "+5 jobs". A contagem real é **+9** (select, proteins, 6 de defesa/AMR, matrix), porque a spec agrupou defesa+AMR como duas linhas de tabela e o AMR são cinco regras encadeadas. O número a usar na verificação final é **+9**; a spec deve ser corrigida na Task 7.

- [ ] **Step 4: Commit**

```bash
git add rules/pangenome.smk
git commit -m "feat(pangenoma): defesa e AMR por membro via use rule

Herdam as regras globais sobre o manifesto dos membros; os corpos
consomem manifesto e nao conhecem wildcard, entao vao inteiros."
```

---

### Task 6: A matriz de três estados

**Files:**
- Create: `scripts/pangenome_matrix.py`
- Modify: `rules/pangenome.smk` (acrescentar a regra ao fim)
- Test: `tests/test_pangenome_matrix.py`

**Interfaces:**
- Consumes: saídas das Tasks 3-5
- Produces:
  - `build_matrix(clusters, members_by_rep, gene_hits, completeness, min_completeness=70.0) -> list[dict]` — uma linha por `(cluster, gene)`, com `states` = `{member_id: 'x'|'.'|'?'}`, `n_present`, `n_evaluable`, `freq`.
  - `summarize_clusters(matrix_rows, members_by_rep, completeness, min_completeness=70.0) -> list[dict]` — uma linha por cluster com `n_members`, `n_members_avaliaveis`, `n_genes_core`, `n_genes_variaveis`, `n_genes_singleton`.
  - `gene_hits` é `{member_id: set[str]}` — nomes de sistema de defesa e de ARG normalizado.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_pangenome_matrix.py
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", "scripts"))

from pangenome_matrix import build_matrix, summarize_clusters

MEMBERS = {"R1": ["m1", "m2", "m3", "m4"]}
COMPLETE = {"m1": 95.0, "m2": 88.0, "m3": 91.0, "m4": 42.0}   # m4 e ruim
HITS = {"m1": {"RM_Type_II", "Gabija"}, "m2": {"RM_Type_II"},
        "m3": {"RM_Type_II", "Gabija"}, "m4": set()}


def _rows():
    return {r["gene"]: r for r in build_matrix(["R1"], MEMBERS, HITS, COMPLETE)}


class TestBuildMatrix:
    def test_incomplete_member_is_unassessable_not_absent(self):
        # O ponto cientifico central: m4 esta 42% completo. Um "." ali
        # afirmaria que o organismo nao tem o gene, quando o que houve foi
        # a regiao nao ter montado.
        assert _rows()["Gabija"]["states"]["m4"] == "?"
        assert _rows()["RM_Type_II"]["states"]["m4"] == "?"

    def test_unassessable_member_is_out_of_the_denominator(self):
        # Gabija em 2 de 3 AVALIAVEIS, nao 2 de 4.
        row = _rows()["Gabija"]
        assert row["n_present"] == 2
        assert row["n_evaluable"] == 3
        assert row["freq"] == "2/3"

    def test_present_and_absent_are_distinguished(self):
        row = _rows()["Gabija"]
        assert row["states"]["m1"] == "x"
        assert row["states"]["m2"] == "."

    def test_gene_in_every_evaluable_member(self):
        row = _rows()["RM_Type_II"]
        assert row["n_present"] == 3
        assert row["freq"] == "3/3"


class TestSummarizeClusters:
    def test_core_uses_90_percent_not_99(self):
        # RM_Type_II esta em 3/3 avaliaveis -> core. Gabija em 2/3 (67%)
        # -> variavel. O limiar de 99% zeraria o core com MAG.
        rows = build_matrix(["R1"], MEMBERS, HITS, COMPLETE)
        summary = {s["representative_id"]: s
                   for s in summarize_clusters(rows, MEMBERS, COMPLETE)}["R1"]
        assert summary["n_genes_core"] == 1
        assert summary["n_genes_variaveis"] == 1

    def test_evaluable_member_count_excludes_the_incomplete_one(self):
        rows = build_matrix(["R1"], MEMBERS, HITS, COMPLETE)
        summary = summarize_clusters(rows, MEMBERS, COMPLETE)[0]
        assert summary["n_members"] == 4
        assert summary["n_members_avaliaveis"] == 3

    def test_singleton_gene_is_counted(self):
        hits = dict(HITS, m1=HITS["m1"] | {"CBASS"})
        rows = build_matrix(["R1"], MEMBERS, hits, COMPLETE)
        summary = summarize_clusters(rows, MEMBERS, COMPLETE)[0]
        assert summary["n_genes_singleton"] == 1

    def test_cluster_with_no_evaluable_member_does_not_divide_by_zero(self):
        rows = build_matrix(["R1"], MEMBERS, HITS,
                            {m: 10.0 for m in MEMBERS["R1"]})
        summary = summarize_clusters(rows, MEMBERS,
                                     {m: 10.0 for m in MEMBERS["R1"]})[0]
        assert summary["n_members_avaliaveis"] == 0
        assert summary["n_genes_core"] == 0
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n snakemake python -m pytest tests/test_pangenome_matrix.py -q`
Expected: FAIL — `ModuleNotFoundError: No module named 'pangenome_matrix'`

- [ ] **Step 3: Write minimal implementation**

```python
#!/usr/bin/env python3
"""Matriz gene x membro dos clusters elegiveis, em TRES estados.

Ausencia em MAG nao e ausencia no organismo. Um "." pode significar "o
organismo nao tem o gene" ou "o MAG esta 74% completo e a regiao nao
montou"; tratar os dois como o mesmo zero produz conclusao errada com
aparencia de dado -- a mesma familia do done.txt vazio lido como zero
biologico (docs/ROADMAP_SIMPLIFICACAO.md).

Por isso: 'x' presente, '.' ausente, '?' NAO AVALIAVEL. E o '?' fica fora
do denominador da frequencia.
"""

MIN_COMPLETENESS = 70.0   # mesmo piso do mag_bakta
CORE_FRACTION = 0.90      # 99% zera o core com MAG (metaFun)


def _evaluable(members, completeness, min_completeness):
    return [m for m in members
            if completeness.get(m, 0.0) >= min_completeness]


def build_matrix(clusters, members_by_rep, gene_hits, completeness,
                 min_completeness=MIN_COMPLETENESS):
    """Uma linha por (cluster, gene)."""
    rows = []
    for rep in clusters:
        members = members_by_rep.get(rep, [])
        evaluable = set(_evaluable(members, completeness, min_completeness))
        genes = sorted({g for m in members for g in gene_hits.get(m, set())})
        for gene in genes:
            states, n_present = {}, 0
            for m in members:
                if m not in evaluable:
                    states[m] = "?"
                elif gene in gene_hits.get(m, set()):
                    states[m] = "x"
                    n_present += 1
                else:
                    states[m] = "."
            rows.append({
                "representative_id": rep,
                "gene": gene,
                "states": states,
                "n_present": n_present,
                "n_evaluable": len(evaluable),
                "freq": f"{n_present}/{len(evaluable)}" if evaluable else "0/0",
            })
    return rows


def summarize_clusters(matrix_rows, members_by_rep, completeness,
                       min_completeness=MIN_COMPLETENESS):
    """Uma linha por cluster: e isto que decide se a fase 2 se justifica."""
    by_rep = {}
    for row in matrix_rows:
        by_rep.setdefault(row["representative_id"], []).append(row)

    summaries = []
    for rep, rows in sorted(by_rep.items()):
        members = members_by_rep.get(rep, [])
        n_eval = len(_evaluable(members, completeness, min_completeness))
        core = variable = singleton = 0
        for row in rows:
            if not n_eval or row["n_present"] == 0:
                continue
            if row["n_present"] >= CORE_FRACTION * n_eval:
                core += 1
            else:
                variable += 1
            if row["n_present"] == 1:
                singleton += 1
        summaries.append({
            "representative_id":    rep,
            "n_members":            len(members),
            "n_members_avaliaveis": n_eval,
            "n_genes_core":         core,
            "n_genes_variaveis":    variable,
            "n_genes_singleton":    singleton,
        })
    return summaries
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n snakemake python -m pytest tests/test_pangenome_matrix.py -q`
Expected: PASS (9 testes)

- [ ] **Step 5: Escrever a regra que alimenta os dois**

Acrescentar a `rules/pangenome.smk`:

```python
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
```

- [ ] **Step 6: Run the full suite**

Run: `conda run -n snakemake python -m pytest tests/ -q`
Expected: PASS, 194 testes (171 + 5 + 9 + 9)

- [ ] **Step 7: Commit**

```bash
git add scripts/pangenome_matrix.py tests/test_pangenome_matrix.py rules/pangenome.smk
git commit -m "feat(pangenoma): matriz gene x membro em tres estados

'?' para membro abaixo de 70% de completude, fora do denominador:
ausencia em MAG incompleto nao e ausencia no organismo."
```

---

### Task 7: Ligar no relatório, verificar e documentar

**Files:**
- Modify: `scripts/report/data_loaders.py:32-58` (`STATUS_TRACKED_GLOBAL_TOOLS`)
- Modify: `docs/superpowers/specs/2026-08-19-pangenoma-clusters-defesa-design.md` (corrigir "+5 jobs" para "+9")
- Modify: `docs/ROADMAP_SIMPLIFICACAO.md`
- Modify: `CLAUDE.md`

**Interfaces:**
- Consumes: todas as saídas anteriores
- Produces: nada consumido por código

- [ ] **Step 1: Rastrear o status das regras novas**

Em `scripts/report/data_loaders.py`, dentro de `STATUS_TRACKED_GLOBAL_TOOLS`:

```python
    # Fase 1 do pangenoma (2026-08-19). "Nenhum cluster elegivel" e
    # resultado legitimo -- tem de aparecer como skipped:, jamais como
    # painel vazio, que e indistinguivel de "a ferramenta quebrou".
    "mag_pangenome_select":       "mag_catalog/pangenome/select_done.txt",
    "mag_pangenome_proteins":     "mag_catalog/pangenome/proteins/done.txt",
    "mag_pangenome_defensefinder": "mag_catalog/pangenome/defensefinder/done.txt",
    "mag_pangenome_amr_consensus": "mag_catalog/pangenome/amr_consensus/done.txt",
    "mag_pangenome_matrix":       "mag_catalog/pangenome/done.txt",
```

- [ ] **Step 2: Verificação final de DAG**

```bash
SP=/tmp/pangenome_dag
mkdir -p $SP
sed 's|^outdir:.*|outdir: "'$SP'/dagout"|' config_amazon_18-08-26.yaml > $SP/config_dagtest.yaml
git stash
conda run -n snakemake snakemake -n --forceall --configfile $SP/config_dagtest.yaml \
  | grep -E "^[a-z_]+ +[0-9]+$" > $SP/dag_antes.txt
git stash pop
conda run -n snakemake snakemake -n --forceall --configfile $SP/config_dagtest.yaml \
  | grep -E "^[a-z_]+ +[0-9]+$" > $SP/dag_depois.txt
diff $SP/dag_antes.txt $SP/dag_depois.txt
```

Expected: exatamente 9 linhas novas, todas com contagem **1**, e o `total` subindo em 9. **Nenhuma regra existente muda de contagem.** Se alguma mudar, a herança com `use rule` roubou input de outra regra — parar e investigar.

- [ ] **Step 3: Verificação de env/container**

Run: `conda run -n snakemake python scripts/check_env_container_sync.py`
Expected: o nº de pares comparados sobe (as regras herdadas declaram conda+container); **nenhuma divergência nova**. As duas do prodigal são anteriores.

- [ ] **Step 4: Corrigir o número na spec**

Em `docs/superpowers/specs/2026-08-19-pangenoma-clusters-defesa-design.md`, seção "6. Verificação": trocar `**+5 jobs globais fixos**` por `**+9 jobs globais fixos**` e explicar que o AMR são cinco regras encadeadas, não uma.

- [ ] **Step 5: Registrar no roadmap**

Acrescentar a `docs/ROADMAP_SIMPLIFICACAO.md`, antes da lista de "Ideias ainda não avaliadas", uma seção com: o portão e seus critérios; o número real de clusters elegíveis medido nos dados da Amazônia; o veredito sobre a fase 2 (com o número de clusters com ≥5 membros avaliáveis); e a nota de que `compute_defense_islands` deixou de ser exclusiva do relatório.

- [ ] **Step 6: Atualizar o CLAUDE.md**

Acrescentar `rules/pangenome.smk` à lista de módulos em `## Project Structure`, com uma frase sobre o portão, o padrão `mag_bakta` (nunca checkpoint) e os três estados da matriz.

- [ ] **Step 7: Commit**

```bash
git add scripts/report/data_loaders.py docs/ CLAUDE.md
git commit -m "docs(pangenoma): registra a fase 1 e rastreia o status das regras"
```

---

## Critério de saída da fase 1

Depois de rodar nos dados reais, `cluster_summary.tsv` responde se a fase 2 se justifica:

- **≥ 1 cluster com ≥ 5 membros avaliáveis E `n_genes_variaveis > 0`** → escrever o plano da fase 2 (PPanGGOLiN, `-K 3`, core 90%).
- **Caso contrário** → a matriz é a resposta, e a fase 2 **não deve ser construída**. Registrar o número medido no roadmap para que a decisão não seja reaberta por intuição.
