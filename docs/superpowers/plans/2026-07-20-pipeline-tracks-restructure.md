# Pipeline Tracks Restructure (Plano 1) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Reestruturar o VAPOR para seleção de tracks (reads/viral/prok) por configuração, com integração vírus↔bactéria condicional, VAMB removido dos binners per-sample, report adaptativo e entry/exit points no CLI — mantendo comportamento idêntico ao atual quando todas as tracks estão ligadas.

**Architecture:** Um Snakefile único com `include:` (sem sistema de módulo). A lógica de resolução/validação de config vira um módulo Python puro testável (`pipeline_config.py`). O `rule all` deixa de ser um bloco chapado e passa a ser composto por funções-alvo (`_t_*()`) gated por flags derivadas. VAMB sai do binning per-sample. O report ganha um `TRACKS` JS constant que esconde abas de tracks desligadas.

**Tech Stack:** Snakemake (Python DSL), Python 3.11, pytest para helpers puros, `snakemake -n` (dry-run) para composição de DAG, ECharts/D3 (report JS).

## Global Constraints

- **Commits NÃO devem incluir linhas de co-autoria do Claude** (`Co-Authored-By: Claude ...`). Regra explícita do usuário ("sem dar autoria pra claude"). Vale para TODO commit deste plano.
- `config.yaml` é a fonte canônica; overrides via `vapor --set KEY=VALUE` (mecanismo `--config` do Snakemake). Não quebrar essa convenção.
- Defaults que preservam o comportamento atual: `tracks.reads=false`, `tracks.viral=true`, `tracks.prok=true`, `use_host_defense=true`, `coassembly.enabled=false`, `cobinning_multisplit=false`.
- Execução de co-assembly/co-binning é do **Plano 2**. Neste plano, config de co-assembly é apenas **parseada e validada**; os alvos de grupo são stubs que retornam `[]`.
- HTML write usa `encoding='utf-8'` (Windows). Não remover.

---

## File Structure

| Arquivo | Responsabilidade | Ação |
|---|---|---|
| `pipeline_config.py` (repo root) | Funções puras: resolver flags derivadas + validar config | **Criar** |
| `Snakefile` | Importar/chamar `pipeline_config`; atribuir vars; `rule all` por funções-alvo | Modificar |
| `rules/prok_binning.smk` | Remover VAMB do binning per-sample (scaffold2bin + regra vamb) | Modificar |
| `vapor.py` | Novos args `--track` / `--until` + helper de tradução | Modificar |
| `rules/report.smk` | Passar `tracks` como param de `generate_report` | Modificar |
| `scripts/generate_report.py` + `scripts/report/renderer.py` | Injetar `const TRACKS` no HTML | Modificar |
| `scripts/report/components/app.js` | Esconder abas de tracks desligadas | Modificar |
| `tests/conftest.py` | Inserir repo root no `sys.path` | **Criar** |
| `tests/test_pipeline_config.py` | Testes de resolução/validação | **Criar** |
| `tests/test_vapor_cli.py` | Testes de tradução `--track`/`--until` | **Criar** |
| `pyproject.toml` | Adicionar `pipeline_config` a `py-modules` | Modificar |

---

### Task 1: Módulo de config testável (`pipeline_config.py`)

**Files:**
- Create: `pipeline_config.py`
- Create: `tests/conftest.py`
- Create: `tests/test_pipeline_config.py`
- Modify: `pyproject.toml` (py-modules)

**Interfaces:**
- Produces:
  - `resolve_pipeline_config(config: dict) -> dict` — retorna dict com chaves derivadas:
    `track_reads, track_viral, track_prok, use_host_defense, integration_enabled,
    coassembly_enabled, coassembly_grouping, coassembly_viral, coassembly_binning,
    cobinning_multisplit` (todas bool, exceto `coassembly_grouping: str`).
  - `validate_pipeline_config(config: dict) -> None` — levanta `ValueError` com mensagem clara em config inválida.

- [ ] **Step 1: Write the failing test**

Create `tests/conftest.py`:

```python
import sys
from pathlib import Path

# Repo root (parent of tests/) on sys.path so `import pipeline_config` / `import vapor` work.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
```

Create `tests/test_pipeline_config.py`:

```python
import pytest
from pipeline_config import resolve_pipeline_config, validate_pipeline_config


def test_defaults_preserve_current_behavior():
    r = resolve_pipeline_config({})
    assert r["track_reads"] is False
    assert r["track_viral"] is True
    assert r["track_prok"] is True
    assert r["use_host_defense"] is True
    assert r["integration_enabled"] is True   # viral and prok and host_defense
    assert r["coassembly_enabled"] is False
    assert r["cobinning_multisplit"] is False


def test_integration_requires_both_tracks_and_flag():
    assert resolve_pipeline_config(
        {"tracks": {"viral": True, "prok": False}})["integration_enabled"] is False
    assert resolve_pipeline_config(
        {"tracks": {"viral": True, "prok": True}, "use_host_defense": False}
    )["integration_enabled"] is False


def test_grouping_none_disables_coassembly():
    r = resolve_pipeline_config(
        {"coassembly": {"enabled": True, "grouping": "none"}})
    assert r["coassembly_enabled"] is False


def test_reads_only_track():
    r = resolve_pipeline_config(
        {"tracks": {"reads": True, "viral": False, "prok": False}})
    assert r["track_reads"] is True
    assert r["track_viral"] is False
    assert r["integration_enabled"] is False


def test_validate_rejects_no_analysis():
    with pytest.raises(ValueError, match="pelo menos uma"):
        validate_pipeline_config(
            {"tracks": {"reads": False, "viral": False, "prok": False}})


def test_validate_rejects_bad_grouping():
    with pytest.raises(ValueError, match="grouping"):
        validate_pipeline_config(
            {"coassembly": {"enabled": True, "grouping": "wrong"}})


def test_validate_accepts_defaults():
    validate_pipeline_config({})   # não levanta
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_pipeline_config.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'pipeline_config'`

- [ ] **Step 3: Write minimal implementation**

Create `pipeline_config.py`:

```python
"""Pure config resolution + validation for the VAPOR Snakefile.

Kept import-free of Snakemake so it can be unit-tested with pytest and imported
from the Snakefile alike.
"""

_VALID_GROUPING = ("metadata", "all", "none")


def _b(d, key, default):
    """Bool coercion tolerant of YAML strings ('false'/'true')."""
    v = d.get(key, default)
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1", "yes", "on")
    return bool(v)


def resolve_pipeline_config(config: dict) -> dict:
    tracks = config.get("tracks", {}) or {}
    coas = config.get("coassembly", {}) or {}

    track_reads = _b(tracks, "reads", False)
    track_viral = _b(tracks, "viral", True)
    track_prok = _b(tracks, "prok", True)
    use_host_defense = _b(config, "use_host_defense", True)

    grouping = str(coas.get("grouping", "metadata")).strip().lower()
    coassembly_enabled = _b(coas, "enabled", False) and grouping != "none"

    return {
        "track_reads": track_reads,
        "track_viral": track_viral,
        "track_prok": track_prok,
        "use_host_defense": use_host_defense,
        "integration_enabled": track_viral and track_prok and use_host_defense,
        "coassembly_enabled": coassembly_enabled,
        "coassembly_grouping": grouping,
        "coassembly_viral": _b(coas, "viral", True),
        "coassembly_binning": _b(coas, "binning", True),
        "cobinning_multisplit": _b(config, "cobinning_multisplit", False),
    }


def validate_pipeline_config(config: dict) -> None:
    r = resolve_pipeline_config(config)

    grouping = r["coassembly_grouping"]
    if grouping not in _VALID_GROUPING:
        raise ValueError(
            f"coassembly.grouping inválido: '{grouping}'. "
            f"Use um de: {', '.join(_VALID_GROUPING)}."
        )

    any_analysis = (
        r["track_reads"] or r["track_viral"] or r["track_prok"]
        or r["coassembly_enabled"] or r["cobinning_multisplit"]
    )
    if not any_analysis:
        raise ValueError(
            "Nenhuma análise habilitada: ligue pelo menos uma de "
            "tracks.reads/viral/prok, coassembly.enabled ou cobinning_multisplit."
        )
```

Modify `pyproject.toml` — add `pipeline_config` to `py-modules`:

```toml
[tool.setuptools]
py-modules = ["vapor", "pipeline_config"]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_pipeline_config.py -v`
Expected: PASS (7 passed)

- [ ] **Step 5: Commit**

```bash
git add pipeline_config.py tests/conftest.py tests/test_pipeline_config.py pyproject.toml
git commit -m "feat: pure config resolution + validation module for pipeline tracks"
```

---

### Task 2: Wire `pipeline_config` into Snakefile

**Files:**
- Modify: `Snakefile` (imports region ~line 51; flags region ~line 222)

**Interfaces:**
- Consumes: `resolve_pipeline_config`, `validate_pipeline_config` from Task 1.
- Produces: módulo-level vars usadas pela `rule all` na Task 3:
  `TRACK_READS, TRACK_VIRAL, TRACK_PROK, USE_HOST_DEFENSE, INTEGRATION_ENABLED,
  COASSEMBLY_ENABLED, COASSEMBLY_GROUPING, COASSEMBLY_VIRAL, COASSEMBLY_BINNING,
  COBINNING_MULTISPLIT` (mantém também `READS_CLASSIFY_ENABLED` existente).

- [ ] **Step 1: Add import at top of Snakefile**

In `Snakefile`, after the existing `from pathlib import Path` (~line 51), add:

```python
import sys as _sys
_sys.path.insert(0, PIPELINE_DIR)
from pipeline_config import resolve_pipeline_config, validate_pipeline_config
```

(`PIPELINE_DIR` is already defined right below as `workflow.basedir`; move these two
lines to *after* the `PIPELINE_DIR = workflow.basedir` assignment — insert them
immediately following the `SCRIPTS_DIR = ...` line.)

- [ ] **Step 2: Assign derived vars + validate**

In `Snakefile`, replace the existing line (~222):

```python
# ── Reads-only classification (Sylph) ─────────────────────────────────
READS_CLASSIFY_ENABLED = config.get("reads_classify", False)
```

with:

```python
# ── Reads-only classification (Sylph) ─────────────────────────────────
READS_CLASSIFY_ENABLED = config.get("reads_classify", False)

# ── Pipeline tracks + co-assembly (resolved + validated) ──────────────
validate_pipeline_config(config)
_PCFG = resolve_pipeline_config(config)
TRACK_READS          = _PCFG["track_reads"]
TRACK_VIRAL          = _PCFG["track_viral"]
TRACK_PROK           = _PCFG["track_prok"]
USE_HOST_DEFENSE     = _PCFG["use_host_defense"]
INTEGRATION_ENABLED  = _PCFG["integration_enabled"]
COASSEMBLY_ENABLED   = _PCFG["coassembly_enabled"]
COASSEMBLY_GROUPING  = _PCFG["coassembly_grouping"]
COASSEMBLY_VIRAL     = _PCFG["coassembly_viral"]
COASSEMBLY_BINNING   = _PCFG["coassembly_binning"]
COBINNING_MULTISPLIT = _PCFG["cobinning_multisplit"]
```

- [ ] **Step 3: Add config block to `config.yaml`**

In `config.yaml`, add near the `reads_classify:` key:

```yaml
# ── Seleção de tracks ──────────────────────────
tracks:
  reads:  false
  viral:  true
  prok:   true

use_host_defense: true

# ── Co-assembly / co-binning (opt-in; execução no Plano 2) ──
coassembly:
  enabled:  false
  grouping: metadata     # metadata | all | none
  viral:    true
  binning:  true
cobinning_multisplit: false
```

- [ ] **Step 4: Verify Snakefile parses (dry-run, default config)**

Run: `conda run -n snakemake snakemake -n --configfile config.yaml --cores 1 2>&1 | tail -20`
Expected: dry-run completes DAG build with no Python error (job list printed).
(If databases aren't present, ignore missing-input errors unrelated to config parsing;
the critical check is NO `ValueError` / `NameError` from config resolution.)

- [ ] **Step 5: Verify validation fires**

Run: `conda run -n snakemake snakemake -n --configfile config.yaml --config tracks='{"reads": false, "viral": false, "prok": false}' --cores 1 2>&1 | grep -i "Nenhuma análise"`
Expected: prints the "Nenhuma análise habilitada" ValueError message.

- [ ] **Step 6: Commit**

```bash
git add Snakefile config.yaml
git commit -m "feat: wire tracks/coassembly config resolution into Snakefile"
```

---

### Task 3: Refatorar `rule all` em funções-alvo

**Files:**
- Modify: `Snakefile` (`rule all`, ~lines 333–445)

**Interfaces:**
- Consumes: vars da Task 2.
- Produces: `rule all` gated por track. Alvos de co-assembly são stub `[]` (Plano 2).

- [ ] **Step 1: Add target-builder functions above `rule all`**

In `Snakefile`, immediately before `rule all:` (~line 333), insert:

```python
# ══════════════════════════════════════════════════════════════════════
#  TARGET BUILDERS — cada track contribui sua fatia de alvos ao rule all.
#  Includes permanecem; só o rule all decide o que efetivamente roda.
# ══════════════════════════════════════════════════════════════════════

def _t_foundation():
    t = []
    if not LONG_READS:
        t += expand(f"{OUTDIR}/{{sample}}/qc_raw/done.txt", sample=SAMPLES)
    else:
        t += expand(f"{OUTDIR}/{{sample}}/qc_lr/nanoplot_done.txt", sample=SAMPLES)
    if TRACK_VIRAL or TRACK_PROK:
        t += expand(f"{OUTDIR}/{{sample}}/mmseqs/{{sample}}_rep_seq.fasta", sample=SAMPLES)
        t += expand(f"{OUTDIR}/{{sample}}/quast/report.tsv", sample=SAMPLES)
        if LONG_READS:
            t += expand(f"{OUTDIR}/{{sample}}/assembly/lr/merged/done.txt", sample=SAMPLES)
        # mapping/depth feeds both viral (vRhyme coverage) and prok (binning depth)
        t += expand(f"{OUTDIR}/{{sample}}/mapping/{{sample}}.sorted.bam", sample=SAMPLES)
        t += expand(f"{OUTDIR}/{{sample}}/mapping/{{sample}}_depth.txt", sample=SAMPLES)
        if not LONG_READS:
            t += expand(f"{OUTDIR}/{{sample}}/mapping/bwa_mem_done.txt", sample=SAMPLES)
    return t


def _t_reads():
    # reads track roda se explicitamente selecionada OU se reads_classify legado ligado
    if not (TRACK_READS or READS_CLASSIFY_ENABLED):
        return []
    return [
        f"{OUTDIR}/reads_classify/reads_classify_done.txt",
        f"{OUTDIR}/final/reads_classify/done.txt",
    ]


def _t_viral():
    if not TRACK_VIRAL:
        return []
    t = []
    t += expand(f"{OUTDIR}/{{sample}}/viral/virsorter2/final-viral-combined.fa", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/genomad/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/vibrant/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_viral_consensus.fasta", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_tool_support.tsv", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/checkv/quality_summary.tsv", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/bins/vrhyme/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/checkv_vrhyme/quality_summary.tsv", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/consensus/{{sample}}_viral_nonredundant.fasta", sample=SAMPLES)
    if VOTU_CLUSTERING_ENABLED:
        t += expand(f"{OUTDIR}/{{sample}}/viral/votu/vOTU_clusters.tsv", sample=SAMPLES)
        t += expand(f"{OUTDIR}/{{sample}}/viral/votu/votu_all_reps.fasta", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/votu/{{sample}}_vOTU_table.tsv", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/votu/{{sample}}_vOTU_abundance.tsv", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/taxonomy/taxonomy_done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/viral/vcontact3/done.txt", sample=SAMPLES)
    if DEFENSE_AMR_VIRAL_ENABLED:
        t += expand(f"{OUTDIR}/{{sample}}/viral/defensefinder/done.txt", sample=SAMPLES)
        t += expand(f"{OUTDIR}/{{sample}}/viral/dbapis/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/abundance/viral_abundance.tsv", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/annotation/pharokka/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/annotation/phold/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/annotation/genome_maps/phage_maps_done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/annotation/genome_maps/virus_maps_done.txt", sample=SAMPLES)
    return t


def _t_prok():
    if not TRACK_PROK:
        return []
    t = []
    t += expand(f"{OUTDIR}/{{sample}}/bins/metabat2/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/bins/semibin2/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/bins/comebin/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/bins/binette/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/bins/checkm2/quality_report.tsv", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/bins/gtdbtk/done.txt", sample=SAMPLES)
    if GUNC_ENABLED:
        t += expand(f"{OUTDIR}/{{sample}}/bins/gunc/GUNC.progenomes_2.1.maxCSS_level.tsv", sample=SAMPLES)
    if MAG_DEREP_ENABLED:
        t += expand(f"{OUTDIR}/{{sample}}/bins/derep/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/bins/mmseqs_taxonomy_prok/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/bins/proteins/done.txt", sample=SAMPLES)
    if DEFENSE_AMR_ENABLED:
        t += expand(f"{OUTDIR}/{{sample}}/bins/defensefinder/done.txt", sample=SAMPLES)
        t += expand(f"{OUTDIR}/{{sample}}/bins/amrfinderplus/done.txt", sample=SAMPLES)
        t += expand(f"{OUTDIR}/{{sample}}/bins/rgi/done.txt", sample=SAMPLES)
        t += expand(f"{OUTDIR}/{{sample}}/bins/deeparg/done.txt", sample=SAMPLES)
        if AMR_CONSENSUS_ENABLED:
            t += expand(f"{OUTDIR}/{{sample}}/bins/amr_consensus/done.txt", sample=SAMPLES)
    if ABRICATE_ENABLED:
        t += expand(f"{OUTDIR}/{{sample}}/bins/abricate/done.txt", sample=SAMPLES)
    if ARGNORM_ENABLED:
        t += expand(f"{OUTDIR}/{{sample}}/bins/argnorm/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/abundance/prok_abundance.tsv", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/annotation/bakta/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/annotation/eggnog/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/annotation/kegg_decoder/done.txt", sample=SAMPLES)
    t += expand(f"{OUTDIR}/{{sample}}/annotation/genome_maps/prok_maps_done.txt", sample=SAMPLES)
    return t


def _t_integration():
    if not INTEGRATION_ENABLED:
        return []
    return expand(f"{OUTDIR}/{{sample}}/viral/phist/done.txt", sample=SAMPLES)


def _t_coassembly():
    # Execução de co-assembly/co-binning é do Plano 2. Stub por ora.
    return []


def _t_report():
    t = [
        f"{OUTDIR}/diversity/diversity_done.txt",
        expand(f"{OUTDIR}/{{sample}}/final/done.txt", sample=SAMPLES),
        f"{OUTDIR}/benchmarks/pipeline_timing_summary.tsv",
        f"{OUTDIR}/report.html",
    ]
    if not LONG_READS:
        t.append(f"{OUTDIR}/multiqc_report/multiqc_report.html")
    # flatten (expand returns lists)
    flat = []
    for x in t:
        flat += x if isinstance(x, list) else [x]
    return flat
```

- [ ] **Step 2: Replace `rule all` input with the builders**

In `Snakefile`, replace the entire `rule all:` `input:` block (from `rule all:` through the
closing of its `input:` list, ~lines 333–445) with:

```python
rule all:
    input:
        *_t_foundation(),
        *_t_reads(),
        *_t_viral(),
        *_t_prok(),
        *_t_integration(),
        *_t_coassembly(),
        *_t_report(),
```

- [ ] **Step 3: Dry-run — full pipeline (all tracks) still resolves**

Run: `conda run -n snakemake snakemake -n --configfile config.yaml --cores 1 2>&1 | grep -E "report.html|phist|checkm2|virsorter2" | head`
Expected: alvos de viral (virsorter2), prok (checkm2) e integração (phist) e report.html aparecem no plano.

- [ ] **Step 4: Dry-run — viral-only omits prok + integration**

Run: `conda run -n snakemake snakemake -n --configfile config.yaml --config tracks='{"reads": false, "viral": true, "prok": false}' --cores 1 2>&1 > /tmp/vapor_viral_only.txt; grep -c "checkm2" /tmp/vapor_viral_only.txt; grep -c "phist" /tmp/vapor_viral_only.txt; grep -c "virsorter2" /tmp/vapor_viral_only.txt`
Expected: `checkm2` → 0, `phist` → 0 (integração desligada sem prok), `virsorter2` → >0.

- [ ] **Step 5: Dry-run — reads-only omits assembly**

Run: `conda run -n snakemake snakemake -n --configfile config.yaml --config tracks='{"reads": true, "viral": false, "prok": false}' reads_classify=true --cores 1 2>&1 > /tmp/vapor_reads_only.txt; grep -c "rep_seq" /tmp/vapor_reads_only.txt; grep -c "reads_classify_done" /tmp/vapor_reads_only.txt`
Expected: `rep_seq` → 0 (sem assembly), `reads_classify_done` → >0.

- [ ] **Step 6: Commit**

```bash
git add Snakefile
git commit -m "refactor: compose rule all from per-track target builders"
```

---

### Task 4: Remover VAMB dos binners per-sample

**Files:**
- Modify: `rules/prok_binning.smk` (regra `vamb` ~192–248; `prepare_scaffold2bin` input ~394–398 e bloco vamb ~427–433; docstrings ~299–300)

**Interfaces:**
- Consumes: nada novo.
- Produces: binning per-sample com 3 binners (MetaBAT2, SemiBin2, COMEBin). A regra
  `vamb` deixa de existir no per-sample (será recriada como co-binning no Plano 2).

- [ ] **Step 1: Remove `vamb` from `prepare_scaffold2bin` input**

In `rules/prok_binning.smk`, in `rule prepare_scaffold2bin` `input:` (~394–398), delete the line:

```python
        vamb    = rules.vamb.output.done,
```

so the input becomes:

```python
    input:
        mb2     = rules.metabat2.output.done,
        sb2     = rules.semibin2.output.done,
        comebin = rules.comebin.output.done,
```

- [ ] **Step 2: Remove the VAMB scaffold2bin conversion block**

In the same rule's `run:` body (~427–433), delete:

```python
        vamb_clusters = f"{params.s}/bins/vamb/clusters.tsv"
        if os.path.exists(vamb_clusters):
            with open(vamb_clusters) as fin, open(f"{outdir}/vamb_s2b.tsv", "w") as fout:
                for line in fin:
                    parts = line.strip().split("\t")
                    if len(parts) == 2 and not parts[0].startswith("clusterid"):
                        fout.write(f"{parts[1]}\t{parts[0]}\n")
```

- [ ] **Step 3: Delete the entire `rule vamb`**

In `rules/prok_binning.smk`, delete `rule vamb:` and its whole body (~192–248, from
`rule vamb:` up to but not including the next `rule ...`). Also update the header
comment listing binners (~9) and the comebin docstring "alongside MetaBAT2 + VAMB +
SemiBin2" (~299) to drop VAMB:

Change (~299–300):

```python
    Added as 4th binner alongside MetaBAT2 + VAMB + SemiBin2.
    Binette uses all 4 scaffold2bin files for refined bin selection.
```

to:

```python
    3rd per-sample binner alongside MetaBAT2 + SemiBin2 (VAMB moved to co-binning).
    Binette uses the 3 scaffold2bin files for refined bin selection.
```

- [ ] **Step 4: Remove the stray `vamb` target from rule all**

The Task 3 refactor already omits `bins/vamb/done.txt` from `_t_prok()`. Confirm no other
reference remains:

Run: `grep -rn "bins/vamb\|rules.vamb\|rule vamb" Snakefile rules/`
Expected: no matches (empty output).

- [ ] **Step 5: Dry-run — binette resolves with 3 binners, no vamb**

Run: `conda run -n snakemake snakemake -n --configfile config.yaml --cores 1 2>&1 > /tmp/vapor_novamb.txt; grep -c "bins/vamb" /tmp/vapor_novamb.txt; grep -c "prepare_scaffold2bin\|scaffold2bin" /tmp/vapor_novamb.txt`
Expected: `bins/vamb` → 0; scaffold2bin still present (>0).

- [ ] **Step 6: Commit**

```bash
git add Snakefile rules/prok_binning.smk
git commit -m "refactor: remove VAMB from per-sample binning (reserved for co-binning)"
```

---

### Task 5: CLI `--track` e `--until` (`vapor.py`)

**Files:**
- Modify: `vapor.py` (argparse ~283; `build_command` ~207–223)
- Create: `tests/test_vapor_cli.py`

**Interfaces:**
- Produces: helper `_track_overrides(track_csv: str) -> list[str]` — traduz `"viral,prok"`
  em `["tracks={...}"]` para `--config`. E `--until STAGE` mapeado via
  `_STAGE_ALIASES` para uma regra representativa passada ao `--until` do Snakemake.

- [ ] **Step 1: Write the failing test**

Create `tests/test_vapor_cli.py`:

```python
import json
import vapor


def test_track_overrides_single():
    out = vapor._track_overrides("viral")
    assert len(out) == 1
    key, _, val = out[0].partition("=")
    assert key == "tracks"
    d = json.loads(val)
    assert d == {"reads": False, "viral": True, "prok": False}


def test_track_overrides_multi():
    d = json.loads(vapor._track_overrides("viral,prok")[0].split("=", 1)[1])
    assert d == {"reads": False, "viral": True, "prok": True}


def test_track_overrides_reads():
    d = json.loads(vapor._track_overrides("reads")[0].split("=", 1)[1])
    assert d["reads"] is True and d["viral"] is False


def test_stage_alias_maps_to_rule():
    assert vapor._STAGE_ALIASES["assembly"] == "mmseqs2"
    assert vapor._STAGE_ALIASES["qc"] == "fastp"
    assert vapor._STAGE_ALIASES["viral"] == "viral_consensus"
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_vapor_cli.py -v`
Expected: FAIL — `AttributeError: module 'vapor' has no attribute '_track_overrides'`

- [ ] **Step 3: Add helpers + args + wiring in `vapor.py`**

In `vapor.py`, add near the top-level helpers (after `_PATH_KEYS`, ~line 76):

```python
import json as _json

_VALID_TRACKS = ("reads", "viral", "prok")

# Friendly stage name → representative rule for `--until`.
_STAGE_ALIASES = {
    "qc":       "fastp",
    "assembly": "mmseqs2",
    "viral":    "viral_consensus",
    "binning":  "binette",
}


def _track_overrides(track_csv: str) -> list[str]:
    """Translate --track viral,prok into a Snakemake --config tracks={...} override."""
    wanted = {t.strip().lower() for t in track_csv.split(",") if t.strip()}
    unknown = wanted - set(_VALID_TRACKS)
    if unknown:
        raise SystemExit(
            f"ERROR: unknown track(s): {', '.join(sorted(unknown))}. "
            f"Valid: {', '.join(_VALID_TRACKS)}."
        )
    d = {t: (t in wanted) for t in _VALID_TRACKS}
    return [f"tracks={_json.dumps(d)}"]
```

Add the arguments in `main()` after the `--singularity-args` block (~line 282):

```python
    parser.add_argument(
        "--track",
        metavar="LIST",
        default=None,
        help="Run only these tracks (comma-separated): reads,viral,prok. "
             "Overrides config.yaml tracks. Example: --track viral,prok",
    )
    parser.add_argument(
        "--until",
        metavar="STAGE",
        default=None,
        help="Stop after a stage: qc | assembly | viral | binning "
             "(maps to Snakemake --until on a representative rule).",
    )
```

In `build_command`, extend the `--set` handling (~207–209). Replace:

```python
    # --set KEY=VALUE overrides passed to Snakemake's --config mechanism.
    if args.set_config:
        cmd += ["--config"] + args.set_config
```

with:

```python
    # --set KEY=VALUE + --track overrides, both via Snakemake's --config mechanism.
    config_overrides = list(args.set_config)
    if getattr(args, "track", None):
        config_overrides += _track_overrides(args.track)
    if config_overrides:
        cmd += ["--config"] + config_overrides

    if getattr(args, "until", None):
        stage = args.until.strip().lower()
        if stage not in _STAGE_ALIASES:
            raise SystemExit(
                f"ERROR: unknown --until stage '{stage}'. "
                f"Valid: {', '.join(_STAGE_ALIASES)}."
            )
        cmd += ["--until", _STAGE_ALIASES[stage]]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_vapor_cli.py -v`
Expected: PASS (4 passed)

- [ ] **Step 5: Verify CLI help lists new flags**

Run: `python vapor.py -h 2>&1 | grep -E "\-\-track|\-\-until"`
Expected: both `--track` and `--until` appear.

- [ ] **Step 6: Commit**

```bash
git add vapor.py tests/test_vapor_cli.py
git commit -m "feat: vapor CLI --track and --until entry points"
```

---

### Task 6: Report adaptativo — esconder abas de tracks desligadas

**Files:**
- Modify: `rules/report.smk` (`generate_report` params ~51–69)
- Modify: `scripts/generate_report.py` (params binding) and/or `scripts/report/renderer.py` (JS constant injection)
- Modify: `scripts/report/components/app.js` (nav tab visibility ~191–210)

**Interfaces:**
- Consumes: vars de track do Snakefile via `params`.
- Produces: `const TRACKS = {reads, viral, prok, integration, coassembly}` no HTML; app.js
  esconde `.nav-tab` cujo track está `false`.

- [ ] **Step 1: Pass tracks as a param to generate_report**

In `rules/report.smk`, in `rule generate_report:` `params:` (~51), add after `low_depth_mode`:

```python
        tracks = {
            "reads":       bool(TRACK_READS or READS_CLASSIFY_ENABLED),
            "viral":       bool(TRACK_VIRAL),
            "prok":        bool(TRACK_PROK),
            "integration": bool(INTEGRATION_ENABLED),
            "coassembly":  bool(COASSEMBLY_ENABLED),
        },
```

- [ ] **Step 2: Inject TRACKS constant in the report data script**

In `scripts/report/renderer.py`, locate the `data_dict` passed to `_build_data_script`
(the dict containing `"READS_CLASSIFY": reads_classify_data,` ~line 432). Add an entry:

```python
        "TRACKS": tracks_param,
```

And near the top of `render_report` where snakemake params are read, add (matching the
existing param-reading style):

```python
    tracks_param = dict(getattr(snakemake.params, "tracks", {}) or {})
    tracks_param.setdefault("reads", True)
    tracks_param.setdefault("viral", True)
    tracks_param.setdefault("prok", True)
    tracks_param.setdefault("integration", True)
    tracks_param.setdefault("coassembly", False)
```

(If `render_report` does not receive `snakemake` directly, thread `tracks_param` from
`scripts/generate_report.py`, which binds `snakemake`, into the render call — mirror how
`cfg_params`/other params already flow.)

- [ ] **Step 3: Hide disabled-track nav tabs in app.js**

In `scripts/report/components/app.js`, after the tabs are initialized (near the
`const tabs = document.querySelectorAll('.nav-tab')` / showTab wiring, ~191), add:

```javascript
  // Hide nav tabs whose track did not run (TRACKS injected from the report data).
  (function hideDisabledTrackTabs() {
    if (typeof TRACKS === 'undefined' || !TRACKS) return;
    const TAB_TRACK = {
      viral: 'viral',
      prokaryotic: 'prok',
      hostdefense: 'integration',
      'reads-classify': 'reads',
    };
    document.querySelectorAll('.nav-tab').forEach(btn => {
      const key = TAB_TRACK[btn.dataset.tab];
      if (key && TRACKS[key] === false) btn.style.display = 'none';
    });
  })();
```

- [ ] **Step 4: Regenerate report locally and verify tabs hidden**

Since a full pipeline run is heavy, verify by generating from existing outputs on the
server, or assert the constant is present in a produced report. Minimal check — grep the
generated HTML for the constant and the hide logic in the bundled JS:

Run (after any report regeneration): `grep -c "const TRACKS" <OUTDIR>/report.html; grep -c "hideDisabledTrackTabs" <OUTDIR>/report.html`
Expected: both → 1.

If no report is regenerated in this environment, at minimum confirm the JS source bundles
correctly:

Run: `grep -c "hideDisabledTrackTabs" scripts/report/components/app.js`
Expected: 1.

- [ ] **Step 5: Commit**

```bash
git add rules/report.smk scripts/report/renderer.py scripts/generate_report.py scripts/report/components/app.js
git commit -m "feat: adaptive report — hide nav tabs for disabled tracks"
```

---

## Self-Review

**Spec coverage:**
- Config `tracks` block + derived vars + validação → Tasks 1–2 ✓
- `rule all` por funções-alvo → Task 3 ✓
- `use_host_defense` + integração condicional → Task 3 (`_t_integration` via `INTEGRATION_ENABLED`) ✓
- VAMB removido do per-sample; Binette consolida 3 → Task 4 ✓
- Report adaptativo (abas condicionais) → Task 6 ✓
- Entry/exit points (`--track`, `--until`) → Task 5 ✓
- COMEBin obrigatório-se-GPU: **não regride** — Task 3 mantém `comebin/done.txt` sempre
  (a regra `comebin` já auto-skipa sem GPU via `COMEBIN_ENABLED`); nenhuma mudança de
  obrigatoriedade é necessária neste plano. Documentado aqui como intencionalmente inalterado.
- Co-assembly config parseada/validada; execução no Plano 2 → Tasks 1–2 (parse/validate),
  `_t_coassembly()` stub em Task 3 ✓

**Placeholder scan:** nenhum "TBD/TODO"; todo passo tem código/comando concreto.

**Type consistency:** `resolve_pipeline_config` keys usadas na Task 2 batem com as
definidas na Task 1. `_track_overrides`/`_STAGE_ALIASES` na Task 5 batem com os testes.
`TRACKS` keys (reads/viral/prok/integration/coassembly) consistentes entre Tasks 1/6.

**Nota de execução:** vários steps de verificação dependem de `conda run -n snakemake`
e de DBs. Onde os DBs não estiverem presentes, o critério é a **ausência de erro de
config/DAG** (ValueError/NameError), não a resolução completa de inputs de ferramentas.
```
