# Co-assembly / Co-binning Track (Plano 2) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Adicionar um track opcional de co-assembly (MEGAHIT por grupo) que alimenta co-binning VAMB, mais VAMB multi-split independente, com grouping por metadata TSV, saídas de grupo paralelas (CheckM2/GTDB) e uma aba única "Co-assembly" no report.

**Architecture:** Novo `rules/coassembly.smk` com um novo namespace de wildcard `{group}` (paralelo ao `{sample}` existente). Grupos vêm de um metadata TSV. Cada grupo: MEGAHIT co-assembly → cada sample do grupo mapeado de volta → matriz de abundância multi-amostra → VAMB co-binning → CheckM2 + GTDB de grupo. VAMB multi-split é um caminho independente (concatena assemblies individuais). Saídas ficam lado a lado com per-sample, **sem derep global**. O `_t_coassembly()` stub do Plano 1 passa a retornar alvos reais.

**Tech Stack:** Snakemake, Python 3.11, MEGAHIT, BWA-MEM2/minimap2, VAMB v5, CheckM2, GTDB-Tk, pytest, `snakemake -n`.

## Global Constraints

- **Commits NÃO incluem co-autoria do Claude** (`Co-Authored-By: Claude ...`). Regra do usuário. Vale para TODO commit.
- **Depende do Plano 1 mergeado** (branch de tracks). Consome: `COASSEMBLY_ENABLED`, `COASSEMBLY_GROUPING`, `COASSEMBLY_BINNING`, `COBINNING_MULTISPLIT` (já resolvidos em `pipeline_config.py`); `_t_coassembly()` stub em `Snakefile`; VAMB já removido do per-sample.
- **Co-assembly e cobinning OFF por padrão** — feature opt-in. Comportamento default do pipeline inalterado.
- **Sem derep global.** Saídas de grupo são paralelas às per-sample.
- **Consumidor viral-na-co-assembly é do Plano 3** — fora de escopo aqui. Este plano cobre só o consumidor de binning (prok) + multi-split + report.
- `config.yaml` canônico; overrides via `vapor --set`.
- HTML write usa `encoding='utf-8'`.

---

## File Structure

| Arquivo | Responsabilidade | Ação |
|---|---|---|
| `coassembly_groups.py` (repo root) | Parser puro do metadata TSV → dict grupo→[samples]; resolução do modo (metadata/all) | **Criar** |
| `rules/coassembly.smk` | Regras: megahit_coassembly, coassembly_mapback, vamb_cobinning, cobinning_multisplit, checkm2_group, gtdbtk_group | **Criar** |
| `Snakefile` | `include: "rules/coassembly.smk"`; construir `GROUPS`; `_t_coassembly()` real | Modificar |
| `config.yaml` | chave `coassembly.metadata` (path do TSV) | Modificar |
| `rules/report.smk` | inputs opcionais de co-assembly em `generate_report` (quando habilitado) | Modificar |
| `scripts/report/data_loaders.py` | `load_coassembly(...)` — CheckM2/GTDB de grupo | Modificar |
| `scripts/report/renderer.py` | injetar `COASSEMBLY_DATA` | Modificar |
| `scripts/report/components/shell.html` | nav-tab + painel "Co-assembly" | Modificar |
| `scripts/report/components/coassembly.js` | render da aba (tabela MAGs de grupo) | **Criar** |
| `tests/test_coassembly_groups.py` | testes do parser de grupos | **Criar** |

---

### Task 1: Parser de grupos do metadata (`coassembly_groups.py`)

**Files:**
- Create: `coassembly_groups.py`
- Create: `tests/test_coassembly_groups.py`
- Modify: `pyproject.toml` (py-modules += "coassembly_groups")

**Interfaces:**
- Produces:
  - `parse_groups(metadata_path: str, samples: list[str], mode: str) -> dict[str, list[str]]`
    — `mode="all"` → `{"all": samples}`; `mode="metadata"` → lê TSV com colunas
    `sample`, `group`, agrupa só samples presentes em `samples`. Ignora grupos vazios.
    Levanta `ValueError` se `mode="metadata"` e o arquivo não existe ou não tem as colunas.

- [ ] **Step 1: Write the failing test**

Create `tests/test_coassembly_groups.py`:

```python
import textwrap
import pytest
from coassembly_groups import parse_groups


def _write(tmp_path, content):
    p = tmp_path / "meta.tsv"
    p.write_text(textwrap.dedent(content))
    return str(p)


def test_mode_all_single_group():
    g = parse_groups("", ["s1", "s2", "s3"], "all")
    assert g == {"all": ["s1", "s2", "s3"]}


def test_metadata_groups(tmp_path):
    meta = _write(tmp_path, """\
        sample\tgroup
        s1\trio
        s2\trio
        s3\tsolo
    """)
    g = parse_groups(meta, ["s1", "s2", "s3"], "metadata")
    assert g == {"rio": ["s1", "s2"], "solo": ["s3"]}


def test_metadata_ignores_unknown_samples(tmp_path):
    meta = _write(tmp_path, """\
        sample\tgroup
        s1\trio
        sX\trio
    """)
    g = parse_groups(meta, ["s1"], "metadata")
    assert g == {"rio": ["s1"]}


def test_metadata_missing_file_raises():
    with pytest.raises(ValueError, match="metadata"):
        parse_groups("/no/such.tsv", ["s1"], "metadata")


def test_metadata_missing_columns_raises(tmp_path):
    meta = _write(tmp_path, "sample\tfoo\ns1\tbar\n")
    with pytest.raises(ValueError, match="colunas"):
        parse_groups(meta, ["s1"], "metadata")
```

- [ ] **Step 2: Run test to verify it fails**

Run: `python -m pytest tests/test_coassembly_groups.py -v`
Expected: FAIL — `ModuleNotFoundError: No module named 'coassembly_groups'`

- [ ] **Step 3: Write minimal implementation**

Create `coassembly_groups.py`:

```python
"""Pure grouping resolver for the co-assembly track. Snakemake-free (pytest-able)."""
import csv
import os


def parse_groups(metadata_path: str, samples: list, mode: str) -> dict:
    mode = (mode or "metadata").strip().lower()
    sample_set = set(samples)

    if mode == "all":
        return {"all": list(samples)}

    # mode == "metadata"
    if not metadata_path or not os.path.exists(metadata_path):
        raise ValueError(
            f"coassembly.grouping=metadata requer um metadata TSV existente "
            f"(coassembly.metadata). Caminho inválido: {metadata_path!r}"
        )
    groups: dict = {}
    with open(metadata_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        cols = {c.lower(): c for c in (reader.fieldnames or [])}
        if "sample" not in cols or "group" not in cols:
            raise ValueError(
                "metadata TSV precisa das colunas 'sample' e 'group'; "
                f"encontrado: {reader.fieldnames}"
            )
        for row in reader:
            s = (row[cols["sample"]] or "").strip()
            g = (row[cols["group"]] or "").strip()
            if not s or not g or s not in sample_set:
                continue
            groups.setdefault(g, []).append(s)
    return groups
```

Modify `pyproject.toml`:

```toml
py-modules = ["vapor", "pipeline_config", "coassembly_groups"]
```

- [ ] **Step 4: Run test to verify it passes**

Run: `python -m pytest tests/test_coassembly_groups.py -v`
Expected: PASS (5 passed)

- [ ] **Step 5: Commit**

```bash
git add coassembly_groups.py tests/test_coassembly_groups.py pyproject.toml
git commit -m "feat: pure metadata grouping resolver for co-assembly track"
```

---

### Task 2: Construir GROUPS no Snakefile + config

**Files:**
- Modify: `Snakefile` (após bloco de resolução do Plano 1; e `include` list)
- Modify: `config.yaml` (chave `coassembly.metadata`)

**Interfaces:**
- Consumes: `parse_groups` (Task 1); `COASSEMBLY_ENABLED`, `COASSEMBLY_GROUPING` (Plano 1).
- Produces: `GROUPS: dict[str,list[str]]` no Snakefile (vazio se co-assembly off).

- [ ] **Step 1: Add metadata key to config**

In `config.yaml`, under the `coassembly:` block from Plano 1, add:

```yaml
coassembly:
  enabled:  false
  grouping: metadata
  metadata: ""          # path do TSV com colunas sample,group (grouping=metadata)
  viral:    true
  binning:  true
```

- [ ] **Step 2: Build GROUPS in Snakefile**

In `Snakefile`, after the Plano-1 track resolution block (after `COBINNING_MULTISPLIT = ...`),
add:

```python
from coassembly_groups import parse_groups

COASSEMBLY_METADATA = _expand(config.get("coassembly", {}).get("metadata", "")) \
    if config.get("coassembly", {}).get("metadata", "") else ""

GROUPS = {}
if COASSEMBLY_ENABLED or COBINNING_MULTISPLIT:
    GROUPS = parse_groups(COASSEMBLY_METADATA, list(SAMPLES.keys()), COASSEMBLY_GROUPING)
    print(f"[Snakemake] Co-assembly groups: {({g: len(s) for g, s in GROUPS.items()})}")
```

- [ ] **Step 3: Add the include**

In `Snakefile`, in the `include:` list, after `include: "rules/reads_classify.smk"`, add:

```python
include: "rules/coassembly.smk"
```

(Creating `rules/coassembly.smk` happens in Task 3; add an empty placeholder now so the
include resolves:)

Run: `printf '# co-assembly / co-binning rules (Plano 2)\n' > rules/coassembly.smk`

- [ ] **Step 4: Dry-run — default (co-assembly off) unaffected**

Run: `conda run -n snakemake snakemake -n --configfile config.yaml --cores 1 2>&1 | grep -ci "error"`
Expected: 0 config/parse errors (GROUPS empty, include resolves).

- [ ] **Step 5: Commit**

```bash
git add Snakefile config.yaml rules/coassembly.smk
git commit -m "feat: build co-assembly GROUPS from metadata + wire include"
```

---

### Task 3: MEGAHIT co-assembly por grupo

**Files:**
- Modify: `rules/coassembly.smk`

**Interfaces:**
- Consumes: `GROUPS`, `_clean_r1`/`_clean_r2` helpers (Snakefile), `MEGAHIT_MEM`, `MEGAHIT_PRESET`.
- Produces: `{OUTDIR}/coassembly/{group}/megahit/final.contigs.fa` per group.
  Wildcard `{group}` constrained to `GROUPS` keys.

- [ ] **Step 1: Write the co-assembly rule**

In `rules/coassembly.smk`, add:

```python
# Group → list of sample R1/R2 (post host-removal), for co-assembly input.
def _group_r1(wc):
    return [_clean_r1(type("W", (), {"sample": s})) for s in GROUPS[wc.group]]

def _group_r2(wc):
    if SINGLE_END:
        return []
    return [_clean_r2(type("W", (), {"sample": s})) for s in GROUPS[wc.group]]


wildcard_constraints:
    group = "|".join(re.escape(g) for g in GROUPS) if GROUPS else "^$",


rule megahit_coassembly:
    """MEGAHIT co-assembly of all samples in a group (comma-separated inputs)."""
    input:
        r1 = _group_r1,
        r2 = _group_r2,
    output:
        contigs = f"{OUTDIR}/coassembly/{{group}}/megahit/final.contigs.fa",
    log:
        f"{OUTDIR}/coassembly/{{group}}/logs/megahit.log"
    benchmark:
        f"{OUTDIR}/coassembly/{{group}}/benchmarks/megahit.tsv"
    conda: "../envs/env_assembly.yaml"
    container: CONTAINERS.get("megahit")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/coassembly/{{group}}/megahit",
        preset = (f"--presets {MEGAHIT_PRESET}"
                  if MEGAHIT_PRESET in ("meta-sensitive", "meta-large") else ""),
    shell:
        """
        rm -rf {params.outdir}
        R1=$(echo {input.r1} | tr ' ' ',')
        if [ -n "{input.r2}" ]; then
            R2=$(echo {input.r2} | tr ' ' ',')
            megahit -1 "$R1" -2 "$R2" -t {threads} -m {MEGAHIT_MEM} \
                {params.preset} -o {params.outdir} > {log} 2>&1
        else
            megahit -r "$R1" -t {threads} -m {MEGAHIT_MEM} \
                {params.preset} -o {params.outdir} > {log} 2>&1
        fi
        """
```

Ensure `import re` is available in `coassembly.smk` (add `import re` at its top if the
Snakefile's import doesn't propagate — Snakemake includes share the global namespace, so
`re` imported in the Snakefile is available; add a local `import re` only if a dry-run
`NameError` appears).

- [ ] **Step 2: Dry-run with a synthetic group config**

Prepare a tiny metadata TSV mapping the discovered samples into one group and dry-run:

Run:
```bash
python - <<'PY'
import yaml
# assumes >=1 sample discoverable; writes a metadata grouping all samples as 'g1'
PY
conda run -n snakemake snakemake -n \
  --configfile config.yaml \
  --config coassembly='{"enabled": true, "grouping": "all", "binning": false, "viral": false}' \
  "$(python -c "import yaml;c=yaml.safe_load(open('config.yaml'));print(c['outdir'])")/coassembly/all/megahit/final.contigs.fa" \
  --cores 1 2>&1 | grep -E "megahit_coassembly|final.contigs" | head
```
Expected: `megahit_coassembly` job for group `all` appears in the plan.

- [ ] **Step 3: Commit**

```bash
git add rules/coassembly.smk
git commit -m "feat: MEGAHIT co-assembly per group"
```

---

### Task 4: Map-back + matriz de abundância multi-amostra

**Files:**
- Modify: `rules/coassembly.smk`

**Interfaces:**
- Consumes: `megahit_coassembly.output.contigs`; `_clean_r1/_clean_r2` per sample; mapping
  tools (bwa-mem2 short / minimap2 long) as used in `rules/mapping.smk`.
- Produces:
  - `{OUTDIR}/coassembly/{group}/mapping/{sample}_depth.txt` per (group, sample)
  - `{OUTDIR}/coassembly/{group}/vamb/abundance.tsv` (multi-sample matrix)

- [ ] **Step 1: Add index + per-sample map-back rules**

In `rules/coassembly.smk`, add (short-read path via bwa-mem2; mirror `rules/mapping.smk`
conventions — the implementer must confirm the exact bwa-mem2/samtools/jgi_summarize
invocation against `rules/mapping.smk`):

```python
rule coassembly_index:
    input:
        contigs = rules.megahit_coassembly.output.contigs,
    output:
        done = f"{OUTDIR}/coassembly/{{group}}/mapping/index.done",
    conda: "../envs/env_mapping.yaml"
    container: CONTAINERS.get("mapping")
    threads: THREADS
    params:
        prefix = f"{OUTDIR}/coassembly/{{group}}/mapping/idx",
    shell:
        """
        bwa-mem2 index -p {params.prefix} {input.contigs}
        touch {output.done}
        """


rule coassembly_mapback:
    """Map one sample's reads back to its group's co-assembly → per-sample depth."""
    input:
        idx = rules.coassembly_index.output.done,
        contigs = rules.megahit_coassembly.output.contigs,
        r1 = lambda wc: _clean_r1(type("W", (), {"sample": wc.sample})),
        r2 = lambda wc: _clean_r2(type("W", (), {"sample": wc.sample})),
    output:
        depth = f"{OUTDIR}/coassembly/{{group}}/mapping/{{sample}}_depth.txt",
    conda: "../envs/env_mapping.yaml"
    container: CONTAINERS.get("mapping")
    threads: THREADS
    params:
        prefix = f"{OUTDIR}/coassembly/{{group}}/mapping/idx",
        bam    = f"{OUTDIR}/coassembly/{{group}}/mapping/{{sample}}.sorted.bam",
    shell:
        """
        if [ -n "{input.r2}" ]; then
            bwa-mem2 mem -t {threads} {params.prefix} {input.r1} {input.r2} \
                | samtools sort -@ {threads} -o {params.bam} -
        else
            bwa-mem2 mem -t {threads} {params.prefix} {input.r1} \
                | samtools sort -@ {threads} -o {params.bam} -
        fi
        samtools index {params.bam}
        jgi_summarize_bam_contig_depths --outputDepth {output.depth} {params.bam}
        """
```

- [ ] **Step 2: Build the multi-sample abundance matrix (VAMB v5 format)**

Add a rule that joins per-sample depths into one abundance TSV keyed by contig, with one
column per sample (VAMB v5 `--abundance_tsv` expects `contigname\t<s1>\t<s2>...`):

```python
rule coassembly_abundance:
    input:
        depths = lambda wc: expand(
            f"{OUTDIR}/coassembly/{wc.group}/mapping/{{sample}}_depth.txt",
            sample=GROUPS[wc.group]),
    output:
        matrix = f"{OUTDIR}/coassembly/{{group}}/vamb/abundance.tsv",
    run:
        import csv, os
        samples = GROUPS[wildcards.group]
        # jgi depth files: columns contigName, contigLen, totalAvgDepth, <bam>-var...
        # Use column 'totalAvgDepth' (index 2) per sample file, keyed by contigName.
        cov = {}
        contig_order = []
        for s, path in zip(samples, input.depths):
            with open(path) as fh:
                r = csv.reader(fh, delimiter="\t")
                header = next(r)
                for row in r:
                    c = row[0]
                    if c not in cov:
                        cov[c] = {}
                        contig_order.append(c)
                    cov[c][s] = row[2]
        os.makedirs(os.path.dirname(output.matrix), exist_ok=True)
        with open(output.matrix, "w") as out:
            out.write("contigname\t" + "\t".join(samples) + "\n")
            for c in contig_order:
                out.write(c + "\t" + "\t".join(cov[c].get(s, "0") for s in samples) + "\n")
```

- [ ] **Step 3: Dry-run resolves map-back + abundance for a group**

Run:
```bash
OUT=$(python -c "import yaml;print(yaml.safe_load(open('config.yaml'))['outdir'])")
conda run -n snakemake snakemake -n --configfile config.yaml \
  --config coassembly='{"enabled": true, "grouping": "all", "binning": true, "viral": false}' \
  "$OUT/coassembly/all/vamb/abundance.tsv" --cores 1 2>&1 | grep -E "coassembly_mapback|coassembly_abundance" | head
```
Expected: map-back (one per sample) + abundance rule appear.

- [ ] **Step 4: Commit**

```bash
git add rules/coassembly.smk
git commit -m "feat: co-assembly map-back + multi-sample abundance matrix"
```

---

### Task 5: VAMB co-binning + multi-split

**Files:**
- Modify: `rules/coassembly.smk`

**Interfaces:**
- Consumes: `coassembly_abundance.output.matrix`, `megahit_coassembly.output.contigs`.
- Produces: `{OUTDIR}/coassembly/{group}/vamb/done.txt` (+ `bins/*.fa`). For multi-split:
  `{OUTDIR}/coassembly/multisplit/vamb/done.txt`.

- [ ] **Step 1: VAMB co-binning rule (on co-assembly)**

In `rules/coassembly.smk`, add:

```python
rule vamb_cobinning:
    """VAMB v5 co-binning on a group's co-assembly using the multi-sample abundance."""
    input:
        contigs   = rules.megahit_coassembly.output.contigs,
        abundance = rules.coassembly_abundance.output.matrix,
    output:
        done = f"{OUTDIR}/coassembly/{{group}}/vamb/done.txt",
    log:
        f"{OUTDIR}/coassembly/{{group}}/logs/vamb.log"
    conda: "../envs/env_vamb.yaml"
    container: CONTAINERS.get("vamb")
    threads: THREADS
    params:
        outdir = f"{OUTDIR}/coassembly/{{group}}/vamb",
    shell:
        """
        rm -rf {params.outdir}/run
        vamb bin default \
            --outdir {params.outdir}/run \
            --fasta {input.contigs} \
            --abundance_tsv {input.abundance} \
            -p {threads} > {log} 2>&1 || echo "[VAMB] failed" | tee -a {log}
        mkdir -p {params.outdir}
        touch {output.done}
        """
```

(The implementer must confirm the VAMB v5 env name `env_vamb.yaml` matches the one the
removed per-sample rule used — reuse it; do not create a new env.)

- [ ] **Step 2: VAMB multi-split rule (independent of co-assembly)**

Add a rule that concatenates per-sample `rep_seq.fasta` (renaming contigs with a sample
prefix `S<name>C...`), maps all samples to the catalog, and runs VAMB, then splits. Since
this mirrors Task 4's mapping against a catalog, scope it to produce a single global
catalog under `coassembly/multisplit/`:

```python
rule multisplit_catalog:
    input:
        reps = expand(f"{OUTDIR}/{{sample}}/mmseqs/{{sample}}_rep_seq.fasta", sample=SAMPLES),
    output:
        catalog = f"{OUTDIR}/coassembly/multisplit/catalog.fasta",
    run:
        import os
        os.makedirs(os.path.dirname(output.catalog), exist_ok=True)
        with open(output.catalog, "w") as out:
            for s, path in zip(list(SAMPLES.keys()), input.reps):
                with open(path) as fh:
                    for line in fh:
                        if line.startswith(">"):
                            out.write(f">S{s}C{line[1:]}")
                        else:
                            out.write(line)
```

Add `multisplit_mapback`, `multisplit_abundance`, and `multisplit_vamb` rules mirroring
Task 4 + Step 1 but against `coassembly/multisplit/catalog.fasta` with all SAMPLES as the
abundance columns, output `{OUTDIR}/coassembly/multisplit/vamb/done.txt`. (The implementer
reuses the same mapping/abundance/vamb shell bodies from Task 4 and Step 1, parameterized
to the multisplit paths.)

- [ ] **Step 3: Dry-run — co-binning and multisplit resolve**

Run:
```bash
OUT=$(python -c "import yaml;print(yaml.safe_load(open('config.yaml'))['outdir'])")
conda run -n snakemake snakemake -n --configfile config.yaml \
  --config coassembly='{"enabled": true, "grouping": "all", "binning": true, "viral": false}' cobinning_multisplit=true \
  "$OUT/coassembly/all/vamb/done.txt" "$OUT/coassembly/multisplit/vamb/done.txt" \
  --cores 1 2>&1 | grep -E "vamb_cobinning|multisplit_vamb" | head
```
Expected: both `vamb_cobinning` and `multisplit_vamb` jobs appear.

- [ ] **Step 4: Commit**

```bash
git add rules/coassembly.smk
git commit -m "feat: VAMB co-binning (co-assembly) + multi-split binning"
```

---

### Task 6: CheckM2 + GTDB-Tk de grupo

**Files:**
- Modify: `rules/coassembly.smk`

**Interfaces:**
- Consumes: `vamb_cobinning.output.done` (group bins dir).
- Produces:
  - `{OUTDIR}/coassembly/{group}/checkm2/quality_report.tsv`
  - `{OUTDIR}/coassembly/{group}/gtdbtk/done.txt`

- [ ] **Step 1: Add group-level CheckM2 + GTDB rules**

In `rules/coassembly.smk`, add rules mirroring the per-sample `checkm2` and `gtdbtk` rules
(the implementer copies the shell bodies from `rules/prok_binning.smk`, pointing the bins
input dir at `{OUTDIR}/coassembly/{group}/vamb/run/bins` and DBs at `CHECKM2_DB`/`GTDBTK_DB`):

```python
rule checkm2_group:
    input:
        done = rules.vamb_cobinning.output.done,
    output:
        report = f"{OUTDIR}/coassembly/{{group}}/checkm2/quality_report.tsv",
    conda: "../envs/env_checkm2.yaml"
    container: CONTAINERS.get("checkm2")
    threads: THREADS
    params:
        bins = f"{OUTDIR}/coassembly/{{group}}/vamb/run/bins",
        outdir = f"{OUTDIR}/coassembly/{{group}}/checkm2",
    shell:
        """
        checkm2 predict --threads {threads} --force \
            --input {params.bins} --extension fa \
            --output-directory {params.outdir} \
            --database_path {CHECKM2_DB} \
            || : ; test -f {output.report} || \
            (mkdir -p {params.outdir} && printf 'Name\tCompleteness\tContamination\n' > {output.report})
        """


rule gtdbtk_group:
    input:
        report = rules.checkm2_group.output.report,
    output:
        done = f"{OUTDIR}/coassembly/{{group}}/gtdbtk/done.txt",
    conda: "../envs/env_gtdbtk.yaml"
    container: CONTAINERS.get("gtdbtk")
    threads: THREADS
    params:
        bins = f"{OUTDIR}/coassembly/{{group}}/vamb/run/bins",
        outdir = f"{OUTDIR}/coassembly/{{group}}/gtdbtk",
    shell:
        """
        export GTDBTK_DATA_PATH={GTDBTK_DB}
        gtdbtk classify_wf --genome_dir {params.bins} --extension fa \
            --out_dir {params.outdir} --cpus {threads} --skip_ani_screen \
            || echo "[GTDB group] no bins / failed"
        touch {output.done}
        """
```

(The implementer must confirm env file names + the VAMB v5 bins output subpath
(`run/bins` vs `run/vae_clusters` etc.) against the real VAMB output layout.)

- [ ] **Step 2: Dry-run — group QC/taxonomy resolve**

Run:
```bash
OUT=$(python -c "import yaml;print(yaml.safe_load(open('config.yaml'))['outdir'])")
conda run -n snakemake snakemake -n --configfile config.yaml \
  --config coassembly='{"enabled": true, "grouping": "all", "binning": true, "viral": false}' \
  "$OUT/coassembly/all/gtdbtk/done.txt" --cores 1 2>&1 | grep -E "checkm2_group|gtdbtk_group" | head
```
Expected: both appear.

- [ ] **Step 3: Commit**

```bash
git add rules/coassembly.smk
git commit -m "feat: group-level CheckM2 + GTDB-Tk for co-binning MAGs"
```

---

### Task 7: `_t_coassembly()` real no Snakefile

**Files:**
- Modify: `Snakefile` (`_t_coassembly` function from Plano 1)

**Interfaces:**
- Consumes: `GROUPS`, `COASSEMBLY_ENABLED`, `COASSEMBLY_BINNING`, `COBINNING_MULTISPLIT`.
- Produces: real target list for the co-assembly track.

- [ ] **Step 1: Replace the stub**

In `Snakefile`, replace the `_t_coassembly()` stub from Plano 1 with:

```python
def _t_coassembly():
    t = []
    if COASSEMBLY_ENABLED and COASSEMBLY_BINNING:
        for g in GROUPS:
            t.append(f"{OUTDIR}/coassembly/{g}/gtdbtk/done.txt")
            t.append(f"{OUTDIR}/coassembly/{g}/checkm2/quality_report.tsv")
    if COBINNING_MULTISPLIT:
        t.append(f"{OUTDIR}/coassembly/multisplit/vamb/done.txt")
    return t
```

(Co-assembly viral targets are added in Plano 3.)

- [ ] **Step 2: Dry-run — enabling co-assembly adds group targets**

Run:
```bash
conda run -n snakemake snakemake -n --configfile config.yaml \
  --config coassembly='{"enabled": true, "grouping": "all", "binning": true, "viral": false}' \
  --cores 1 2>&1 | grep -c "coassembly/all/gtdbtk"
```
Expected: >0. And with defaults (co-assembly off):
Run: `conda run -n snakemake snakemake -n --configfile config.yaml --cores 1 2>&1 | grep -c "coassembly/"`
Expected: 0.

- [ ] **Step 3: Commit**

```bash
git add Snakefile
git commit -m "feat: wire co-assembly targets into rule all"
```

---

### Task 8: Aba "Co-assembly" no report

**Files:**
- Modify: `rules/report.smk` (optional inputs when COASSEMBLY_ENABLED)
- Modify: `scripts/report/data_loaders.py` (`load_coassembly`)
- Modify: `scripts/report/renderer.py` (inject `COASSEMBLY_DATA`)
- Modify: `scripts/report/components/shell.html` (nav-tab + panel)
- Create: `scripts/report/components/coassembly.js`

**Interfaces:**
- Consumes: group CheckM2/GTDB TSVs.
- Produces: `COASSEMBLY_DATA` JS constant `{groups: [{group, mags: [{bin, completeness, contamination, classification}]}]}`; single nav tab `data-tab="coassembly"`.

- [ ] **Step 1: Loader**

In `scripts/report/data_loaders.py`, add:

```python
def load_coassembly(outdir, groups):
    """Collect group-level CheckM2 + GTDB into report records (empty if none)."""
    import os, glob
    out = []
    for g in groups:
        checkm2 = os.path.join(outdir, "coassembly", g, "checkm2", "quality_report.tsv")
        mags = []
        if os.path.exists(checkm2):
            for row in load_tsv(checkm2):
                mags.append({
                    "bin": row.get("Name", ""),
                    "completeness": safe_float(row.get("Completeness", 0)),
                    "contamination": safe_float(row.get("Contamination", 0)),
                    "classification": "",
                })
        # join GTDB classification if present
        for summ in glob.glob(os.path.join(outdir, "coassembly", g, "gtdbtk", "**", "*summary.tsv"), recursive=True):
            for row in load_tsv(summ):
                name = row.get("user_genome", "")
                for m in mags:
                    if m["bin"] == name:
                        m["classification"] = row.get("classification", "")
        if mags:
            out.append({"group": g, "mags": mags})
    return {"groups": out, "has_data": bool(out)}
```

- [ ] **Step 2: report.smk optional inputs + renderer injection**

In `rules/report.smk`, add to `generate_report` `input:` (mirroring the reads_classify
optional-input pattern):

```python
        **({
            "coassembly_sentinel": expand(
                f"{OUTDIR}/coassembly/{{group}}/gtdbtk/done.txt", group=list(GROUPS.keys()))
        } if (COASSEMBLY_ENABLED and COASSEMBLY_BINNING and GROUPS) else {}),
```

and to `params:` add `coassembly_groups = list(GROUPS.keys())`.

In `scripts/report/renderer.py`, call the loader and inject the constant:

```python
    coassembly_data = load_coassembly(outdir, list(getattr(snakemake.params, "coassembly_groups", []) or []))
```

and add `"COASSEMBLY_DATA": coassembly_data,` to the data_dict. Import `load_coassembly`
in the `from .data_loaders import (...)` block.

- [ ] **Step 3: shell.html tab + panel + JS include**

In `scripts/report/components/shell.html`, add a nav-tab after the "Reads Survey" button:

```html
      <button role="tab" class="nav-tab" data-tab="coassembly" aria-selected="false" tabindex="-1">Co-assembly</button>
```

and a panel section (mirroring an existing `<section class="tab-panel">`):

```html
    <section id="tab-coassembly" class="tab-panel" role="tabpanel" hidden>
      <h2>Co-assembly MAGs</h2>
      <div id="coassembly-empty" class="muted" style="display:none">Co-assembly not run.</div>
      <div id="coassembly-table" class="table-wrap"></div>
    </section>
```

Register `coassembly.js` in the renderer's JS include list (`renderer.py`, the
`for js_file in [...]` list): add `"coassembly.js"`.

Also extend the Plano-1 `hideDisabledTrackTabs` `TAB_TRACK` map in `app.js` with
`coassembly: 'coassembly'` so the tab hides when the track is off.

- [ ] **Step 4: coassembly.js render**

Create `scripts/report/components/coassembly.js`:

```javascript
/* coassembly.js — Co-assembly MAGs tab */
(function () {
  'use strict';
  const CA = typeof COASSEMBLY_DATA !== 'undefined' ? COASSEMBLY_DATA : null;

  window.renderCoassembly = function () {
    const empty = document.getElementById('coassembly-empty');
    const tbl = document.getElementById('coassembly-table');
    if (!CA || !CA.has_data) { if (empty) empty.style.display = ''; if (tbl) tbl.innerHTML = ''; return; }
    const rows = [];
    for (const g of CA.groups) for (const m of g.mags) {
      rows.push(`<tr><td>${g.group}</td><td>${m.bin}</td>` +
        `<td>${(+m.completeness).toFixed(1)}%</td><td>${(+m.contamination).toFixed(1)}%</td>` +
        `<td>${m.classification || '—'}</td></tr>`);
    }
    tbl.innerHTML = `<table class="vapor-table"><thead><tr>` +
      `<th>Group</th><th>MAG</th><th>Completeness</th><th>Contamination</th><th>GTDB</th>` +
      `</tr></thead><tbody>${rows.join('')}</tbody></table>`;
  };
})();
```

Add `renderCoassembly` to the render dispatch list in `app.js` (the `forEach` over render
functions, ~line 338 in Plano 1's app.js).

- [ ] **Step 5: Verify JS bundles + loader importable**

Run: `python -c "import sys; sys.path.insert(0,'scripts'); from report.data_loaders import load_coassembly; print('ok')" 2>&1 | tail -1`
Expected: `ok` (or, if relative-import package, run via the report package entry — at
minimum `grep -c "def load_coassembly" scripts/report/data_loaders.py` → 1).

Run: `grep -c "renderCoassembly" scripts/report/components/coassembly.js scripts/report/components/app.js`
Expected: coassembly.js → ≥1, app.js → ≥1.

- [ ] **Step 6: Commit**

```bash
git add rules/report.smk scripts/report/data_loaders.py scripts/report/renderer.py scripts/report/components/shell.html scripts/report/components/coassembly.js scripts/report/components/app.js
git commit -m "feat: Co-assembly report tab (group MAGs table)"
```

---

## Self-Review

**Spec coverage:**
- Grouping por metadata TSV (+ modo `all`) → Tasks 1–2 ✓
- MEGAHIT co-assembly por grupo → Task 3 ✓
- Map-back + abundância multi-amostra → Task 4 ✓
- VAMB co-binning (co-assembly) + multi-split independente → Task 5 ✓
- CheckM2 + GTDB de grupo → Task 6 ✓
- Alvos no rule all → Task 7 ✓
- Aba única "Co-assembly" no report → Task 8 ✓
- **Sem derep global** → nenhuma task adiciona derep cross-source ✓
- Consumidor viral-na-co-assembly → **Plano 3** (fora de escopo, documentado) ✓

**Placeholder scan:** dois pontos deixam explícito "o implementer confirma contra o
código real" (interface VAMB v5 bins subpath; nomes de env; corpo shell de mapping) — não
são placeholders de conteúdo, são checagens de integração necessárias porque espelham
regras existentes. Todo passo tem código concreto.

**Type consistency:** `parse_groups` retorno (dict grupo→samples) usado consistentemente
em GROUPS/`_t_coassembly`/rules. `COASSEMBLY_DATA` shape (groups/mags/has_data) consistente
entre loader e coassembly.js.

**Dependência de integração:** Tasks 4–6 espelham regras existentes (`mapping.smk`,
`prok_binning.smk`) — o implementer DEVE ler essas regras e copiar os corpos shell reais
(bwa-mem2/minimap2/jgi_summarize/checkm2/gtdbtk) em vez de assumir os comandos aqui, que
são aproximações estruturais. Isto está sinalizado em cada task.

## Plano 3 (deferido, não escrito): consumidor viral na co-assembly

Escopo: rodar detecção viral (VS2 + GeNomad + VIBRANT) → consenso → CheckV → vRhyme →
vOTU → taxonomia sobre os contigs de `megahit_coassembly`, com map-back per-sample para
abundância, produzindo vOTUs de grupo, integrados na aba "Co-assembly". Espelha o track
viral inteiro; será planejado como Plano 3 após o Plano 2 mergear.
```
