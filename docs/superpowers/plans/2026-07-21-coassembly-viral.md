# Co-assembly Viral Consumer (Plano 3) — Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Run the CORE viral pipeline (detection → consensus → CheckV → vOTU clustering → taxonomy) on each group's co-assembled contigs, producing per-group vOTUs surfaced in the report's "Co-assembly" tab. Works for both short reads and long reads.

**Architecture:** New rules in `rules/coassembly.smk` that mirror the per-sample viral rules (`virsorter2`/`genomad`/`vibrant`/`viral_consensus` in viral_detection.smk; `checkv`/`skani_votu`/`skani_cluster` in viral_binning.smk; `mmseqs_taxonomy_viral`/`viral_taxonomy` in taxonomy.smk) but read the canonical group assembly `{OUTDIR}/coassembly/{group}/contigs.fa` and write under `{OUTDIR}/coassembly/{group}/viral/…`. All these tools operate on a contig FASTA, so the consumer is **mode-agnostic** (no read mapping) — it runs the same for SR (megahit) and LR (flye). Gated on `COASSEMBLY_ENABLED and COASSEMBLY_VIRAL`.

**Tech Stack:** Snakemake, VirSorter2, GeNomad, VIBRANT, CheckV, skani, MMseqs2, (optional vConTACT3), Python 3.11, `snakemake -n`.

## Global Constraints

- **Commits MUST NOT include any Claude co-authorship line** (`Co-Authored-By: Claude …`). Vale para TODO commit.
- **Depends on Plano 2 merged** (co-assembly track on master). Consumes: `COASSEMBLY_ENABLED`, `COASSEMBLY_VIRAL`, `GROUPS`, canonical `{OUTDIR}/coassembly/{group}/contigs.fa`.
- **Mode-agnostic:** these rules run for SR *and* LR — do NOT wrap them in `if not LONG_READS:` (that guard is only for co-binning). The co-assembly contigs exist for both modes.
- **Off by default:** `coassembly.viral` only matters when `coassembly.enabled=true`. With co-assembly off, none of these rules trigger.
- **Mirror faithfully:** each rule copies the real per-sample rule's env/container/flags — the implementer MUST read the referenced rule and reuse its shell body, changing only the input (→ group contigs), output paths (→ `coassembly/{group}/viral/…`), and wildcard (`{sample}`→`{group}`).
- **No read mapping / abundance in this plan** — vOTU abundance and vRhyme vMAGs are out of scope (they need coverage). Core = detection → consensus → CheckV → vOTU → taxonomy.
- Dry-run harness: env `smk8`, synthetic fastqs + `syn_meta.tsv` under the scratchpad; verify jobs resolve per group for SR and LR.

---

## File Structure

| Arquivo | Responsabilidade | Ação |
|---|---|---|
| `rules/coassembly.smk` | New block `if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:` with the mirrored viral rules (detection, consensus, CheckV, vOTU, taxonomy) | Modificar |
| `Snakefile` | `_t_coassembly()` — request group viral targets when `COASSEMBLY_VIRAL` | Modificar |
| `rules/report.smk` | optional viral-group inputs on `generate_report` | Modificar |
| `scripts/report/data_loaders.py` | extend `load_coassembly` (or add loader) with group vOTU/taxonomy | Modificar |
| `scripts/report/components/coassembly.js` | render group vOTUs section | Modificar |

---

### Task 1: Group viral detection + consensus

**Files:** Modify `rules/coassembly.smk` (new `if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:` block).

**Interfaces:**
- Consumes: `{OUTDIR}/coassembly/{group}/contigs.fa` (Plano 2), `MIN_VIRAL_TOOLS`, `VIRAL_CONSENSUS_MODE`, `SCORE_VS2_MIN`, `SCORE_GENOMAD_MIN`, DB vars (`VS2_DB`, `GENOMAD_DB`, `_VIBRANT_BASE`).
- Produces:
  - `{OUTDIR}/coassembly/{group}/viral/virsorter2/final-viral-combined.fa`
  - `{OUTDIR}/coassembly/{group}/viral/genomad/done.txt`
  - `{OUTDIR}/coassembly/{group}/viral/vibrant/done.txt`
  - `{OUTDIR}/coassembly/{group}/viral/consensus/{group}_viral_consensus.fasta`
  - `{OUTDIR}/coassembly/{group}/viral/consensus/{group}_tool_support.tsv`

- [ ] **Step 1: Add `coassembly_virsorter2`, `coassembly_genomad`, `coassembly_vibrant`**

Read `rule virsorter2`/`rule genomad`/`rule vibrant` in `rules/viral_detection.smk`. Copy each verbatim into a new block at the end of `rules/coassembly.smk`:
```python
if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:
    rule coassembly_virsorter2:
        # <copy of rule virsorter2, with:>
        #   input.contigs = f"{OUTDIR}/coassembly/{{group}}/contigs.fa"
        #   all output/log/benchmark/params paths: {OUTDIR}/{{sample}}/... → {OUTDIR}/coassembly/{{group}}/...
        #   same conda/container/threads/shell (env_viral / CONTAINERS.get("virsorter"), VS2_DB, etc.)
    ...  # genomad + vibrant likewise
```
Keep env/container/DB flags identical to the originals. Only the `{sample}`→`{group}` wildcard and the `coassembly/{group}` path prefix change.

- [ ] **Step 2: Add `coassembly_viral_consensus`**

Read `rule viral_consensus` (viral_detection.smk:125). Copy it, changing:
- `input.contigs` → `f"{OUTDIR}/coassembly/{{group}}/contigs.fa"`
- `input.vs2_done`/`genomad_done`/`vibrant_done` → `rules.coassembly_virsorter2.output...` / `rules.coassembly_genomad.output.done` / `rules.coassembly_vibrant.output.done`
- outputs → `{OUTDIR}/coassembly/{{group}}/viral/consensus/{{group}}_viral_consensus.fasta` and `..._tool_support.tsv`
- `{wildcards.sample}` → `{wildcards.group}` in the shell/run body.
Keep `MIN_VIRAL_TOOLS`/`VIRAL_CONSENSUS_MODE`/score logic identical.

- [ ] **Step 3: Dry-run — detection+consensus resolve per group (SR)**

Run:
```bash
SC=$(ls -d /tmp/claude-*/-home-lucas-metagen-pipe-final-claude/*/scratchpad | head -1)
FQ="$SC/synfq"; OUT="$SC/synout"; META="$SC/syn_meta.tsv"
CA='coassembly={"enabled":true,"grouping":"metadata","metadata":"'"$META"'","binning":false,"viral":true}'
conda run -n smk8 snakemake -n "$OUT/coassembly/g1/viral/consensus/g1_viral_consensus.fasta" \
  --configfile config.yaml --config fastq_dir="$FQ" outdir="$OUT" "$CA" --cores 1 2>&1 \
  | awk '/Job stats:/,/^$/' | grep -E "coassembly_virsorter2|coassembly_genomad|coassembly_vibrant|coassembly_viral_consensus|megahit_coassembly"
```
Expected: all five appear once for g1.

- [ ] **Step 4: Dry-run — LR mode (flye contigs → same viral rules)**

Run the same targeting `.../contigs.fa` with `fastq_dir=$SC/synlr long_reads=true`; expected: `flye_coassembly` + the viral rules resolve (mode-agnostic).

- [ ] **Step 5: Commit**

```bash
git add rules/coassembly.smk
git commit -m "feat: co-assembly viral detection + consensus (per-group, SR+LR)"
```

---

### Task 2: Group CheckV

**Files:** Modify `rules/coassembly.smk`.

**Interfaces:**
- Consumes: `coassembly_viral_consensus` output fasta; `CHECKV_DB`.
- Produces: `{OUTDIR}/coassembly/{group}/viral/checkv/quality_summary.tsv`.

- [ ] **Step 1: Add `coassembly_checkv`**

Read `rule checkv` (viral_binning.smk:18). Copy it, changing input to the group consensus fasta (`rules.coassembly_viral_consensus.output.fasta`), outputs to `{OUTDIR}/coassembly/{{group}}/viral/checkv/…`, `{sample}`→`{group}`. Keep `env`/container/`CHECKV_DB` identical.

- [ ] **Step 2: Dry-run**

Run (SR harness from Task 1) targeting `$OUT/coassembly/g1/viral/checkv/quality_summary.tsv`; expected `coassembly_checkv` + upstream resolve.

- [ ] **Step 3: Commit**

```bash
git add rules/coassembly.smk
git commit -m "feat: co-assembly CheckV on group viral genomes"
```

---

### Task 3: Group vOTU clustering

**Files:** Modify `rules/coassembly.smk`.

**Interfaces:**
- Consumes: group viral consensus fasta (or CheckV-filtered set — mirror what per-sample `skani_votu` inputs); `VOTU_ANI`, `VOTU_AF`.
- Produces: `{OUTDIR}/coassembly/{group}/viral/votu/vOTU_clusters.tsv` and `votu_all_reps.fasta`.

- [ ] **Step 1: Add `coassembly_skani_votu` + `coassembly_skani_cluster`**

Read `rule skani_votu` (viral_binning.smk:244) and its downstream `skani_cluster` (pure-Python clustering) + the vOTU representative rule. Copy them, changing input to the group viral genomes, outputs to `{OUTDIR}/coassembly/{{group}}/viral/votu/…`, `{sample}`→`{group}`. Keep skani ANI/AF params identical. Include the "no clusters / skani empty" guards from the originals.

- [ ] **Step 2: Dry-run**

Target `$OUT/coassembly/g1/viral/votu/votu_all_reps.fasta`; expected the skani rules + upstream resolve.

- [ ] **Step 3: Commit**

```bash
git add rules/coassembly.smk
git commit -m "feat: co-assembly vOTU clustering (skani) per group"
```

---

### Task 4: Group viral taxonomy

**Files:** Modify `rules/coassembly.smk`.

**Interfaces:**
- Consumes: group vOTU representative fasta (`votu_all_reps.fasta`), GeNomad taxonomy from `coassembly_genomad`, `INPHARED_DB` / MMseqs viral DB.
- Produces: `{OUTDIR}/coassembly/{group}/viral/taxonomy/viral_taxonomy_merged.tsv` + `{OUTDIR}/coassembly/{group}/viral/taxonomy/taxonomy_done.txt`.

- [ ] **Step 1: Add `coassembly_mmseqs_taxonomy_viral` + `coassembly_viral_taxonomy`**

Read `rule mmseqs_taxonomy_viral` (taxonomy.smk:100) and `rule viral_taxonomy` (the merge). Copy them, input = group vOTU reps + GeNomad output, outputs → `{OUTDIR}/coassembly/{{group}}/viral/taxonomy/…`, `{sample}`→`{group}`. Keep the merge logic (deepest-rank-wins) identical.

**Scope note:** vConTACT3 is the heaviest tool and is a taxonomy source in the per-sample track. For the Plano-3 core, taxonomy is built from **GeNomad + MMseqs2/INPHARED** (both already available for the group). vConTACT3 on group vOTUs is a **deferred add-on** — do NOT wire it here unless explicitly requested; note it in the rule docstring as a future source.

- [ ] **Step 2: Dry-run**

Target `$OUT/coassembly/g1/viral/taxonomy/taxonomy_done.txt`; expected the taxonomy rules + upstream resolve.

- [ ] **Step 3: Commit**

```bash
git add rules/coassembly.smk
git commit -m "feat: co-assembly viral taxonomy (GeNomad + MMseqs2) per group"
```

---

### Task 5: Wire into rule all + report

**Files:** Modify `Snakefile` (`_t_coassembly`), `rules/report.smk`, `scripts/report/data_loaders.py`, `scripts/report/components/coassembly.js`.

**Interfaces:**
- Consumes: group viral taxonomy_done + vOTU reps + CheckV.
- Produces: group vOTU targets in `rule all`; a "Group vOTUs" section in the Co-assembly report tab.

- [ ] **Step 1: `_t_coassembly()` viral targets**

In `Snakefile` `_t_coassembly()`, add — when `COASSEMBLY_ENABLED and COASSEMBLY_VIRAL` (for ALL modes, no LONG_READS guard) — per group:
```python
t.append(f"{OUTDIR}/coassembly/{g}/viral/taxonomy/taxonomy_done.txt")
t.append(f"{OUTDIR}/coassembly/{g}/viral/checkv/quality_summary.tsv")
t.append(f"{OUTDIR}/coassembly/{g}/viral/votu/votu_all_reps.fasta")
```
(Keep the existing binning targets from Plano 2 as-is.)

- [ ] **Step 2: report.smk optional inputs + loader + JS**

- `rules/report.smk`: add an optional viral-group input on `generate_report` gated on `COASSEMBLY_ENABLED and COASSEMBLY_VIRAL and GROUPS`, e.g. `expand(f"{OUTDIR}/coassembly/{{group}}/viral/taxonomy/taxonomy_done.txt", group=list(GROUPS.keys()))`.
- `scripts/report/data_loaders.py`: extend `load_coassembly` (or add `load_coassembly_viral`) to read each group's `viral/checkv/quality_summary.tsv` + `viral/taxonomy/viral_taxonomy_merged.tsv` → a `vOTUs` list per group `{count, taxonomy summary}`. Tolerant of missing files.
- `scripts/report/components/coassembly.js`: render a "Group vOTUs" sub-section (count + top families) alongside the existing group-MAGs table.

- [ ] **Step 3: Dry-run — viral targets appear; default unaffected**

Run co-assembly-on (viral:true) dry-run → the three viral targets per group present; default (co-assembly off) → 0 coassembly rules.

- [ ] **Step 4: Commit**

```bash
git add Snakefile rules/report.smk scripts/report/data_loaders.py scripts/report/components/coassembly.js
git commit -m "feat: wire co-assembly group vOTUs into rule all + report tab"
```

---

## Self-Review

**Spec coverage (core: detection → consensus → CheckV → vOTU → taxonomy, SR+LR):**
- Detection + consensus → Task 1 ✓
- CheckV → Task 2 ✓
- vOTU clustering → Task 3 ✓
- Taxonomy (GeNomad + MMseqs2; vConTACT3 deferred) → Task 4 ✓
- rule all + report → Task 5 ✓
- Mode-agnostic (SR+LR, no read mapping) → all tasks are contig-based, gated on `COASSEMBLY_VIRAL` without a LONG_READS guard ✓

**Placeholder scan:** the rules are specified as "mirror rule X with these exact input/output/wildcard changes" (the proven Plano-2 pattern) rather than inlined shell — each references the concrete source rule and the exact path/wildcard changes. No vague "handle edge cases".

**Type consistency:** output paths are consistent across tasks (`coassembly/{group}/viral/{virsorter2,genomad,vibrant,consensus,checkv,votu,taxonomy}/…`); Task 5 targets exactly the paths Tasks 2–4 produce.

**Scope guards:** vConTACT3, vOTU abundance, vRhyme vMAGs, and annotation/host/defense on group vOTUs are explicitly OUT of scope (deferred), matching the user's "core" decision.

**Integration note:** Task 1's detection rules are heavy (VS2+GeNomad+VIBRANT per group). On the server this roughly doubles viral-detection cost for grouped samples — it is opt-in via `coassembly.viral`. Real-tool execution is the server gate; the sandbox verifies only DAG resolution.
