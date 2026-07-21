# Co-assembly Parity — Viral Filter Fix + Provirus Trim + vRhyme (Plan 4)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax.

**Goal:** Correct the co-assembly track to match per-sample architecture: (1) filter free-living viruses out of the co-binning VAMB input (Plan 2 gap), (2) CheckV-trim proviruses in the vOTU set, (3) add group vRhyme (viral vMAGs) using the shared co-assembly coverage, independent of the prok co-binning flag.

**Architecture:** Per-sample maps reads to the FULL assembly (shared BAM/depth); viral/prok separation is at the BINNING INPUT (`filter_viral_for_prok` for prok; vRhyme on viral genomes). Mirror this for the co-assembly: keep the shared map-back against the full `contigs.fa`, add a group viral filter feeding VAMB nonviral contigs, add provirus trimming for the vOTU chain, and add a multi-sample group vRhyme consuming the group's per-sample BAMs.

**Tech Stack:** Snakemake, CheckV, VAMB, vRhyme, bwa-mem2, Python 3.11, `snakemake -n`.

## Global Constraints
- **Commits MUST NOT include any Claude co-authorship line.** Vale para TODO commit.
- **Depends on Plans 2+3** (co-assembly track + viral consumer). Builds on branch `feat/coassembly-viral`.
- **Short-read only** for the coverage-dependent parts (VAMB, vRhyme). The map-back rules are already defined under `if not LONG_READS:` in `rules/coassembly.smk` and run when their output is requested — no rule-definition change needed for "decoupling".
- **Mirror faithfully:** copy the real per-sample rule's env/container/flags; change only inputs/outputs/wildcard.
- Off by default. Dry-run harness: env `smk8`, `synfq`/`synlr`/`syn_meta.tsv` in scratchpad.

---

### Task 1: Filter free-living viruses out of the co-binning VAMB input (fix Plan 2 gap)

**Files:** Modify `rules/coassembly.smk`.

**Problem:** `vamb_cobinning` currently runs on the FULL `coassembly/{group}/contigs.fa` — it does NOT exclude free-living viral contigs, unlike per-sample (`filter_viral_for_prok`). So prok group MAGs can be virus-contaminated.

**Interfaces:**
- Consumes: `coassembly/{group}/contigs.fa`, `coassembly_viral_consensus` fasta, `coassembly_checkv` summary, `coassembly_genomad` done.
- Produces: `coassembly/{group}/prok_input/{group}_contigs_nonviral.fasta`; VAMB now bins this.

- [ ] **Step 1: Add `coassembly_filter_viral_for_prok`**

Read `rule filter_viral_for_prok` (rules/prok_binning.smk:40). Mirror it inside the `if not LONG_READS:` block of coassembly.smk: inputs = group `contigs.fa` + `rules.coassembly_viral_consensus.output.fasta` + `rules.coassembly_checkv.output.summary` + `rules.coassembly_genomad.output.done`; output `coassembly/{{group}}/prok_input/{{group}}_contigs_nonviral.fasta`; `{sample}`→`{group}`. Keep the "viral MINUS provirus, provirus-bearing stays" run logic IDENTICAL. NOTE: this rule requires `coassembly_viral_consensus`/`coassembly_checkv`/`coassembly_genomad`, which live in the `if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL:` block — so the filter can only run when viral is on. Guard: if `COASSEMBLY_VIRAL` is false, fall back to the full contigs (no filter possible without viral detection). Implement `coassembly_filter_viral_for_prok` only when `COASSEMBLY_VIRAL`; otherwise VAMB uses `contigs.fa` directly (document this).

- [ ] **Step 2: Repoint `vamb_cobinning` to the nonviral fasta + filter the abundance matrix**

`vamb_cobinning.input.contigs` → the nonviral fasta (when COASSEMBLY_VIRAL) else `contigs.fa`. CRITICAL: VAMB requires the `--abundance_tsv` rows to exactly match the `--fasta` contigs. The `coassembly_abundance` matrix covers ALL contigs; when binning the nonviral subset, filter the abundance to the nonviral contig set in the vamb shell (mirror the old per-sample vamb's `grep "^>" {fasta} | ... awk 'keep[$1]... '` filtering — recover it from git `git show 275340e:rules/prok_binning.smk` rule vamb). Build a filtered abundance TSV inside `vamb_cobinning` before calling vamb.

- [ ] **Step 3: Dry-run — VAMB now depends on the nonviral filter (viral on)**

```bash
SC=$(ls -d /tmp/claude-*/-home-lucas-metagen-pipe-final-claude/*/scratchpad | head -1)
FQ="$SC/synfq"; OUT="$SC/synout"; META="$SC/syn_meta.tsv"
CA='coassembly={"enabled":true,"grouping":"metadata","metadata":"'"$META"'","binning":true,"viral":true}'
conda run -n smk8 snakemake -n "$OUT/coassembly/g1/vamb/done.txt" --configfile config.yaml --config fastq_dir="$FQ" outdir="$OUT" "$CA" --cores 1 2>&1 | awk '/Job stats:/,/^$/' | grep -E "coassembly_filter_viral_for_prok|vamb_cobinning|coassembly_checkv"
```
Expect: `coassembly_filter_viral_for_prok` + `vamb_cobinning` resolve; VAMB depends on the filter. Also verify default off → 0. Paste into report.

- [ ] **Step 4: Commit** — `git add rules/coassembly.smk && git commit -m "fix: filter free-living viruses out of co-assembly co-binning input (match per-sample)"`

---

### Task 2: Provirus trim for the vOTU chain

**Files:** Modify `rules/coassembly.smk`.

**Interfaces:**
- Consumes: `coassembly_checkv` `viruses`/`proviruses` fna outputs (already produced, currently unused).
- Produces: `coassembly/{group}/viral/checkv/{group}_viral_trimmed.fasta`; skani/vOTU now cluster this instead of the raw consensus.

- [ ] **Step 1: Add `coassembly_viral_trimmed`**

Read the provirus-trim logic in `rule viral_nonredundant` (rules/viral_binning.smk:125-205) — ONLY the CheckV-trim part (viruses.fna + proviruses.fna → trimmed set; provirus header `orig|start_end` = trimmed region). Add a rule `coassembly_viral_trimmed` (in the COASSEMBLY_VIRAL block) that builds `coassembly/{{group}}/viral/checkv/{{group}}_viral_trimmed.fasta` from `rules.coassembly_checkv.output.viruses` + `.proviruses`. Do NOT include the multi-assembler dedup / bins-first parts (single assembler here). Keep it a pure-Python `run:` block.

- [ ] **Step 2: Repoint skani/vOTU to the trimmed set**

Change `coassembly_skani_votu`, `coassembly_skani_cluster`, `coassembly_viral_votu_reps` fasta inputs from `rules.coassembly_viral_consensus.output.fasta` → `rules.coassembly_viral_trimmed.output.fasta`.

- [ ] **Step 3: Dry-run** — target `votu_all_reps.fasta`; expect `coassembly_viral_trimmed` in the chain before skani. Paste into report.

- [ ] **Step 4: Commit** — `git add rules/coassembly.smk && git commit -m "feat: CheckV-trim proviruses before co-assembly vOTU clustering"`

---

### Task 3: Group vRhyme (viral vMAGs) + CheckV on vMAGs

**Files:** Modify `rules/coassembly.smk`.

**Interfaces:**
- Consumes: `coassembly_viral_trimmed` fasta, the group's per-sample BAMs `coassembly/{group}/mapping/{sample}.sorted.bam` (Plan 2 map-back — defined under `if not LONG_READS:`), `GROUPS`.
- Produces: `coassembly/{group}/bins/vrhyme/done.txt`, `coassembly/{group}/viral/checkv_vrhyme/quality_summary.tsv`.

- [ ] **Step 1: Add `coassembly_vrhyme` (multi-sample)**

Read `rule vrhyme` (rules/viral_binning.smk). Mirror it in a `if COASSEMBLY_ENABLED and COASSEMBLY_VIRAL and not LONG_READS:` guard (vRhyme needs coverage → short reads). Inputs: `viral = rules.coassembly_viral_trimmed.output.fasta`; `bams = lambda wc: expand("{OUTDIR}/coassembly/{group}/mapping/{sample}.sorted.bam", OUTDIR=OUTDIR, group=wc.group, sample=GROUPS[wc.group])` (ALL group samples' BAMs — multi-sample differential coverage). Adapt the vrhyme shell to pass multiple BAMs (`-b bam1 bam2 ...`). Output `coassembly/{{group}}/bins/vrhyme/done.txt`. Keep env_vrhyme/container.

- [ ] **Step 2: Add `coassembly_checkv_vrhyme`**

Read `rule checkv_vrhyme` (rules/viral_binning.smk:93). Mirror it: input `rules.coassembly_vrhyme.output.done`; output `coassembly/{{group}}/viral/checkv_vrhyme/quality_summary.tsv`; keep the "empty summary if no bins" guard.

- [ ] **Step 3: Dry-run** — target `coassembly/g1/viral/checkv_vrhyme/quality_summary.tsv`; expect `coassembly_vrhyme` (with 2 BAMs for g1) + `coassembly_checkv_vrhyme` + the group map-back rules resolve. Verify LR mode does NOT define these (short-read guard). Paste into report.

- [ ] **Step 4: Commit** — `git add rules/coassembly.smk && git commit -m "feat: group vRhyme vMAGs + CheckV (multi-sample coverage, short reads)"`

---

### Task 4: Wire into rule all + report

**Files:** Modify `Snakefile` (`_t_coassembly`), `rules/report.smk`, `scripts/report/data_loaders.py`, `scripts/report/components/coassembly.js`.

- [ ] **Step 1: `_t_coassembly()` — add vRhyme target**

When `COASSEMBLY_ENABLED and COASSEMBLY_VIRAL and not LONG_READS`, add per group `coassembly/{g}/viral/checkv_vrhyme/quality_summary.tsv`. (Viral taxonomy/checkv/votu targets from Plan 3 stay.)

- [ ] **Step 2: report.smk optional input + loader + JS**

- report.smk: optional `coassembly_vrhyme_sentinel` gated on `COASSEMBLY_ENABLED and COASSEMBLY_VIRAL and not LONG_READS and GROUPS`.
- data_loaders.py `load_coassembly`: read `viral/checkv_vrhyme/quality_summary.tsv` per group → `n_vmags` (count) + quality tiers; add to group record. Tolerant of missing.
- coassembly.js: show `n_vmags` in the group vOTUs section.

- [ ] **Step 3: Dry-run + tests** — viral+SR → vrhyme target present; LR → absent; default → 0; `python3 -m pytest tests/ -q`. Paste into report.

- [ ] **Step 4: Commit** — `git add -A && git commit -m "feat: wire group vRhyme vMAGs into rule all + report"`

---

## Self-Review

- Plan 2 gap (no viral filter before co-binning) → Task 1 ✓ (VAMB now bins nonviral; abundance filtered to match).
- Provirus trim → Task 2 ✓ (uses CheckV's own trimmed outputs; no multi-assembler dedup).
- vRhyme (viral vMAGs) → Task 3 ✓ (multi-sample coverage via the shared group BAMs; short-read only; independent of the `binning` flag — it requests the BAMs directly).
- Report → Task 4 ✓.
- Coverage decoupling: confirmed a NON-issue — map-back rules are defined under `if not LONG_READS:` and run when requested; vRhyme requesting the BAMs triggers them even with `binning: false`.
- Placeholders: rules specified as "mirror rule X with these input/output/wildcard changes"; the abundance-filter (Task 1 Step 2) and the multi-BAM vrhyme (Task 3 Step 1) are the two spots needing real adaptation, both flagged explicitly.
