# Co-assembly Prok Functional — Annotation + AMR + Defense (Plan 5)

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development. Steps use `- [ ]`.

**Goal:** Give the co-assembly group MAGs the same functional layer as per-sample prok bins: protein prediction, GUNC, MAG dereplication, AMR (AMRFinderPlus, RGI/CARD, DeepARG, argNorm, consensus), defense (DefenseFinder), ABRicate (VFDB/PlasmidFinder), and annotation (Bakta, eggNOG/KEGG) — all on `coassembly/{group}/vamb/run/bins/*.fna`.

**Architecture:** Mirror the per-sample rules into `rules/coassembly.smk` (inside `if COASSEMBLY_ENABLED and COASSEMBLY_BINNING and not LONG_READS:`), changing only the bins dir (→ VAMB group bins, **`.fna` extension** not `.fa`), the proteins/output paths (→ `coassembly/{group}/bins/...`), and the `{sample}`→`{group}` wildcard. Reuse the module-level helpers `_read_manifest`/`_concat_proteins` (prok_binning.smk, in scope since it's included before coassembly.smk).

## Global Constraints
- **Commits MUST NOT include any Claude co-authorship line.**
- Depends on Plans 2–4 (co-assembly + co-binning + parity). Branch `feat/coassembly-parity`.
- Short-read only (co-binning is SR-only). Off by default.
- **Bins are `.fna`** (VAMB v5), unlike per-sample Binette `.fa` — every glob/`-x`/extension must be `fna`.
- Mirror faithfully: copy env/container/DB flags from the real rule; change only input bins dir + output paths + wildcard.

---

### Task 1: `coassembly_prok_bin_proteins` (foundation — every AMR/defense tool needs it)
Mirror `rule prok_bin_proteins` (prok_binning.smk:456). bins_dir = `coassembly/{group}/vamb/run/bins`, glob `*.fna` (not `*.fa`); low_depth path uses the group nonviral fasta (`coassembly/{group}/prok_input/{group}_contigs_nonviral.fasta`); outputs `coassembly/{group}/bins/proteins/{manifest.txt,done.txt}`. env_viral/prodigal. Reuse `_read_manifest`/`_concat_proteins`.
Verify: dry-run target `coassembly/g1/bins/proteins/done.txt` resolves after vamb_cobinning. Commit `feat: co-assembly group protein prediction (prodigal per MAG)`.

### Task 2: GUNC + MAG derep (galah)
Mirror `rule gunc` (prok_binning.smk:638) and `rule galah_derep` (700), gated on `GUNC_ENABLED`/`MAG_DEREP_ENABLED` respectively. bins_dir → group VAMB bins (`.fna`); outputs `coassembly/{group}/bins/{gunc,derep}/...`. Commit `feat: co-assembly group GUNC + MAG dereplication`.

### Task 3: AMR tools (AMRFinderPlus, RGI/CARD, DeepARG)
Mirror `amrfinderplus` (defense_amr.smk:374), `rgi_card` (432), `deeparg` (517). Input = group proteins manifest (via `_concat_proteins`) + group bins. Outputs `coassembly/{group}/bins/{amrfinderplus,rgi,deeparg}/...`. Keep CARD_DB/DEEPARG_DB. Commit `feat: co-assembly group AMR (AMRFinderPlus, RGI, DeepARG)`.

### Task 4: Defense + ABRicate + argNorm + AMR consensus
Mirror `defensefinder` (47), `abricate` (606), `argnorm_normalize` (690), `amr_consensus` (767) on the group MAGs/proteins. Gated on `DEFENSE_AMR_ENABLED`/`ABRICATE_ENABLED`/`ARGNORM_ENABLED`/`AMR_CONSENSUS_ENABLED`. Commit `feat: co-assembly group defense + ABRicate + argNorm + AMR consensus`.

### Task 5: Annotation (Bakta, eggNOG, KEGG)
Mirror `bakta` (annotation.smk:213), `eggnog_prok` (308), `extract_kegg_kos` (375) on the group MAGs. Keep BAKTA_DB/EGGNOG_DB + completeness/contamination thresholds. Commit `feat: co-assembly group annotation (Bakta, eggNOG, KEGG)`.

### Task 6: Wire into rule all + report
`_t_coassembly()`: add the group prok-functional sentinels when `COASSEMBLY_BINNING and not LONG_READS`. report.smk optional input + load_coassembly (AMR/defense counts per group MAG set) + coassembly.js (a group AMR/defense summary). Commit `feat: wire co-assembly group prok-functional into rule all + report`.

---

## Self-Review
- All per-sample prok-functional rules covered (Tasks 1–5); wired (Task 6).
- `.fna` extension for VAMB bins flagged in Global Constraints + Task 1.
- Each task mirrors named real rules with input/output/wildcard changes (the proven pattern); DB paths kept.
- Short-read only (co-binning gate). Off by default. Each verified by dry-run.
