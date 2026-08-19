# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

A comprehensive Snakemake-based metagenomics and virome pipeline for analyzing viral and prokaryotic communities from environmental samples. It supports both short-read (Illumina paired-end) and long-read (Nanopore ONT / PacBio HiFi) sequencing data. The pipeline is implemented in Snakemake DSL with Python scripts, utilizing isolated conda environments for all dependencies, without requiring an installable Python package.

The pipeline processes raw sequencing data through multiple stages: quality control, assembly, deduplication, viral detection, binning (both viral and prokaryotic), taxonomy assignment, host prediction, and final reporting. It integrates numerous bioinformatics tools to provide a complete analysis workflow for metagenomic studies.

## Project Structure

The repository is organized as follows:

- **Snakefile**: Main workflow file containing configuration, sample discovery, and rule inclusions.
- **rules/**: Directory containing modular Snakemake rules, each handling a specific functional block:
  - `qc.smk`: Quality control and trimming (fastp for short reads; NanoPlot, porechop_abi, Filtlong for long reads). There is no FastQC or Trim Galore rule.
  - `host_removal.smk`: Optional host decontamination (bwa-mem2 for short reads, minimap2 for long reads, then samtools filtering). Active only when `host_genome` is set.
  - `assembly.smk`: Genome assembly. One assembler per track, both chosen on 2026-08-18: MEGAHIT for short reads (metaSPAdes/metaviralSPAdes removed), and for long reads a single assembler picked by `lr_tech` — Flye + Medaka polishing for ONT, metaMDBG for HiFi (hifiasm and the `merge_lr` step were removed). **Flye has no `--min-contig-len`**, so `scripts/filter_min_length.py` applies `min_contig` after it (both `flye_lr` and `flye_coassembly`) — until 2026-08-19 ONT was the one track where that config key silently did nothing. metaMDBG's output name varies by version and recent ones write `contigs.fasta.gz`; the rule covers both and records a real `failed:` status instead of promoting an empty FASTA.
  - `cobra.smk`: Optional COBRA contig extension — a single run with the MEGAHIT k-mer params (the per-assembler `cobra_spades` twin went out with metaSPAdes).
  - `quast.smk`: Assembly quality assessment with QUAST.
  - `viral_detection.smk`: Viral sequence detection with VirSorter2 and geNomad, plus consensus generation. DeepVirFinder and CenoteTaker3 are NOT wired in; VIBRANT was removed on 2026-08-18.
  - `mapping.smk`: Read mapping to contigs (BWA-MEM2 for short reads; minimap2 for long reads) and coverage via `jgi_summarize_bam_contig_depths`.
  - `viral_binning.smk`: Viral binning with CheckV quality assessment and vRhyme clustering, followed by a second CheckV pass on the bins. `rule viral_nonredundant` (bins-first) applies, since 2026-08-18, the composite length/quality/bin gate of item (e) in `docs/ROADMAP_SIMPLIFICACAO.md` to unbinned sequences: a short contig with no vRhyme bin is kept only if CheckV calls it Complete/High-quality/Medium-quality (the fixed `MQ_TIERS` in `scripts/viral_length_gate.py` — deliberately NOT the configurable `VIRAL_KEEP_TIERS`, which ships expanded to all five tiers and would make the arm a no-op) or is `>= VIRAL_MIN_CONTIG` bp — the group co-assembly equivalent is `coassembly_viral_nonredundant` in `coassembly.smk`. Dropped sequences are not discarded silently: both rules also write `{sample|group}_viral_discarded.fasta` (bare `contig_id` headers) and a `_viral_discarded.tsv` sidecar with the CheckV evidence per contig, copied to `final/viral/`. **CheckV provirus headers are `{contig}_{n}`** (e.g. `>k141_219139_1 1-13933/18998`, the range being a *description* after the space) — NOT `orig_id|start_end`, which a comment here claimed until 2026-08-19 and which made the trimming silently never happen, emitting host flanks as viral. Resolve them only through `scripts/checkv_provirus.py`, which matches candidates against known contig IDs instead of guessing a delimiter. geNomad uses a different convention, `contig|provirus_START_END` — do not mix the two.
  - `votu_catalog.smk`: Global vOTU catalog — pools all viral sets with source-prefixed IDs, a single-pass `skani triangle --sparse`, ICTV-standard clustering (95% ANI + 85% AF), two representative tiers (`all`, `mq`), and vOTU x sample presence/abundance matrices from read recruitment. Replaces the former per-sample clustering. Also runs, once for the whole catalog (principle "(h)" in `docs/ROADMAP_SIMPLIFICACAO.md`): prodigal-gv gene calling (`votu_prodigal` — it really does call `prodigal-gv` since 2026-08-19; before that it called plain `prodigal` while three docs claimed otherwise, truncating ORFs in TAG-recoding phages), MMseqs2 taxonomy, pharokka + phold annotation, and — since 2026-08-18 — `votu_defensefinder_viral`/`votu_dbapis_viral` (anti-defense on viral ORFs, moved off per-sample/per-group `defense_amr.smk` rules because both already consumed this file's global `.faa`; the per-group run had been silently writing catalog-wide systems into a group-scoped path).
  - `prok_binning.smk`: Prokaryotic binning with MetaBAT2 + SemiBin2, consolidated by Binette, then CheckM2 (quality) and GUNC (chimerism). MaxBin2 is not used; VAMB appears only in `coassembly.smk`. **Dereplication, GTDB-Tk and `prok_bin_proteins` left this file on 2026-08-19** — see `mag_catalog.smk`. The `_read_manifest`/`_concat_proteins` helpers stayed: the global defense/AMR rules use them.
  - `mag_catalog.smk`: Global prokaryotic MAG catalog — the prokaryotic analogue of the vOTU catalog, same principle as item "(h)": compute on the representative, inherit on the member. Pools every sample's Binette bins and every group's VAMB bins under namespaced filenames (`{source}__{bin}.fa`, since Binette emits `binette_bin1` in *every* sample and VAMB emits bare integers — pooling without a prefix would silently overwrite different organisms), re-keys the CheckM2 reports onto those IDs (galah matches `--checkm2-quality-report` by genome name; concatenating the raw reports would leave galah with quality for *no* genome and it would pick representatives by some other criterion without saying so), runs **one** global `galah cluster`, and writes `provenance.tsv` + `mag_membership.tsv`. `mag_catalog_gtdbtk` runs `classify_wf` once over the representatives — 39 runs (32 samples + 7 groups) became 1. `rule gtdbtk` and `gtdbtk_group` are now **views** over the global table, joined through membership and writing the bin's ORIGINAL name, so `phist`, `finalize.smk` and the report keep their contract. Pure helpers live in `scripts/mag_catalog.py` (21 tests). `low_depth_mode` — and with it the `contigs_pseudogenome` path of `prok_bin_proteins` — **was removed from the tool on 2026-08-19**; the last version that has it is the git tag `v-lowdepth`. The catalog is unconditional. On the same day the defense/AMR/prok-taxonomy analyses migrated onto the representatives: `prok_bin_proteins`, `defensefinder`, `amrfinderplus`, `rgi_card`, `deeparg`, `abricate`, `argnorm_normalize`, `amr_consensus` and `mmseqs_taxonomy_prok` — plus their seven `coassembly_*` twins — were **deleted**, and each now runs once as `mag_catalog_proteins` / `mag_*` over the catalog (`defense_amr.smk`, `taxonomy.smk`). `bakta`, `eggnog_prok` and `extract_kegg_kos` followed on the same day (`mag_bakta`, `mag_eggnog_prok`, `mag_extract_kegg_kos` in `annotation.smk`), so **nothing prokaryotic downstream of binning is per-sample any more** — only binning, CheckM2 and GUNC are, and they must be (the coverage that separates bins is that sample's). The `multisplit` track stays out of the catalog.
  - **The views** (`mag_views_sample` / `mag_views_group`, at the end of `defense_amr.smk` — they live there only because they reference `rules.mag_*`, which the include order defines after `mag_catalog.smk`) are the ONLY thing that writes `{sample}/bins/...` and `coassembly/{group}/bins/...`. Same paths as before, so `finalize.smk` and the report keep their contract; no `ruleorder` is needed because there is no second producer. A view is INHERITANCE, not identity: every bin gets its representative's row (filtering the global table by source prefix would return almost nothing — that was the `viral_taxonomy` bug of 2026-08-18). Two variants, because the tools disagree on where the genome lives: a `genome` column (DefenseFinder, ABRicate) and a prefix on the ID (AMRFinderPlus `Protein identifier`, RGI `ORF_ID`, DeepARG `#ARG`, MMseqs2 `qseqid`, argNorm, `locus` in the consensus). **Never cut the ID at the first `__`**: a catalog protein is `S1__binette_bin1__k141_1_5`, so that cut yields `S1` and attributes every AMR hit to the SAMPLE instead of the MAG — `resolve_prefixed_id` matches against known representatives instead. After the view rewrites the prefix to the original bin name, the report's first-`__` cut is correct again. The per-source protein manifest points at the REPRESENTATIVE's `.faa` on purpose: the defense tables carry its protein IDs, and that is what `compute_defense_islands` matches against.
  - `taxonomy.smk`: MMseqs2 taxonomy for prokaryotic MAGs (`mag_mmseqs_taxonomy_prok` — global over the catalog representatives since 2026-08-19, distributed by `mag_views_sample`) and `viral_taxonomy` — which since 2026-08-18 is no longer a computation but a per-sample/per-group **view** over the global catalog table, joining through `provenance.tsv` and writing bare member IDs. Viral prodigal-gv gene calling and the MMseqs2 viral tiers (INPHARED + optional custom DB) live in `votu_catalog.smk`. vConTACT3 was removed from the pipeline on 2026-08-17. Diamond is not used here — it appears in `votu_catalog.smk` (dbAPIS).
  - `host_prediction.smk`: Phage-host prediction using PHIST — over the **global vOTU catalog** (`votu_catalog_reps.mq_fasta`) against the sample's own MAGs. `coassembly_phist` does the same over the group's VAMB MAGs — it used the group's own `votu_mq_reps.fasta` until 2026-08-19, when the group-local skani chain was removed; both tracks now share one viral ID space. PHIST is really `kmer-db build` + `kmer-db new2all` + `phist`, and kmer-db emits one result row per input **file**, which is why `split_viral_fastas.py` exists: both sides need one FASTA per genome plus a list file. Note these results describe assembled contigs (`k141_…`); they share no key with the sylph track's IMG/VR reference genomes (`t__IMGVR_UViG_…`).
  - `annotation.smk`: bakta + eggNOG-mapper + KO extraction on prokaryotic MAGs — **global over the catalog representatives since 2026-08-19**. `mag_eggnog_prok` prefixes every protein with `{genome}__` when concatenating the bakta FAAs (bakta's own IDs are locus tags like `LLOGBO_00001`, which carry no link to the MAG) and writes `eggnog_annotations.tsv` with emapper's `##` preamble stripped — a plain DictReader took `## <date>` as the header and returned None for every column, which is why the report's COG chart counted every protein as "Function unknown". `ko_per_mag.tsv`'s `mag` column is the real genome now, resolved against known genomes, not the regex-derived locus-tag prefix it used to be. pharokka and phold moved to `votu_catalog.smk` on 2026-08-18 (principle "(h)"). **The genome maps were removed entirely on 2026-08-19** (`genome_map_prok`, `votu_genome_map_phage`, `votu_genome_map_virus`, `scripts/genome_map*.py`, the report panel, the `genome_map_*` config keys) — they took with them `votu_catalog_quality_summary` and `votu_catalog_genomad_genes`, which existed only to feed them.
  - `defense_amr.smk`: Defense systems (DefenseFinder with its built-in `--antidefensefinder` pass) and AMR (AMRFinderPlus, RGI/CARD, DeepARG) on prokaryotic MAGs, plus ABRicate (VFDB+PlasmidFinder) and argNorm. **Every rule here is global since 2026-08-19** (`mag_*`, one run over the catalog representatives), and the file ends with the two per-source view rules — see `mag_catalog.smk`. There is no PADLOC rule. Anti-defense detection on viral ORFs (DefenseFinder `--antidefensefinder` + dbAPIS/Diamond) moved OUT of this file on 2026-08-18 — it now runs once for the whole vOTU catalog as `votu_defensefinder_viral`/`votu_dbapis_viral` in `votu_catalog.smk` (second half of principle "(h)"), not per sample/group here.
  - `abundance.smk`: CoverM abundance (viral and prokaryotic) and diversity metrics. `votu_abundance`'s `representative` column is **always namespaced** (`{source}|{contig}`) and the representative's length comes from the catalog's reps FASTA — until 2026-08-19 it was half-stripped, which zeroed the metric for every vOTU represented by another sample and broke the report's join for the local ones. In `compute_diversity.py`, **Simpson and Chao1 are computed on read COUNTS, never on RPKM** (both are count estimators: f1/f2 and `a*(a-1)`); without counts they are written empty and the report skips them instead of plotting 0.
  - `reads_classify.smk`: Assembly-independent read profiling with sylph + sylph-tax, prevalence filtering, OTU table and host collapsing. **Neither BACPHLIP nor PHIST run in this track** (removed 2026-08-19): both need the reference genome *sequences*, and only sylph's `.syldb` k-mer sketches exist on disk — a sketch cannot yield sequence. They would also be redundant: `bacphlip_votu` (`votu_catalog.smk`) already predicts lifestyle over the user's own vOTUs and `rule phist` already predicts hosts for them. Host annotation here comes from the **database** instead, via `rule reads_host_map` (`scripts/reads_classify/build_host_map.py`), which recovers the `Virus_host (if viral)` column that `sylph-tax merge` discards. `collapse_by_host.py` selects viral rows by realm prefix (`r__`, since IMG/VR lineages have no `d__Viruses` and never reach `s__`) and keeps hierarchy leaves, not a fixed rank. Its output carries `host_source` (`db`/`none`) plus a `viral_host_assignments.tsv` sidecar. **The merged table's column headers are the reads FILE PATHS sylph was given** (`Sample_file`), which in paired-end runs is the `-1` file — `{sample}_R1.fastq.gz`, not `{sample}.fastq.gz`. Resolve them to sample names only through `load_reads_classify`'s `_fastq_to_sample(col, samples)`, which matches against the known sample list; stripping the extension alone silently zeroed the entire reads track on every PE run until 2026-08-19. `sylph_merge`/`reads_host_map` list this run's sample directories instead of globbing `{OUTDIR}/*/` (a stale `.sylphmpa` from a renamed sample used to join the merge). **sylph profiles the TRIMMED reads, and the host-removed ones when `host_genome` is set** (same `_clean_r1`/`_clean_r2`/`_clean_lr` helpers as the assembly) — so a reads-only run still triggers fastp/filtlong and, if configured, host removal. It read the raw FASTQs until 2026-08-19, which meant adapters and host reads entered the profile and `host_genome` had no effect here.
  - `coassembly.smk`: Mirrors most of the per-sample workflow for co-assembly groups, with rules prefixed `coassembly_`/`multisplit_`. Adds VAMB co-binning. `rule coassembly_viral_nonredundant` (short-read groups only — long-read co-assembly has no group-level vRhyme) mirrors `viral_nonredundant`'s bins-first + composite length/quality/bin gate over the already CheckV-trimmed, already vRhyme-binned group set, and is what `_catalog_sources()` (`votu_catalog.smk`) actually pools for short-read groups since 2026-08-18 — before that, group vMAGs from `coassembly_vrhyme` never reached the global vOTU catalog. Long-read groups still source from `coassembly_viral_trimmed`'s pre-binning output. **The group's own skani vOTU chain was REMOVED on 2026-08-19** (`coassembly_skani_votu`/`coassembly_skani_cluster`/`coassembly_viral_votu_reps`, −279 lines): it clustered the very fasta the global catalog already pools, so the run carried two vOTU definitions with different ID spaces and possibly different representatives. Group vOTU counts, representative lengths and the accumulation curve now come from the catalog filtered to that source (`load_catalog_clusters_by_source`, which returns the BARE member — it joins the group's VAMB abundance matrix — and the NAMESPACED representative). The `use_votu` flag went with it: nothing read it, and the catalog is not an optional step.
  - `finalize.smk`: Organizing and finalizing output files. Per-sample `rule organize_outputs` (and its group twin `coassembly_organize_outputs`) copy the item-(e) discard sidecar into `final/viral/viral_discarded.{fasta,tsv}`, alongside `viral_nonredundant.fasta`.
  - `report.smk`: Generating interactive HTML reports (ECharts + D3 + Plotly) and MultiQC aggregation.
- **scripts/**: Auxiliary Python scripts. Notable ones:
  - `generate_report.py`: Creates the standalone interactive HTML report.
  - `votu_catalog.py`, `make_votu_table.py`: vOTU catalog clustering and tables.
  - `prepare_mmseqs_taxdb.py`, `prepare_diamond_db.py`: Custom taxonomy database preparation.
  - `consolidate_amr.py`, `compute_diversity.py`, `genome_map*.py`, `split_viral_fastas.py`, `pin_containers.py`.
  - `reads_classify/`: helpers for the sylph track (`make_otu.py`, `filter_by_prevalence.py`, `collapse_by_host.py`, `build_host_map.py`, `build_imgvr_taxonomy.py`). `bacphlip_lifestyle.py` was removed on 2026-08-19 with the rule that called it.
- **config.yaml**: All runtime parameters (see Configuration below).
- **containers.yaml**: Per-tool container images consumed via `CONTAINERS.get(<name>)`.
- **INSTALL.md**, **README.md**: Installation guide and project documentation.

For a full per-stage inventory of every tool the pipeline actually runs, see `docs/VAPOR_TOOLS_MAP.md` (derived from the rules, not from prose).

## Configuration

All runtime parameters live in **`config.yaml`**, loaded by `configfile: "config.yaml"` (`Snakefile:82`). The block that follows in the `Snakefile` only unpacks `config[...]` into module-level constants — edit `config.yaml`, not the `Snakefile`. Key variables include:

| Variable | Purpose |
|---|---|
| `FASTQ_DIR` | Input directory for raw FASTQ files |
| `OUTDIR` | Output root directory |
| `THREADS` | Default CPU threads per rule |
| `LONG_READS` | Boolean: `True` for ONT/HiFi long reads, `False` for Illumina PE |
| `MIN_CONTIG` | Minimum contig length post-assembly (bp). Default 1000 — this is no longer viral detection's length gate (item (e) below); binners clamp their own floor when it's set below it (vRhyme 2000, MetaBAT2 1500). |
| `VIRAL_MIN_CONTIG` | Length-only arm of the composite viral length/quality/bin gate applied AFTER vRhyme (`rule viral_nonredundant` / `coassembly_viral_nonredundant`). Default 5000 bp — see item (e) of `docs/ROADMAP_SIMPLIFICACAO.md`. |
| `VIRAL_CONSENSUS_MODE` | Consensus strategy for viral detection: `"count"`, `"score"`, or `"hybrid"` |
| `MEGAHIT_MEM` | RAM limit for MEGAHIT (bytes) |
| Database paths | Paths to required databases (CheckV, VirSorter2, GeNomad, etc.) |
| Custom Diamond DBs | Optional custom databases for enhanced taxonomy |

Additional parameters control tool-specific settings, such as minimum scores for viral tools, binning environments, and long-read specific options.

## Sample Discovery

`find_samples()` in the `Snakefile` detects samples from FASTQ file naming:
- **Short reads, paired**: `*_R1*/*_R2*` or `*_1.*/*_2.*` — both conventions are supported.
- **Short reads, single**: `*.fq.gz` / `*.fastq.gz` with no `_R1`/`_R2` suffix (`single_end: true`).
- **Long reads**: single `*.fastq.gz` / `*.fq.gz`.

Samples are stored in the `SAMPLES` dictionary. For co-assembly, `GROUPS` is built separately from `coassembly_metadata` / `coassembly_grouping`.

## Running the Pipeline

Execute from the project directory containing `Snakefile` and input data.

```bash
# Activate Snakemake environment
conda activate snakemake

# Dry-run to validate workflow
snakemake -n --use-conda --cores 32

# Full execution
snakemake --use-conda --cores 32

# Visualize directed acyclic graph (DAG)
snakemake --dag | dot -Tsvg > dag.svg

# Force re-run a specific rule
snakemake --use-conda --cores 32 --forcerun <rule_name>

# Run up to a specific target
snakemake --use-conda --cores 32 results/<sample>/viral/taxonomy/taxonomy_done.txt
```

## Central Hub: the assembly

Each sample's assembly IS the central reference for every downstream analysis
— viral detection, mapping, binning, taxonomy, host prediction. Resolved by the
`_sample_contigs()` helper in the `Snakefile`: MEGAHIT's `final.contigs.fa` for
short reads, the Flye+Medaka or metaMDBG output for long reads.

There is no deduplication step. `merge_contigs` + MMseqs2 `easy-linclust`
(`rep_seq.fasta`) were removed on 2026-08-18 (item "(d)" in
`docs/ROADMAP_SIMPLIFICACAO.md`): `easy-linclust` existed to collapse redundancy
**between assemblers**, and items (b)/(c) left one assembler per track, so that
source of redundancy no longer exists. The co-assembly track had always run this
way (`coassembly/{group}/contigs.fa` straight into everything), so this aligned
the two tracks rather than inventing a new pattern.

Contig IDs are now the assembler's own (`k141_10`, no `MEGAHIT_` prefix).
Uniqueness across samples is guaranteed where it matters — the global vOTU
catalog namespaces every contig by its source (`votu_catalog_pool`).

## Dependencies and Environments

The pipeline uses 25 isolated conda environment files in `envs/`, all named `env_*.yaml`. Rules also carry a per-rule `container:` resolved from `containers.yaml`; there is no global `containerized:` directive. Key environments include:
- `env_qc`: Quality control tools.
- `env_assembly`: Assemblers + MMseqs2 (used only by the MMseqs2 taxonomy rules — the dedup step is gone).
- `env_mapping`: Mapping tools.
- `env_viral`: Viral detection and taxonomy.
- And more for binning, reporting, etc.

All environments are created via `snakemake --use-conda --cores 1 --create-envs-only`.

## Databases

Requires ~500 GB of pre-downloaded databases:
- CheckV, VirSorter2, geNomad, INPHARED, CheckM2, GTDB-Tk, pharokka, phold, bakta, eggNOG, CARD, DeepARG, DefenseFinder models, dbAPIS, sylph-tax.
- Optional custom MMseqs2 databases for improved taxonomy (`custom_prok_mmseqs_db`, `custom_viral_mmseqs_db`).

Paths must be set in `config.yaml`.

## Setup and Installation

Refer to `INSTALL.md` for step-by-step instructions:
1. Install Miniforge3.
2. Create Snakemake environment.
3. Build all conda environments.
4. Download and configure databases.
5. Edit `config.yaml`.

System requirements: Linux, 16+ cores, 64+ GB RAM, 600+ GB disk.

## Key Features

- **Modular Design**: Each step is a separate Snakemake rule for flexibility.
- **Multi-Platform Support**: Handles Illumina and long-read data seamlessly.
- **Comprehensive Analysis**: Covers QC, assembly, viral/prokaryotic binning, taxonomy, host prediction.
- **Consensus Approaches**: Viral detection uses multi-tool consensus.
- **Interactive Reporting**: Generates detailed HTML reports with ECharts/D3 visualizations.
- **Reproducibility**: Environment locking and version pinning.

## Troubleshooting

- Ensure all database paths are correctly set.
- Use dry-run to check for errors before full execution.
- Monitor RAM usage for memory-intensive steps.
- For long reads, set `long_reads: true` in `config.yaml` and specify `lr_tech`.

This pipeline is designed for high-throughput metagenomic analysis, providing end-to-end processing from raw reads to curated results and reports.

## HTML Report — `scripts/report/` package

**The old monolithic `generate_report.py` (~2974 lines) no longer exists.** It was refactored into a package; `scripts/generate_report.py` is now a 19-line entry point that does `from report.renderer import build_report`.

### Layout

| Path | Lines | Role |
|---|---|---|
| `scripts/generate_report.py` | 19 | Snakemake entry point only. Adds `scripts/` to `sys.path`, calls `build_report(snakemake)`. |
| `scripts/report/data_loaders.py` | ~1919 | All parsing/loading. ~65 `load_*` / `parse_*` / `collect_*` functions, one per data source. |
| `scripts/report/renderer.py` | ~496 | Orchestration: calls the loaders, assembles the data dict, injects into the HTML shell. |
| `scripts/report/components/shell.html` | ~611 | Page skeleton with `{{CSS}}`, `{{ECHARTS_JS}}`, `{{D3_JS}}`, `{{DATA_JSON}}`, `{{APP_JS}}` placeholders. |
| `scripts/report/components/base.css` | ~674 | Design system. |
| `scripts/report/components/app.js` | ~783 | Tab routing, shared helpers, bootstrap. |
| `scripts/report/components/*.js` | — | One file per tab: `overview`, `sequencing`, `viral`, `prokaryotic`, `annotation`, `hostdefense`, `diversity`, `reads_classify`, `coassembly`, `about`, `export`. |
| `scripts/report/assets/` | — | Vendored `echarts.min.js` and `d3.min.js`, inlined at build time for a standalone HTML. |

### Critical coding rules

- **There are no f-strings in `renderer.py`.** Substitution is explicit `str.replace()` against `{{PLACEHOLDER}}` tokens (`renderer.py:478-482`). The old rule about doubling every `{`/`}` in the HTML/JS block is obsolete — **do not** double braces in the component files.
- **Edit JS and CSS in `components/*.js` and `components/base.css`**, as plain files. No Python escaping applies to them, so `\t` and `\n` are normal JS escapes again.
- **Charting is ECharts + D3**, both vendored in `assets/`. Plotly is no longer used by the report.
- `_jsstr()` escapes `</` as `<\/` when serializing JSON into a `<script>` tag — keep that when adding new data payloads.
- To add a data source: write a `load_*` in `data_loaders.py`, wire it in `_build()` in `renderer.py`, add it to the data dict, then consume it from the relevant `components/*.js`.
- `load_tool_status()` reads per-rule `done.txt` status lines (`ok` / `skipped: <reason>` / `failed: <reason>`). A failed tool must render as a gap, never as a count of 0 — an empty `done.txt` once made a disk-full AMRFinderPlus run read as a biological zero.
