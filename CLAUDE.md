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
  - `assembly.smk`: Genome assembly. One assembler per track, both chosen on 2026-08-18: MEGAHIT for short reads (metaSPAdes/metaviralSPAdes removed), and for long reads a single assembler picked by `lr_tech` — Flye + Medaka polishing for ONT, metaMDBG for HiFi (hifiasm and the `merge_lr` step were removed).
  - `cobra.smk`: Optional COBRA contig extension — a single run with the MEGAHIT k-mer params (the per-assembler `cobra_spades` twin went out with metaSPAdes).
  - `merge_dedup.smk`: Prefixing/length-filtering the MEGAHIT contigs (`merge_contigs`, kept as a stage even with one assembler) and deduplication using MMseqs2 `easy-linclust` to create a central reference FASTA.
  - `quast.smk`: Assembly quality assessment with QUAST.
  - `viral_detection.smk`: Viral sequence detection with VirSorter2 and geNomad, plus consensus generation. DeepVirFinder and CenoteTaker3 are NOT wired in; VIBRANT was removed on 2026-08-18.
  - `mapping.smk`: Read mapping to contigs (BWA-MEM2 for short reads; minimap2 for long reads) and coverage via `jgi_summarize_bam_contig_depths`.
  - `viral_binning.smk`: Viral binning with CheckV quality assessment and vRhyme clustering, followed by a second CheckV pass on the bins.
  - `votu_catalog.smk`: Global vOTU catalog — pools all viral sets with source-prefixed IDs, a single-pass `skani triangle --sparse`, ICTV-standard clustering (95% ANI + 85% AF), two representative tiers (`all`, `mq`), and vOTU x sample presence/abundance matrices from read recruitment. Replaces the former per-sample clustering. Also runs, once for the whole catalog (principle "(h)" in `docs/ROADMAP_SIMPLIFICACAO.md`): prodigal-gv gene calling, MMseqs2 taxonomy, pharokka + phold annotation, genome maps, and — since 2026-08-18 — `votu_defensefinder_viral`/`votu_dbapis_viral` (anti-defense on viral ORFs, moved off per-sample/per-group `defense_amr.smk` rules because both already consumed this file's global `.faa`; the per-group run had been silently writing catalog-wide systems into a group-scoped path).
  - `prok_binning.smk`: Prokaryotic binning with MetaBAT2 + SemiBin2, consolidated by Binette, then CheckM2 (quality), GUNC (chimerism), galah (dereplication) and GTDB-Tk (taxonomy). MaxBin2 is not used; VAMB appears only in `coassembly.smk`.
  - `taxonomy.smk`: MMseqs2 taxonomy for prokaryotic bins (`mmseqs_taxonomy_prok`) and `viral_taxonomy` — which since 2026-08-18 is no longer a computation but a per-sample/per-group **view** over the global catalog table, joining through `provenance.tsv` and writing bare member IDs. Viral prodigal-gv gene calling and the MMseqs2 viral tiers (INPHARED + optional custom DB) live in `votu_catalog.smk`. vConTACT3 was removed from the pipeline on 2026-08-17. Diamond is not used here — it appears in `votu_catalog.smk` (dbAPIS).
  - `host_prediction.smk`: Phage-host prediction using PHIST.
  - `annotation.smk`: bakta + eggNOG-mapper on prokaryotic MAGs and prokaryotic genome maps. pharokka, phold and the phage/virus genome maps moved to `votu_catalog.smk` on 2026-08-18 (principle "(h)") — they run once over the global catalog, not per sample.
  - `defense_amr.smk`: Defense systems (DefenseFinder with its built-in `--antidefensefinder` pass) and AMR (AMRFinderPlus, RGI/CARD, DeepARG) on prokaryotic bins, plus ABRicate (VFDB+PlasmidFinder) and argNorm, and computes defense islands. There is no PADLOC rule. Anti-defense detection on viral ORFs (DefenseFinder `--antidefensefinder` + dbAPIS/Diamond) moved OUT of this file on 2026-08-18 — it now runs once for the whole vOTU catalog as `votu_defensefinder_viral`/`votu_dbapis_viral` in `votu_catalog.smk` (second half of principle "(h)"), not per sample/group here.
  - `abundance.smk`: CoverM abundance (viral and prokaryotic) and diversity metrics.
  - `reads_classify.smk`: Assembly-independent read profiling with sylph + sylph-tax, prevalence filtering, OTU table and host collapsing. Optional BACPHLIP lifestyle prediction, which runs only over sylph-detected reference genomes (requires `reads_classify_genome_fasta`) — not over assembled vOTUs.
  - `coassembly.smk`: Mirrors most of the per-sample workflow for co-assembly groups, with rules prefixed `coassembly_`/`multisplit_`. Adds VAMB co-binning.
  - `finalize.smk`: Organizing and finalizing output files.
  - `report.smk`: Generating interactive HTML reports (ECharts + D3 + Plotly) and MultiQC aggregation.
- **scripts/**: Auxiliary Python scripts. Notable ones:
  - `generate_report.py`: Creates the standalone interactive HTML report.
  - `votu_catalog.py`, `make_votu_table.py`: vOTU catalog clustering and tables.
  - `prepare_mmseqs_taxdb.py`, `prepare_diamond_db.py`: Custom taxonomy database preparation.
  - `consolidate_amr.py`, `compute_diversity.py`, `genome_map*.py`, `split_viral_fastas.py`, `pin_containers.py`.
  - `reads_classify/`: helpers for the sylph track (`bacphlip_lifestyle.py`, `make_otu.py`, `filter_by_prevalence.py`, `collapse_by_host.py`, `build_imgvr_taxonomy.py`).
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
| `MIN_CONTIG` | Minimum contig length post-assembly (bp) |
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

## Central Hub: `rep_seq.fasta`

The deduplicated representative sequences (`{sample}_rep_seq.fasta`, generated via MMseqs2 at 95% identity) serve as the central reference for all downstream analyses. This ensures consistency across viral detection, mapping, binning, taxonomy, and host prediction.

## Dependencies and Environments

The pipeline uses 25 isolated conda environment files in `envs/`, all named `env_*.yaml`. Rules also carry a per-rule `container:` resolved from `containers.yaml`; there is no global `containerized:` directive. Key environments include:
- `env_qc`: Quality control tools.
- `env_assembly`: Assemblers + MMseqs2 (dedup, plus the MMseqs2 taxonomy rules).
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
