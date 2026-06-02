# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Project Is

A comprehensive Snakemake-based metagenomics and virome pipeline for analyzing viral and prokaryotic communities from environmental samples. It supports both short-read (Illumina paired-end) and long-read (Nanopore ONT / PacBio HiFi) sequencing data. The pipeline is implemented in Snakemake DSL with Python scripts, utilizing isolated conda environments for all dependencies, without requiring an installable Python package.

The pipeline processes raw sequencing data through multiple stages: quality control, assembly, deduplication, viral detection, binning (both viral and prokaryotic), taxonomy assignment, host prediction, and final reporting. It integrates numerous bioinformatics tools to provide a complete analysis workflow for metagenomic studies.

## Project Structure

The repository is organized as follows:

- **Snakefile**: Main workflow file containing configuration, sample discovery, and rule inclusions.
- **rules/**: Directory containing modular Snakemake rules, each handling a specific functional block:
  - `qc.smk`: Quality control and trimming (FastQC, Trim Galore for short reads; NanoPlot, Porechop, Filtlong for long reads).
  - `assembly.smk`: Genome assembly (MEGAHIT + metaSPAdes for short reads; Flye + hifiasm + Medaka for long reads).
  - `merge_dedup.smk`: Merging assemblies and deduplication using MMseqs2 to create a central reference FASTA.
  - `quast.smk`: Assembly quality assessment with QUAST.
  - `viral_detection.smk`: Viral sequence detection using multiple tools (VirSorter2, GeNomad, DeepVirFinder, CenoteTaker3, VIBRANT) and consensus generation.
  - `mapping.smk`: Read mapping to contigs (BWA-MEM2 for short reads; minimap2 for long reads) and coverage calculation.
  - `viral_binning.smk`: Viral binning with CheckV quality assessment and vRhyme clustering.
  - `prok_binning.smk`: Prokaryotic binning using MetaBAT2, MaxBin2, VAMB, SemiBin2, followed by Binette consolidation, CheckM2 quality check, and GTDB-Tk taxonomy.
  - `taxonomy.smk`: Taxonomy assignment via Prodigal gene prediction and Diamond BLAST against databases, including viral taxonomy.
  - `host_prediction.smk`: Phage-host prediction using PHIST and iPHoP.
  - `finalize.smk`: Organizing and finalizing output files.
  - `report.smk`: Generating interactive HTML reports with Plotly and MultiQC aggregation.
- **scripts/**: Auxiliary Python scripts:
  - `filter_checkv_hq.py`: Filters viral FASTAs based on CheckV quality tiers.
  - `generate_report.py`: Creates the standalone interactive HTML report.
  - `merge_lr_assemblies.py`: Merges long-read assemblies from Flye and hifiasm.
  - `prepare_diamond_db.py`: Prepares custom Diamond databases for taxonomy.
  - `split_viral_fastas.py`: Splits viral FASTAs for per-genome analysis.
- **INSTALL.md**: Detailed installation guide for dependencies and databases.
- **README-new.md** and **README-old.md**: General project documentation.
- **skills-lock.json**: Environment lock file for reproducibility.

## Configuration

All runtime parameters are defined in the configuration block at the top of `Snakefile` (lines 46–121). Key variables include:

| Variable | Purpose |
|---|---|
| `FASTQ_DIR` | Input directory for raw FASTQ files |
| `OUTDIR` | Output root directory |
| `THREADS` | Default CPU threads per rule |
| `LONG_READS` | Boolean: `True` for ONT/HiFi long reads, `False` for Illumina PE |
| `MIN_CONTIG` | Minimum contig length post-assembly (bp) |
| `VIRAL_CONSENSUS_MODE` | Consensus strategy for viral detection: `"count"`, `"score"`, or `"hybrid"` |
| `SPADES_MEM` / `MEGAHIT_MEM` | RAM limits for assemblers (GB or bytes) |
| Database paths | Paths to required databases (CheckV, VirSorter2, GeNomad, etc.) |
| Custom Diamond DBs | Optional custom databases for enhanced taxonomy |

Additional parameters control tool-specific settings, such as minimum scores for viral tools, binning environments, and long-read specific options.

## Sample Discovery

The pipeline automatically detects samples based on FASTQ file naming:
- **Short reads**: Paired files matching `*_R1*.fq.gz` and `*_R2*.fq.gz`.
- **Long reads**: Single files matching `*.fastq.gz` or `*.fq.gz`.

Samples are stored in the `SAMPLES` dictionary, enabling dynamic workflow execution.

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

The pipeline uses 23 isolated conda environments for tool isolation. Key environments include:
- `env_qc`: Quality control tools.
- `env_assembly`: Assemblers.
- `env_mmseqs`: Deduplication.
- `env_mapping`: Mapping tools.
- `env_viral`: Viral detection and taxonomy.
- And more for binning, reporting, etc.

All environments are created via `snakemake --use-conda --cores 1 --create-envs-only`.

## Databases

Requires ~500 GB of pre-downloaded databases:
- CheckV, VirSorter2, GeNomad, VIBRANT, CenoteTaker3, INPHARED, vConTACT3, iPHoP, CheckM2, GTDB-Tk.
- Optional custom Diamond databases for improved taxonomy.

Paths must be set in `Snakefile`.

## Setup and Installation

Refer to `INSTALL.md` for step-by-step instructions:
1. Install Miniforge3.
2. Create Snakemake environment.
3. Build all conda environments.
4. Download and configure databases.
5. Edit `Snakefile` configuration.

System requirements: Linux, 16+ cores, 64+ GB RAM, 600+ GB disk.

## Key Features

- **Modular Design**: Each step is a separate Snakemake rule for flexibility.
- **Multi-Platform Support**: Handles Illumina and long-read data seamlessly.
- **Comprehensive Analysis**: Covers QC, assembly, viral/prokaryotic binning, taxonomy, host prediction.
- **Consensus Approaches**: Viral detection uses multi-tool consensus.
- **Interactive Reporting**: Generates detailed HTML reports with Plotly visualizations.
- **Reproducibility**: Environment locking and version pinning.

## Troubleshooting

- Ensure all database paths are correctly set.
- Use dry-run to check for errors before full execution.
- Monitor RAM usage for memory-intensive steps.
- For long reads, set `LONG_READS = True` and specify `LR_TECH`.

This pipeline is designed for high-throughput metagenomic analysis, providing end-to-end processing from raw reads to curated results and reports.

## generate_report.py — Line Map (~2974 lines total)

Use this map to call `Read` with `offset`+`limit` instead of reading the whole file.

### Python blocks (data loading + figures)

| Lines | Content |
|---|---|
| 1–22 | Imports, snakemake input/output/params binding |
| 23–52 | `_inp()`, `cfg_params` dict (pipeline config for About tab) |
| 53–125 | Helpers: `safe_float`, `safe_int`, `parse_tsv`, `parse_quast_all`, `parse_support`, `compute_bin_abundance` |
| 127–210 | FastQC + trimming: `parse_fastqc_zip`, `collect_fastqc`, `parse_trim_report`, `parse_mapping_rate` |
| 212–244 | Depth/contig helpers: `collect_depth_data`, `parse_fasta_lengths` |
| 246–410 | Viral data: `collect_viral_tool_counts`, `collect_viral_scores`, `collect_vrhyme_stats`, `collect_binner_counts`, `parse_genomad_taxonomy`, `parse_ct3_taxonomy`, `parse_checkm2_phyla` |
| 420–500 | Tool version detection: `_get_ver`, `_dvf_version`, `_get_ver2`, `tool_versions` dict |
| 503–555 | Data loading loop setup: `_path_dict`, per-sample dict init, main `for sample in samples` loop |
| 556–595 | Design system: `TEAL/AMBER/GREEN/RED/GRAY/PURPLE/CYAN/BLUE/ORANGE/PALETTE`, Plotly template `mite_light` registration |
| 596–630 | FastQC figures: `fig_fq_reads`, `fig_fq_qual`, `fig_fq_gc`, `fig_trim`, `fig_mapping` |
| 631–730 | Assembly figures: `fig_asm_prog`, `fig_assembly`, `fig_contig_len`, `fig_cov` |
| 731–870 | Viral figures: `fig_tool_counts`, `fig_viral`, `fig_scores`, `fig_checkv`, `fig_checkv_vrh`, `fig_checkv_scatter`, `fig_viral_len`, `fig_vrhyme`, `fig_vrhyme_detail` |
| 871–1000 | Bin figures: `fig_binner_total`, `fig_das_tax`, `fig_checkm2` (with HQ/MQ zones + shapes), `fig_cm2_hist`, `fig_bin_size`, `fig_abundance` |
| 1001–1028 | Preliminary taxonomy figures: `fig_tax_viral` (sunburst), `fig_tax_ct3`, `fig_tax_bac` |
| 1029–1048 | Overview dict build; `figs_json` dict (first serialization pass) |
| 1029–1095 | Low-level load helpers: `load_tsv`, `load_csv` |
| 1052–1100 | `load_vibrant()` — VIBRANT scaffold summary + AMG individuals/pathways |
| 1101–1135 | `load_vcontact3()` — vConTACT3 v3 CSV (filters Reference==True) |
| 1136–1188 | `load_viral_taxonomy()` — reads `viral_taxonomy_merged.tsv`, builds unified records |
| 1189–1248 | `load_gtdbtk()` — bac120 + ar53 summary TSVs, GTDB classification parsing |
| 1214–1248 | `load_custom_prok()` — Diamond custom prok hits, majority-vote per bin |
| 1249–1295 | `load_phist()`, `load_iphop()` — host prediction loaders |
| 1296–1385 | Path building for taxonomy/vcontact3/gtdbtk/phist/iphop; data load calls; novelty + MIMAG metrics; overview back-fill |
| 1109–1130 | **vConTACT3 network loading** — `vc3_network_data` dict; `_path_dict(getattr(snakemake.input,'network_json',...))` |
| 1385–1398 | `enrich_taxonomy_with_checkv()` — joins CheckV completeness into tax_data |
| 1399–1432 | `merge_prok_taxonomy()` — GTDB-Tk priority + Diamond fallback + CheckM2 merge |
| 1433–1452 | `load_tool_support_matrix()` — tool agreement matrix for heatmap |
| 1453–1541 | Task 5 figures: `fig_read_funnel` (go.Funnel), `fig_tool_heatmap` (go.Heatmap), `fig_viral_depth` (go.Histogram) |
| 1685–1710 | Final JSON serialization: `tax_json`, `vc3_network_json`, `figs_json_str`; HTML tables |

### HTML/CSS/JS block

| Lines | Content |
|---|---|
| 1624–1755 | `_CSS` — plain string (NOT f-string), full design system. Edit CSS here. No `{{}}` escaping needed. |
| 1756–1765 | `html = (...)` — string concatenation building `<head>` + `<style>` + `_CSS` |
| 1766–1810 | f-string start: `<nav id="sidebar">` — sidebar tabs, theme button |
| 1811–1960 | HTML panels: `panel-overview`, `panel-readqc`, `panel-assembly`, `panel-viral` |
| 1961–2060 | HTML panel: `panel-taxonomy` — Viral Taxonomy (source pie, family bar, sunburst, master table), Novelty section, **vConTACT3 Network section** (`p-vc3-network`), Prokaryotic Taxonomy |
| 2060–2115 | HTML panels: `panel-hostpred`, `panel-bins`, `panel-abundance`, `panel-annotation`, `panel-about` |
| 2220–2240 | JS constants: `SAMPLES`, `OVERVIEW`, `FIGS`, `TAX_DATA`, `VC3_NETWORK`, `GTDB_DATA`, etc. |
| 2240–2360 | JS helpers: `toggleTheme`, `showTab`, `makeSampleDropdown`, `rf`, `rfFiltered`, `sourceBadge`, `qualBadge`, `makeTable`, `filterTable` |
| 2360–2410 | `rf(...)` calls for all pre-built Plotly figures; `makeSampleDropdown` wiring |
| 2389–2470 | **`renderVC3Network(sample)`** — vConTACT3 network graph (spectral scatter+edges, coloured by cluster) |
| 2473–2560 | `makeSampleDropdown('sample-ctrl-taxonomy', ...)` — calls `renderVC3Network` + tax functions |
| 2560–2700 | Viral taxonomy JS: `renderTaxSourcePie`, `renderTaxFamilyBar`, `renderTaxSunburst`; prok: `renderProkDomainBar`, `renderProkTopPhyla`, `renderProkMasterTable`; GTDB sunburst |
| 2700–2800 | Host prediction JS: `renderPhistSection`; virus-host network |
| 2800–2900 | Overview grid render: `OV_GROUPS`, `scCard`, run summary, per-sample metric cards |
| ~2960 | `</script></div></div></body></html>"""` — f-string ends here |
| ~2960–2974 | `os.makedirs`, `open(out_html, 'w', encoding='utf-8')`, error handler with traceback log |

### Critical coding rules for this file

- **CSS is in `_CSS` (~lines 1624–1755): plain string, no f-string.** `{` and `}` are literal — no doubling needed.
- **Everything from ~line 1766 to ~2960 is one f-string.** All `{` and `}` in JS/HTML must be doubled: `{{` `}}`. Exception: `{python_var}` for actual Python interpolation.
- **`\t` and `\n` inside JS strings within the f-string must be `\\t` and `\\n`** — otherwise Python interpolates them as real tab/newline characters, breaking the JS (SyntaxError in browser).
- **Plotly 6.x:** `go.Sunburst` does not accept `showlegend` parameter — remove it if present.
- **HTML write uses `encoding='utf-8'`** (line ~2970) — required on Windows; do not remove.
- **`VC3_NETWORK`** JS constant: dict keyed by sample name, each value `{nodes:[...], edges:[...]}` loaded from `network_layout.json` per sample.

