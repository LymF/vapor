# VAPOR Pipeline — Snakemake to Nextflow DSL2 Migration Plan

## Executive Summary

This document outlines the strategy for migrating the VAPOR (Viral And Prokaryotic Omics Resource) metagenomics pipeline from Snakemake to Nextflow DSL2. The pipeline currently spans 17 rule files, ~60+ rules, 23 conda environments, and processes both short-read (Illumina PE/SE) and long-read (ONT/HiFi) data through QC, assembly, viral detection, binning, taxonomy, and reporting stages.

**Key goals:**
- Preserve all analytical functionality and results reproducibility
- Leverage Nextflow's native cloud execution (AWS Batch, Google Cloud, Azure)
- Enable seamless deployment via Seqera Platform
- Maintain the modular architecture with nf-core conventions
- Support Wave containers for reproducible, lock-file-based environments

---

## Table of Contents

1. [Architecture Overview](#1-architecture-overview)
2. [Module Mapping: Snakemake Rules → Nextflow Processes](#2-module-mapping)
3. [Subworkflow Design](#3-subworkflow-design)
4. [Channel Architecture](#4-channel-architecture)
5. [Configuration Strategy](#5-configuration-strategy)
6. [Container Strategy](#6-container-strategy)
7. [Migration Phases](#7-migration-phases)
8. [Testing Strategy](#8-testing-strategy)
9. [Platform Integration](#9-platform-integration)
10. [Risk Assessment](#10-risk-assessment)

---

## 1. Architecture Overview

### Current Snakemake Structure

```
Snakefile (main)
├── config.yaml          (all parameters)
├── containers.lock.yaml (pinned container URIs)
├── rules/
│   ├── qc.smk              (BLOCK 1)   — fastp, NanoPlot, Porechop, Filtlong
│   ├── host_removal.smk    (BLOCK 1.5) — bwa-mem2/minimap2 host decontamination
│   ├── assembly.smk        (BLOCK 2)   — MEGAHIT, metaSPAdes, metaviralSPAdes, Flye, hifiasm, metaMDBG, Medaka
│   ├── merge_dedup.smk     (BLOCK 3)   — contig merging + MMseqs2 deduplication
│   ├── quast.smk           (BLOCK 4)   — assembly QC
│   ├── viral_detection.smk (BLOCK 5)   — VirSorter2, GeNomad, VIBRANT, consensus
│   ├── cobra.smk           (BLOCK 5.5) — COBRA contig extension (optional, SR PE only)
│   ├── mapping.smk         (BLOCK 6)   — BWA-MEM2/minimap2, depth calculation
│   ├── viral_binning.smk   (BLOCK 7)   — CheckV, vRhyme, vOTU clustering (skani)
│   ├── prok_binning.smk    (BLOCK 8)   — MetaBAT2, VAMB, SemiBin2, COMEBin, Binette, CheckM2, GTDB-Tk, GUNC, galah
│   ├── taxonomy.smk        (BLOCK 9)   — Prodigal, MMseqs2 taxonomy, vConTACT3, viral taxonomy
│   ├── host_prediction.smk (BLOCK 10)  — PHIST
│   ├── abundance.smk       (BLOCK 10.2)— CoverM, diversity metrics
│   ├── annotation.smk      (BLOCK 10.3)— Pharokka, phold, Bakta, eggNOG, KEGG, genome maps
│   ├── defense_amr.smk     (BLOCK 10.5)— DefenseFinder, AMRFinderPlus, RGI, DeepARG, ABRicate, dbAPIS
│   ├── finalize.smk        (BLOCK 11)  — output organization
│   └── report.smk          (BLOCK 12)  — HTML report generation, MultiQC
└── scripts/
    ├── generate_report.py
    ├── merge_lr_assemblies.py
    ├── filter_checkv_hq.py
    ├── split_viral_fastas.py
    └── prepare_diamond_db.py
```

### Target Nextflow Structure

```
main.nf
├── nextflow.config
├── nextflow_schema.json
├── modules/
│   ├── local/
│   │   ├── fastp/main.nf
│   │   ├── megahit/main.nf
│   │   ├── metaspades/main.nf
│   │   ├── flye/main.nf
│   │   ├── hifiasm_meta/main.nf
│   │   ├── metamdbg/main.nf
│   │   ├── mmseqs2_dedup/main.nf
│   │   ├── virsorter2/main.nf
│   │   ├── genomad/main.nf
│   │   ├── vibrant/main.nf
│   │   ├── viral_consensus/main.nf
│   │   ├── checkv/main.nf
│   │   ├── vrhyme/main.nf
│   │   ├── vcontact3/main.nf
│   │   ├── metabat2/main.nf
│   │   ├── vamb/main.nf
│   │   ├── semibin2/main.nf
│   │   ├── comebin/main.nf
│   │   ├── binette/main.nf
│   │   ├── checkm2/main.nf
│   │   ├── gtdbtk/main.nf
│   │   ├── gunc/main.nf
│   │   ├── galah/main.nf
│   │   ├── phist/main.nf
│   │   ├── cobra/main.nf
│   │   ├── pharokka/main.nf
│   │   ├── phold/main.nf
│   │   ├── bakta/main.nf
│   │   ├── eggnog_mapper/main.nf
│   │   ├── kegg_decoder/main.nf
│   │   ├── coverm/main.nf
│   │   ├── defensefinder/main.nf
│   │   ├── amrfinderplus/main.nf
│   │   ├── rgi/main.nf
│   │   ├── deeparg/main.nf
│   │   ├── abricate/main.nf
│   │   ├── argnorm/main.nf
│   │   ├── dbapis/main.nf
│   │   ├── skani/main.nf
│   │   └── generate_report/main.nf
│   └── nf-core/       (reusable nf-core modules)
│       ├── fastqc/
│       ├── nanoplot/
│       ├── porechop_abi/
│       ├── filtlong/
│       ├── bwa_mem2/
│       ├── minimap2/
│       ├── samtools/
│       ├── quast/
│       ├── prodigal/
│       └── multiqc/
├── subworkflows/
│   └── local/
│       ├── qc_short_reads.nf
│       ├── qc_long_reads.nf
│       ├── host_removal.nf
│       ├── assembly_short_reads.nf
│       ├── assembly_long_reads.nf
│       ├── merge_dedup.nf
│       ├── viral_detection.nf
│       ├── viral_binning.nf
│       ├── prok_binning.nf
│       ├── taxonomy.nf
│       ├── annotation_viral.nf
│       ├── annotation_prok.nf
│       ├── defense_amr.nf
│       ├── abundance.nf
│       └── reporting.nf
├── conf/
│   ├── base.config
│   ├── modules.config
│   ├── test.config
│   └── profiles/
│       ├── docker.config
│       ├── singularity.config
│       ├── conda.config
│       └── wave.config
├── assets/
│   └── samplesheet_schema.json
├── bin/
│   ├── generate_report.py
│   ├── merge_lr_assemblies.py
│   ├── filter_checkv_hq.py
│   └── split_viral_fastas.py
└── tests/
    ├── main.nf.test
    └── modules/
```


---

## 2. Module Mapping: Snakemake Rules → Nextflow Processes

### BLOCK 1 — Quality Control

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `fastqc_raw` | `FASTQC` | ✅ `nf-core/fastqc` | Standard QC |
| `fastp_trim` | `FASTP` | ✅ `nf-core/fastp` | Adapter trimming + filtering |
| `nanoplot_raw` | `NANOPLOT` | ✅ `nf-core/nanoplot` | Long-read QC |
| `porechop_abi` | `PORECHOP_ABI` | ✅ `nf-core/porechop_abi` | Adapter removal |
| `filtlong` | `FILTLONG` | ✅ `nf-core/filtlong` | Length/quality filtering |

### BLOCK 1.5 — Host Removal (Optional)

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `bwa_index_host` | `BWA_MEM2_INDEX` | ✅ `nf-core/bwamem2/index` | Index host genome |
| `host_remove_sr` | `BWA_MEM2_MEM` + `SAMTOOLS_FASTQ` | ✅ | Map + extract unmapped |
| `host_remove_lr` | `MINIMAP2_ALIGN` + `SAMTOOLS_FASTQ` | ✅ | Long-read host removal |

### BLOCK 2 — Assembly

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `megahit` | `MEGAHIT` | ✅ `nf-core/megahit` | Primary SR assembler |
| `metaspades` | `METASPADES` | Local module | metaSPAdes + metaviralSPAdes |
| `flye` | `FLYE` | Local module | LR metagenome assembler |
| `hifiasm_meta` | `HIFIASM_META` | Local module | HiFi assembler |
| `metamdbg` | `METAMDBG` | Local module | HiFi/ONT R10 assembler |
| `medaka` | `MEDAKA` | Local module | ONT polishing |

### BLOCK 3 — Merge & Deduplicate

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `merge_contigs` | `MERGE_CONTIGS` | Local (script) | Concatenate + length filter |
| `mmseqs2_dedup` | `MMSEQS2_CLUSTER` | Local module | 95% ANI deduplication → rep_seq.fasta |

### BLOCK 4 — Assembly QC

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `quast` | `QUAST` | ✅ `nf-core/quast` | Assembly statistics |

### BLOCK 5 — Viral Detection

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `virsorter2` | `VIRSORTER2` | Local module | Viral prediction |
| `genomad` | `GENOMAD` | Local module | geNomad end-to-end |
| `vibrant` | `VIBRANT` | Local module | VIBRANT viral detection |
| `viral_consensus` | `VIRAL_CONSENSUS` | Local (script) | Multi-tool consensus logic |

### BLOCK 5.5 — COBRA (Optional)

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `cobra_extend` | `COBRA` | Local module | Contig extension (SR PE only) |

### BLOCK 6 — Read Mapping

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `bwa_index` | `BWA_MEM2_INDEX` | ✅ | Index rep_seq.fasta |
| `bwa_mem` | `BWA_MEM2_MEM` | ✅ | Map reads to contigs |
| `minimap2_lr` | `MINIMAP2_ALIGN` | ✅ | Long-read mapping |
| `samtools_sort` | `SAMTOOLS_SORT` | ✅ | Sort + index BAM |
| `calc_depth` | `CALC_DEPTH` | Local (script) | MetaBAT2 jgi_summarize_bam |

### BLOCK 7 — Viral Binning

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `checkv_initial` | `CHECKV` | Local module | Quality assessment |
| `vrhyme` | `VRHYME` | Local module | Viral binning |
| `checkv_vrhyme` | `CHECKV` (2nd pass) | Reuse module | Post-binning QC |
| `votu_cluster` | `SKANI_CLUSTER` | Local module | vOTU at 95% ANI / 85% AF |
| `votu_table` | `VOTU_TABLE` | Local (script) | vOTU membership + abundance |

### BLOCK 8 — Prokaryotic Binning

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `filter_viral_contigs` | `FILTER_VIRAL` | Local (script) | Remove viral from binning input |
| `metabat2` | `METABAT2` | ✅ `nf-core/metabat2` | Depth-based binning |
| `vamb` | `VAMB` | Local module | VAE-based binning |
| `semibin2` | `SEMIBIN2` | Local module | Semi-supervised binning |
| `comebin` | `COMEBIN` | Local module | Contrastive learning binning |
| `binette` | `BINETTE` | Local module | Bin refinement / consolidation |
| `checkm2` | `CHECKM2` | Local module | MAG quality assessment |
| `gtdbtk` | `GTDBTK_CLASSIFYWF` | ✅ (partial) | Genome taxonomy |
| `gunc` | `GUNC` | Local module | Chimera detection |
| `galah_derep` | `GALAH` | Local module | MAG dereplication |

### BLOCK 9 — Taxonomy

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `prodigal_orfs` | `PRODIGAL` | ✅ `nf-core/prodigal` | ORF prediction |
| `mmseqs_taxonomy_prok` | `MMSEQS2_TAXONOMY` | Local module | Prokaryotic LCA taxonomy |
| `mmseqs_taxonomy_viral` | `MMSEQS2_TAXONOMY` | Reuse module | Custom viral taxonomy |
| `genomad_taxonomy` | Part of `GENOMAD` | — | Viral taxonomy from geNomad |
| `vcontact3` | `VCONTACT3` | Local module | Viral clustering + taxonomy |
| `merge_viral_taxonomy` | `MERGE_VIRAL_TAX` | Local (script) | Consolidate viral taxonomy |

### BLOCK 10 — Host Prediction

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `phist` | `PHIST` | Local module | Phage-host prediction |

### BLOCK 10.2 — Abundance & Diversity

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `coverm_viral` | `COVERM` | Local module | Viral genome coverage |
| `coverm_prok` | `COVERM` | Reuse module | MAG coverage |
| `diversity_metrics` | `DIVERSITY` | Local (script) | Alpha/beta diversity |

### BLOCK 10.3 — Annotation

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `pharokka` | `PHAROKKA` | Local module | Phage annotation |
| `phold` | `PHOLD` | Local module | Phage structural annotation |
| `bakta` | `BAKTA` | ✅ `nf-core/bakta` | Prokaryotic annotation |
| `eggnog_mapper` | `EGGNOG_MAPPER` | Local module | Functional annotation |
| `kegg_decoder` | `KEGG_DECODER` | Local (script) | Metabolic pathway decoding |
| `genome_maps` | `GENOME_MAPS` | Local (script) | Circular genome visualizations |

### BLOCK 10.5 — Defense & AMR

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `defensefinder_prok` | `DEFENSEFINDER` | Local module | Defense systems in MAGs |
| `amrfinderplus` | `AMRFINDERPLUS` | ✅ `nf-core/amrfinderplus` | AMR gene detection |
| `rgi` | `RGI` | Local module | CARD-based AMR |
| `deeparg` | `DEEPARG` | Local module | Deep learning AMR |
| `abricate` | `ABRICATE` | ✅ `nf-core/abricate` | Multi-DB screening |
| `argnorm` | `ARGNORM` | Local module | ARG normalization |
| `defensefinder_viral` | `DEFENSEFINDER` | Reuse module | Defense in viral genomes |
| `dbapis` | `DBAPIS` | Local module | Anti-defense systems |
| `defense_islands` | `DEFENSE_ISLANDS` | Local (script) | Co-localization analysis |

### BLOCK 12 — Reporting

| Snakemake Rule | Nextflow Process | nf-core Module? | Notes |
|---|---|---|---|
| `generate_report` | `GENERATE_REPORT` | Local (script) | Interactive HTML report |
| `multiqc` | `MULTIQC` | ✅ `nf-core/multiqc` | Aggregated QC report |


---

## 3. Subworkflow Design

The pipeline will be organized into **15 subworkflows** that mirror the logical blocks. Each subworkflow is self-contained and communicates via typed channels.

### 3.1 Main Workflow (main.nf)

```groovy
workflow {
    // Parse samplesheet
    ch_input = channel.fromPath(params.input)
        .splitCsv(header: true)
        .map { row -> create_meta_channel(row) }

    // Route by sequencing type
    if (params.long_reads) {
        QC_LONG_READS(ch_input)
        ch_clean = QC_LONG_READS.out.reads
    } else {
        QC_SHORT_READS(ch_input)
        ch_clean = QC_SHORT_READS.out.reads
    }

    // Optional host removal
    if (params.host_genome) {
        HOST_REMOVAL(ch_clean, params.host_genome)
        ch_clean = HOST_REMOVAL.out.reads
    }

    // Assembly (platform-aware branching)
    if (params.long_reads) {
        ASSEMBLY_LONG_READS(ch_clean)
        ch_contigs = ASSEMBLY_LONG_READS.out.contigs
    } else {
        ASSEMBLY_SHORT_READS(ch_clean)
        ch_contigs = ASSEMBLY_SHORT_READS.out.contigs
    }

    // Central hub: deduplicated representative sequences
    MERGE_DEDUP(ch_contigs)
    ch_rep_seq = MERGE_DEDUP.out.rep_seq

    // Mapping (needed for binning depth)
    MAPPING(ch_clean, ch_rep_seq)

    // Viral branch
    VIRAL_DETECTION(ch_rep_seq)
    VIRAL_BINNING(VIRAL_DETECTION.out.viral_seqs, MAPPING.out.bam, MAPPING.out.depth)

    // Prokaryotic branch
    PROK_BINNING(ch_rep_seq, MAPPING.out.depth, VIRAL_DETECTION.out.viral_ids)

    // Taxonomy & host prediction
    TAXONOMY(VIRAL_BINNING.out.viral_bins, PROK_BINNING.out.mags, ch_rep_seq)

    // Annotation
    ANNOTATION_VIRAL(VIRAL_BINNING.out.viral_bins)
    ANNOTATION_PROK(PROK_BINNING.out.mags)

    // Abundance & diversity
    ABUNDANCE(MAPPING.out.bam, VIRAL_BINNING.out.viral_bins, PROK_BINNING.out.mags)

    // Defense & AMR
    DEFENSE_AMR(PROK_BINNING.out.mags, VIRAL_BINNING.out.viral_orfs)

    // Reporting
    REPORTING(/* collect all QC + results channels */)
}
```

### 3.2 Subworkflow Catalog

| Subworkflow | Processes | Key Inputs | Key Outputs |
|---|---|---|---|
| `QC_SHORT_READS` | FASTQC, FASTP | raw reads | trimmed reads, QC reports |
| `QC_LONG_READS` | NANOPLOT, PORECHOP_ABI, FILTLONG | raw LR | filtered reads, QC reports |
| `HOST_REMOVAL` | BWA_MEM2/MINIMAP2, SAMTOOLS | clean reads, host ref | decontaminated reads |
| `ASSEMBLY_SHORT_READS` | MEGAHIT, METASPADES | trimmed reads | raw contigs per assembler |
| `ASSEMBLY_LONG_READS` | FLYE, HIFIASM_META, METAMDBG, MEDAKA | filtered LR | raw contigs per assembler |
| `MERGE_DEDUP` | MERGE_CONTIGS, MMSEQS2_CLUSTER, QUAST | multiple contig sets | rep_seq.fasta |
| `VIRAL_DETECTION` | VIRSORTER2, GENOMAD, VIBRANT, VIRAL_CONSENSUS | rep_seq | viral contigs + tool support TSV |
| `VIRAL_BINNING` | CHECKV, VRHYME, SKANI_CLUSTER, VOTU_TABLE | viral seqs, BAM, depth | viral bins, vOTU table |
| `PROK_BINNING` | METABAT2, VAMB, SEMIBIN2, COMEBIN, BINETTE, CHECKM2, GTDBTK, GUNC, GALAH | rep_seq, depth, viral IDs | quality-filtered MAGs |
| `TAXONOMY` | PRODIGAL, MMSEQS2_TAX, VCONTACT3, MERGE_TAX | viral bins, MAGs | taxonomy tables |
| `ANNOTATION_VIRAL` | PHAROKKA, PHOLD, GENOME_MAPS | viral bins | annotated viral genomes |
| `ANNOTATION_PROK` | BAKTA, EGGNOG_MAPPER, KEGG_DECODER, GENOME_MAPS | MAGs | annotated MAGs |
| `ABUNDANCE` | COVERM, DIVERSITY | BAM, bins | abundance tables, diversity metrics |
| `DEFENSE_AMR` | DEFENSEFINDER, AMRFINDERPLUS, RGI, DEEPARG, ABRICATE, ARGNORM, DBAPIS | MAG proteins, viral ORFs | defense/AMR tables |
| `REPORTING` | GENERATE_REPORT, MULTIQC | all results | HTML reports |


---

## 4. Channel Architecture

### 4.1 Meta Map Convention

All channels carry a `[meta, ...]` tuple where `meta` is a Groovy map:

```groovy
// meta structure
[
    id:          "sample_name",
    single_end:  false,        // or true for SE/LR
    long_reads:  false,        // true for ONT/HiFi
    lr_tech:     "ont"         // "ont", "hifi", or null
]
```

### 4.2 Key Channel Shapes

```
ch_reads:      [ meta, [ reads ] ]           // PE: [R1, R2], SE/LR: [reads]
ch_contigs:    [ meta, contigs.fasta ]        // per-assembler
ch_rep_seq:    [ meta, rep_seq.fasta ]        // central hub — deduplicated
ch_bam:        [ meta, sorted.bam, sorted.bam.bai ]
ch_depth:      [ meta, depth.txt ]
ch_viral_seqs: [ meta, viral_consensus.fasta, tool_support.tsv ]
ch_viral_bins: [ meta, bins_dir/ ]
ch_mags:       [ meta, mags_dir/, quality_report.tsv ]
ch_taxonomy:   [ meta, taxonomy.tsv ]
```

### 4.3 Branching Strategy

The pipeline uses conditional logic at the workflow level (not process level) to route between SR/LR paths:

```groovy
// In main.nf — branch by sequencing type
ch_input.branch { meta, reads ->
    long_reads: meta.long_reads
    short_reads: !meta.long_reads
}
.set { ch_branched }

// Each branch feeds its own subworkflow
QC_SHORT_READS(ch_branched.short_reads)
QC_LONG_READS(ch_branched.long_reads)

// Merge back after assembly
ch_contigs = ASSEMBLY_SHORT_READS.out.contigs
    .mix(ASSEMBLY_LONG_READS.out.contigs)
```

### 4.4 Hub-and-Spoke Pattern

The `rep_seq.fasta` acts as the central hub. Multiple downstream subworkflows consume it independently:

```
                    ┌─► VIRAL_DETECTION
                    │
rep_seq.fasta ──────┼─► MAPPING ─────┬─► VIRAL_BINNING
                    │                 └─► PROK_BINNING
                    └─► TAXONOMY
```

This maps naturally to Nextflow's dataflow model — a single channel can be consumed by multiple processes without explicit fan-out.

### 4.5 Optional Steps as Channel Gates

Optional pipeline steps (host removal, COBRA, GUNC, galah) use Nextflow's channel emptiness pattern:

```groovy
// If host_genome is provided, run host removal; otherwise pass reads through
if (params.host_genome) {
    HOST_REMOVAL(ch_reads, file(params.host_genome))
    ch_clean = HOST_REMOVAL.out.reads
} else {
    ch_clean = ch_reads
}
```


---

## 5. Configuration Strategy

### 5.1 Parameter Hierarchy

```
nextflow.config          (base defaults, process resources, profiles)
├── conf/base.config     (CPU/memory/time defaults per process label)
├── conf/modules.config  (publishDir, ext.args per process)
├── conf/test.config     (minimal test dataset)
└── params via CLI or Seqera Platform Launchpad
```

### 5.2 Parameter Surface

Parameters map directly from current `config.yaml`:

```groovy
params {
    // Input/Output
    input               = null          // samplesheet CSV
    outdir              = './results'

    // Sequencing mode
    long_reads          = false
    single_end          = false
    lr_tech             = 'ont'         // 'ont' or 'hifi'

    // Assembly
    use_spades          = true
    spades_mem          = 250
    spades_kmers        = 'auto'
    megahit_preset      = 'meta-large'
    min_contig          = 1000

    // Deduplication
    min_seq_id          = 0.95

    // Viral detection
    viral_consensus_mode = 'hybrid'
    min_viral_tools      = 2
    score_vs2_min        = 0.5
    score_genomad_min    = 0.7

    // Host removal (optional)
    host_genome          = null
    host_index           = null

    // Databases (required)
    checkv_db           = null
    vs2_db              = null
    genomad_db          = null
    vibrant_base        = null
    checkm2_db          = null
    gtdbtk_db           = null
    inphared_db         = null
    vcontact3_db        = null

    // Databases (optional)
    custom_prok_mmseqs_db  = null
    custom_viral_mmseqs_db = null
    pharokka_db            = null
    phold_db               = null
    bakta_db               = null
    eggnog_db              = null
    card_db                = null
    deeparg_db             = null
    apis_db                = null
    gunc_db                = null

    // Feature toggles
    cobra_enabled       = false
    comebin_enabled     = true
    gunc_enabled        = true
    mag_derep_enabled   = true
    votu_clustering_enabled = true
    defense_amr_enabled = true
    use_gpu             = false

    // Long-read specific
    lr_min_len          = 1000
    lr_min_mean_q       = 10
    lr_flye_overlap     = 3000
    lr_medaka_model     = 'r1041_e82_400bps_sup_v5.0.0'
    lr_metaMDBG         = true
}
```

### 5.3 Resource Configuration (conf/base.config)

```groovy
process {
    // Default resources
    cpus   = { 4 * task.attempt }
    memory = { 8.GB * task.attempt }
    time   = { 4.h * task.attempt }

    // Labels for resource tiers
    withLabel: process_low {
        cpus   = { 2 * task.attempt }
        memory = { 4.GB * task.attempt }
        time   = { 2.h * task.attempt }
    }
    withLabel: process_medium {
        cpus   = { 8 * task.attempt }
        memory = { 32.GB * task.attempt }
        time   = { 8.h * task.attempt }
    }
    withLabel: process_high {
        cpus   = { 16 * task.attempt }
        memory = { 64.GB * task.attempt }
        time   = { 16.h * task.attempt }
    }
    withLabel: process_high_memory {
        cpus   = { 16 * task.attempt }
        memory = { 200.GB * task.attempt }
        time   = { 24.h * task.attempt }
    }

    // Tool-specific overrides
    withName: 'METASPADES' {
        cpus   = { 16 * task.attempt }
        memory = { params.spades_mem.GB * task.attempt }
        time   = { 72.h * task.attempt }
    }
    withName: 'MEGAHIT' {
        cpus   = { 16 * task.attempt }
        memory = { 64.GB * task.attempt }
    }
    withName: 'GTDBTK_CLASSIFYWF' {
        cpus   = { 32 }
        memory = { 220.GB * task.attempt }
        time   = { 48.h * task.attempt }
    }
    withName: 'VAMB' {
        cpus   = { 16 * task.attempt }
        memory = { 64.GB * task.attempt }
        accelerator = params.use_gpu ? 1 : 0
    }
}
```

### 5.4 Samplesheet Format

Replace the auto-discovery pattern with explicit samplesheet CSV:

```csv
sample,fastq_1,fastq_2,long_reads
sample_A,/data/sample_A_R1.fq.gz,/data/sample_A_R2.fq.gz,false
sample_B,/data/sample_B_R1.fq.gz,/data/sample_B_R2.fq.gz,false
sample_C,/data/sample_C.fastq.gz,,true
```

Validated via `nextflow_schema.json` with nf-schema plugin.


---

## 6. Container Strategy

### 6.1 Approach

Use **Wave containers** built from conda lock files for maximum reproducibility:

- Each process gets a dedicated container built from a `conda.yml` spec
- Wave builds containers on-demand with exact package versions
- The existing `containers.lock.yaml` maps directly to Nextflow `container` directives
- Fall back to Docker Hub / Biocontainers for standard tools with existing images

### 6.2 Container Mapping

```groovy
// In conf/modules.config or per-process
process {
    withName: 'VIRSORTER2' {
        conda    = "${projectDir}/envs/env_viral.yml"
        container = 'community.wave.seqera.io/library/virsorter2:2.2.4--abc123'
    }
    withName: 'GENOMAD' {
        conda    = "${projectDir}/envs/env_viral.yml"
        container = 'community.wave.seqera.io/library/genomad:1.8.0--def456'
    }
    withName: 'CHECKM2' {
        conda    = "${projectDir}/envs/env_checkm2.yml"
        container = 'community.wave.seqera.io/library/checkm2:1.0.2--ghi789'
    }
}
```

### 6.3 Environment Consolidation

Current pipeline has 23 conda environments. For Nextflow, consolidate where tools share dependencies:

| Current Env | Tools | Nextflow Strategy |
|---|---|---|
| `env_qc` | fastp, FastQC | Separate nf-core containers (already available) |
| `env_assembly` | MEGAHIT, SPAdes, MMseqs2 | Split: MEGAHIT separate, SPAdes separate, MMseqs2 separate |
| `env_viral` | VirSorter2, GeNomad, VIBRANT | Keep separate — conflicting deps |
| `env_mapping` | BWA-MEM2, samtools, minimap2 | Use nf-core biocontainers |
| `env_binning` | MetaBAT2, VAMB, SemiBin2 | Split per tool — VAMB needs GPU |
| `env_checkm2` | CheckM2 | Dedicated container |
| `env_gtdbtk` | GTDB-Tk | Dedicated container (large) |
| `env_taxonomy` | Prodigal, Diamond, MMseqs2 | Split per tool |
| `env_annotation` | Pharokka, Bakta, eggNOG | Split per tool |
| `env_defense` | DefenseFinder, AMRFinderPlus, RGI | Split per tool |
| `env_report` | Python (pandas, jinja2, echarts) | Dedicated lightweight container |

**Rule:** One tool = one container in Nextflow. This avoids dependency conflicts and enables independent version bumps.

### 6.4 Wave Integration

```groovy
// nextflow.config
wave {
    enabled = true
    strategy = ['conda']
}

docker {
    enabled = true
}
```

This allows Nextflow to automatically build containers from conda specs when no pre-built image is specified.


---

## 7. Migration Phases

### Phase 1: Foundation (Weeks 1–2)

**Goal:** Scaffold project, migrate QC + Assembly core

- [ ] Initialize Nextflow project structure (`main.nf`, `nextflow.config`, module dirs)
- [ ] Create `nextflow_schema.json` with all parameters
- [ ] Create samplesheet schema + validation
- [ ] Implement `QC_SHORT_READS` subworkflow (FASTQC + FASTP)
- [ ] Implement `QC_LONG_READS` subworkflow (NANOPLOT + PORECHOP_ABI + FILTLONG)
- [ ] Implement `HOST_REMOVAL` subworkflow (optional)
- [ ] Implement `ASSEMBLY_SHORT_READS` subworkflow (MEGAHIT + METASPADES)
- [ ] Implement `ASSEMBLY_LONG_READS` subworkflow (FLYE + HIFIASM + METAMDBG + MEDAKA)
- [ ] Implement `MERGE_DEDUP` subworkflow (merge + MMseqs2 + QUAST)
- [ ] Create test profile with minimal dataset
- [ ] Validate: QC → Assembly → rep_seq.fasta output matches Snakemake

**Deliverable:** Pipeline runs QC through deduplication, producing identical rep_seq.fasta

### Phase 2: Viral Analysis (Weeks 3–4)

**Goal:** Complete viral detection, binning, and taxonomy branch

- [ ] Implement `VIRAL_DETECTION` subworkflow (VS2 + GeNomad + VIBRANT + consensus)
- [ ] Implement `MAPPING` processes (BWA-MEM2/minimap2 + depth)
- [ ] Implement `VIRAL_BINNING` subworkflow (CheckV + vRhyme + skani vOTU)
- [ ] Implement COBRA optional contig extension
- [ ] Implement viral taxonomy (geNomad taxonomy + vConTACT3 + merge)
- [ ] Implement `PHIST` host prediction
- [ ] Validate: viral consensus FASTA, CheckV quality, vOTU tables match

**Deliverable:** Full viral branch operational

### Phase 3: Prokaryotic Analysis (Weeks 5–6)

**Goal:** Complete prokaryotic binning and taxonomy branch

- [ ] Implement `PROK_BINNING` subworkflow
  - [ ] Viral contig filtering
  - [ ] MetaBAT2, VAMB, SemiBin2, COMEBin
  - [ ] Binette consolidation
  - [ ] CheckM2 quality assessment
  - [ ] GUNC chimera detection
  - [ ] galah MAG dereplication
  - [ ] GTDB-Tk classification
- [ ] Implement prokaryotic MMseqs2 taxonomy
- [ ] Validate: MAG quality reports, GTDB taxonomy match

**Deliverable:** Full prokaryotic branch operational

### Phase 4: Annotation & Defense (Weeks 7–8)

**Goal:** Complete functional annotation and defense/AMR analysis

- [ ] Implement `ANNOTATION_VIRAL` (Pharokka + phold + genome maps)
- [ ] Implement `ANNOTATION_PROK` (Bakta + eggNOG + KEGG decoder + genome maps)
- [ ] Implement `ABUNDANCE` (CoverM + diversity metrics)
- [ ] Implement `DEFENSE_AMR` subworkflow
  - [ ] Prokaryotic: DefenseFinder, AMRFinderPlus, RGI, DeepARG, ABRicate, argNorm
  - [ ] Viral: DefenseFinder, dbAPIS, defense islands
- [ ] Validate: annotation outputs, AMR tables match

**Deliverable:** All analysis branches complete

### Phase 5: Reporting & Platform (Weeks 9–10)

**Goal:** Reporting, Seqera Platform integration, final validation

- [ ] Implement `REPORTING` subworkflow (HTML report + MultiQC)
- [ ] Port `generate_report.py` to `bin/` directory
- [ ] Configure Seqera Platform pipeline
- [ ] Set up compute environments (AWS Batch / local)
- [ ] Full end-to-end validation against Snakemake reference outputs
- [ ] Performance benchmarking (runtime, memory, cost)
- [ ] Documentation (README, usage docs)
- [ ] Release v1.0.0

**Deliverable:** Production-ready pipeline on Seqera Platform


---

## 8. Testing Strategy

### 8.1 Test Levels

| Level | What | How | When |
|---|---|---|---|
| Unit | Individual processes | nf-test per module | Every module created |
| Integration | Subworkflow chains | nf-test with small dataset | Per phase completion |
| End-to-end | Full pipeline | Reference dataset comparison | Phase 5 |
| Regression | Output reproducibility | Snapshot tests vs. Snakemake | Final validation |

### 8.2 Test Dataset

Create a minimal test dataset (subset of real data):
- 2 samples × 10,000 read pairs (SR) or 5,000 reads (LR)
- Pre-computed database subsets (CheckV, GeNomad mini DBs)
- Expected outputs stored in `tests/expected/`

### 8.3 nf-test Example

```groovy
// tests/modules/local/virsorter2/main.nf.test
nextflow_process {
    name "Test VIRSORTER2"
    script "../../../modules/local/virsorter2/main.nf"
    process "VIRSORTER2"

    test("Should detect viral contigs") {
        when {
            params {
                vs2_db = "${projectDir}/tests/data/vs2_db_mini"
            }
            process {
                """
                input[0] = [
                    [id: 'test_sample'],
                    file("${projectDir}/tests/data/test_contigs.fasta")
                ]
                input[1] = file(params.vs2_db)
                """
            }
        }
        then {
            assert process.success
            assert path(process.out.viral_seqs[0][1]).exists()
        }
    }
}
```

### 8.4 Validation Criteria

For each migrated block, verify:
1. **Identical outputs** — same contigs, same bins, same taxonomy assignments
2. **Equivalent performance** — runtime within 20% of Snakemake
3. **Resource usage** — no memory/CPU regressions
4. **Error handling** — graceful failures with informative messages
5. **Resume capability** — `-resume` skips completed tasks correctly


---

## 9. Platform Integration

### 9.1 Seqera Platform Configuration

The migrated pipeline will be deployed on Seqera Platform with:

- **Launchpad** — one-click execution with pre-filled parameters
- **Parameter schema** — `nextflow_schema.json` renders a form UI
- **Compute environments** — AWS Batch (primary), local (dev/test)
- **Data Studios** — post-run exploration of results
- **Reports** — custom HTML report + MultiQC via `reports` config block

### 9.2 Reports Configuration

```groovy
// nextflow.config
reports {
    "results/report.html" {
        display = "VAPOR Interactive Report"
    }
    "results/multiqc_report.html" {
        display = "MultiQC Summary"
    }
}
```

### 9.3 Database Management

Large databases (CheckV, GTDB-Tk, GeNomad etc.) will be:
1. **Pre-staged** on shared storage (S3/EFS) accessible to compute environments
2. **Referenced** via absolute paths in Platform-level params
3. **Versioned** in a `databases.config` profile for reproducibility

```groovy
// conf/databases_aws.config
params {
    checkv_db    = 's3://vapor-databases/checkv/v1.5'
    vs2_db       = 's3://vapor-databases/virsorter2/v2'
    genomad_db   = 's3://vapor-databases/genomad/v1.7'
    gtdbtk_db    = 's3://vapor-databases/gtdbtk/r220'
    checkm2_db   = 's3://vapor-databases/checkm2/v2'
}
```

### 9.4 Compute Environment Requirements

| Process Group | vCPUs | Memory | Storage | GPU | Time |
|---|---|---|---|---|---|
| QC / Mapping | 4–8 | 8–16 GB | 50 GB | ❌ | 1–2h |
| Assembly (MEGAHIT) | 16 | 64 GB | 100 GB | ❌ | 4–8h |
| Assembly (SPAdes) | 16 | 250 GB | 200 GB | ❌ | 12–72h |
| Viral Detection | 8–16 | 32–64 GB | 50 GB | ❌ | 2–8h |
| Binning (MetaBAT2/SemiBin2) | 8–16 | 32 GB | 50 GB | ❌ | 2–4h |
| Binning (VAMB/COMEBin) | 16 | 64 GB | 50 GB | ✅ (optional) | 4–8h |
| GTDB-Tk | 32 | 220 GB | 100 GB | ❌ | 12–48h |
| CheckM2 | 8 | 16 GB | 20 GB | ❌ | 1–2h |
| Annotation (Bakta) | 8 | 32 GB | 50 GB | ❌ | 2–8h |
| Report generation | 2 | 4 GB | 10 GB | ❌ | <30min |

### 9.5 AWS Batch Configuration

Recommended Forge settings:
- **maxCpus:** 500 (allows ~30 samples in parallel)
- **Instance types:** m5, r5 (high-memory for GTDB-Tk/SPAdes), g4dn (GPU for VAMB)
- **EBS autoscale:** enabled, 100 GB initial + 500 GB max
- **Spot instances:** enabled for fault-tolerant steps (QC, mapping)
- **Fusion filesystem:** enabled for S3 direct access


---

## 10. Risk Assessment

### 10.1 Technical Risks

| Risk | Impact | Likelihood | Mitigation |
|---|---|---|---|
| Conda environment conflicts in containers | Tool fails at runtime | Medium | One tool per container; use Wave with pinned versions |
| SPAdes memory overflow on cloud | Task killed, retry cost | High | Dynamic memory with `task.attempt` multiplier; cap at r5.8xlarge |
| GTDB-Tk database too large for EBS | Task fails on disk | Medium | Use Fusion/S3 direct access or pre-stage on EFS |
| vRhyme non-deterministic output | Snapshot test failures | Low | Use fixed seed; tolerance-based assertions |
| VAMB GPU requirement | Can't run on CPU-only CEs | Low | Conditional GPU allocation; CPU fallback |
| Channel cardinality mismatch | Silent data loss | High | Explicit `.count()` assertions in tests; join discipline |
| Database version drift | Non-reproducible results | Medium | Pin DB versions in config; checksum validation |

### 10.2 Operational Risks

| Risk | Impact | Mitigation |
|---|---|---|
| Learning curve for NF DSL2 | Slower development | Use nf-core modules where possible; reference this plan |
| Parallel migration (maintaining both) | Double maintenance burden | Hard cutover after Phase 5 validation |
| Cloud cost overruns during testing | Budget concerns | Use test profile with minimal data; spot instances |
| Resume cache invalidation | Wasted compute | Consistent container tags; avoid changing publishDir |

### 10.3 Migration-Specific Risks

| Snakemake Pattern | Nextflow Equivalent | Risk |
|---|---|---|
| Wildcard expansion | Channel operators (map, combine) | Medium — requires channel shape discipline |
| `config["key"]` access | `params.key` | Low — direct mapping |
| `--use-conda` environment isolation | Wave containers | Low — equivalent functionality |
| Checkpoint rules (dynamic DAG) | Not applicable here | N/A — this pipeline is static DAG |
| `ruleorder` disambiguation | Process selection via channels | Medium — requires explicit routing |
| Touch sentinel files (`done.txt`) | Process `output: path "done.txt"` or implicit | Low — Nextflow tracks task completion natively |
| Script-relative paths (`scripts/`) | `bin/` directory or `template` directive | Low — `bin/` is auto-added to PATH |

---

## Appendix A: nf-core Module Availability Check

Modules confirmed available in nf-core/modules (as of 2026-07):

- ✅ fastqc, fastp, nanoplot, porechop_abi, filtlong
- ✅ bwamem2 (index + mem), minimap2 (index + align)
- ✅ samtools (sort, index, view, fastq, stats)
- ✅ megahit, quast, prodigal
- ✅ metabat2 (jgi_summarize_bam_contig_depths + metabat2)
- ✅ multiqc, bakta, amrfinderplus, abricate
- ⚠️ gtdbtk (classify_wf — check version compatibility)
- ❌ virsorter2, genomad, vibrant, vrhyme, vcontact3 (need local modules)
- ❌ checkv, checkm2, semibin2, comebin, binette (need local modules)
- ❌ pharokka, phold, defensefinder, rgi, deeparg (need local modules)
- ❌ cobra, phist, coverm, skani, galah, metamdbg (need local modules)

**Estimate:** ~15 nf-core modules reusable, ~35 local modules to author.

---

## Appendix B: Key Differences Summary

| Aspect | Snakemake | Nextflow |
|---|---|---|
| Execution model | Pull-based (targets) | Push-based (dataflow) |
| Parallelism | Per-rule threads | Per-process resources + executor |
| Environment | `--use-conda` per rule | Container per process (Docker/Singularity) |
| Configuration | `config.yaml` + CLI | `nextflow.config` + params + profiles |
| Sample handling | Wildcard expansion | Channel from samplesheet |
| Dependencies | Implicit (file-based) | Explicit (channel connections) |
| Cloud execution | Requires plugins (e.g., snakemake-executor-plugin-*) | Native (AWS Batch, Google Batch, Azure) |
| Resume | Per-rule via timestamps | Per-task via content hash |
| Reporting | Custom scripts | Seqera Platform + custom reports block |
| Monitoring | Log files | Seqera Platform real-time UI |

---

## Appendix C: Estimated Effort

| Phase | Weeks | Processes | Complexity |
|---|---|---|---|
| Phase 1: Foundation | 2 | ~15 | Medium (scaffolding + assembly) |
| Phase 2: Viral | 2 | ~12 | High (multi-tool consensus logic) |
| Phase 3: Prokaryotic | 2 | ~12 | High (multi-binner + refinement) |
| Phase 4: Annotation & Defense | 2 | ~15 | Medium (many tools, simple wiring) |
| Phase 5: Reporting & Platform | 2 | ~5 | Low (integration + polish) |
| **Total** | **10 weeks** | **~60 processes** | |

---

*Document generated: 2026-07-17*
*Pipeline version: VAPOR v1.0 (Snakemake)*
*Target: Nextflow DSL2 strict syntax (26.04+)*
