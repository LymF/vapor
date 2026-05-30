# VAPOR — Viral And Prokaryotic mOdular pipelineR

A modular Snakemake pipeline for metagenomics and viromics.
Supports **short reads (Illumina PE)** and **long reads (ONT / PacBio HiFi)**.

---

## Project structure

```
vapor/
│
├── Snakefile                    ← workflow entry point + rule includes
├── config.yaml                  ← all user-editable parameters
│
├── rules/                       ← modular rule files (one per functional block)
│   ├── qc.smk                   # BLOCK 1  — fastp, NanoPlot, Porechop_ABI, Filtlong
│   ├── assembly.smk             # BLOCK 2  — MEGAHIT, metaSPAdes, metaviralSPAdes,
│   │                            #            Flye, hifiasm, Medaka, metaMDBG
│   ├── merge_dedup.smk          # BLOCK 3  — merge_contigs, MMseqs2
│   ├── quast.smk                # BLOCK 4  — QUAST
│   ├── viral_detection.smk      # BLOCK 5  — VirSorter2, GeNomad, VIBRANT, consensus
│   ├── mapping.smk              # BLOCK 6  — BWA-MEM2 / minimap2, calc_depth
│   ├── viral_binning.smk        # BLOCK 7  — CheckV, vRhyme, CheckV(vRhyme)
│   ├── prok_binning.smk         # BLOCK 8  — MetaBAT2, VAMB, SemiBin2, COMEBin,
│   │                            #            Binette, CheckM2, GTDB-Tk
│   ├── taxonomy.smk             # BLOCK 9  — Prodigal, Diamond, vConTACT3
│   ├── host_prediction.smk      # BLOCK 10 — PHIST
│   ├── annotation.smk           # BLOCK 11 — Pharokka, Phold, Bakta, EggNOG, genome maps
│   ├── abundance.smk            # BLOCK 12 — CoverM (viral + prok), diversity
│   ├── finalize.smk             # BLOCK 13 — organize_outputs, benchmarks
│   └── report.smk               # BLOCK 14 — interactive HTML report, MultiQC
│
├── scripts/                     ← auxiliary Python scripts
│   ├── filter_checkv_hq.py      # filter viral FASTA for vConTACT3 input
│   ├── merge_lr_assemblies.py   # merge Flye + hifiasm + metaMDBG assemblies
│   ├── prepare_diamond_db.py    # build Diamond DB + metadata TSV
│   ├── split_viral_fastas.py    # split viral FASTA into per-genome files for PHIST
│   ├── genome_map_universal.py  # circular genome maps (phage / virus / prok)
│   ├── compute_diversity.py     # alpha, beta diversity, Procrustes
│   ├── make_votu_table.py       # consolidated vOTU summary table
│   └── generate_report.py       # standalone Plotly HTML report
│
└── envs/                        ← reproducible conda environment definitions
    ├── env_qc.yaml
    ├── env_assembly.yaml        # includes MMseqs2 + metaMDBG
    ├── env_flye.yaml            # includes hifiasm + hifiasm_meta
    ├── env_medaka.yaml
    ├── env_lr_utils.yaml
    ├── env_mapping.yaml
    ├── env_viral.yaml
    ├── env_genomad.yaml
    ├── phage_vibrant.yaml
    ├── env_vrhyme.yaml
    ├── env_binning.yaml
    ├── env_comebin.yaml
    ├── env_binette.yaml
    ├── env_checkm2.yaml
    ├── env_gtdbtk.yaml
    ├── env_phist.yaml
    ├── env_annotation.yaml      # includes circular genome maps (pycirclize)
    ├── env_vcontact3.yaml
    └── env_coverm.yaml          # includes diversity (numpy + scipy)
```

---

## Running the pipeline

```bash
# Activate the Snakemake environment (required)
conda activate snakemake

# Dry-run — validate the workflow without executing any jobs
vapor --dry-run

# Full execution with 32 cores
vapor --threads 32

# Use a custom config file
vapor --threads 32 --config /path/to/config.yaml

# Visualize the DAG (requires graphviz)
vapor --dag

# Force re-run specific rules
vapor --threads 32 --forcerun viral_consensus

# Resume an interrupted run
vapor --threads 32 --rerun-incomplete

# Unlock directory after a crash
vapor --unlock

# Run up to a specific output file
vapor --threads 32 --target results/sample1/viral/taxonomy/taxonomy_done.txt
```

---

## Configuration

All parameters are defined in **`config.yaml`** — no `.smk` files need to be edited.

### Key parameters

| Variable | Description |
|---|---|
| `fastq_dir` | Input directory containing FASTQ files |
| `outdir` | Output root directory |
| `threads` | Default CPU threads per rule |
| `min_contig` | Minimum contig length after assembly (bp) |
| `long_reads` | `true` for Nanopore / PacBio data |
| `lr_tech` | `"ont"` or `"hifi"` |
| `viral_consensus_mode` | `"count"` / `"score"` / `"hybrid"` |
| `min_viral_tools` | Minimum tools agreeing for viral consensus |
| `use_gpu` | `true` to enable GPU in VAMB, SemiBin2, COMEBin, GeNomad |

### Database paths

```yaml
checkv_db:    "/path/to/checkv-db-v1.5"
vs2_db:       "/path/to/virsorter2"
genomad_db:   "/path/to/genomad_db"
vibrant_base: "/path/to/vibrant-1.0.1"
checkm2_db:   "/path/to/uniref100.KO.1.dmnd"
inphared_db:  "/path/to/inphared"
vcontact3_db: "/path/to/vcontact3"
gtdbtk_db:    "/path/to/gtdbtk/release226"
pharokka_db:  "/path/to/pharokka"
phold_db:     "/path/to/phold_db"
bakta_db:     "/path/to/bakta/db"
eggnog_db:    "/path/to/eggnog"

# Optional — leave "" to skip
custom_viral_dmnd: ""
custom_viral_meta: ""
custom_prok_dmnd:  ""
custom_prok_meta:  ""
```

---

## Conda environments

All environments are defined as YAML files in `envs/` for full reproducibility.
Create all environments at once (one-time setup):

```bash
snakemake --snakefile Snakefile --use-conda --cores 1 --create-envs-only
```

| Environment | Main tools |
|---|---|
| `env_qc` | fastp, quast, multiqc |
| `env_assembly` | megahit, spades, metaMDBG, mmseqs2 |
| `env_flye` | flye, hifiasm, hifiasm_meta |
| `env_medaka` | medaka (ONT polishing) |
| `env_lr_utils` | nanoplot, porechop_abi, filtlong |
| `env_mapping` | bwa-mem2, minimap2, samtools |
| `env_viral` | virsorter2, checkv, prodigal, diamond |
| `env_genomad` | genomad |
| `phage_vibrant` | vibrant |
| `env_vrhyme` | vrhyme |
| `env_binning` | metabat2, vamb, semibin2 |
| `env_comebin` | comebin |
| `env_binette` | binette |
| `env_checkm2` | checkm2 |
| `env_gtdbtk` | gtdb-tk |
| `env_phist` | phist, kmer-db |
| `env_annotation` | pharokka, phold, bakta, eggnog-mapper, pycirclize |
| `env_vcontact3` | vcontact3 |
| `env_coverm` | coverm, numpy, scipy |

---

## Central hub: `rep_seq.fasta`

Deduplicated representative sequences (`{sample}_rep_seq.fasta`, generated by MMseqs2 at 95% identity) serve as the central reference for all downstream analyses: viral detection, read mapping, binning, taxonomy assignment, and host prediction.

---

## Auxiliary scripts

### `prepare_diamond_db.py`
```bash
# IMG NR — viral subset
python3 scripts/prepare_diamond_db.py \
    --faa img_nr.faa --format img \
    --img-tax taxonOId2Taxonomy.tsv \
    --filter-domain Viruses \
    --out /path/to/img_viral --threads 32

# NCBI RefSeq viral
python3 scripts/prepare_diamond_db.py \
    --faa viral.protein.faa --format ncbi \
    --out /path/to/refseq_viral --threads 32
```

### `filter_checkv_hq.py`
```bash
python3 scripts/filter_checkv_hq.py \
    quality_summary.tsv viral_consensus.fasta hq_viral.fasta
```

### `split_viral_fastas.py`
```bash
python3 scripts/split_viral_fastas.py \
    viral_consensus.fasta vrhyme_dir/ output_fastas_dir/
```

---

## System requirements

| Component | Minimum | Recommended |
|---|---|---|
| OS | Linux (Ubuntu 20.04+) | Ubuntu 22.04 LTS |
| CPU | 16 cores | 32–64 cores |
| RAM | 64 GB | 128–256 GB |
| Disk (databases) | 400 GB | 600 GB |
| Disk (results) | 200 GB / project | SSD preferred |
