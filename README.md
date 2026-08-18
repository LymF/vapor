# VAPOR — Viral And Prokaryotic mOdular pipelineR

A modular Snakemake pipeline for metagenomics and viromics.
Supports **short reads (Illumina PE and SE)** and **long reads (ONT / PacBio HiFi)**.

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
│   ├── assembly.smk             # BLOCK 2  — MEGAHIT (SR); Flye+Medaka (ONT) / metaMDBG (HiFi)
│   ├── quast.smk                # BLOCK 4  — QUAST
│   ├── viral_detection.smk      # BLOCK 5  — VirSorter2, GeNomad, consensus
│   ├── mapping.smk              # BLOCK 6  — BWA-MEM2 / minimap2, calc_depth
│   ├── cobra.smk                # BLOCK 5.5— COBRA contig extension (optional, SR PE)
│   ├── viral_binning.smk        # BLOCK 7  — CheckV, provirus trim, vRhyme, CheckV(vRhyme)
│   ├── votu_catalog.smk         # BLOCK 7.5— global vOTU catalog: pool (namespaced IDs) →
│   │                            #            skani triangle --sparse → cluster (95% ANI +
│   │                            #            85% AF) → 3-tier reps → read recruitment →
│   │                            #            presence/abundance matrices. Replaces the
│   │                            #            old per-sample skani clustering.
│   ├── prok_binning.smk         # BLOCK 8  — viral→prok filter, MetaBAT2, VAMB,
│   │                            #            SemiBin2, Binette, GUNC,
│   │                            #            CheckM2, galah derep, GTDB-Tk
│   ├── taxonomy.smk             # BLOCK 9  — Prodigal, MMseqs2/INPHARED (LCA),
│   │                            #            MMseqs2/INPHARED, MMseqs2/Custom
│   ├── host_prediction.smk      # BLOCK 10 — PHIST
│   ├── defense_amr.smk          # BLOCK 10.5 — DefenseFinder/AMRFinderPlus/RGI/DeepARG/
│   │                            #              ABRicate/argNorm (bins) + DefenseFinder/
│   │                            #              dbAPIS on viral ORFs, defense islands
│   ├── annotation.smk           # BLOCK 11 — Pharokka, Phold, Bakta, EggNOG, genome maps
│   ├── abundance.smk            # BLOCK 12 — CoverM (viral + prok), diversity
│   ├── finalize.smk             # BLOCK 13 — organize_outputs, benchmarks
│   ├── report.smk               # BLOCK 14 — interactive HTML report, MultiQC
│   └── reads_classify.smk       # BLOCK 15 — reads-only classification (Sylph + sylph-tax)
│
├── scripts/                     ← auxiliary Python scripts
│   ├── prepare_diamond_db.py    # build Diamond DB + metadata TSV
│   ├── prepare_mmseqs_taxdb.py   # build MMseqs2 seqTaxDB (--format img/ncbi/inphared)
│   ├── split_viral_fastas.py    # split viral FASTA into per-genome files for PHIST
│   ├── genome_map_universal.py  # circular genome maps (phage / virus / prok)
│   ├── compute_diversity.py     # alpha, beta diversity, Procrustes
│   ├── make_votu_table.py       # consolidated vOTU summary table (joins the global catalog)
│   ├── votu_catalog.py          # pure logic for the global vOTU catalog: pooling,
│   │                            #   sparse-skani parsing, single-linkage clustering (ICTV)
│   ├── generate_report.py           # wrapper → scripts/report/ (ECharts + D3 HTML report)
│   ├── pin_containers.py            # resolve quay.io tags → containers.lock.yaml
│   └── reads_classify/              # reads-only classification helpers
│       ├── filter_by_prevalence.py  #   filter taxa by minimum prevalence across samples
│       ├── make_otu.py              #   reformat merged table to #OTU_ID (QIIME2/phyloseq)
│       ├── collapse_by_host.py      #   aggregate viral abundance by predicted host genus
│       ├── bacphlip_lifestyle.py    #   extract detected genomes + run BACPHLIP lifestyle
│       └── build_imgvr_taxonomy.py  #   convert IMG/VR _meta.tsv to sylph-tax custom taxonomy TSV
│
├── containers.yaml              ← container version definitions (edit to update)
├── containers.lock.yaml         ← resolved quay.io URIs (generated, commit this)
│
├── docker/                      ← custom Dockerfiles for tools without official images
│   ├── Dockerfile.genome-map    # pycirclize + matplotlib + biopython
│   ├── Dockerfile.medaka-gpu    # medaka with CUDA (optional GPU mode)
│
└── envs/                        ← reproducible conda environment definitions
    ├── env_qc.yaml
    ├── env_assembly.yaml        # includes MMseqs2 + metaMDBG
    ├── env_flye.yaml            # metaFlye (ONT assembly)
    ├── env_medaka.yaml
    ├── env_lr_utils.yaml
    ├── env_mapping.yaml
    ├── env_viral.yaml
    ├── env_genomad.yaml
    ├── env_vrhyme.yaml
    ├── env_cobra.yaml           # COBRA contig extension (optional)
    ├── env_binning.yaml
    ├── env_binette.yaml
    ├── env_checkm2.yaml
    ├── env_gtdbtk.yaml
    ├── env_phist.yaml
    ├── env_annotation.yaml      # includes circular genome maps (pycirclize)
    ├── env_coverm.yaml          # includes diversity (numpy + scipy)
    ├── env_gunc.yaml            # GUNC (MAG chimera detection)
    ├── env_derep.yaml           # skani + galah (vOTU clustering + MAG dereplication)
    └── env_reads_classify.yaml  # sylph + sylph-tax + BACPHLIP (reads-only classification)
```

---

## Running the pipeline

Activate the Snakemake environment first:

```bash
conda activate snakemake
```

`vapor` auto-detects the available runtime: **apptainer → singularity → conda**.

Bind mounts for Apptainer/Singularity are derived **automatically** from `config.yaml`
(`fastq_dir`, `outdir`, all database paths). No manual `--bind` needed.

```bash
# Dry-run (validate without executing)
vapor --dry-run

# Full execution (auto-detects runtime)
vapor --threads 32

# Force a specific executor
vapor --executor apptainer --threads 32
vapor --executor singularity --threads 32
vapor --executor conda --threads 32

# GPU pass-through for Apptainer (medaka, VAMB, SemiBin2)
vapor --executor apptainer --singularity-args '--nv' --threads 32

# Extra bind mounts beyond config paths (e.g. scratch on NAS)
vapor --executor apptainer --singularity-args '--bind /mnt/nas' --threads 32
```

### Container mode setup (one-time)

```bash
# Resolve exact image tags → containers.lock.yaml
python3 scripts/pin_containers.py
```

### Conda mode setup (one-time, offline-capable)

```bash
# Create all environments
snakemake --snakefile Snakefile --use-conda --cores 1 --create-envs-only
```

### Other commands

```bash
# Custom config file
vapor --threads 32 --config /path/to/config.yaml

# Visualize DAG (requires graphviz)
vapor --dag

# Force re-run specific rules
vapor --threads 32 --forcerun viral_consensus

# Resume interrupted run
vapor --threads 32 --rerun-incomplete

# Unlock directory after crash
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
| `use_gpu` | `true` to enable GPU in VAMB, SemiBin2, GeNomad |
| `single_end` | `true` for SE short reads (Illumina SE or Ion Torrent) |
| `prok_filter_viral` | `true` to remove free-living viral contigs from prok binners |
| `prok_filter_keep_provirus` | `true` to preserve provirus-containing contigs (default) |
| `gunc_enabled` | `true` to run GUNC chimera detection on Binette final bins |
| `gunc_db` | Path to `gunc_db_progenomes2.1.dmnd` |
| `votu_presence_min_coverage` | % of a vOTU representative a sample's reads must cover to count as "recruited" presence (default `75.0`, Roux et al. 2017) |
| `votu_recruit_min_identity` | Min. read identity (%) for catalog recruitment; `null` = 95 (short/HiFi) or 85 (ONT) |
| `mag_derep_enabled` | `true` to dereplicte MAGs with galah before GTDB-Tk |
| `mag_derep_ani` | ANI threshold for MAG dereplication (default `95.0`) |
| `cobra_enabled` | `true` to extend viral contigs with COBRA (SR PE only, default `false`) |
| `cobra_megahit_maxk` | COBRA k-mer max for MEGAHIT contigs (default `141`) |
| `reads_classify` | `true` to enable reads-only classification module (Sylph, default `false`) |
| `reads_classify_dbs` | Dict with up to 4 optional DB paths: `imgvr`, `uhgv`, `gtdb`, `custom` |
| `reads_classify_tax_dir` | Directory for sylph-tax taxonomy files (auto-populated if empty) |
| `reads_classify_min_prevalence` | Minimum fraction of samples a taxon must appear in (default `0.0`) |
| `reads_classify_min_kmers` | Minimum k-mer matches for sylph profiling (default `20`) |
| `reads_classify_genome_fasta` | Path to reference genome FASTA for BACPHLIP lifestyle prediction (optional) |
| `reads_classify_virulence_threshold` | BACPHLIP score threshold for virulent classification (default `0.5`) |

### Database paths

```yaml
checkv_db:    "/path/to/checkv-db-v1.5"
vs2_db:       "/path/to/virsorter2"
genomad_db:   "/path/to/genomad_db"
checkm2_db:   "/path/to/uniref100.KO.1.dmnd"
inphared_db:  "/path/to/inphared"
gtdbtk_db:    "/path/to/gtdbtk/release226"
pharokka_db:  "/path/to/pharokka"
phold_db:     "/path/to/phold_db"
bakta_db:     "/path/to/bakta/db"
eggnog_db:    "/path/to/eggnog"

# Optional MMseqs2 seqTaxDBs — leave "" to skip (no Diamond option for either)
custom_prok_mmseqs_db:  ""   # e.g. IMG NR
custom_viral_mmseqs_db: ""   # e.g. IMG/VR
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
| `env_assembly` | megahit, metaMDBG, mmseqs2 (also reused by mmseqs_taxonomy_viral/prok) |
| `env_flye` | flye |
| `env_medaka` | medaka (ONT polishing) |
| `env_lr_utils` | nanoplot, porechop_abi, filtlong |
| `env_mapping` | bwa-mem2, minimap2, samtools |
| `env_viral` | virsorter2, checkv, prodigal, diamond |
| `env_genomad` | genomad |
| `env_vrhyme` | vrhyme |
| `env_cobra` | cobra-meta, blast (optional) |
| `env_binning` | metabat2, vamb, semibin2 |
| `env_binette` | binette |
| `env_checkm2` | checkm2 |
| `env_gtdbtk` | gtdb-tk |
| `env_phist` | phist, kmer-db |
| `env_annotation` | pharokka, phold, bakta, eggnog-mapper, pycirclize |
| `env_coverm` | coverm, numpy, scipy |
| `env_gunc` | gunc, diamond, prodigal (MAG chimera detection) |
| `env_derep` | skani, galah (vOTU clustering + MAG dereplication) |
| `env_reads_classify` | sylph, sylph-tax, pandas, biopython, bacphlip (reads-only classification) |

---

## Central hub: the assembly

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

---

## Genome quality improvements

VAPOR includes four optional quality-enhancement steps, all enabled by default and individually configurable in `config.yaml`:

| Feature | Rule | Benefit |
|---|---|---|
| **Viral → prok filter** | `filter_viral_for_prok` | Removes free-living viral contigs from prokaryotic binner input; provirus-containing contigs are preserved (detected via CheckV + GeNomad metadata) |
| **GUNC** | `gunc` | Detects chimeric MAGs by checking taxon consistency across Diamond-annotated genes; report appears in final summary |
| **skani vOTU catalog** | `votu_catalog_skani`, `votu_catalog_cluster` | Clusters the pooled viral genomes of every sample and co-assembly group at ICTV standard (95% ANI + 85% AF) using skani pairwise ANI; replaces the simpler MMseqs2 identity grouping, and clusters once globally instead of per sample |
| **galah MAG derep** | `galah_derep` | Dereplicates prokaryotic bins using CheckM2 quality scores; selects highest-quality representative per cluster before GTDB-Tk |
| **COBRA** | `cobra_megahit`, `cobra_merge` | Extends fragmented viral contigs by traversing the MEGAHIT assembly graph k-mer overlap; longest extension wins; disabled by default — beneficial for low-diversity viromes |

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

### `prepare_mmseqs_taxdb.py`
Builds an MMseqs2 seqTaxDB (real per-query LCA via `mmseqs taxonomy`,
avoids the "spurious specificity" of best-hit-only methods). Same
multi-format pattern as `prepare_diamond_db.py` above, but for MMseqs2's
`taxonomy` module — `--format img/ncbi/inphared/imgvr`, each hardcoding the
header-parsing + lineage-loading logic for its own source (not
interchangeable). Not auto-built by the pipeline for any format — run once
manually before processing samples:
```bash
conda activate env_assembly  # already has mmseqs2, no dedicated env needed

# INPHARED — used by rule mmseqs_taxonomy_viral (replaces an earlier
# Diamond/INPHARED best-hit tier; this is now the only INPHARED-based source)
python3 scripts/prepare_mmseqs_taxdb.py \
    --faa  14Apr2025_vConTACT2_proteins.faa --format inphared \
    --inphared-tax 14Apr2025_data.tsv \
    --out  /path/to/inphared/inphared_mmseqs_taxdb --threads 32

# IMG NR — used by rule mmseqs_taxonomy_prok (custom_prok_mmseqs_db)
python3 scripts/prepare_mmseqs_taxdb.py \
    --faa img_unrestricted_isolates_nr.faa --format img \
    --img-tax taxonOId2Taxonomy.tsv \
    --out /path/to/img_nr_mmseqs --threads 32

# NCBI RefSeq/GenBank — for a custom seqTaxDB from a plain protein set
python3 scripts/prepare_mmseqs_taxdb.py \
    --faa viral.1.protein.faa --format ncbi \
    --ncbi-tax taxonomy.tsv \
    --out /path/to/refseq_viral_mmseqs --threads 32

# IMG/VR — used by rule mmseqs_taxonomy_custom_viral (custom_viral_mmseqs_db;
# replaces an earlier diamond_custom_viral best-hit tier, removed 2026-06-23)
python3 scripts/prepare_mmseqs_taxdb.py \
    --faa IMGVR_all_proteins-high_confidence.faa --format imgvr \
    --imgvr-tax IMGVR_all_Sequence_information-high_confidence.tsv \
    --out /path/to/imgvr_mmseqs --threads 32
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

### `votu_catalog.py`
Pure-logic module for the global vOTU catalog (`rules/votu_catalog.smk`): pooling every
sample's/group's viral set with namespaced IDs, parsing `skani triangle --sparse` output,
single-linkage clustering at the ICTV threshold (95% ANI + 85% AF), and picking
representatives. No CLI — it is imported directly by the `votu_catalog_*` rules and
covered by `tests/test_votu_catalog.py`.

### `reads_classify/build_imgvr_taxonomy.py`

Converts an IMG/VR `_meta.tsv` file (with columns: `accession`, `phylum`, `class`, `order`,
`family`, `genus`, `organism`, `domain`) to the two-column TSV format required by
`sylph-tax taxprof -t` for custom database annotation.

```bash
python3 scripts/reads_classify/build_imgvr_taxonomy.py \
    IMGVR_v7.1_meta.tsv \
    imgvr_v71_taxonomy.tsv
```

See [Sylph pre-built databases](https://sylph-docs.github.io/pre%E2%80%90built-databases/) for the
list of pre-built sylph DBs (IMG/VR 4.1, UHGV, GTDB r232) and their download links.

---

## System requirements

| Component | Minimum | Recommended |
|---|---|---|
| OS | Linux (Ubuntu 20.04+) | Ubuntu 22.04 LTS |
| CPU | 16 cores | 32–64 cores |
| RAM | 64 GB | 128–256 GB |
| Disk (databases) | 400 GB | 600 GB |
| Disk (results) | 200 GB / project | SSD preferred |
