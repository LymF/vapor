# Installation Guide — VAPOR

This guide covers installing all software dependencies and databases required by the pipeline.
Follow the steps in order — database downloads are the most time-consuming part (~400 GB total).

---

## Table of Contents

1. [System requirements](#1-system-requirements)
2. [Install Miniforge3](#2-install-miniforge3)
3. [Install Snakemake](#3-install-snakemake)
4. [Create conda environments](#4-create-conda-environments)
5. [Install Apptainer (container mode)](#5-install-apptainer-container-mode)
6. [Install databases](#6-install-databases)
7. [Optional: custom Diamond databases](#7-optional-custom-diamond-databases)
8. [Configure config.yaml](#8-configure-configyaml)
9. [Verify installation](#9-verify-installation)

---

## 1. System requirements

| Component | Minimum | Recommended |
|---|---|---|
| OS | Linux (Ubuntu 20.04+) | Ubuntu 22.04 LTS |
| CPU | 16 cores | 32–64 cores |
| RAM | 64 GB | 128–256 GB |
| Disk (databases) | 400 GB | 600 GB |
| Disk (results) | 200 GB / project | SSD preferred |

> **Memory note:** metaSPAdes routinely uses 80–150 GB RAM for large metagenomes. vConTACT3 can use 30–60 GB. Plan accordingly.

---

## 2. Install Miniforge3

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda init bash
source ~/.bashrc
```

---

## 3. Install Snakemake

Snakemake must run in its own dedicated environment.

```bash
mamba create -n snakemake -c conda-forge -c bioconda \
    "snakemake>=8.0" python=3.11 -y
conda activate snakemake
snakemake --version   # should print 8.x.x
```

---

## 4. Create conda environments

Environments are defined in `envs/*.yaml`. Snakemake can create them automatically:

```bash
conda activate snakemake
snakemake --use-conda --cores 1 --create-envs-only
```

Or create manually (useful for testing individual environments):

```bash
# Short-read QC
mamba create -n env_qc -c conda-forge -c bioconda \
    fastp quast multiqc -y

# Assembly (MEGAHIT + SPAdes + metaMDBG + MMseqs2)
mamba create -n env_assembly -c conda-forge -c bioconda \
    megahit spades metamdbg mmseqs2 -y

# Long-read assembly: Flye + hifiasm (merged)
mamba create -n env_flye -c conda-forge -c bioconda \
    flye hifiasm hifiasm_meta -y

# ONT polishing (GPU via CUDA)
mamba create -n env_medaka -c conda-forge -c bioconda \
    medaka -y

# Long-read QC
mamba create -n env_lr_utils -c conda-forge -c bioconda \
    nanoplot filtlong porechop_abi -y

# Read mapping
mamba create -n env_mapping -c conda-forge -c bioconda \
    bwa-mem2 minimap2 samtools -y

# Viral detection + taxonomy
mamba create -n env_viral -c conda-forge -c bioconda \
    "virsorter=2" diamond prodigal checkv -y

# GeNomad
mamba create -n env_genomad -c conda-forge -c bioconda \
    genomad -y

# VIBRANT 1.2.1
mamba create -n phage_vibrant -c conda-forge -c bioconda \
    "vibrant=1.2.1" hmmer prodigal -y

# vRhyme
mamba create -n env_vrhyme -c conda-forge -c bioconda \
    vrhyme -y

# COBRA (viral contig extension — optional, short-read PE only)
mamba create -n env_cobra -c conda-forge -c bioconda \
    "cobra-meta=1.3.0" blast -y

# Prokaryotic binning
mamba create -n env_binning -c conda-forge -c bioconda \
    metabat2 vamb semibin2 -y

# COMEBin (transformer-based, rank 1 in 2025 benchmark)
mamba create -n env_comebin -c conda-forge -c bioconda -c pytorch \
    comebin pytorch -y

# Binette
mamba create -n env_binette -c conda-forge -c bioconda \
    binette -y

# CheckM2
mamba create -n env_checkm2 -c conda-forge -c bioconda \
    checkm2 -y

# GTDB-Tk
mamba create -n env_gtdbtk -c conda-forge -c bioconda \
    "gtdbtk=2.4" -y

# PHIST
mamba create -n env_phist -c conda-forge -c bioconda \
    phist -y

# Annotation (Pharokka + Phold + Bakta + EggNOG + circular genome maps)
mamba create -n env_annotation -c conda-forge -c bioconda \
    pharokka phold bakta eggnog-mapper pycirclize matplotlib biopython -y

# vConTACT3
mamba create -n env_vcontact3 -c conda-forge -c bioconda \
    vcontact3 -y

# CoverM + diversity (alpha, beta, Procrustes)
mamba create -n env_coverm -c conda-forge -c bioconda \
    coverm numpy scipy -y

# GUNC (MAG chimera detection)
mamba create -n env_gunc -c conda-forge -c bioconda \
    "gunc=1.1.1" diamond prodigal -y

# skani + galah (vOTU clustering + MAG dereplication)
mamba create -n env_derep -c conda-forge -c bioconda \
    "skani=0.3.1" "galah=0.4.2" -y
```

---

## 5. Install Apptainer (container mode)

Container mode is the recommended way to run VAPOR for reproducibility and portability.
Skip this section if you plan to use conda mode only (`vapor --executor conda`).

### Install Apptainer

```bash
# Ubuntu 22.04 / Debian 12
sudo apt-get update
sudo apt-get install -y software-properties-common
sudo add-apt-repository -y ppa:apptainer/ppa
sudo apt-get update
sudo apt-get install -y apptainer

apptainer --version   # should print 1.1+

# Create singularity → apptainer symlink (for Snakemake compatibility)
sudo ln -sf $(which apptainer) /usr/local/bin/singularity
```

### Populate the container lock file

Run once from the VAPOR directory. The script queries quay.io/biocontainers to
resolve the exact image tag for each tool and writes `containers.lock.yaml`.
All images — including custom ones published on GHCR — are already built and
available; this step only downloads the lock file metadata, not the images themselves
(those are pulled on first use by Apptainer).

```bash
conda activate snakemake
python3 scripts/pin_containers.py
```

> **Note:** `containers.lock.yaml` is committed to the repository and updated
> with each VAPOR release, so in most cases you can skip this step and use the
> lock file that ships with the code.

---

## 6. Install databases

Set a base directory for all databases:

```bash
DB_BASE="/path/to/your/databases"
mkdir -p "$DB_BASE"
```

> **Container mode (no conda required)**
> Every download step below has a Docker alternative — no conda environment is needed.
> The general pattern is:
> ```bash
> docker run --rm -v "$DB_BASE:/dbs" quay.io/biocontainers/<tool>:<tag> <command /dbs/...>
> ```
> Exact image tags (with build hash) can be obtained by running `python3 scripts/pin_containers.py`
> which generates `containers.lock.yaml`, or looked up directly at
> `https://quay.io/repository/biocontainers/<tool>?tab=tags`.
> The tool versions in `containers.yaml` are the reference versions used by VAPOR.

---

### CheckV

**Conda:**
```bash
conda activate env_checkm2
checkv download_database "$DB_BASE/checkv"
# Output: $DB_BASE/checkv/checkv-db-v1.5/
conda deactivate
```

**Docker:**
```bash
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/checkv:1.0.3--pyhdfd78af_0 \
    checkv download_database /dbs/checkv
```

---

### VirSorter2

**Conda:**
```bash
conda activate env_viral
virsorter setup -d "$DB_BASE/virsorter2" -j 4
conda deactivate
```

**Docker:**
```bash
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/virsorter:2.2.4--pyhdfd78af_1 \
    virsorter setup -d /dbs/virsorter2 -j 4
```

---

### GeNomad

**Conda:**
```bash
conda activate env_genomad
genomad download-database "$DB_BASE/genomad"
# Output: $DB_BASE/genomad/genomad_db/
conda deactivate
```

**Docker:**
```bash
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/genomad:1.11.2--pyhdfd78af_0 \
    genomad download-database /dbs/genomad
```

---

### VIBRANT

The databases are version-agnostic (same files for 1.0.1 and 1.2.1).
Set `vibrant_base` in `config.yaml` to whichever directory you create.

**Conda:**
```bash
mkdir -p "$DB_BASE/vibrant/databases"
mkdir -p "$DB_BASE/vibrant/files"
conda activate phage_vibrant

cd "$DB_BASE/vibrant/databases"
wget https://zenodo.org/record/4543735/files/VIBRANT_v1.0.1.tar.gz
tar -xzf VIBRANT_v1.0.1.tar.gz --strip-components=2 "VIBRANT_v1.0.1/databases/"

cd "$DB_BASE/vibrant/files"
git clone --depth 1 https://github.com/AnantharamanLab/VIBRANT /tmp/vibrant_src
cp /tmp/vibrant_src/files/* "$DB_BASE/vibrant/files/"
rm -rf /tmp/vibrant_src

# Press-index HMM databases
for hmm in "$DB_BASE/vibrant/databases/"*.HMM; do
    hmmpress "$hmm"
done
conda deactivate
```

**Container:** (wget outside container; hmmpress via hmmer image)
```bash
mkdir -p "$DB_BASE/vibrant/databases" "$DB_BASE/vibrant/files"
cd "$DB_BASE/vibrant/databases"
wget https://zenodo.org/record/4543735/files/VIBRANT_v1.0.1.tar.gz
tar -xzf VIBRANT_v1.0.1.tar.gz --strip-components=2 "VIBRANT_v1.0.1/databases/"
cd "$DB_BASE/vibrant/files"
git clone --depth 1 https://github.com/AnantharamanLab/VIBRANT /tmp/vibrant_src
cp /tmp/vibrant_src/files/* ./ && rm -rf /tmp/vibrant_src

# hmmpress via apptainer/singularity
apptainer exec docker://quay.io/biocontainers/hmmer:3.4--hdbdd923_1 \
    bash -c 'for hmm in '"$DB_BASE"'/vibrant/databases/*.HMM; do hmmpress "$hmm"; done'
```

> **Existing installs:** if your databases are already at `vibrant-1.0.1/`, you can keep
> that path and just point `vibrant_base` there — no need to move files.

---

### INPHARED

**Conda:**
```bash
mkdir -p "$DB_BASE/inphared"
cd "$DB_BASE/inphared"

# Check https://github.com/RyanCook94/inphared for the latest date tag
DATE="1Feb2024"
wget "https://millardlab-inphared.s3.climb.ac.uk/${DATE}_genomes.fa"
wget "https://millardlab-inphared.s3.climb.ac.uk/${DATE}_data_excluding_refseq.tsv"
wget "https://millardlab-inphared.s3.climb.ac.uk/${DATE}_vConTACT2_proteins.faa"
wget "https://millardlab-inphared.s3.climb.ac.uk/${DATE}_vConTACT2_gene2genome.csv"

conda activate env_viral
diamond makedb \
    --in "${DATE}_vConTACT2_proteins.faa" \
    --db "inphared_proteins" \
    --threads 32
conda deactivate
```

**Docker:** (wget needs no container; diamond via Docker)
```bash
mkdir -p "$DB_BASE/inphared" && cd "$DB_BASE/inphared"
DATE="1Feb2024"
wget "https://millardlab-inphared.s3.climb.ac.uk/${DATE}_genomes.fa"
wget "https://millardlab-inphared.s3.climb.ac.uk/${DATE}_data_excluding_refseq.tsv"
wget "https://millardlab-inphared.s3.climb.ac.uk/${DATE}_vConTACT2_proteins.faa"
wget "https://millardlab-inphared.s3.climb.ac.uk/${DATE}_vConTACT2_gene2genome.csv"

docker run --rm -v "$DB_BASE/inphared:/dbs" \
    quay.io/biocontainers/diamond:2.1.8--h43eeafb_0 \
    diamond makedb \
        --in "/dbs/${DATE}_vConTACT2_proteins.faa" \
        --db "/dbs/inphared_proteins" \
        --threads 32
```

---

### vConTACT3

**Conda:**
```bash
mkdir -p "$DB_BASE/vcontact3"
conda activate env_vcontact3
vcontact3 prepare_database \
    --output "$DB_BASE/vcontact3" \
    --threads 32
conda deactivate
# Note the version string printed — set vcontact3_ver in config.yaml accordingly
```

**Docker:**
```bash
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/vcontact3:3.1.6--pyhdfd78af_0 \
    vcontact3 prepare_database --output /dbs/vcontact3 --threads 32
```

---

### CheckM2

**Conda:**
```bash
conda activate env_checkm2
checkm2 database --download --path "$DB_BASE/checkm2"
# Output: $DB_BASE/checkm2/CheckM2_database/uniref100.KO.1.dmnd
conda deactivate
```

**Docker:**
```bash
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/checkm2:1.1.0--pyh7e72e81_0 \
    checkm2 database --download --path /dbs/checkm2
```

---

### GTDB-Tk

**Conda:**
```bash
mkdir -p "$DB_BASE/gtdbtk"
conda activate env_gtdbtk
download-db.sh "$DB_BASE/gtdbtk"
conda deactivate
```

**Docker:**
```bash
mkdir -p "$DB_BASE/gtdbtk"
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/gtdbtk:2.6.1--pyhdfd78af_0 \
    download-db.sh /dbs/gtdbtk
```

---

### GUNC (chimera detection)

Required only when `gunc_enabled: true` in `config.yaml` (default).
The progenomes_2.1 DIAMOND database is ~13 GB.

**Conda:**
```bash
mkdir -p "$DB_BASE/gunc"
conda activate env_gunc
gunc download_db -db progenomes_2.1 "$DB_BASE/gunc"
# Output: $DB_BASE/gunc/gunc_db_progenomes2.1.dmnd
conda deactivate
```

**Docker:**
```bash
mkdir -p "$DB_BASE/gunc"
docker run --rm -v "$DB_BASE/gunc:/dbs" \
    quay.io/biocontainers/gunc:1.1.1--pyhdfd78af_0 \
    gunc download_db -db progenomes_2.1 /dbs
# Output: $DB_BASE/gunc/gunc_db_progenomes2.1.dmnd
```

Then set in `config.yaml`:

```yaml
gunc_db: "/path/to/your/databases/gunc/gunc_db_progenomes2.1.dmnd"
```

To disable GUNC entirely set `gunc_enabled: false`.

---

### Pharokka

**Conda:**
```bash
conda activate env_annotation
pharokka_install_databases.py -o "$DB_BASE/pharokka"
conda deactivate
```

**Docker:**
```bash
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/pharokka:1.9.1--pyhdfd78af_0 \
    pharokka_install_databases.py -o /dbs/pharokka
```

---

### Phold

**Conda:**
```bash
conda activate env_annotation
phold install -d "$DB_BASE/phold_db"
conda deactivate
```

**Docker:**
```bash
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/phold:1.2.5--pyhdfd78af_0 \
    phold install -d /dbs/phold_db
```

---

### Bakta

**Conda:**
```bash
conda activate env_annotation
bakta_db download --output "$DB_BASE/bakta" --type full
conda deactivate
```

**Docker:**
```bash
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/bakta:1.12.0--pyhdfd78af_0 \
    bakta_db download --output /dbs/bakta --type full
```

---

### EggNOG-mapper

**Conda:**
```bash
conda activate env_annotation
export EGGNOG_DATA_DIR="$DB_BASE/eggnog"
download_eggnog_data.py -y -f --data_dir "$DB_BASE/eggnog"
conda deactivate
```

**Docker:**
```bash
docker run --rm -v "$DB_BASE:/dbs" \
    quay.io/biocontainers/eggnog-mapper:2.1.13--pyhdfd78af_0 \
    download_eggnog_data.py -y -f --data_dir /dbs/eggnog
```

---

## 7. Optional: custom Diamond databases

Custom databases allow classifying contigs and bins not covered by primary databases.

```bash
# IMG NR — viral subset
conda activate env_viral
python3 scripts/prepare_diamond_db.py \
    --faa img_nr.faa --format img \
    --img-tax taxonOId2Taxonomy.tsv \
    --filter-domain Viruses \
    --out "$DB_BASE/img/img_viral" --threads 32
# Output: img_viral.dmnd + img_viral_meta.tsv

# IMG NR — prokaryote subset
python3 scripts/prepare_diamond_db.py \
    --faa img_nr.faa --format img \
    --img-tax taxonOId2Taxonomy.tsv \
    --filter-domain Bacteria,Archaea \
    --out "$DB_BASE/img/img_prok" --threads 32
conda deactivate
```

---

## 8. Configure config.yaml

```yaml
# Input / output
fastq_dir: "fastqs"
outdir:    "results"
threads:   32

# Data type
long_reads: false
lr_tech:    "ont"   # "ont" | "hifi"

# Database paths — update to match your installation
checkv_db:     "/path/to/checkv/checkv-db-v1.5"
vs2_db:        "/path/to/virsorter2"
genomad_db:    "/path/to/genomad/genomad_db"
vibrant_base:  "/path/to/vibrant-1.0.1"
checkm2_db:    "/path/to/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
gtdbtk_db:     "/path/to/gtdbtk/release226"
inphared_db:   "/path/to/inphared"
vcontact3_db:  "/path/to/vcontact3"
vcontact3_ver: "230"
pharokka_db:   "/path/to/pharokka"
phold_db:      "/path/to/phold_db"
bakta_db:      "/path/to/bakta/db"
eggnog_db:     "/path/to/eggnog"

# Custom databases — leave "" to skip
custom_viral_dmnd: ""
custom_viral_meta: ""
custom_prok_dmnd:  ""
custom_prok_meta:  ""
```

---

## 9. Verify installation

```bash
conda activate snakemake

# Create dummy FASTQs to test sample detection
mkdir -p fastqs
touch fastqs/TEST_R1.fastq.gz fastqs/TEST_R2.fastq.gz

# Dry run via VAPOR CLI
vapor --dry-run 2>&1 | tail -20

# Clean up
rm fastqs/TEST_R1.fastq.gz fastqs/TEST_R2.fastq.gz
```

### Check individual environments

```bash
for env in env_qc env_assembly env_flye env_medaka env_lr_utils \
           env_mapping env_viral env_genomad phage_vibrant env_vrhyme \
           env_cobra env_binning env_comebin env_binette env_checkm2 \
           env_gtdbtk env_phist env_annotation env_vcontact3 env_coverm \
           env_gunc env_derep; do
    echo -n "$env: "
    conda run -n "$env" python --version 2>/dev/null || echo "MISSING"
done
```

---

## Disk space summary

| Database | Approximate size |
|---|---|
| CheckV v1.5 | 3 GB |
| VirSorter2 | 4 GB |
| GeNomad | 3 GB |
| VIBRANT v1.0.1 | 12 GB |
| INPHARED | 2 GB |
| vConTACT3 | 15 GB |
| CheckM2 | 3 GB |
| GTDB-Tk release 226 | 85 GB |
| Pharokka (PHROGS) | 1 GB |
| Phold | 10 GB |
| Bakta (full) | 30 GB |
| EggNOG-mapper | 50 GB |
| GUNC progenomes_2.1 | 13 GB |
| IMG viral (optional) | 5–20 GB |
| IMG prokaryote (optional) | 50–200 GB |
| **Total (without custom)** | **~231 GB** |

---

## Troubleshooting

**metaSPAdes runs out of memory**
Reduce `spades_mem` in `config.yaml`. If RAM is the bottleneck, use MEGAHIT only.

**VIBRANT: `-f older` error**
This pipeline uses `cd` to the output directory before running VIBRANT, passing only `-f nucl` explicitly. Ensure you are using the current version of the `Snakefile`.

**vConTACT3 takes too long**
vConTACT3 can take several hours for large viral datasets. Run with `--cores 32` and confirm that `vcontact3_ver` in `config.yaml` matches your installed database version.

**GTDB-Tk: pplacer error**
Confirm that `gtdbtk_db` points to the correct directory and that the database version is compatible with your installed GTDB-Tk (`gtdbtk check_install`).

**GUNC: `gunc_db` not found**
Set the full path to `gunc_db_progenomes2.1.dmnd` in `config.yaml` (key `gunc_db`). To disable GUNC entirely set `gunc_enabled: false`.

**galah: no bins to dereplicute**
galah requires CheckM2 quality scores. If CheckM2 failed or produced no output, galah will exit with an error. Set `mag_derep_enabled: false` to skip dereplication and feed all Binette bins directly to GTDB-Tk.

**COBRA: nearly all contigs report `orphan_end`**
This is expected for highly diverse environmental communities (seawater, soil, wastewater) where the assembly graph is fragmented. COBRA extends contigs using paired-end k-mer overlap; if the paired reads do not cover both ends of a contig, extension is not possible. Leave `cobra_enabled: false` (default) for such samples. Enable it primarily for low-diversity viromes where circular genomes are likely.

**COBRA: not applicable for long reads**
COBRA was designed for paired-end short-read assemblies. Set `cobra_enabled: false` whenever `long_reads: true`.
