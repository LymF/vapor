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
7. [Reads-only classification databases (reads_classify module, optional)](#7-reads-only-classification-databases-reads_classify-module-optional)
8. [Optional: custom MMseqs2 seqTaxDBs](#8-optional-custom-mmseqs2-seqtaxdbs)
9. [Configure config.yaml](#9-configure-configyaml)
10. [Verify installation](#10-verify-installation)

---

## 1. System requirements

| Component | Minimum | Recommended |
|---|---|---|
| OS | Linux (Ubuntu 20.04+) | Ubuntu 22.04 LTS |
| CPU | 16 cores | 32–64 cores |
| RAM | 64 GB | 128–256 GB |
| Disk (databases) | 400 GB | 600 GB |
| Disk (results) | 200 GB / project | SSD preferred |

> **Memory note:** metaSPAdes routinely uses 80–150 GB RAM for large metagenomes. Plan accordingly.

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

# Assembly (MEGAHIT + SPAdes + metaMDBG + MMseqs2). MMseqs2 here is also
# reused by rules mmseqs_taxonomy_viral/mmseqs_taxonomy_prok -- no separate
# env needed for those.
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

# vRhyme
mamba create -n env_vrhyme -c conda-forge -c bioconda \
    vrhyme -y

# COBRA (viral contig extension — optional, short-read PE only)
mamba create -n env_cobra -c conda-forge -c bioconda \
    "cobra-meta=1.3.0" blast -y

# Prokaryotic binning
mamba create -n env_binning -c conda-forge -c bioconda \
    metabat2 vamb semibin2 -y


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

# Defense systems: DefenseFinder (+ built-in AntiDefenseFinder)
# PADLOC was evaluated as a 2nd complementary detector but dropped -- its
# biocontainers image ships a broken BusyBox `rm` that can never build the
# database inside a container (confirmed on litrp4), and every tool in this
# pipeline must work via Docker/Apptainer, not conda-only.
mamba create -n env_defense -c conda-forge -c bioconda \
    "defense-finder=3.0.0" hmmer -y

# RGI / CARD (curated AMR) — isolated env, RGI's pins are finicky alongside others
mamba create -n env_rgi -c conda-forge -c bioconda \
    "rgi=6.0.5" -y

# DeepARG (exploratory/deep-learning AMR) — isolated env, legacy Python2/Theano stack
mamba create -n env_deeparg -c conda-forge -c bioconda \
    "deeparg=1.0.4" -y

# Annotation (Pharokka + Phold + Bakta + EggNOG + circular genome maps)
# AMRFinderPlus (ncbi-amrfinderplus) is already pinned here too — curated AMR
# reuses this env instead of a dedicated one.
mamba create -n env_annotation -c conda-forge -c bioconda \
    pharokka phold bakta eggnog-mapper ncbi-amrfinderplus pycirclize matplotlib biopython -y

# CoverM + diversity (alpha, beta, Procrustes)
mamba create -n env_coverm -c conda-forge -c bioconda \
    coverm numpy scipy -y

# GUNC (MAG chimera detection)
mamba create -n env_gunc -c conda-forge -c bioconda \
    "gunc=1.1.1" diamond prodigal -y

# skani + galah (vOTU clustering + MAG dereplication)
mamba create -n env_derep -c conda-forge -c bioconda \
    "skani=0.3.1" "galah=0.4.2" -y

# Reads-only classification (Sylph + sylph-tax + BACPHLIP)
# Required only when reads_classify: true in config.yaml
mamba create -n env_reads_classify -c conda-forge -c bioconda \
    "sylph>=0.6" "sylph-tax>=1.7" python pandas biopython -y
conda activate env_reads_classify
pip install bacphlip
conda deactivate
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

### INPHARED

Files are available at **https://github.com/RyanCook94/inphared** — go to the repository
and download the two required files manually (click the filename → "Raw" → Save as):

| File | Purpose |
|------|---------|
| `14Apr2025_data.tsv` | Phage metadata (taxonomy: Realm/Kingdom/Phylum/Class/Order/Family/Sub-family/Genus) |
| `14Apr2025_vConTACT2_proteins.faa` | Protein sequences for MMseqs2 taxonomy |

> Check the repository for a newer date prefix (e.g. `Jan2026`) and use that instead of `14Apr2025`.

After downloading both files to `$DB_BASE/inphared`, build the MMseqs2 seqTaxDB that
`rule mmseqs_taxonomy_viral` queries against (real per-query LCA, see
`scripts/prepare_mmseqs_taxdb.py --format inphared` -- the same script also builds
the IMG_NR seqTaxDB for prokaryote bins via `--format img`, and supports
`--format ncbi` for plain NCBI RefSeq/GenBank protein sets). It is **not**
auto-built on first run -- the rule just skips gracefully if it's missing (same
pattern as `mmseqs_taxonomy_prok`/`custom_prok_mmseqs_db`) -- so build it once
manually, before running the pipeline:

```bash
mkdir -p "$DB_BASE/inphared"
# Copy downloaded files to $DB_BASE/inphared, then:
cd "$DB_BASE/inphared"

conda activate env_assembly
python3 /path/to/pipe-final-claude/scripts/prepare_mmseqs_taxdb.py \
    --faa  $(ls *_vConTACT2_proteins.faa | sort | tail -1) \
    --format inphared \
    --inphared-tax $(ls *_data.tsv | sort | tail -1) \
    --out  ./inphared_mmseqs_taxdb --threads 32
conda deactivate
```

Output: `$DB_BASE/inphared/inphared_mmseqs_taxdb/seqTaxDB` (the exact path the rule
expects).

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

### DefenseFinder (defense systems / anti-defense)

Required only when `defense_amr_enabled: true` in `config.yaml` (default).
The database is small (HMM models, a few hundred MB) and auto-downloaded by
`rules/defense_amr.smk` on first run — **no `config.yaml` path is needed**.
If your compute nodes have no internet access (common on HPC clusters), pre-fetch
it once on a login/head node.

> PADLOC was evaluated as a complementary 2nd defense-system detector (catches
> some systems DefenseFinder misses) but was dropped: its biocontainers image
> ships a BusyBox `rm` that lacks the `-d` flag PADLOC's own `--db-update`
> script needs, so the database can never be built inside that container
> (confirmed broken on litrp4 2026-06-19). It runs fine via conda, but every
> tool in this pipeline must work via Docker/Apptainer, not conda-only.

`rules/defense_amr.smk` passes an explicit `--models-dir` to both `defense-finder
update` and `defense-finder run`, set via `defense_finder_models_db` in
`config.yaml` — **do not rely on the tool's own default**
(`$HOME/.macsyfinder`). Which `$HOME` that resolves to is not stable across
conda/Apptainer/cwd, and a batch run can silently end up split across two
different, each only-partially-populated caches: confirmed on litrp4
2026-06-20, where a multi-sample run left every sample with 0 systems even
though MacSyFinder itself finished cleanly (its own tmp output was fully
populated) — the run was reading/writing `$HOME/vapor/.macsyfinder` in one
context and `$HOME/.macsyfinder` in another. Pre-fetch into the *same* path
you put in `defense_finder_models_db`:

**Conda:**
```bash
mkdir -p "$DB_BASE/defense_finder_models"
conda activate env_defense
defense-finder update --models-dir "$DB_BASE/defense_finder_models"
conda deactivate
```
> **Pin `defense-finder=3.0.0`, not `2.0.0`.** Earlier `defense-finder` releases
> (2.0.0/2.0.1) fail against current CasFinder releases with
> `macsypy.error.MacsypyError: ... has not the right version. version
> supported is '2.0'` — reported in
> [mdmparis/defense-finder#95](https://github.com/mdmparis/defense-finder/issues/95)
> (same macsyfinder 2.1.4 + CasFinder 3.1.1 combo), with a maintainer attempt to
> pin CasFinder to a compatible version abandoned in
> [#101](https://github.com/mdmparis/defense-finder/issues/101). `3.0.0` resolves
> a compatible CasFinder (3.1.0) automatically and runs clean — confirmed live
> on litrp4 2026-06-19. The env/container pins in this repo already use `3.0.0`.

**Docker:**
```bash
mkdir -p "$DB_BASE/defense_finder_models"
docker run --rm -v "$DB_BASE/defense_finder_models:/models" \
    quay.io/biocontainers/defense-finder:3.0.0--pyhdfd78af_0 \
    defense-finder update --models-dir /models
```

Then set in `config.yaml`:
```yaml
defense_finder_models_db: "/path/to/your/databases/defense_finder_models"
```

> **GitHub rate limit.** `defense-finder update` calls GitHub's API even when
> models are already cached, which exhausts the unauthenticated 60-req/hour
> limit fast across a multi-sample batch (confirmed on litrp4 2026-06-20: every
> sample after the first few failed with `"maximum number of request per
> hour"`). `rules/defense_amr.smk` already skips the call once
> `defense_finder_models_db` is non-empty on disk, so this only matters for the
> very first run — set `GITHUB_TOKEN` in the environment beforehand if you hit
> it anyway.

`rules/defense_amr.smk` also reuses this exact DefenseFinder/AntiDefenseFinder
setup on **viral** proteins (`rule defensefinder_viral`, gated by
`defense_amr_viral_enabled` in `config.yaml`) — same models, same
`defense_finder_models_db`, nothing extra to configure.

---

### dbAPIS (viral-side anti-defense, complementary to DefenseFinder)

Yan, Y., Zheng, J., Zhang, X., & Yin, Y. (2023). *dbAPIS: a database of
anti-prokaryotic immune system genes*. Nucleic Acids Research.
https://doi.org/10.1093/nar/gkad932 — https://pro.unl.edu/dbAPIS
(moved from the old bcb.unl.edu host; bcb.unl.edu now 302-redirects to the
bare /dbAPIS homepage instead of the requested file, so any old direct
bcb.unl.edu/dbAPIS/downloads/<file> URL silently fetches an HTML page
instead of the real file — confirmed 2026-06-23. Use pro.unl.edu's
`download_file.php?file=<name>` endpoint instead.)

`rule dbapis_viral` (`rules/defense_amr.smk`) runs DIAMOND blastp (already a
pipeline dependency, no new tool) against dbAPIS's curated anti-defense
protein set on viral ORFs (`rules.prodigal_viral.output.faa`). dbAPIS is
sequence-similarity-based (no genetic-architecture rule like
MacSyFinder/DefenseFinder), better suited to the single scattered
anti-defense genes typically found in small phage genomes — kept as a
separate, complementary detector, never merged with DefenseFinder's calls
(same "never merge tiers" rule as AMR curated/exploratory). Reported
separately in the report's "Viral Anti-Defense" tab.

The database is tiny (~4,400 curated proteins, a few MB) — auto-downloaded
and cached on first run, same pattern as `card_db`/`deeparg_db`:

```bash
mkdir -p "$DB_BASE/dbapis"
wget -O "$DB_BASE/dbapis/anti_defense.pep" \
    "https://pro.unl.edu/dbAPIS/download_file.php?file=anti_defense.pep"
# Family -> (gene name, inhibited defense-system) mapping, one row per
# APIS family -- used by the report to translate bare family IDs like
# "APIS331" into a readable label (see load_dbapis_viral in
# scripts/report/data_loaders.py).
wget -O "$DB_BASE/dbapis/seed_and_familyrep_all_infor.tsv" \
    "https://pro.unl.edu/dbAPIS/download_file.php?file=seed_and_familyrep_all_infor.tsv"
conda activate env_viral   # has diamond
diamond makedb --in "$DB_BASE/dbapis/anti_defense.pep" -d "$DB_BASE/dbapis/APIS_db"
conda deactivate
```

Then set in `config.yaml`:
```yaml
apis_db: "/path/to/your/databases/dbapis"
```

If `apis_db` is left empty, the rule auto-populates `{outdir}/dbapis_db` on
first run instead — fine for a single-machine setup, but pre-fetching into a
shared path (like the example above) avoids re-downloading per output
directory. Gated by `defense_amr_viral_enabled` in `config.yaml` (default:
enabled).

---

### ABRicate (VFDB + PlasmidFinder screening)

Self-contained: every database ABRicate supports ships bundled with the
conda/container package itself — no separate download step. `rules/defense_amr.smk`
only screens **VFDB** (virulence factors) and **PlasmidFinder** (plasmid
replicons), the two databases not already covered by AMRFinderPlus/RGI/DeepARG
above; it is not used for AMR calling (ABRicate's bundled AMR databases are
flat BLASTN screens with no point-mutation or SNP/variant models, a downgrade
vs. the curated/ML tools already in the pipeline). Nothing to configure beyond
`abricate_enabled: true` in `config.yaml` (default).

---

### argNorm (AMR → ARO normalization)

Also self-contained — no database to pre-fetch, it ships its own ARO mapping
tables. `rules/defense_amr.smk` normalizes AMRFinderPlus and DeepARG gene
calls onto the Antibiotic Resistance Ontology (ARO) so they can be compared
with each other and with RGI (RGI already speaks ARO natively via CARD, so
it is **not** routed through argNorm — a hAMRonization bridge for this was
tried and dropped: argNorm 1.1.0 has no working RGI support despite the docs
implying otherwise, confirmed on litrp4 2026-06-20). Nothing to configure
beyond `argnorm_enabled: true` in `config.yaml` (default).

---

### AMRFinderPlus (curated AMR)

Self-managed: `rules/defense_amr.smk` runs `amrfinder -u` automatically on first
use. AMRFinderPlus installs its database under the conda environment's own
install prefix (e.g. `$CONDA_PREFIX/share/amrfinderplus/data/`), not `$HOME` —
so, like DefenseFinder above, a bare `docker run --rm <image> amrfinder -u`
with no volume mount downloads it into a throwaway container and loses it the
instant the container exits. Use conda to pre-fetch on an offline-compute node:

**Conda (recommended):**
```bash
conda activate env_annotation
amrfinder -u
conda deactivate
```

**Docker:** not recommended for pre-fetching — the exact internal install path
varies by image build, so there is no `-v` mount that's reliably correct across
versions. If you must use Docker, inspect the image first
(`docker run --rm --entrypoint sh quay.io/biocontainers/ncbi-amrfinderplus:4.2.7--hf69ffd2_0 -c 'amrfinder --database_path'`
or check its Dockerfile) to find the real path, then mount that directory.
Exact tag confirmed via `https://quay.io/api/v1/repository/biocontainers/ncbi-amrfinderplus/tag/` —
re-check there if this stops resolving (BioContainers rebuilds bump the hash suffix).

---

### RGI / CARD (curated AMR)

Unlike the tools above, CARD has no built-in auto-update — `rules/defense_amr.smk`
downloads and loads it automatically into `card_db` on first run, but you can
pre-fetch it the same way:

**Conda:**
```bash
mkdir -p "$DB_BASE/card"
curl -sL https://card.mcmaster.ca/latest/data -o "$DB_BASE/card/card_data.tar.bz2"
tar -xjf "$DB_BASE/card/card_data.tar.bz2" -C "$DB_BASE/card"   # archive members are "./card.json" etc — extract everything
conda activate env_rgi
rgi load --card_json "$DB_BASE/card/card.json" --local
conda deactivate
```

**Docker:**
```bash
mkdir -p "$DB_BASE/card"
curl -sL https://card.mcmaster.ca/latest/data -o "$DB_BASE/card/card_data.tar.bz2"
tar -xjf "$DB_BASE/card/card_data.tar.bz2" -C "$DB_BASE/card"   # archive members are "./card.json" etc — extract everything
docker run --rm -v "$DB_BASE/card:/dbs" quay.io/biocontainers/rgi:6.0.5--pyh05cac1d_0 \
    rgi load --card_json /dbs/card.json --local
```
Exact tag confirmed via `https://quay.io/api/v1/repository/biocontainers/rgi/tag/` —
re-check there if this stops resolving (BioContainers rebuilds bump the hash suffix).

Then set in `config.yaml`:

```yaml
card_db: "/path/to/your/databases/card"
```

---

### DeepARG (exploratory/deep-learning AMR)

> bioconda's `deeparg=1.0.4` is the classic Python2/Theano codebase, not the
> newer PyTorch/HuggingFace rewrite some docs describe — confirmed live on
> litrp4 (`--threads` doesn't exist on `deeparg predict`, and `-d/--data-path`
> is a required argument, not auto-managed). Data must be fetched once via
> `deeparg download_data` into a directory set by `deeparg_db` in `config.yaml`.

**Conda:**
```bash
mkdir -p "$DB_BASE/deeparg"
conda activate env_deeparg
deeparg download_data -o "$DB_BASE/deeparg"
conda deactivate
```

**Docker:**
```bash
mkdir -p "$DB_BASE/deeparg"
docker run --rm -v "$DB_BASE/deeparg:/dbs" \
    quay.io/biocontainers/deeparg:1.0.4--pyhdfd78af_0 \
    deeparg download_data -o /dbs
```

Then set in `config.yaml`:
```yaml
deeparg_db: "/path/to/your/databases/deeparg"
```

---

## 7. Reads-only classification databases (reads_classify module, optional)

This module (BLOCK 15) uses **Sylph** for direct taxonomic profiling from raw reads.
It runs independently of the assembly pipeline and is disabled by default (`reads_classify: false` in `config.yaml`).

Up to four databases can be enabled simultaneously; configure each in `reads_classify_dbs` in `config.yaml`.

> **Pre-built databases**: Sylph maintains official pre-built `.syldb` files for the most common databases.
> For the current list and download links, see:
> **https://sylph-docs.github.io/pre%E2%80%90built-databases/**

---

### UHGV (Unified Human Gut Virome, ~0.4 GB)

Pre-built for gut virome studies — 171 k vOTUs at 100% ANI resolution.

```bash
mkdir -p "$DB_BASE/sylph"
wget -O "$DB_BASE/sylph/uhgv_c100_dbv1.syldb" \
    "https://zenodo.org/records/14884392/files/uhgv_c100_dbv1.syldb"
```

Set in `config.yaml`:
```yaml
reads_classify_dbs:
  uhgv: "/path/to/databases/sylph/uhgv_c100_dbv1.syldb"
```

---

### IMG/VR (pre-built v4.1, ~2 GB)

2.9M viral genomes at 200-mer resolution. Note: your server may have IMG/VR 7.1 raw files — the pre-built DB covers v4.1. To use a local IMG/VR export instead, see [Custom syldb from own FASTA](#custom-syldb-from-own-fasta) below.

```bash
wget -O "$DB_BASE/sylph/imgvr_c200_v0.3.0.syldb" \
    "https://zenodo.org/records/14884392/files/imgvr_c200_v0.3.0.syldb"
```

Set in `config.yaml`:
```yaml
reads_classify_dbs:
  imgvr: "/path/to/databases/sylph/imgvr_c200_v0.3.0.syldb"
```

---

### GTDB r232 (prokaryotes)

All GTDB r232 representative genomes, for prokaryotic profiling alongside viral DBs.
Two compression levels are available — choose based on your RAM/disk budget:

| Variant | Disk | RAM usage | Sensitivity |
|---|---|---|---|
| `c1000` (recommended) | ~5 GB | lower | standard |
| `c200` (high sensitivity) | ~24 GB | higher | more hits at low ANI |

```bash
# Standard (c1000, ~5 GB) — recommended for most uses
wget -O "$DB_BASE/sylph/gtdb-r232-c1000-dbv1.syldb" \
    "https://zenodo.org/records/14884392/files/gtdb-r232-c1000-dbv1.syldb"

# High-sensitivity (c200, ~24 GB) — for low-abundance or divergent taxa
wget -O "$DB_BASE/sylph/gtdb-r232-c200-dbv1.syldb" \
    "https://faust.compbio.cs.cmu.edu/sylph-stuff/gtdb-r232-c200-dbv1.syldb"
```

Set in `config.yaml`:
```yaml
reads_classify_dbs:
  gtdb: "/path/to/databases/sylph/gtdb-r232-c1000-dbv1.syldb"   # or c200 variant
  gtdb_version: "r232"
```

---

### Custom syldb from own FASTA

You can build a Sylph database from any collection of genome sequences (FASTA).
For example, to use IMG/VR 7.1 nucleotide sequences you have locally:

**Step 1 — Build the `.syldb`:**

```bash
# Requires sylph ≥0.6 (conda activate env_reads_classify)
conda activate env_reads_classify
sylph sketch -c 100 -i /path/to/IMGVR_v7.1_genomes/*.fna \
    -o "$DB_BASE/sylph/imgvr_v71" -t 32
# Output: $DB_BASE/sylph/imgvr_v71.syldb
conda deactivate
```

> `-c 100` sets the compression factor (lower = more sensitive, larger file). Use `-c 200` for very large collections.

**Step 2 — Build a taxonomy TSV for sylph-tax:**

For IMG/VR exports that include a metadata TSV (`_meta.tsv`), use the helper script:

```bash
python3 scripts/reads_classify/build_imgvr_taxonomy.py \
    IMGVR_v7.1_meta.tsv \
    "$DB_BASE/sylph/imgvr_v71_taxonomy.tsv"
```

For any other source, build a two-column TSV (no header) where each row is:
```
genome_id    d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species
```
Empty ranks are written as `p__`, `c__`, etc.

**Step 3 — Set paths in `config.yaml`:**
```yaml
reads_classify_dbs:
  custom: "/path/to/databases/sylph/imgvr_v71.syldb"
  custom_taxonomy: "/path/to/databases/sylph/imgvr_v71_taxonomy.tsv"
```

---

### Taxonomy files (one-time download, all databases)

The first pipeline run automatically downloads sylph-tax taxonomy annotation files (~50 MB)
via `rule sylph_tax_download`. To pre-fetch manually:

```bash
mkdir -p "$DB_BASE/sylph/taxonomy"
conda activate env_reads_classify
export SYLPH_TAXONOMY_CONFIG="$DB_BASE/sylph/taxonomy/config.json"
sylph-tax download --download-to "$DB_BASE/sylph/taxonomy"
conda deactivate
```

Then set in `config.yaml`:
```yaml
reads_classify_tax_dir: "/path/to/databases/sylph/taxonomy"
```

---

## 8. Optional: custom MMseqs2 seqTaxDBs

Custom databases allow classifying contigs and bins not covered by primary
databases. Both the prokaryote and viral custom paths run exclusively
through MMseqs2 `taxonomy` (real per-genome LCA) now -- there is no Diamond
best-hit/majority-vote option for either (`diamond_custom_prok` and
`diamond_custom_viral` were both removed entirely, the latter 2026-06-23).
MMseqs2 avoids "spurious specificity" (von Meijenfeldt et al. 2019, CAT/BAT)
on the divergent/environmental genomes these custom DBs exist to cover in
the first place.

Same script for both (`scripts/prepare_mmseqs_taxdb.py`), different
`--format`: `img`/`ncbi` for prokaryotes, `inphared`/`imgvr` for viral
(plus `ncbi` works for either). Each format is its own header-parsing +
lineage-loading pair (see the script's own docstring); a source database
that doesn't match one of these needs its own small format added there,
not a new separate script.

### Prokaryotes (e.g. IMG NR)

`rule mmseqs_taxonomy_prok` is the only source for custom prokaryote
taxonomy. The report aggregates per-protein LCA to genome level with a
second LCA pass (`load_mmseqs_taxonomy_prok` in
`scripts/report/data_loaders.py`), not a vote, for the same reason.

```bash
conda activate env_assembly
python3 scripts/prepare_mmseqs_taxdb.py \
    --faa img_unrestricted_isolates_nr.faa --format img \
    --img-tax taxonOId2Taxonomy.tsv \
    --out "$DB_BASE/img/img_nr_mmseqs" --threads 32
conda deactivate
# Output: $DB_BASE/img/img_nr_mmseqs/seqTaxDB -> set as custom_prok_mmseqs_db
```

> `mmseqs createdb` on the full IMG NR protein set (hundreds of millions of
> sequences) can take several hours — the script skips it on a re-run if
> `seqTaxDB.dbtype` already exists, so it's safe to re-run after an
> interrupted `createtaxdb` step.

### Viral (e.g. IMG/VR)

`rule mmseqs_taxonomy_custom_viral` is the only source for custom viral
taxonomy beyond INPHARED (`mmseqs_taxonomy_viral`, always-on) and GeNomad
(fallback). Both run alongside `mmseqs_taxonomy_viral` and are compared by
resolved rank depth in `viral_taxonomy`, not a fixed priority.

IMG/VR's high-confidence export ships
`IMGVR_all_proteins-high_confidence.faa` (protein headers
`>{UVIG}|{Taxon_oid}|{rest}`) and
`IMGVR_all_Sequence_information-high_confidence.tsv` (`UVIG` join key +
`Taxonomic classification` column, geNomad-generated rank-prefixed lineage
e.g. `r__Realm;k__Kingdom;p__Phylum;...`) -- confirmed against a real
IMG/VR 2022-12-19 export.

```bash
conda activate env_assembly
python3 scripts/prepare_mmseqs_taxdb.py \
    --faa IMGVR_all_proteins-high_confidence.faa --format imgvr \
    --imgvr-tax IMGVR_all_Sequence_information-high_confidence.tsv \
    --out "$DB_BASE/img/IMG_VR/imgvr_mmseqs" --threads 32
conda deactivate
# Output: $DB_BASE/img/IMG_VR/imgvr_mmseqs/seqTaxDB -> set as custom_viral_mmseqs_db
```

---

## 9. Configure config.yaml

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
checkm2_db:    "/path/to/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
gtdbtk_db:     "/path/to/gtdbtk/release226"
inphared_db:   "/path/to/inphared"
pharokka_db:   "/path/to/pharokka"
phold_db:      "/path/to/phold_db"
bakta_db:      "/path/to/bakta/db"
eggnog_db:     "/path/to/eggnog"

# Defense systems + AMR (AMRFinderPlus self-manages its own small DB —
# CARD, DeepARG and DefenseFinder's models need an explicit path each)
defense_amr_enabled:         true
card_db:                   "/path/to/card"
deeparg_db:                "/path/to/deeparg"
defense_finder_models_db:  "/path/to/defense_finder_models"
abricate_enabled:            true   # VFDB + PlasmidFinder only, self-contained
argnorm_enabled:              true   # AMRFinderPlus + DeepARG -> ARO, self-contained

# Viral-side defense/anti-defense (Han et al. 2026 cold seep paper) —
# DefenseFinder reused on viral ORFs (no extra DB) + dbAPIS (small, own DB)
defense_amr_viral_enabled:   true
apis_db:                     "/path/to/dbapis"

# Custom MMseqs2 seqTaxDBs — leave "" to skip (no Diamond option for either)
custom_prok_mmseqs_db:  ""   # e.g. IMG NR
custom_viral_mmseqs_db: ""   # e.g. IMG/VR

# Reads-only classification (BLOCK 15 — disabled by default)
reads_classify: false
reads_classify_dbs:
  imgvr:    ""   # path to imgvr_c200_v0.3.0.syldb or custom IMG/VR syldb
  uhgv:     ""   # path to uhgv_c100_dbv1.syldb
  gtdb:     ""   # path to gtdb-r232-c1000-dbv1.syldb
  gtdb_version: "r232"
  custom:   ""   # any user-built .syldb
  custom_taxonomy: ""   # two-column TSV (genome_id TAB lineage)
reads_classify_tax_dir: ""          # leave "" to auto-populate under outdir
reads_classify_min_prevalence: 0.0  # 0 = no filter; 0.1 = ≥10% of samples
reads_classify_min_kmers: 20        # minimum k-mer matches per genome
reads_classify_genome_fasta: ""     # optional: reference FASTA for BACPHLIP
reads_classify_virulence_threshold: 0.5
```

---

## 10. Verify installation

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
           env_mapping env_viral env_genomad env_vrhyme \
           env_cobra env_binning env_binette env_checkm2 \
           env_gtdbtk env_phist env_annotation env_coverm \
           env_gunc env_derep env_defense env_rgi env_deeparg \
           env_reads_classify; do
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
| INPHARED | 2 GB |
| CheckM2 | 3 GB |
| GTDB-Tk release 226 | 85 GB |
| Pharokka (PHROGS) | 1 GB |
| Phold | 10 GB |
| Bakta (full) | 30 GB |
| EggNOG-mapper | 50 GB |
| GUNC progenomes_2.1 | 13 GB |
| DefenseFinder models (auto) | <1 GB |
| AMRFinderPlus DB (auto) | <1 GB |
| CARD (RGI) | ~1 GB |
| DeepARG model + DB (auto) | ~2 GB |
| IMG viral (optional) | 5–20 GB |
| IMG prokaryote (optional) | 50–200 GB |
| Sylph UHGV (optional, reads_classify) | 0.4 GB |
| Sylph IMG/VR 4.1 (optional, reads_classify) | 2 GB |
| Sylph GTDB r232 c1000 (optional, reads_classify) | 5 GB |
| Sylph GTDB r232 c200 (optional, reads_classify) | 24 GB |
| Sylph taxonomy files (optional, reads_classify) | <0.1 GB |
| **Total (without custom)** | **~242 GB** |

---

## Troubleshooting

**metaSPAdes runs out of memory**
Reduce `spades_mem` in `config.yaml`. If RAM is the bottleneck, use MEGAHIT only.

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
