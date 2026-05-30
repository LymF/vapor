# Guia de Instalação — VAPOR

Este guia cobre a instalação de todas as dependências e bancos de dados necessários.
Siga os passos em ordem — o download dos bancos é a etapa mais demorada (~400 GB no total).

---

## Índice

1. [Requisitos de sistema](#1-requisitos-de-sistema)
2. [Instalar Miniforge3](#2-instalar-miniforge3)
3. [Instalar Snakemake](#3-instalar-snakemake)
4. [Criar ambientes conda](#4-criar-ambientes-conda)
5. [Instalar bancos de dados](#5-instalar-bancos-de-dados)
6. [Bancos Diamond customizados (opcional)](#6-bancos-diamond-customizados-opcional)
7. [Configurar config.yaml](#7-configurar-configyaml)
8. [Verificar instalação](#8-verificar-instalação)

---

## 1. Requisitos de sistema

| Componente | Mínimo | Recomendado |
|---|---|---|
| OS | Linux (Ubuntu 20.04+) | Ubuntu 22.04 LTS |
| CPU | 16 cores | 32–64 cores |
| RAM | 64 GB | 128–256 GB |
| Disco (bancos) | 400 GB | 600 GB |
| Disco (resultados) | 200 GB / projeto | SSD preferencial |

> **Nota de memória:** metaSPAdes usa 80–150 GB RAM para metagenomas grandes. vConTACT3 pode usar 30–60 GB. Planeje de acordo.

---

## 2. Instalar Miniforge3

```bash
curl -L -O https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh
bash Miniforge3-Linux-x86_64.sh -b -p "$HOME/miniforge3"
source "$HOME/miniforge3/etc/profile.d/conda.sh"
conda init bash
source ~/.bashrc
```

---

## 3. Instalar Snakemake

```bash
mamba create -n snakemake -c conda-forge -c bioconda \
    "snakemake>=8.0" python=3.11 -y
conda activate snakemake
snakemake --version   # deve imprimir 8.x.x
```

---

## 4. Criar ambientes conda

Os ambientes são definidos nos YAMLs em `envs/`. O Snakemake pode criá-los automaticamente:

```bash
conda activate snakemake
snakemake --use-conda --cores 1 --create-envs-only
```

Ou criar manualmente (útil para testar individualmente):

```bash
# QC (short reads: fastp, quast, multiqc)
mamba create -n env_qc -c conda-forge -c bioconda \
    fastp quast multiqc -y

# Assembly (MEGAHIT + SPAdes + metaMDBG + MMseqs2)
mamba create -n env_assembly -c conda-forge -c bioconda \
    megahit spades metamdbg mmseqs2 -y

# Long reads: Flye + hifiasm (merged)
mamba create -n env_flye -c conda-forge -c bioconda \
    flye hifiasm hifiasm_meta -y

# Medaka — polimento ONT (GPU via CUDA)
mamba create -n env_medaka -c conda-forge -c bioconda \
    medaka -y

# Long-read QC
mamba create -n env_lr_utils -c conda-forge -c bioconda \
    nanoplot filtlong porechop_abi -y

# Mapeamento
mamba create -n env_mapping -c conda-forge -c bioconda \
    bwa-mem2 minimap2 samtools -y

# Detecção viral + taxonomia
mamba create -n env_viral -c conda-forge -c bioconda \
    "virsorter=2" diamond prodigal checkv -y

# GeNomad
mamba create -n env_genomad -c conda-forge -c bioconda \
    genomad -y

# VIBRANT
mamba create -n phage_vibrant -c conda-forge -c bioconda \
    "python=3.9" "scikit-learn=0.23" numpy pandas matplotlib \
    hmmer prodigal biopython -y
conda activate phage_vibrant
pip install vibrant
conda deactivate

# vRhyme
mamba create -n env_vrhyme -c conda-forge -c bioconda \
    vrhyme -y

# Binning prokarioto
mamba create -n env_binning -c conda-forge -c bioconda \
    metabat2 vamb semibin2 -y

# COMEBin (transformer binner, rank 1 em 2025)
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

# Anotação (Pharokka + Phold + Bakta + EggNOG + mapas circulares)
mamba create -n env_annotation -c conda-forge -c bioconda \
    pharokka phold bakta eggnog-mapper pycirclize matplotlib biopython -y

# vConTACT3
mamba create -n env_vcontact3 -c conda-forge -c bioconda \
    vcontact3 -y

# CoverM + diversidade (alpha, beta, Procrustes)
mamba create -n env_coverm -c conda-forge -c bioconda \
    coverm numpy scipy -y
```

---

## 5. Instalar bancos de dados

Defina um diretório base para todos os bancos:

```bash
DB_BASE="/path/to/your/databases"
mkdir -p "$DB_BASE"
```

### CheckV

```bash
conda activate env_checkm2
checkv download_database "$DB_BASE/checkv"
# Resultado: $DB_BASE/checkv/checkv-db-v1.5/
conda deactivate
```

### VirSorter2

```bash
conda activate env_viral
virsorter setup -d "$DB_BASE/virsorter2" -j 4
conda deactivate
```

### GeNomad

```bash
conda activate env_genomad
genomad download-database "$DB_BASE/genomad"
# Resultado: $DB_BASE/genomad/genomad_db/
conda deactivate
```

### VIBRANT

```bash
mkdir -p "$DB_BASE/vibrant-1.0.1/databases"
mkdir -p "$DB_BASE/vibrant-1.0.1/files"
conda activate phage_vibrant

cd "$DB_BASE/vibrant-1.0.1/databases"
wget https://zenodo.org/record/4543735/files/VIBRANT_v1.0.1.tar.gz
tar -xzf VIBRANT_v1.0.1.tar.gz --strip-components=2 "VIBRANT_v1.0.1/databases/"

cd "$DB_BASE/vibrant-1.0.1/files"
git clone --depth 1 https://github.com/AnantharamanLab/VIBRANT /tmp/vibrant_src
cp /tmp/vibrant_src/files/* "$DB_BASE/vibrant-1.0.1/files/"
rm -rf /tmp/vibrant_src

for hmm in "$DB_BASE/vibrant-1.0.1/databases/"*.HMM; do
    hmmpress "$hmm"
done
conda deactivate
```

### INPHARED

```bash
mkdir -p "$DB_BASE/inphared"
cd "$DB_BASE/inphared"

# Verificar a data mais recente em: https://github.com/RyanCook94/inphared
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

### vConTACT3

```bash
mkdir -p "$DB_BASE/vcontact3"
conda activate env_vcontact3
vcontact3 prepare_database \
    --output "$DB_BASE/vcontact3" \
    --threads 32
conda deactivate
# Anotar a versão impressa — configurar VCONTACT3_VER no config.yaml
```

### CheckM2

```bash
conda activate env_checkm2
checkm2 database --download --path "$DB_BASE/checkm2"
# Resultado: $DB_BASE/checkm2/CheckM2_database/uniref100.KO.1.dmnd
conda deactivate
```

### GTDB-Tk

```bash
mkdir -p "$DB_BASE/gtdbtk"
conda activate env_gtdbtk
download-db.sh "$DB_BASE/gtdbtk"
conda deactivate
```

### Pharokka

```bash
conda activate env_annotation
pharokka_install_databases.py -o "$DB_BASE/pharokka"
conda deactivate
```

### Phold

```bash
conda activate env_annotation
phold install -d "$DB_BASE/phold_db"
conda deactivate
```

### Bakta

```bash
conda activate env_annotation
bakta_db download --output "$DB_BASE/bakta" --type full
conda deactivate
```

### EggNOG-mapper

```bash
conda activate env_annotation
export EGGNOG_DATA_DIR="$DB_BASE/eggnog"
download_eggnog_data.py -y -f --data_dir "$DB_BASE/eggnog"
conda deactivate
```

---

## 6. Bancos Diamond customizados (opcional)

Permitem classificar contigs e bins não cobertos pelos bancos primários.

```bash
# IMG NR — subconjunto viral
conda activate env_viral
python3 scripts/prepare_diamond_db.py \
    --faa img_nr.faa --format img \
    --img-tax taxonOId2Taxonomy.tsv \
    --filter-domain Viruses \
    --out "$DB_BASE/img/img_viral" --threads 32
# Produz: img_viral.dmnd + img_viral_meta.tsv

# IMG NR — subconjunto prokarioto
python3 scripts/prepare_diamond_db.py \
    --faa img_nr.faa --format img \
    --img-tax taxonOId2Taxonomy.tsv \
    --filter-domain Bacteria,Archaea \
    --out "$DB_BASE/img/img_prok" --threads 32
conda deactivate
```

---

## 7. Configurar config.yaml

```yaml
# Entrada/saída
fastq_dir: "fastqs"
outdir:    "results"
threads:   32

# Tipo de dado
long_reads: false
lr_tech:    "ont"   # "ont" | "hifi"

# Bancos de dados — substituir pelos seus caminhos
checkv_db:    "/path/to/checkv/checkv-db-v1.5"
vs2_db:       "/path/to/virsorter2"
genomad_db:   "/path/to/genomad/genomad_db"
vibrant_base: "/path/to/vibrant-1.0.1"
checkm2_db:   "/path/to/checkm2/CheckM2_database/uniref100.KO.1.dmnd"
gtdbtk_db:    "/path/to/gtdbtk/release226"
inphared_db:  "/path/to/inphared"
vcontact3_db: "/path/to/vcontact3"
vcontact3_ver: "230"
pharokka_db:  "/path/to/pharokka"
phold_db:     "/path/to/phold_db"
bakta_db:     "/path/to/bakta/db"
eggnog_db:    "/path/to/eggnog"

# Bancos customizados (opcional — deixar "" para pular)
custom_viral_dmnd: ""
custom_viral_meta: ""
custom_prok_dmnd:  ""
custom_prok_meta:  ""
```

---

## 8. Verificar instalação

```bash
conda activate snakemake

# Criar FASTQs de teste
mkdir -p fastqs
touch fastqs/TEST_R1.fastq.gz fastqs/TEST_R2.fastq.gz

# Dry-run
snakemake --use-conda --cores 4 --dry-run 2>&1 | tail -20

# Limpar arquivos de teste
rm fastqs/TEST_R1.fastq.gz fastqs/TEST_R2.fastq.gz
```

### Verificar ambientes individualmente

```bash
for env in env_qc env_assembly env_flye env_medaka env_lr_utils \
           env_mapping env_viral env_genomad phage_vibrant env_vrhyme \
           env_binning env_comebin env_binette env_checkm2 env_gtdbtk \
           env_phist env_annotation env_vcontact3 env_coverm; do
    echo -n "$env: "
    conda run -n "$env" python --version 2>/dev/null || echo "MISSING"
done
```

---

## Espaço em disco

| Banco | Tamanho aproximado |
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
| IMG viral (opcional) | 5–20 GB |
| IMG prokarioto (opcional) | 50–200 GB |
| **Total (sem customizados)** | **~218 GB** |

---

## Troubleshooting

**metaSPAdes sem memória**
Reduza `spades_mem` no `config.yaml`. Se RAM for o gargalo, use apenas MEGAHIT.

**VIBRANT: erro `-f older`**
Este pipeline usa `cd` para o diretório de saída antes de rodar VIBRANT, passando apenas `-f nucl`. Confirme que está usando a versão atual do `Snakefile`.

**vConTACT3 muito lento**
Pode levar horas para datasets virais grandes. Rode com `--cores 32` e confirme que `vcontact3_ver` no config.yaml bate com sua versão do banco.

**GTDB-Tk: erro pplacer**
Confirme que `GTDBTK_DATA_PATH` aponta para o diretório correto e que a versão do banco é compatível com o gtdbtk instalado (`gtdbtk check_install`).
