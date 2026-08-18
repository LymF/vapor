#!/usr/bin/env bash
# install_kraken_phanta.sh — instala Kraken2 + Bracken + baixa bancos
# Uso: bash install_kraken_phanta.sh [--db-dir /path/to/dbs] [--threads 16]
# Requer: conda/mamba, ~300GB disco, ~120GB RAM para banco Standard

set -euo pipefail

# ── Defaults ──────────────────────────────────────────────────────────────────
DB_ROOT="${HOME}/databases/kraken2"
THREADS=16
SKIP_STANDARD=0
SKIP_PHANTA=0

# ── Args ──────────────────────────────────────────────────────────────────────
while [[ $# -gt 0 ]]; do
    case "$1" in
        --db-dir)    DB_ROOT="$2"; shift 2 ;;
        --threads)   THREADS="$2"; shift 2 ;;
        --skip-standard) SKIP_STANDARD=1; shift ;;
        --skip-phanta)   SKIP_PHANTA=1; shift ;;
        *) echo "Argumento desconhecido: $1"; exit 1 ;;
    esac
done

STANDARD_DB="${DB_ROOT}/standard"
PHANTA_DB="${DB_ROOT}/phanta"

echo "============================================================"
echo " Kraken2 + Bracken + Phanta — instalação"
echo " DB root : ${DB_ROOT}"
echo " Threads : ${THREADS}"
echo "============================================================"

# ── 1. Conda env: kraken2 + bracken ──────────────────────────────────────────
echo ""
echo "[1/4] Criando ambiente conda env_kraken ..."

if conda env list | grep -q "^env_kraken "; then
    echo "      env_kraken já existe — pulando criação."
else
    conda create -y -n env_kraken -c conda-forge -c bioconda \
        kraken2 bracken krakentools \
        python=3.11 pandas numpy
    echo "      env_kraken criado."
fi

# Verifica instalação
KRAKEN2=$(conda run -n env_kraken which kraken2)
BRACKEN=$(conda run -n env_kraken which bracken)
echo "      kraken2: ${KRAKEN2}"
echo "      bracken: ${BRACKEN}"

# ── 2. Banco Standard (bacteria + archaea + virus + human) ───────────────────
if [[ $SKIP_STANDARD -eq 0 ]]; then
    echo ""
    echo "[2/4] Baixando banco Kraken2 Standard (~65GB comprimido, ~70GB descomprimido) ..."
    mkdir -p "${STANDARD_DB}"

    # Pre-built index — mais rápido que build from scratch (que leva dias)
    # URL mantida em: https://benlangmead.github.io/aws-indexes/k2
    K2_URL="https://genome-idx.s3.amazonaws.com/kraken/k2_standard_20240605.tar.gz"
    K2_TAR="${DB_ROOT}/k2_standard.tar.gz"

    if [[ -f "${STANDARD_DB}/hash.k2d" ]]; then
        echo "      Banco Standard já existe — pulando download."
    else
        echo "      Baixando: ${K2_URL}"
        echo "      Destino : ${K2_TAR}"
        wget -c --show-progress -O "${K2_TAR}" "${K2_URL}"
        echo "      Extraindo ..."
        tar -xzf "${K2_TAR}" -C "${STANDARD_DB}"
        rm -f "${K2_TAR}"
        echo "      Banco Standard extraído em ${STANDARD_DB}"
    fi

    # Build Bracken para o Standard (necessário para Bracken funcionar)
    BRACKEN_DB_S="${STANDARD_DB}/database150mers.kmer_distrib"
    if [[ -f "${BRACKEN_DB_S}" ]]; then
        echo "      Bracken files já existem — pulando build."
    else
        echo "      Construindo Bracken (read length 150bp, pode levar ~1-2h) ..."
        conda run -n env_kraken bracken-build \
            -d "${STANDARD_DB}" \
            -t "${THREADS}" \
            -l 150
        echo "      Bracken Standard pronto."
    fi
else
    echo "[2/4] Banco Standard — PULADO (--skip-standard)"
fi

# ── 3. Banco Phanta (viral-focused: NCBI viral + INPHARED) ──────────────────
if [[ $SKIP_PHANTA -eq 0 ]]; then
    echo ""
    echo "[3/4] Baixando banco Phanta (~20GB comprimido) ..."
    mkdir -p "${PHANTA_DB}"

    # Phanta pre-built database (Zenodo)
    # Repositório: https://github.com/bhatt-lab/phanta
    # Banco:       https://zenodo.org/record/7661573
    PHANTA_URL="https://zenodo.org/record/7661573/files/phanta_db.tar.gz"
    PHANTA_TAR="${DB_ROOT}/phanta_db.tar.gz"

    if [[ -f "${PHANTA_DB}/hash.k2d" ]]; then
        echo "      Banco Phanta já existe — pulando download."
    else
        echo "      Baixando: ${PHANTA_URL}"
        wget -c --show-progress -O "${PHANTA_TAR}" "${PHANTA_URL}"
        echo "      Extraindo ..."
        tar -xzf "${PHANTA_TAR}" -C "${PHANTA_DB}" --strip-components=1
        rm -f "${PHANTA_TAR}"
        echo "      Banco Phanta extraído em ${PHANTA_DB}"
    fi

    # Build Bracken para Phanta
    BRACKEN_DB_P="${PHANTA_DB}/database150mers.kmer_distrib"
    if [[ -f "${BRACKEN_DB_P}" ]]; then
        echo "      Bracken Phanta files já existem — pulando build."
    else
        echo "      Construindo Bracken para Phanta (~30min) ..."
        conda run -n env_kraken bracken-build \
            -d "${PHANTA_DB}" \
            -t "${THREADS}" \
            -l 150
        echo "      Bracken Phanta pronto."
    fi
else
    echo "[3/4] Banco Phanta — PULADO (--skip-phanta)"
fi

# ── 4. Clone scripts do Phanta (normalização) ────────────────────────────────
echo ""
echo "[4/4] Clonando scripts do Phanta ..."
PHANTA_SCRIPTS="${DB_ROOT}/phanta_scripts"

if [[ -d "${PHANTA_SCRIPTS}/.git" ]]; then
    echo "      Já existe — atualizando."
    git -C "${PHANTA_SCRIPTS}" pull --ff-only
else
    git clone https://github.com/bhatt-lab/phanta.git "${PHANTA_SCRIPTS}"
fi

# ── Resumo ───────────────────────────────────────────────────────────────────
echo ""
echo "============================================================"
echo " Instalação concluída!"
echo ""
echo " Bancos:"
echo "   Standard : ${STANDARD_DB}"
echo "   Phanta   : ${PHANTA_DB}"
echo "   Scripts  : ${PHANTA_SCRIPTS}"
echo ""
echo " Próximo passo:"
echo "   bash scripts/run_kraken_classify.sh \\"
echo "     --fastq-dir /path/to/fastqs \\"
echo "     --outdir    /path/to/output \\"
echo "     --standard-db ${STANDARD_DB} \\"
echo "     --phanta-db   ${PHANTA_DB} \\"
echo "     --phanta-scripts ${PHANTA_SCRIPTS} \\"
echo "     --threads   ${THREADS}"
echo "============================================================"
