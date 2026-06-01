#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────
# run_cobra_standalone.sh — test COBRA-meta outside the Snakemake pipeline
#
# Uses VAPOR artifacts for a sample. IMPORTANT: VAPOR maps reads against
# rep_seq.fasta (MMseqs2 deduplicated), NOT against contigs.fasta.
# COBRA-meta needs coverage for every contig in the full metaSPAdes
# contigs.fasta — so this script remaps reads against contigs.fasta first.
#
# Required VAPOR outputs:
#   - viral/consensus/{sample}_viral_consensus.fasta   (SPADES_-prefixed)
#   - assembly/metaspades/contigs.fasta
#   - trimmed/{sample}_R1_fastp.fq.gz + _R2_fastp.fq.gz
#
# Output:
#   {results-dir}/{sample}/viral/cobra_standalone/cobra_out/
#
# Usage:
#   bash run_cobra_standalone.sh --sample SAMPLE [--results-dir DIR] [--threads 16]
#   bash run_cobra_standalone.sh --sample SAMPLE --use-conda
# ──────────────────────────────────────────────────────────────────────
set -euo pipefail

SAMPLE=""
RESULTS_DIR="results-teste"
THREADS=16
USE_CONDA=0
MINK=21
MAXK=127

IMG_BWA="quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0"
IMG_SAMTOOLS="quay.io/biocontainers/samtools:1.23--h96c455f_0"
IMG_METABAT="quay.io/biocontainers/metabat2:2.18_23_gc869c52--h61f4f8f_0"
IMG_COBRA="quay.io/biocontainers/cobra-meta:1.2.3--pyhdfd78af_0"

while [ $# -gt 0 ]; do
    case "$1" in
        --sample)       SAMPLE="$2"; shift 2 ;;
        --results-dir)  RESULTS_DIR="$2"; shift 2 ;;
        --threads)      THREADS="$2"; shift 2 ;;
        --mink)         MINK="$2"; shift 2 ;;
        --maxk)         MAXK="$2"; shift 2 ;;
        --use-conda)    USE_CONDA=1; shift ;;
        -h|--help)
            sed -n '2,22p' "$0"; exit 0 ;;
        *) echo "Unknown argument: $1" >&2; exit 2 ;;
    esac
done

[ -n "$SAMPLE" ] || { echo "ERROR: --sample required" >&2; exit 2; }

RESULTS_DIR="$(realpath "$RESULTS_DIR")"
SAMPLE_DIR="$RESULTS_DIR/$SAMPLE"
VIRAL="$SAMPLE_DIR/viral/consensus/${SAMPLE}_viral_consensus.fasta"
MSPADES="$SAMPLE_DIR/assembly/metaspades/contigs.fasta"
R1="$SAMPLE_DIR/trimmed/${SAMPLE}_R1_fastp.fq.gz"
R2="$SAMPLE_DIR/trimmed/${SAMPLE}_R2_fastp.fq.gz"
WORK="$SAMPLE_DIR/viral/cobra_standalone"
COBRA_BAM="$WORK/${SAMPLE}_cobra.sorted.bam"
COBRA_DEPTH="$WORK/${SAMPLE}_cobra_depth.txt"

for f in "$VIRAL" "$MSPADES" "$R1" "$R2"; do
    [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

mkdir -p "$WORK"

# ── helpers ────────────────────────────────────────────────────────────
run_bwa_index() {
    if [ "$USE_CONDA" -eq 1 ]; then
        bwa-mem2 index "$@"
    else
        docker run --rm \
            -v "$RESULTS_DIR:$RESULTS_DIR" \
            "$IMG_BWA" bwa-mem2 index "$@"
    fi
}

run_bwa_and_sort() {
    local threads=$1 idx=$2 r1=$3 r2=$4 bam=$5 log=$6
    if [ "$USE_CONDA" -eq 1 ]; then
        bwa-mem2 mem -t "$threads" "$idx" "$r1" "$r2" 2>"$log" \
            | samtools sort -@ "$threads" -o "$bam" -
    else
        docker run --rm \
            -v "$RESULTS_DIR:$RESULTS_DIR" \
            "$IMG_BWA" \
            bwa-mem2 mem -t "$threads" "$idx" "$r1" "$r2" \
            2>"$log" \
        | docker run --rm -i \
            -v "$RESULTS_DIR:$RESULTS_DIR" \
            "$IMG_SAMTOOLS" \
            samtools sort -@ "$threads" -o "$bam" -
    fi
}

run_samtools_index() {
    if [ "$USE_CONDA" -eq 1 ]; then
        samtools index "$@"
    else
        docker run --rm \
            -v "$RESULTS_DIR:$RESULTS_DIR" \
            "$IMG_SAMTOOLS" samtools index "$@"
    fi
}

run_jgi() {
    if [ "$USE_CONDA" -eq 1 ]; then
        jgi_summarize_bam_contig_depths "$@"
    else
        docker run --rm \
            -v "$RESULTS_DIR:$RESULTS_DIR" \
            "$IMG_METABAT" jgi_summarize_bam_contig_depths "$@"
    fi
}

run_cobra() {
    if [ "$USE_CONDA" -eq 1 ]; then
        eval "$(conda shell.bash hook)"
        conda activate env_cobra
        cobra-meta "$@"
    else
        docker run --rm \
            -v "$RESULTS_DIR:$RESULTS_DIR" \
            "$IMG_COBRA" cobra-meta "$@"
    fi
}

# ── Step 1: remap reads against contigs.fasta ──────────────────────────
# VAPOR's depth file is from mapping against rep_seq.fasta (deduplicated).
# COBRA needs coverage for ALL contigs in contigs.fasta, including those
# merged away by MMseqs2. We remap once and cache the result.
IDX="$WORK/contigs_index.bwt.2bit.64"
IDX_PREFIX="$WORK/contigs_index"

if [ ! -s "$COBRA_BAM" ]; then
    echo "[1/4] Indexing metaSPAdes contigs.fasta for BWA-MEM2..."
    if [ ! -s "$IDX" ]; then
        run_bwa_index -p "$IDX_PREFIX" "$MSPADES" \
            > "$WORK/bwa_index.log" 2>&1
    fi

    echo "[1/4] Mapping reads against contigs.fasta (COBRA-specific BAM)..."
    run_bwa_and_sort \
        "$THREADS" "$IDX_PREFIX" "$R1" "$R2" \
        "$COBRA_BAM" "$WORK/bwa_mem.log"
    run_samtools_index "$COBRA_BAM"
    echo "      OK"
else
    echo "[1/4] COBRA BAM exists: $COBRA_BAM"
fi

# ── Step 2: per-contig depth from COBRA BAM ────────────────────────────
if [ ! -s "$COBRA_DEPTH" ]; then
    echo "[2/4] Computing per-contig depth (contigs.fasta)..."
    run_jgi --outputDepth "$COBRA_DEPTH" "$COBRA_BAM" \
        > "$WORK/jgi_depth.log" 2>&1
    echo "      OK"
else
    echo "[2/4] Depth file exists: $COBRA_DEPTH"
fi

# ── Step 3: build query.fa and coverage.tsv ────────────────────────────
echo "[3/4] Building query.fa and coverage.tsv..."

# SPADES_-prefixed viral contigs → strip prefix
awk 'BEGIN{w=0} /^>/{w=0; if($0 ~ /^>SPADES_/){sub("^>SPADES_",">"); w=1; print} next} w{print}' \
    "$VIRAL" > "$WORK/query.fa"

N_QUERY=$(grep -c '^>' "$WORK/query.fa" || echo 0)
echo "      SPADES viral query contigs: $N_QUERY"
if [ "$N_QUERY" -eq 0 ]; then
    echo "      No SPADES viral contigs — nothing to extend."
    exit 0
fi

# Full coverage from contigs.fasta mapping (col 1 = contig, col 3 = avgDepth)
awk 'BEGIN{FS=OFS="\t"} NR==1{next} {print $1, $3}' \
    "$COBRA_DEPTH" > "$WORK/coverage.tsv"

# ── Step 4: run COBRA-meta ─────────────────────────────────────────────
OUT="$WORK/cobra_out"
rm -rf "$OUT"
echo "[4/4] Running COBRA-meta..."

run_cobra \
    -q "$WORK/query.fa" \
    -f "$MSPADES" \
    -m "$COBRA_BAM" \
    -c "$WORK/coverage.tsv" \
    -a metaspades \
    -mink "$MINK" -maxk "$MAXK" \
    -o "$OUT" \
    -t "$THREADS"

# ── Report ─────────────────────────────────────────────────────────────
echo "──────────────────────────────────────────────────────────"
echo " COBRA results for $SAMPLE"
echo "──────────────────────────────────────────────────────────"
echo "Query contigs (SPADES viral):  $N_QUERY"

_count() { [ -s "$1" ] && grep -c '^>' "$1" || echo 0; }

N_CIRC=$(_count  "$OUT/COBRA_category_i_self_circular.fasta")
N_EXCA=$(_count  "$OUT/COBRA_category_ii-a_extended_circular_unique.fasta")
N_EXPA=$(_count  "$OUT/COBRA_category_ii-b_extended_partial_unique.fasta")
N_FAIL=$(_count  "$OUT/COBRA_category_ii-c_extended_failed.fasta")
N_ORPH=$(_count  "$OUT/COBRA_category_iii_orphan_end.fasta")

echo "  cat-i   self-circular:          $N_CIRC"
echo "  cat-ii-a extended circular:     $N_EXCA"
echo "  cat-ii-b extended partial:      $N_EXPA"
echo "  cat-ii-c extended failed:       $N_FAIL"
echo "  cat-iii  orphan end:            $N_ORPH"

if [ -s "$OUT/COBRA_joining_summary.txt" ]; then
    echo ""
    echo "Joining summary:"
    cat "$OUT/COBRA_joining_summary.txt"
fi

echo ""
echo "Full updated assembly: $OUT/contigs.new.fa"
echo "  size: $(du -sh "$OUT/contigs.new.fa" 2>/dev/null | cut -f1)"
echo "──────────────────────────────────────────────────────────"

# non-zero result = at least one category has sequences
TOTAL=$((N_CIRC + N_EXCA + N_EXPA))
if [ "$TOTAL" -eq 0 ] && [ "$N_ORPH" -eq 0 ]; then
    echo "[standalone] COBRA produced no output. Inspect $OUT/log and $OUT/debug.txt"
    exit 1
fi
