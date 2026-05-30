#!/usr/bin/env bash
# ──────────────────────────────────────────────────────────────────────
# run_cobra_standalone.sh — test COBRA-meta outside the Snakemake pipeline
#
# Uses the artifacts that VAPOR has already produced for a sample:
#   - viral/consensus/{sample}_viral_consensus.fasta   (query, SPADES_-prefixed)
#   - assembly/metaspades/contigs.fasta                (all metaspades contigs)
#   - assembly/metaspades/assembly_graph_with_scaffolds.gfa (informational)
#   - mapping/{sample}_depth.txt                       (jgi_summarize output)
#
# Output:
#   results-teste/{sample}/viral/cobra_standalone/cobra_out/
#
# Quick check against the integrated COBRA rule before turning it on at
# scale. Compares pre/post contig counts and total length.
#
# Usage:
#   bash scripts/test/run_cobra_standalone.sh --sample SAMPLE [--results-dir results-teste] [--threads 16]
#   bash scripts/test/run_cobra_standalone.sh --sample SAMPLE --use-conda    # conda env instead of docker
# ──────────────────────────────────────────────────────────────────────
set -euo pipefail

SAMPLE=""
RESULTS_DIR="results-teste"
THREADS=16
USE_CONDA=0
MINK=21
MAXK=127

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

if [ -z "$SAMPLE" ]; then
    echo "ERROR: --sample required" >&2; exit 2
fi

SAMPLE_DIR="$RESULTS_DIR/$SAMPLE"
VIRAL="$SAMPLE_DIR/viral/consensus/${SAMPLE}_viral_consensus.fasta"
MSPADES="$SAMPLE_DIR/assembly/metaspades/contigs.fasta"
DEPTH="$SAMPLE_DIR/mapping/${SAMPLE}_depth.txt"
WORK="$SAMPLE_DIR/viral/cobra_standalone"

for f in "$VIRAL" "$MSPADES" "$DEPTH"; do
    [ -s "$f" ] || { echo "ERROR: missing $f" >&2; exit 1; }
done

mkdir -p "$WORK"

# 1) SPADES_-prefixed viral contigs → strip prefix → query.fa
awk 'BEGIN{w=0} /^>/{w=0; if($0 ~ /^>SPADES_/){sub("^>SPADES_",">"); w=1; print} next} w{print}' \
    "$VIRAL" > "$WORK/query.fa"

N_QUERY=$(grep -c '^>' "$WORK/query.fa" || echo 0)
echo "[standalone] SPADES viral query contigs: $N_QUERY"
if [ "$N_QUERY" -eq 0 ]; then
    echo "[standalone] No SPADES viral contigs — nothing to extend."
    exit 0
fi

# 2) jgi depth → 2-col, SPADES_ only with prefix stripped
awk 'BEGIN{FS=OFS="\t"} NR==1{next} /^SPADES_/{sub("^SPADES_",""); print $1, $3}' \
    "$DEPTH" > "$WORK/coverage.tsv"

# 3) Run COBRA-meta
OUT="$WORK/cobra_out"
rm -rf "$OUT"

if [ "$USE_CONDA" -eq 1 ]; then
    eval "$(conda shell.bash hook)"
    conda activate env_cobra
    cobra-meta \
        -q "$WORK/query.fa" \
        -f "$MSPADES" \
        -c "$WORK/coverage.tsv" \
        -a metaspades \
        -mink "$MINK" -maxk "$MAXK" \
        -o "$OUT" \
        -t "$THREADS"
else
    docker run --rm \
        -v "$PWD:/data" \
        -w /data \
        quay.io/biocontainers/cobra-meta:1.2.3--pyhdfd78af_0 \
        cobra-meta \
            -q "$WORK/query.fa" \
            -f "$MSPADES" \
            -c "$WORK/coverage.tsv" \
            -a metaspades \
            -mink "$MINK" -maxk "$MAXK" \
            -o "$OUT" \
            -t "$THREADS"
fi

# 4) Report
EXT=""
for cand in "$OUT"/COBRA_all_assemblies.fasta "$OUT"/all.cobra.fasta "$OUT"/*_extended_*.fasta "$OUT"/*.new.fasta; do
    for f in $cand; do [ -f "$f" ] && [ -s "$f" ] && EXT="$f" && break 2; done
done

if [ -z "$EXT" ]; then
    echo "[standalone] No extended FASTA produced. Inspect $OUT for logs."
    exit 1
fi

N_OUT=$(grep -c '^>' "$EXT")
echo "──────────────────────────────────────────"
echo "Query contigs (SPADES viral):  $N_QUERY"
echo "COBRA extended contigs:        $N_OUT"
echo "Extended FASTA:                $EXT"
echo "──────────────────────────────────────────"
