#!/bin/bash
# ══════════════════════════════════════════════════════════════════════
# run_iontorrent_assembly.sh
#
# Prepares Ion Torrent (SE, pre-trimmed) samples for the VAPOR pipeline:
#   1. Groups and concatenates FASTQs by sample (GOS[0-9]+)
#   2. metaSPAdes --iontorrent  (primary assembler)
#   3. MEGAHIT -r               (secondary assembler, SE mode)
#   4. Merge assemblies + filter by minimum contig length
#   5. MMseqs2 dedup → rep_seq.fasta  (VAPOR's central hub)
#   6. BWA-MEM2 SE mapping → sorted BAM + depth.txt
#   7. Creates sentinel files so VAPOR skips these steps and continues downstream
#
# Tool resolution order (per tool): conda envs → Docker → Singularity → fail
# Docker image tags are resolved automatically from quay.io/biocontainers API.
#
# After this script completes:
#   conda activate snakemake
#   cd <PIPE_DIR>
#   snakemake --use-conda --cores $THREADS
# ══════════════════════════════════════════════════════════════════════
set -euo pipefail

# ══════════════════════════════════════════════════════════════════════
# CONFIG — edit these variables before running
# ══════════════════════════════════════════════════════════════════════
FASTQ_DIR="$HOME/amazon"               # directory with Ion Torrent .fastq.gz files
OUTDIR="$HOME/amazon/results"          # VAPOR output root (must match config.yaml outdir)
PIPE_DIR="/home/alumnos/lmelo/vapor"   # VAPOR pipeline directory

THREADS=32
MIN_CONTIG=3000                        # must match config.yaml min_contig
MIN_SEQ_ID=0.95                        # must match config.yaml min_seq_id
SPADES_MEM=120                         # SPAdes RAM limit (GB)
MEGAHIT_MEM=120                        # MEGAHIT RAM limit (GB)

# Regex to extract sample name from filename
# e.g. GOS947_08_BC006.2.btrim_l100_a12_w10.fastq.gz → GOS947_08
# The _08 / _3 / _01 markers are part of the sample identity — do NOT strip them
SAMPLE_REGEX='GOS[0-9]+_[0-9]+'

# ── Tool versions + verified quay.io/biocontainers tags ───────────────
# Tags verified against quay.io registry (format: version--hash_build).
# The resolve_image() function tries to find the latest tag dynamically;
# these are used as hardcoded fallbacks if the API is unreachable.
VER_SPADES="4.2.0";   TAG_SPADES="quay.io/biocontainers/spades:4.2.0--h8d6e82b_2"
VER_MEGAHIT="1.2.9";  TAG_MEGAHIT="quay.io/biocontainers/megahit:1.2.9--haf24da9_8"
VER_MMSEQS="18.8cc5c"; TAG_MMSEQS="quay.io/biocontainers/mmseqs2:18.8cc5c--hd6d6fdc_0"
VER_BWA="2.3";        TAG_BWA="quay.io/biocontainers/bwa-mem2:2.3--he70b90d_0"
VER_SAMTOOLS="1.23";  TAG_SAMTOOLS="quay.io/biocontainers/samtools:1.23--h96c455f_0"
VER_METABAT="2.18";   TAG_METABAT="quay.io/biocontainers/metabat2:2.18_23_gc869c52--h61f4f8f_0"
# ══════════════════════════════════════════════════════════════════════


# ══════════════════════════════════════════════════════════════════════
# SETUP
# ══════════════════════════════════════════════════════════════════════

# ── Detect container runtime ──────────────────────────────────────────
DOCKER_CMD=""
SINGULARITY_CMD=""
if command -v docker &>/dev/null && docker info &>/dev/null 2>&1; then
    DOCKER_CMD="docker"
    echo "[setup] Docker: $(docker --version)"
elif command -v singularity &>/dev/null; then
    SINGULARITY_CMD="singularity"
    echo "[setup] Singularity: $(singularity --version)"
else
    echo "[setup] No container runtime — relying on conda envs only"
fi

# ── Resolve correct quay.io/biocontainers image tag ───────────────────
# Tries two APIs in order:
#   1. Docker Registry v2  — GET /v2/.../tags/list  (no auth required for public images)
#   2. quay.io custom API  — backup
# Verifies the winning tag exists with `docker manifest inspect`.
resolve_image() {
    local tool=$1      # e.g. "spades"
    local version=$2   # e.g. "4.2.0"
    local base="quay.io/biocontainers/${tool}"
    local tag=""

    # ── Strategy 1: Docker Registry v2 tags/list ──────────────────────
    if command -v curl &>/dev/null; then
        tag=$(curl -sf --max-time 10 \
            "https://quay.io/v2/biocontainers/${tool}/tags/list" \
            2>/dev/null | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    ver = '${version}'
    tags = [t for t in data.get('tags', []) if t.startswith(ver + '--')]
    tags.sort(key=lambda t: int(t.rsplit('_', 1)[-1]) if t.rsplit('_', 1)[-1].isdigit() else 0)
    print(tags[-1] if tags else '')
except Exception:
    print('')
" 2>/dev/null || true)
    fi

    # ── Strategy 2: quay.io custom API ────────────────────────────────
    if [ -z "$tag" ] && command -v curl &>/dev/null; then
        tag=$(curl -sf --max-time 10 \
            "https://quay.io/api/v1/repository/biocontainers/${tool}/tag/?limit=100" \
            2>/dev/null | python3 -c "
import sys, json
try:
    data = json.load(sys.stdin)
    ver = '${version}'
    tags = [
        t['name'] for t in data.get('tags', [])
        if t['name'].startswith(ver + '--') and '--' in t['name']
    ]
    tags.sort(key=lambda t: int(t.rsplit('_', 1)[-1]) if t.rsplit('_', 1)[-1].isdigit() else 0)
    print(tags[-1] if tags else '')
except Exception:
    print('')
" 2>/dev/null || true)
    fi

    if [ -z "$tag" ]; then
        echo "[setup] WARNING: could not resolve tag for ${tool}:${version} — will try docker pull at runtime" >&2
        # Return version-only; docker will error with a clear message if it doesn't exist
        echo "${base}:${version}"
        return
    fi

    local full_image="${base}:${tag}"

    # ── Verify tag actually exists (manifest inspect, no download) ─────
    if [ -n "$DOCKER_CMD" ]; then
        if ! $DOCKER_CMD manifest inspect "$full_image" &>/dev/null; then
            echo "[setup] WARNING: manifest check failed for $full_image" >&2
            echo "${base}:${version}"
            return
        fi
    fi

    echo "$full_image"
}

# ── Resolve image tags at startup ─────────────────────────────────────
if [ -n "$DOCKER_CMD" ] || [ -n "$SINGULARITY_CMD" ]; then
    echo "[setup] Resolving container image tags from quay.io ..."
    IMG_SPADES=$(resolve_image  "spades"   "$VER_SPADES")
    IMG_MEGAHIT=$(resolve_image "megahit"  "$VER_MEGAHIT")
    IMG_MMSEQS=$(resolve_image  "mmseqs2"  "$VER_MMSEQS")
    IMG_BWA=$(resolve_image     "bwa-mem2" "$VER_BWA")
    IMG_SAMTOOLS=$(resolve_image "samtools" "$VER_SAMTOOLS")
    IMG_METABAT=$(resolve_image "metabat2" "$VER_METABAT")
    echo "  spades   → $IMG_SPADES"
    echo "  megahit  → $IMG_MEGAHIT"
    echo "  mmseqs2  → $IMG_MMSEQS"
    echo "  bwa-mem2 → $IMG_BWA"
    echo "  samtools → $IMG_SAMTOOLS"
    echo "  metabat2 → $IMG_METABAT"
else
    IMG_SPADES="" IMG_MEGAHIT="" IMG_MMSEQS=""
    IMG_BWA="" IMG_SAMTOOLS="" IMG_METABAT=""
fi

# ── Find tools in conda envs ──────────────────────────────────────────
find_in_conda() {
    local name=$1
    local hit

    # 1. snakemake cache inside the pipeline dir (where snakemake normally writes envs)
    hit=$(find "$PIPE_DIR/.snakemake/conda" -name "$name" -type f 2>/dev/null \
          | grep -v __pycache__ | head -1)

    # 2. snakemake cache in home dir (fallback / older runs)
    if [ -z "$hit" ]; then
        hit=$(find ~/.snakemake/conda -name "$name" -type f 2>/dev/null \
              | grep -v __pycache__ | head -1)
    fi

    # 3. named conda envs (conda info --base)/envs/
    if [ -z "$hit" ]; then
        hit=$(find "$(conda info --base 2>/dev/null)/envs" -name "$name" -type f 2>/dev/null \
              | grep -v __pycache__ | head -1)
    fi

    # 4. tool is already in PATH
    if [ -z "$hit" ]; then
        hit=$(command -v "$name" 2>/dev/null || true)
    fi

    echo "$hit"
}

echo "[setup] Locating tools in conda envs ..."
SPADES_BIN=$(find_in_conda "spades.py")
MEGAHIT_BIN=$(find_in_conda "megahit")
MMSEQS_BIN=$(find_in_conda "mmseqs")
BWA_BIN=$(find_in_conda "bwa-mem2")
SAMTOOLS_BIN=$(find_in_conda "samtools")
JGI_BIN=$(find_in_conda "jgi_summarize_bam_contig_depths")

# Report and validate
MISSING=0
for entry in \
    "SPADES_BIN:spades:$IMG_SPADES" \
    "MEGAHIT_BIN:megahit:$IMG_MEGAHIT" \
    "MMSEQS_BIN:mmseqs2:$IMG_MMSEQS" \
    "BWA_BIN:bwa-mem2:$IMG_BWA" \
    "SAMTOOLS_BIN:samtools:$IMG_SAMTOOLS" \
    "JGI_BIN:metabat2(jgi):$IMG_METABAT"
do
    varname="${entry%%:*}"; rest="${entry#*:}"
    toolname="${rest%%:*}"; image="${rest#*:}"
    val="${!varname}"
    if [ -n "$val" ]; then
        echo "  $toolname → conda: $val"
    elif [ -n "$image" ]; then
        echo "  $toolname → container: $image"
    else
        echo "  [ERROR] $toolname: not found in conda and no container runtime available"
        echo "          Fix: snakemake --use-conda --cores 1 --create-envs-only"
        MISSING=1
    fi
done
[ "$MISSING" -eq 0 ] || exit 1
echo ""

# ── Per-tool runner functions ─────────────────────────────────────────
# Each function knows the correct invocation for both conda and container.

_container_run() {
    local image=$1; shift
    if [ -n "$DOCKER_CMD" ]; then
        $DOCKER_CMD run --rm -u "$(id -u):$(id -g)" \
            -v "$FASTQ_DIR:$FASTQ_DIR" -v "$OUTDIR:$OUTDIR" \
            "$image" "$@"
    elif [ -n "$SINGULARITY_CMD" ]; then
        $SINGULARITY_CMD exec --bind "$FASTQ_DIR,$OUTDIR" \
            "docker://$image" "$@"
    else
        echo "[ERROR] No container runtime for: $image"; exit 1
    fi
}

run_spades() {
    # conda: python3 /full/path/to/spades.py
    # container: spades.py is in PATH
    if [ -n "$SPADES_BIN" ]; then
        python3 "$SPADES_BIN" "$@"
    else
        _container_run "$IMG_SPADES" spades.py "$@"
    fi
}

run_megahit() {
    if [ -n "$MEGAHIT_BIN" ]; then
        "$MEGAHIT_BIN" "$@"
    else
        _container_run "$IMG_MEGAHIT" megahit "$@"
    fi
}

run_mmseqs() {
    if [ -n "$MMSEQS_BIN" ]; then
        "$MMSEQS_BIN" "$@"
    else
        _container_run "$IMG_MMSEQS" mmseqs "$@"
    fi
}

run_bwa_index() {
    if [ -n "$BWA_BIN" ]; then
        "$BWA_BIN" index "$@"
    else
        _container_run "$IMG_BWA" bwa-mem2 index "$@"
    fi
}

# BWA-MEM2 | samtools sort — piped, needs both in same execution context
run_bwa_map_and_sort() {
    local threads=$1 idx=$2 reads=$3 bam=$4 log=$5
    if [ -n "$BWA_BIN" ] && [ -n "$SAMTOOLS_BIN" ]; then
        "$BWA_BIN" mem -t "$threads" "$idx" "$reads" 2>"$log" \
            | "$SAMTOOLS_BIN" sort -@ "$threads" -o "$bam"
    elif [ -n "$DOCKER_CMD" ]; then
        $DOCKER_CMD run --rm -u "$(id -u):$(id -g)" \
            -v "$FASTQ_DIR:$FASTQ_DIR" -v "$OUTDIR:$OUTDIR" \
            "$IMG_BWA" \
            bwa-mem2 mem -t "$threads" "$idx" "$reads" \
            2>"$log" \
        | $DOCKER_CMD run --rm -i \
            -u "$(id -u):$(id -g)" \
            -v "$OUTDIR:$OUTDIR" \
            "$IMG_SAMTOOLS" \
            samtools sort -@ "$threads" -o "$bam" -
    elif [ -n "$SINGULARITY_CMD" ]; then
        $SINGULARITY_CMD exec --bind "$FASTQ_DIR,$OUTDIR" \
            "docker://$IMG_BWA" \
            bwa-mem2 mem -t "$threads" "$idx" "$reads" \
            2>"$log" \
        | $SINGULARITY_CMD exec --bind "$OUTDIR" \
            "docker://$IMG_SAMTOOLS" \
            samtools sort -@ "$threads" -o "$bam" -
    fi
}

run_samtools() {
    if [ -n "$SAMTOOLS_BIN" ]; then
        "$SAMTOOLS_BIN" "$@"
    else
        _container_run "$IMG_SAMTOOLS" samtools "$@"
    fi
}

run_jgi() {
    if [ -n "$JGI_BIN" ]; then
        "$JGI_BIN" "$@"
    else
        _container_run "$IMG_METABAT" jgi_summarize_bam_contig_depths "$@"
    fi
}


# ══════════════════════════════════════════════════════════════════════
# Sample discovery
# ══════════════════════════════════════════════════════════════════════
echo "[samples] Scanning $FASTQ_DIR ..."
declare -A SAMPLE_FILES
for f in "$FASTQ_DIR"/*.fastq.gz "$FASTQ_DIR"/*.fq.gz; do
    [ -f "$f" ] || continue
    sample=$(basename "$f" | grep -oP "$SAMPLE_REGEX" | head -1 || true)
    [ -n "$sample" ] || continue
    SAMPLE_FILES[$sample]+="$f "
done

if [ ${#SAMPLE_FILES[@]} -eq 0 ]; then
    echo "[ERROR] No samples found in $FASTQ_DIR matching pattern '$SAMPLE_REGEX'"
    exit 1
fi
echo "[samples] Found: $(echo "${!SAMPLE_FILES[@]}" | tr ' ' '\n' | sort | tr '\n' ' ')"
echo ""


# ══════════════════════════════════════════════════════════════════════
# Per-sample loop
# ══════════════════════════════════════════════════════════════════════
for SAMPLE in $(echo "${!SAMPLE_FILES[@]}" | tr ' ' '\n' | sort); do

    echo "════════════════════════════════════════════════════════"
    echo " Sample: $SAMPLE"
    echo "════════════════════════════════════════════════════════"

    S="$OUTDIR/$SAMPLE"
    mkdir -p \
        "$S/assembly/metaspades" \
        "$S/assembly/megahit" \
        "$S/assembly/metaviral" \
        "$S/mmseqs/tmp" \
        "$S/mapping" \
        "$S/trimmed" \
        "$S/qc_raw" \
        "$S/logs" \
        "$S/benchmarks"

    # ── Step 1: Concatenate FASTQs ────────────────────────────────────
    CAT_FQ="$S/${SAMPLE}_all.fq.gz"
    if [ ! -s "$CAT_FQ" ]; then
        echo "[1/6] Concatenating FASTQs for $SAMPLE..."
        files=( ${SAMPLE_FILES[$SAMPLE]} )
        echo "      Files: ${#files[@]}"
        cat "${files[@]}" > "$CAT_FQ"
        echo "      → $CAT_FQ ($(du -sh "$CAT_FQ" | cut -f1))"
    else
        echo "[1/6] Concat exists: $CAT_FQ"
    fi

    # Symlink as trimmed R1 (pipeline expects this path)
    ln -sf "$(realpath "$CAT_FQ")" "$S/trimmed/${SAMPLE}_R1_fastp.fq.gz" 2>/dev/null || true
    [ -f "$S/trimmed/${SAMPLE}_R2_fastp.fq.gz" ] || touch "$S/trimmed/${SAMPLE}_R2_fastp.fq.gz"
    touch "$S/qc_raw/done.txt"

    # ── Step 2: metaSPAdes --iontorrent ──────────────────────────────
    SPADES_CONTIGS="$S/assembly/metaspades/contigs.fasta"
    if [ ! -s "$SPADES_CONTIGS" ]; then
        echo "[2/6] metaSPAdes --iontorrent (SE)..."
        run_spades --iontorrent \
            -s "$CAT_FQ" \
            -t "$THREADS" -m "$SPADES_MEM" \
            -o "$S/assembly/metaspades" \
            > "$S/logs/metaspades.log" 2>&1 \
        && echo "      OK: $(grep -c '^>' "$SPADES_CONTIGS" 2>/dev/null || echo 0) contigs" \
        || echo "      [WARNING] SPAdes failed — check $S/logs/metaspades.log"
    else
        echo "[2/6] SPAdes output exists: $SPADES_CONTIGS"
    fi
    touch "$S/assembly/metaspades/assembly_graph_with_scaffolds.gfa"
    touch "$S/assembly/metaviral/contigs.fasta"

    # ── Step 3: MEGAHIT SE ────────────────────────────────────────────
    MEGAHIT_CONTIGS="$S/assembly/megahit/final.contigs.fa"
    if [ ! -s "$MEGAHIT_CONTIGS" ]; then
        echo "[3/6] MEGAHIT (SE mode)..."
        rm -rf "$S/assembly/megahit"
        run_megahit \
            -r "$CAT_FQ" \
            -t "$THREADS" \
            -m "${MEGAHIT_MEM}000000000" \
            --min-contig-len "$MIN_CONTIG" \
            -o "$S/assembly/megahit" \
            > "$S/logs/megahit.log" 2>&1 \
        && echo "      OK: $(grep -c '^>' "$MEGAHIT_CONTIGS" 2>/dev/null || echo 0) contigs" \
        || echo "      [WARNING] MEGAHIT failed — check $S/logs/megahit.log"
    else
        echo "[3/6] MEGAHIT output exists: $MEGAHIT_CONTIGS"
    fi

    # ── Step 4: Merge + length filter ─────────────────────────────────
    MERGED_FA="$S/mmseqs/${SAMPLE}_merged_filtered.fasta"
    if [ ! -s "$MERGED_FA" ]; then
        echo "[4/6] Merging assemblies (>= ${MIN_CONTIG} bp)..."
        python3 - \
            "$SPADES_CONTIGS" \
            "$MEGAHIT_CONTIGS" \
            "$S/assembly/metaviral/contigs.fasta" \
            "$MERGED_FA" \
            "$MIN_CONTIG" <<'PYEOF'
import sys, os

def iter_fasta(path, prefix, min_len):
    if not os.path.isfile(path) or os.path.getsize(path) == 0:
        return
    with open(path) as fh:
        header, seq = "", []
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header and len("".join(seq)) >= min_len:
                    yield header, "".join(seq)
                header = line.replace(">", f">{prefix}_", 1)
                seq = []
            else:
                seq.append(line)
        if header and len("".join(seq)) >= min_len:
            yield header, "".join(seq)

spades_fa, megahit_fa, metaviral_fa, out_fa, min_len = \
    sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5])
sources = [(spades_fa, "SPADES"), (megahit_fa, "MEGAHIT"), (metaviral_fa, "METAVIRAL")]
kept = 0
with open(out_fa, "w") as fh:
    for path, pfx in sources:
        for h, s in iter_fasta(path, pfx, min_len):
            fh.write(h + "\n" + s + "\n")
            kept += 1
print(f"      {kept} contigs >= {min_len} bp → {out_fa}")
PYEOF
    else
        echo "[4/6] Merged FASTA exists: $MERGED_FA"
    fi

    # ── Step 5: MMseqs2 dedup ─────────────────────────────────────────
    REP_SEQ="$S/mmseqs/${SAMPLE}_rep_seq.fasta"
    if [ ! -s "$REP_SEQ" ]; then
        echo "[5/6] MMseqs2 easy-linclust (${MIN_SEQ_ID} identity)..."
        run_mmseqs easy-linclust \
            "$MERGED_FA" \
            "$S/mmseqs/$SAMPLE" \
            "$S/mmseqs/tmp" \
            --min-seq-id "$MIN_SEQ_ID" \
            -c 0.9 --cov-mode 2 \
            --threads "$THREADS" \
            > "$S/logs/mmseqs2.log" 2>&1
        echo "      OK: $(grep -c '^>' "$REP_SEQ") representative sequences"
    else
        echo "[5/6] rep_seq.fasta exists: $REP_SEQ"
    fi

    # ── Step 6: BWA-MEM2 SE mapping + depth ──────────────────────────
    BAM="$S/mapping/${SAMPLE}.sorted.bam"
    DEPTH="$S/mapping/${SAMPLE}_depth.txt"
    IDX="$S/mapping/contigs_index"

    if [ ! -f "${IDX}.bwt.2bit.64" ]; then
        echo "[6/6] BWA-MEM2 index..."
        run_bwa_index -p "$IDX" "$REP_SEQ" \
            > "$S/logs/bwa_index.log" 2>&1
    fi

    if [ ! -s "$BAM" ]; then
        echo "[6/6] BWA-MEM2 SE mapping..."
        run_bwa_map_and_sort "$THREADS" "$IDX" "$CAT_FQ" "$BAM" \
            "$S/logs/bwa_mem.log"
        run_samtools index "$BAM"
        run_samtools flagstat "$BAM" \
            > "$S/mapping/flagstat.txt" 2>> "$S/logs/bwa_mem.log"
        touch "$S/mapping/bwa_mem_done.txt"
        echo "      OK: BAM created"
    else
        echo "[6/6] BAM exists: $BAM"
        touch "$S/mapping/bwa_mem_done.txt"
    fi

    if [ ! -s "$DEPTH" ] && [ -s "$BAM" ]; then
        echo "[6/6] Calculating depth..."
        run_jgi --outputDepth "$DEPTH" "$BAM" \
            > "$S/logs/calc_depth.log" 2>&1
        echo "      OK: $DEPTH"
    elif [ -s "$DEPTH" ]; then
        echo "[6/6] Depth file exists: $DEPTH"
    fi

    echo ""
    echo "  [done] $SAMPLE"
    echo ""

done

# ══════════════════════════════════════════════════════════════════════
echo "════════════════════════════════════════════════════════"
echo " All samples processed."
echo ""
echo " Next: edit config.yaml, then run VAPOR:"
echo "   fastq_dir: \"$FASTQ_DIR\""
echo "   outdir:    \"$OUTDIR\""
echo ""
echo "   conda activate snakemake"
echo "   cd $PIPE_DIR"
echo "   snakemake --use-conda --cores $THREADS"
echo ""
echo " VAPOR will skip assembly/mapping (outputs already exist)"
echo " and continue from viral detection, binning, taxonomy, etc."
echo "════════════════════════════════════════════════════════"
