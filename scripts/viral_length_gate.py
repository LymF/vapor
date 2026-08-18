"""viral_length_gate.py — the composite length/quality/bin gate for short viral
contigs, shared by `rule viral_nonredundant` (rules/viral_binning.smk) and its
group equivalent `rule coassembly_viral_nonredundant` (rules/coassembly.smk).

Decision (docs/ROADMAP_SIMPLIFICACAO.md, item (e), 2026-08-18): the length
floor only applies where there is no independent evidence the sequence is a
real genome. A sequence is kept if ANY of:
  1. it sits in a vRhyme bin (the bin itself is the evidence), or
  2. CheckV calls it Complete / High-quality / Medium-quality, or
     completeness >= 50%, or
  3. it is >= VIRAL_MIN_CONTIG bp (default 5000, the IMG/VR / Earth's Virome
     Protocol / MVP cutoff — Roux et al. 2021 NAR 49:D764;
     Coclet, Camargo & Roux 2024 mSystems 9:e00888-24).

Both callers run the gate AFTER vRhyme, so all three arms are always
available — sample/group parity (docs/AUDITORIA_COASSEMBLY_PARES.md). The
one place the `binned` arm never fires is long-read co-assembly groups: that
track has no group-level vRhyme at all (rules/coassembly.smk gates it on
`not LONG_READS`), so `_catalog_sources()` (rules/votu_catalog.smk) sources
those groups straight from the pre-binning CheckV-trimmed set instead of
routing them through this gate.

Kept import-free of Snakemake so it is unit-testable with plain pytest/unittest.
"""

from collections import namedtuple

# The quality arm is a FIXED set, deliberately NOT the configurable
# VIRAL_KEEP_TIERS. Those two look interchangeable and are not:
# VIRAL_KEEP_TIERS is driven by `viral_min_quality`, which ships as
# "not_determined" and therefore expands to ALL FIVE CheckV tiers --
# Low-quality and Not-determined included. Wiring the gate to it made arm 2
# true for every contig CheckV had ever touched, i.e. the whole gate was a
# no-op under the shipped config. This set mirrors MVP's published rule
# (Coclet, Camargo & Roux 2024, mSystems 9:e00888-24): the length floor is
# waived only for genomes with real quality evidence.
MQ_TIERS = frozenset({"Complete", "High-quality", "Medium-quality"})

GateResult = namedtuple("GateResult", ["kept", "reason"])

# reason values, used both as the per-sequence tag and as counter keys.
REASON_BINNED = "binned"
REASON_QUALITY = "quality"
REASON_LENGTH = "length"
REASON_DROPPED = "dropped"


def passes_gate(*, binned, quality_tier, completeness, length,
                 min_length, keep_tiers=MQ_TIERS):
    """Evaluate the composite gate for one sequence.

    binned         : bool — True if the sequence belongs to a vRhyme bin.
    quality_tier   : str  — CheckV `checkv_quality` value, "" if unknown.
    completeness   : float — CheckV `completeness` (%), 0.0 if unknown.
    length         : int  — sequence length in bp.
    min_length     : int  — the length-only floor (VIRAL_MIN_CONTIG).
    keep_tiers     : container of str — accepted CheckV tiers. Defaults to
                     MQ_TIERS; do NOT pass VIRAL_KEEP_TIERS here (see above).

    Returns a GateResult(kept: bool, reason: one of the REASON_* constants).
    The reason always names the FIRST arm satisfied, in the priority order
    bin > quality > length, purely for reporting -- all satisfied arms are
    equally sufficient for keeping the sequence.
    """
    if binned:
        return GateResult(True, REASON_BINNED)
    is_mq = (quality_tier or "") in keep_tiers or completeness >= 50.0
    if is_mq:
        return GateResult(True, REASON_QUALITY)
    if length >= min_length:
        return GateResult(True, REASON_LENGTH)
    return GateResult(False, REASON_DROPPED)


def summarize(results):
    """results: iterable of GateResult -> dict of counts per reason + total."""
    counts = {"total": 0, REASON_BINNED: 0, REASON_QUALITY: 0,
              REASON_LENGTH: 0, REASON_DROPPED: 0}
    for r in results:
        counts["total"] += 1
        counts[r.reason] += 1
    return counts


# ── Discard audit sidecar (item (e), 2026-08-18 follow-up) ─────────────────
# Sequences the gate drops used to simply vanish. `rule viral_nonredundant`
# and `rule coassembly_viral_nonredundant` now also write a
# `{sample|group}_viral_discarded.fasta` of every dropped sequence, plus a
# `{sample|group}_viral_discarded.tsv` sidecar carrying the evidence the gate
# saw for each one -- a short "Not-determined" contig is exactly what a
# genuinely novel virus looks like to CheckV (no close reference to score
# against), so the discard is a precision bet, not a verdict, and has to stay
# inspectable. The FASTA header is deliberately the BARE contig_id -- the
# same id used everywhere else in the pipeline -- with NO reason encoded in
# it, so nothing downstream that keys off the header (BLAST, a second CheckV
# pass, geNomad at a looser threshold, ...) breaks on a decorated id. All
# per-sequence metadata lives in the TSV instead, joined back on contig_id.
DISCARD_TSV_COLUMNS = [
    "contig_id", "length", "checkv_quality", "checkv_completeness",
    "in_vrhyme_bin", "source_id",
]


def format_discard_row(contig_id, length, quality_tier, completeness, source_id):
    """Build one DISCARD_TSV_COLUMNS-ordered row (list of str) for a sequence
    the composite gate dropped.

    contig_id    : str — MUST equal the discard FASTA header for this
                   sequence exactly (the TSV's join key back to the FASTA).
    length       : int — post-trim sequence length in bp (what the gate's
                   length arm actually evaluated).
    quality_tier : str or None — CheckV `checkv_quality`; None (-> "") if
                   CheckV never scored this contig_id at all. Do NOT collapse
                   "never scored" and "scored, tier empty" to the same value
                   here -- that distinction is the point of the sidecar.
    completeness : float or None — CheckV `completeness` (%); None (-> "")
                   under the same never-scored condition as quality_tier.
    source_id    : str — sample or group name, so per-run TSVs concatenate
                   into one auditable table across a whole project.

    `in_vrhyme_bin` is always the literal "False": every row in this table
    is by construction a sequence gate arm 1 did NOT keep (binned sequences
    are never dropped), so the column has one value here -- kept explicit
    rather than omitted so the TSV is self-describing without this function's
    source, and so its schema matches a future run where that stops being
    universally true.
    """
    return [
        str(contig_id),
        str(length),
        quality_tier or "",
        "" if completeness is None else f"{completeness:.2f}",
        "False",
        str(source_id),
    ]
