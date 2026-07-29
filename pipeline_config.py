"""Pure config resolution + validation for the VAPOR Snakefile.

Kept import-free of Snakemake so it can be unit-tested with pytest and imported
from the Snakefile alike.
"""

_VALID_GROUPING = ("metadata", "all", "none")

# CheckV quality tiers ordered worst → best. Used to build the minimum-quality
# gate for the viral subset that feeds taxonomy / host prediction / annotation.
_CHECKV_TIERS = (
    "Not-determined", "Low-quality", "Medium-quality", "High-quality", "Complete",
)
_CHECKV_TIER_RANK = {t: i for i, t in enumerate(_CHECKV_TIERS)}

# Accepted `viral_min_quality` config values → minimum CheckV tier rank kept.
# Aliases (mq/hq/all/hyphens) are tolerated; the canonical set is validated.
_VIRAL_QUALITY_RANK = {
    "not_determined": 0, "notdetermined": 0, "all": 0, "none": 0,
    "low": 1,
    "medium": 2, "mq": 2,
    "high": 3, "hq": 3,
    "complete": 4,
}
_VALID_VIRAL_QUALITY = ("not_determined", "low", "medium", "high", "complete")


def _norm_quality(v) -> str:
    return str(v).strip().lower().replace("-", "_").replace(" ", "_")


def viral_min_quality_rank(min_quality) -> int:
    """CheckV tier rank at/above which contigs are kept (default: medium=2)."""
    return _VIRAL_QUALITY_RANK.get(_norm_quality(min_quality), 2)


def viral_keep_tiers(min_quality) -> frozenset:
    """Set of CheckV tier labels kept at/above the given minimum quality."""
    rank = viral_min_quality_rank(min_quality)
    return frozenset(t for t, r in _CHECKV_TIER_RANK.items() if r >= rank)


def _b(d, key, default):
    """Bool coercion tolerant of YAML strings ('false'/'true')."""
    v = d.get(key, default)
    if isinstance(v, str):
        return v.strip().lower() in ("true", "1", "yes", "on")
    return bool(v)


def resolve_pipeline_config(config: dict) -> dict:
    tracks = config.get("tracks", {}) or {}
    coas = config.get("coassembly", {}) or {}

    track_reads = _b(tracks, "reads", False)
    track_viral = _b(tracks, "viral", True)
    track_prok = _b(tracks, "prok", True)
    use_host_defense = _b(config, "use_host_defense", True)

    grouping = str(coas.get("grouping") or "metadata").strip().lower()
    coassembly_enabled = _b(coas, "enabled", False) and grouping != "none"

    viral_min_quality = _norm_quality(config.get("viral_min_quality", "medium"))

    return {
        "viral_min_quality": viral_min_quality,
        "viral_min_quality_rank": viral_min_quality_rank(viral_min_quality),
        "viral_keep_tiers": viral_keep_tiers(viral_min_quality),
        "track_reads": track_reads,
        "track_viral": track_viral,
        "track_prok": track_prok,
        "use_host_defense": use_host_defense,
        "integration_enabled": track_viral and track_prok and use_host_defense,
        "coassembly_enabled": coassembly_enabled,
        "coassembly_grouping": grouping,
        "coassembly_viral": _b(coas, "viral", True),
        "coassembly_binning": _b(coas, "binning", True),
        "cobinning_multisplit": _b(config, "cobinning_multisplit", False),
    }


def validate_pipeline_config(config: dict) -> None:
    r = resolve_pipeline_config(config)

    grouping = r["coassembly_grouping"]
    if grouping not in _VALID_GROUPING:
        raise ValueError(
            f"coassembly.grouping inválido: '{grouping}'. "
            f"Use um de: {', '.join(_VALID_GROUPING)}."
        )

    vmq = _norm_quality(config.get("viral_min_quality", "medium"))
    if vmq not in _VIRAL_QUALITY_RANK:
        raise ValueError(
            f"viral_min_quality inválido: '{config.get('viral_min_quality')}'. "
            f"Use um de: {', '.join(_VALID_VIRAL_QUALITY)}."
        )

    any_analysis = (
        r["track_reads"] or r["track_viral"] or r["track_prok"]
        or r["coassembly_enabled"] or r["cobinning_multisplit"]
    )
    if not any_analysis:
        raise ValueError(
            "Nenhuma análise habilitada: ligue pelo menos uma de "
            "tracks.reads/viral/prok, coassembly.enabled ou cobinning_multisplit."
        )
