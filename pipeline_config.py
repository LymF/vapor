"""Pure config resolution + validation for the VAPOR Snakefile.

Kept import-free of Snakemake so it can be unit-tested with pytest and imported
from the Snakefile alike.
"""

_VALID_GROUPING = ("metadata", "all", "none")


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

    grouping = str(coas.get("grouping", "metadata")).strip().lower()
    coassembly_enabled = _b(coas, "enabled", False) and grouping != "none"

    return {
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

    any_analysis = (
        r["track_reads"] or r["track_viral"] or r["track_prok"]
        or r["coassembly_enabled"] or r["cobinning_multisplit"]
    )
    if not any_analysis:
        raise ValueError(
            "Nenhuma análise habilitada: ligue pelo menos uma de "
            "tracks.reads/viral/prok, coassembly.enabled ou cobinning_multisplit."
        )
