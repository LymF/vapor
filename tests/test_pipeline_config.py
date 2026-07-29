import pytest
from pipeline_config import (
    resolve_pipeline_config,
    validate_pipeline_config,
    viral_keep_tiers,
    viral_min_quality_rank,
)


def test_defaults_preserve_current_behavior():
    r = resolve_pipeline_config({})
    assert r["track_reads"] is False
    assert r["track_viral"] is True
    assert r["track_prok"] is True
    assert r["use_host_defense"] is True
    assert r["integration_enabled"] is True   # viral and prok and host_defense
    assert r["coassembly_enabled"] is False
    assert r["cobinning_multisplit"] is False


def test_integration_requires_both_tracks_and_flag():
    assert resolve_pipeline_config(
        {"tracks": {"viral": True, "prok": False}})["integration_enabled"] is False
    assert resolve_pipeline_config(
        {"tracks": {"viral": True, "prok": True}, "use_host_defense": False}
    )["integration_enabled"] is False


def test_grouping_none_disables_coassembly():
    r = resolve_pipeline_config(
        {"coassembly": {"enabled": True, "grouping": "none"}})
    assert r["coassembly_enabled"] is False


def test_grouping_explicit_null_falls_back_to_metadata():
    r = resolve_pipeline_config({"coassembly": {"enabled": True, "grouping": None}})
    assert r["coassembly_grouping"] == "metadata"
    assert r["coassembly_enabled"] is True


def test_reads_only_track():
    r = resolve_pipeline_config(
        {"tracks": {"reads": True, "viral": False, "prok": False}})
    assert r["track_reads"] is True
    assert r["track_viral"] is False
    assert r["integration_enabled"] is False


def test_validate_rejects_no_analysis():
    with pytest.raises(ValueError, match="pelo menos uma"):
        validate_pipeline_config(
            {"tracks": {"reads": False, "viral": False, "prok": False}})


def test_validate_rejects_bad_grouping():
    with pytest.raises(ValueError, match="grouping"):
        validate_pipeline_config(
            {"coassembly": {"enabled": True, "grouping": "wrong"}})


def test_validate_accepts_defaults():
    validate_pipeline_config({})   # não levanta


def test_viral_min_quality_defaults_to_medium():
    r = resolve_pipeline_config({})
    assert r["viral_min_quality"] == "medium"
    assert r["viral_min_quality_rank"] == 2
    assert r["viral_keep_tiers"] == frozenset(
        {"Medium-quality", "High-quality", "Complete"})


@pytest.mark.parametrize("value,expected", [
    ("not_determined", {"Not-determined", "Low-quality", "Medium-quality",
                        "High-quality", "Complete"}),
    ("low",     {"Low-quality", "Medium-quality", "High-quality", "Complete"}),
    ("medium",  {"Medium-quality", "High-quality", "Complete"}),
    ("high",    {"High-quality", "Complete"}),
    ("complete", {"Complete"}),
])
def test_viral_keep_tiers_levels(value, expected):
    assert viral_keep_tiers(value) == frozenset(expected)
    # Lower thresholds are strict supersets of higher ones.
    assert viral_keep_tiers("medium") <= viral_keep_tiers(value) or \
        viral_min_quality_rank(value) > 2


def test_viral_min_quality_case_insensitive():
    r = resolve_pipeline_config({"viral_min_quality": "LOW"})
    assert r["viral_min_quality"] == "low"
    assert r["viral_min_quality_rank"] == 1


def test_validate_rejects_bad_viral_min_quality():
    with pytest.raises(ValueError, match="viral_min_quality"):
        validate_pipeline_config({"viral_min_quality": "garbage"})
