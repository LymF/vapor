import pytest
from pipeline_config import resolve_pipeline_config, validate_pipeline_config


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
