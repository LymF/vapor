import json
import vapor


def test_track_overrides_single():
    out = vapor._track_overrides("viral")
    assert len(out) == 1
    key, _, val = out[0].partition("=")
    assert key == "tracks"
    d = json.loads(val)
    assert d == {"reads": False, "viral": True, "prok": False}


def test_track_overrides_multi():
    d = json.loads(vapor._track_overrides("viral,prok")[0].split("=", 1)[1])
    assert d == {"reads": False, "viral": True, "prok": True}


def test_track_overrides_reads():
    d = json.loads(vapor._track_overrides("reads")[0].split("=", 1)[1])
    assert d["reads"] is True and d["viral"] is False


def test_stage_alias_maps_to_rule():
    assert vapor._STAGE_ALIASES["assembly"] == "mmseqs2"
    assert vapor._STAGE_ALIASES["qc"] == "fastp"
    assert vapor._STAGE_ALIASES["viral"] == "viral_consensus"
    assert vapor._STAGE_ALIASES["binning"] == "binette"
