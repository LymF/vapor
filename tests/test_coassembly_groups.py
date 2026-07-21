import textwrap
import pytest
from coassembly_groups import parse_groups


def _write(tmp_path, content):
    p = tmp_path / "meta.tsv"
    p.write_text(textwrap.dedent(content))
    return str(p)


def test_mode_all_single_group():
    g = parse_groups("", ["s1", "s2", "s3"], "all")
    assert g == {"all": ["s1", "s2", "s3"]}


def test_metadata_groups(tmp_path):
    meta = _write(tmp_path, """\
        sample\tgroup
        s1\trio
        s2\trio
        s3\tsolo
    """)
    g = parse_groups(meta, ["s1", "s2", "s3"], "metadata")
    assert g == {"rio": ["s1", "s2"], "solo": ["s3"]}


def test_metadata_ignores_unknown_samples(tmp_path):
    meta = _write(tmp_path, """\
        sample\tgroup
        s1\trio
        sX\trio
    """)
    g = parse_groups(meta, ["s1"], "metadata")
    assert g == {"rio": ["s1"]}


def test_metadata_missing_file_raises():
    with pytest.raises(ValueError, match="metadata"):
        parse_groups("/no/such.tsv", ["s1"], "metadata")


def test_metadata_missing_columns_raises(tmp_path):
    meta = _write(tmp_path, "sample\tfoo\ns1\tbar\n")
    with pytest.raises(ValueError, match="colunas"):
        parse_groups(meta, ["s1"], "metadata")


def test_metadata_reserved_group_name_raises(tmp_path):
    meta = _write(tmp_path, """\
        sample\tgroup
        s1\tmultisplit
        s2\tmultisplit
    """)
    with pytest.raises(ValueError, match="reserved"):
        parse_groups(meta, ["s1", "s2"], "metadata")
