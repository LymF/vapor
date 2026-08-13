"""Regression test against a real 6-sample slice of the June/2026 run.

Skipped when the reference data is absent. See
tests/data/votu_regression/README.md for how to regenerate it.
"""
import os
import pytest
from votu_catalog import parse_skani_sparse, cluster_votus

DATA_DIR = os.path.join(os.path.dirname(__file__), "data", "votu_regression")
POOL = os.path.join(DATA_DIR, "pool.fasta")
ANI = os.path.join(DATA_DIR, "pool_sparse.tsv")

pytestmark = pytest.mark.skipif(
    not (os.path.exists(POOL) and os.path.exists(ANI)),
    reason="reference data absent -- see tests/data/votu_regression/README.md",
)


@pytest.fixture(scope="module")
def clusters():
    ids = [l[1:].split()[0] for l in open(POOL) if l.startswith(">")]
    edges = parse_skani_sparse(ANI, 95.0, 85.0, set(ids))
    return ids, cluster_votus(ids, edges, {})


def test_pool_size_matches_baseline(clusters):
    ids, _ = clusters
    assert len(ids) == 9653


def test_clustering_collapses_to_baseline(clusters):
    _, cl = clusters
    # Exact count is pinned: any drift means the criterion or parser changed.
    assert len(cl) == 5524


def test_reduction_is_substantial(clusters):
    ids, cl = clusters
    reduction = 1 - len(cl) / len(ids)
    assert reduction > 0.40


def test_shared_votus_are_detected(clusters):
    """~10% of vOTUs span more than one sample -- the whole point of pooling."""
    _, cl = clusters
    multi = sum(1 for c in cl if len({m.split("|")[0] for m in c}) > 1)
    assert multi == 577
