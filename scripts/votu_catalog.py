"""votu_catalog.py -- pure logic for the global vOTU catalog.

Kept out of Snakemake `run:` blocks on purpose: the per-sample clustering
this module replaces was broken for its entire lifetime precisely because
it lived in an untested `run:` block. Everything here is importable and
covered by tests/test_votu_catalog.py.
"""
import os


def prefixed_id(source_id, contig_id):
    """Namespace a contig ID by its source.

    Contig IDs are only unique within an assembly -- 'MEGAHIT_k141_10006'
    exists in many samples. Pooling without a prefix silently merges
    unrelated contigs into the same catalog entry.
    """
    return f"{source_id}|{contig_id}"


def build_pool(sources, pool_path, provenance_path):
    """Concatenate viral FASTAs into one pool with namespaced IDs.

    sources: list of (source_type, source_id, fasta_path); source_type is
             'sample' or 'group'.
    Returns {'n_sequences', 'n_sources', 'n_skipped'}.
    Missing or empty sources are skipped and counted, not fatal -- a sample
    with no viral contigs is a real outcome.
    """
    os.makedirs(os.path.dirname(pool_path) or ".", exist_ok=True)
    n_sequences = n_sources = n_skipped = 0

    with open(pool_path, "w") as out, open(provenance_path, "w") as prov:
        prov.write("member_id\tsource_type\tsource_id\toriginal_contig_id\n")
        for source_type, source_id, fasta_path in sources:
            if not fasta_path or not os.path.exists(fasta_path) \
                    or os.path.getsize(fasta_path) == 0:
                n_skipped += 1
                continue
            wrote_any = False
            with open(fasta_path) as fh:
                for line in fh:
                    if line.startswith(">"):
                        contig_id = line[1:].strip().split()[0]
                        member_id = prefixed_id(source_id, contig_id)
                        out.write(f">{member_id}\n")
                        prov.write(
                            f"{member_id}\t{source_type}\t{source_id}\t{contig_id}\n")
                        n_sequences += 1
                        wrote_any = True
                    else:
                        out.write(line)
            if wrote_any:
                n_sources += 1
            else:
                n_skipped += 1

    return {"n_sequences": n_sequences, "n_sources": n_sources,
            "n_skipped": n_skipped}


# skani triangle --sparse column layout (verified against skani 0.3.2):
#   0 Ref_file  1 Query_file  2 ANI  3 Align_fraction_ref
#   4 Align_fraction_query  5 Ref_name  6 Query_name
# The genome names are columns 5 and 6 -- columns 0 and 1 are FILE PATHS.
# Reading names from 0/1 (as the removed per-sample rule did) produces a
# parser that silently never matches anything.
_COL_ANI = 2
_COL_AF_REF = 3
_COL_AF_QUERY = 4
_COL_REF_NAME = 5
_COL_QUERY_NAME = 6
_N_SPARSE_COLS = 7


def parse_skani_sparse(path, ani_min, af_min, valid_ids):
    """Read `skani triangle --sparse` output into vOTU edges.

    An edge is kept when ANI >= ani_min AND max(af_ref, af_query) >= af_min
    (ICTV / Roux et al. 2019). Self-comparisons and names absent from
    valid_ids are dropped.

    Returns a list of (name_a, name_b) tuples.
    """
    edges = []
    if not path or not os.path.exists(path):
        return edges

    with open(path) as fh:
        # Skip the header row. A dense PHYLIP matrix opens with a bare
        # genome count instead, and yields no parseable rows below -- its
        # lines never reach the 7 columns the sparse format guarantees.
        fh.readline()
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < _N_SPARSE_COLS:
                continue
            try:
                ani = float(parts[_COL_ANI])
                af_ref = float(parts[_COL_AF_REF])
                af_query = float(parts[_COL_AF_QUERY])
            except ValueError:
                continue
            a = parts[_COL_REF_NAME]
            b = parts[_COL_QUERY_NAME]
            if a == b:
                continue
            if a not in valid_ids or b not in valid_ids:
                continue
            if ani >= ani_min and max(af_ref, af_query) >= af_min:
                edges.append((a, b))
    return edges
