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


class ClusteringCollapseError(RuntimeError):
    """Raised when clustering produced exactly one cluster per input sequence.

    That outcome is the signature of the defect this module replaces: the
    old rule fed skani's dense matrix to an edge-list parser, so no edge was
    ever built and every contig stayed a singleton for the pipeline's entire
    history. Failing loudly here is the point -- a silent N-in-N-out catalog
    reports assembly redundancy as viral richness.
    """


def pick_representative(members, completeness):
    """Cluster representative: highest CheckV completeness, ties by member order."""
    return max(
        members,
        key=lambda m: (completeness.get(m, 0.0), -members.index(m)),
    )


def cluster_votus(ids, edges, completeness):
    """Single-linkage connected components over the kept edges.

    Returns clusters sorted by size (descending) then by representative ID,
    so that a re-run over the same pool yields the same vOTU labels.
    """
    neigh = {i: set() for i in ids}
    for a, b in edges:
        if a in neigh and b in neigh:
            neigh[a].add(b)
            neigh[b].add(a)

    seen = set()
    clusters = []
    for start in ids:
        if start in seen:
            continue
        comp = []
        stack = [start]
        while stack:
            node = stack.pop()
            if node in seen:
                continue
            seen.add(node)
            comp.append(node)
            # Sort neighbors before pushing: `neigh[node]` is a set, and
            # Python randomizes string hashing per process (PYTHONHASHSEED),
            # so an unsorted traversal visits neighbors in a different order
            # on every run. That changes cluster member order and, on a
            # completeness tie, the representative pick_representative()
            # returns -- silently reassigning vOTU_* labels across reruns
            # of the same pool.
            stack.extend(sorted(n for n in neigh[node] if n not in seen))
        clusters.append(comp)

    clusters.sort(key=lambda c: (-len(c), pick_representative(c, completeness)))
    return clusters


def assign_votu_ids(clusters, completeness=None):
    """Turn clusters into (votu_id, representative, member) rows."""
    completeness = completeness or {}
    rows = []
    for idx, members in enumerate(clusters, start=1):
        votu_id = f"vOTU_{idx:05d}"
        rep = pick_representative(members, completeness)
        for member in members:
            rows.append((votu_id, rep, member))
    return rows


def write_clusters(clusters, n_input, path, completeness=None):
    """Write vOTU_clusters.tsv, refusing to emit a catalog that collapsed nothing."""
    if n_input > 1 and len(clusters) == n_input:
        raise ClusteringCollapseError(
            f"no clusters formed: {n_input} input sequences produced "
            f"{len(clusters)} clusters. This is the signature of a skani "
            f"output-format or parser mismatch -- check that "
            f"`skani triangle --sparse` was used and that genome names are "
            f"read from columns 5 and 6."
        )
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w") as fh:
        fh.write("votu_id\trepresentative\tmember\n")
        for votu_id, rep, member in assign_votu_ids(clusters, completeness):
            fh.write(f"{votu_id}\t{rep}\t{member}\n")
    return len(clusters)


# ── geNomad x IDs do catalogo ─────────────────────────────────────────────
#
# Tres convencoes de nome se encontram aqui, e nenhuma e igual a outra:
#
#   geNomad (provirus)  "{contig}|provirus_{start}_{end}"
#   CheckV  (aparado)   "{contig}_{n}"
#   catalogo            "{source}|{id_que_chegou_ao_pool}"
#
# O que chega ao pool e a sequencia APARADA pelo CheckV, entao o ID do
# catalogo de um profago e "{source}|{contig}_{n}" -- enquanto a linha do
# geNomad para a mesma sequencia esta sob "{source}|{contig}|provirus_X_Y".
# Casar as duas exige normalizar as DUAS pontas: e exatamente o que
# `_load_catalog_completeness` ja faz do lado do CheckV, e o que faltava do
# lado do geNomad, deixando todo profago sem taxonomia geNomad em silencio.


def genomad_base_contig(seq_name):
    """'k141_10|provirus_1_5000' -> 'k141_10'. Sem sufixo, devolve como veio."""
    return (seq_name or "").split("|", 1)[0]


def resolve_genomad_key(catalog_id, genomad_index, known_by_source):
    """Chave do geNomad para um ID do catalogo, ou None.

    `genomad_index` e chaveado por "{source}|{contig_base}"; `known_by_source`
    da, por fonte, os contigs base que o geNomad conhece -- e contra esse
    conjunto que o sufixo de provirus do CheckV e desfeito, em vez de se
    adivinhar um delimitador (mesma disciplina de checkv_provirus.py).
    """
    if catalog_id in genomad_index:
        return catalog_id
    source, sep, bare = (catalog_id or "").partition("|")
    if not sep:
        return None
    known = known_by_source.get(source) or set()
    base, suffix_sep, suffix = bare.rpartition("_")
    if suffix_sep and suffix.isdigit() and base in known:
        return f"{source}|{base}"
    return None
