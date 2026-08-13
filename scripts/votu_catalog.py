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
