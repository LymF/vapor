#!/usr/bin/env python3
"""
vcontact3_network_layout.py — Pre-compute 2D spectral layout for the
vConTACT3 protein-sharing network (HMMprofile_mLayer_shared_ntw.h5).

Usage:
    python3 vcontact3_network_layout.py <h5_file> <final_assignments_csv> <output_json>

Output JSON schema:
    {"nodes": [{"id", "x", "y", "cluster", "class", "is_novel", "family", "genus"}],
     "edges": [{"source", "target", "weight"}]}
"""
import sys, os, re, json, csv
import numpy as np
import h5py
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import laplacian as sp_laplacian
from scipy.sparse.linalg import eigsh


def main(h5_path, csv_path, out_json, min_edge_weight=1):
    # ── Load sparse matrix ────────────────────────────────────────────────────
    print(f"[vc3_layout] Loading {h5_path}", flush=True)
    with h5py.File(h5_path, 'r') as f:
        keys = list(f.keys())
        print(f"[vc3_layout] h5 fields: {keys}", flush=True)
        data    = f['data'][:]
        indices = f['indices'][:]
        indptr  = f['indptr'][:]
        shape   = tuple(int(x) for x in f['shape'][:])
        # Genome names may be stored under various field names
        genome_names = None
        for fn in ('index', 'Index', 'genomes', 'names', 'node_labels',
                   'genome_names', 'labels'):
            if fn in f:
                raw = f[fn][:]
                genome_names = [n.decode('utf-8') if isinstance(n, bytes) else str(n)
                                for n in raw]
                print(f"[vc3_layout] Genome names from field '{fn}': {len(genome_names)}", flush=True)
                break

    mat = csr_matrix((data.astype(np.float64), indices, indptr), shape=shape)
    print(f"[vc3_layout] Matrix {shape}, nnz={mat.nnz}", flush=True)

    if genome_names is None:
        print("[vc3_layout] WARNING: genome names not found in h5; cannot filter to query genomes", flush=True)
        json.dump({"nodes": [], "edges": []}, open(out_json, 'w'))
        return

    # ── Load query genome assignments from final_assignments.csv ─────────────
    query_info = {}   # genome_id → {cluster, class, is_novel, family, genus}
    if csv_path and os.path.exists(csv_path):
        with open(csv_path, newline='', encoding='utf-8') as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                genome = row.get('Genome', row.get('genome', '')).strip()
                if not genome:
                    continue
                ref_val = row.get('Reference', row.get('reference', ''))
                if str(ref_val).lower() in ('true', '1', 'yes'):
                    continue
                fam_pred = (row.get('family_prediction', '') or '').strip()
                gen_pred = (row.get('genus_prediction',  '') or '').strip()
                cls_pred = (row.get('class_prediction',  '') or '').strip()
                # Determine cluster: "novel_family_7_of_Caudoviricetes" → "family_7"
                is_novel = fam_pred.lower().startswith('novel_')
                if is_novel:
                    m = re.match(r'novel_(family_\d+)', fam_pred, re.I)
                    cluster = m.group(1) if m else 'novel'
                elif fam_pred and fam_pred.lower() not in ('', 'singleton', 'nd', 'none',
                                                            'unclassified'):
                    cluster = fam_pred
                else:
                    cluster = 'singleton'
                query_info[genome] = {
                    'cluster':  cluster,
                    'class':    cls_pred or 'Unknown',
                    'is_novel': is_novel,
                    'family':   fam_pred,
                    'genus':    gen_pred,
                }
    print(f"[vc3_layout] Query genomes (from CSV): {len(query_info)}", flush=True)

    # ── Map genome names → matrix indices ────────────────────────────────────
    name_to_idx = {n: i for i, n in enumerate(genome_names)}
    q_genomes   = [g for g in query_info if g in name_to_idx]
    q_idx       = [name_to_idx[g] for g in q_genomes]

    print(f"[vc3_layout] Query genomes found in matrix: {len(q_idx)}", flush=True)

    if len(q_idx) < 3:
        print("[vc3_layout] Too few query genomes for layout — saving empty JSON", flush=True)
        json.dump({"nodes": [], "edges": []},
                  open(out_json, 'w', encoding='utf-8'))
        return

    # ── Extract query–query submatrix ─────────────────────────────────────────
    q_sub = mat[q_idx, :][:, q_idx].astype(np.float64)
    # Symmetrize (take mean of both directions — the matrix may be directed)
    q_sub = (q_sub + q_sub.T) * 0.5
    n = len(q_idx)
    print(f"[vc3_layout] Submatrix {n}×{n}, nnz={q_sub.nnz}", flush=True)

    # ── Compute 2D spectral layout ─────────────────────────────────────────────
    # Use only edges above threshold to get a cleaner layout signal
    q_filt = q_sub.copy()
    if min_edge_weight > 1:
        q_filt.data[q_filt.data < min_edge_weight] = 0.0
        q_filt.eliminate_zeros()
        if q_filt.nnz < n:
            q_filt = q_sub  # fall back to unfiltered if too sparse

    L = sp_laplacian(q_filt, normed=True)
    k = min(4, n - 1)
    try:
        eigenvalues, eigenvectors = eigsh(L, k=k, which='SM', tol=1e-4, maxiter=2000)
        order = np.argsort(eigenvalues)
        eigenvectors = eigenvectors[:, order]
        # Skip the trivial constant eigenvector (eigenvalue ≈ 0); use cols 1 and 2
        xy = eigenvectors[:, 1:3] if k >= 3 else eigenvectors[:, :2]
    except Exception as exc:
        print(f"[vc3_layout] eigsh failed ({exc}) — using random placement", flush=True)
        rng = np.random.default_rng(42)
        xy  = rng.standard_normal((n, 2))

    # Normalise coordinates to [-1, 1]
    for col in range(xy.shape[1]):
        mn, mx = xy[:, col].min(), xy[:, col].max()
        if mx - mn > 0:
            xy[:, col] = (xy[:, col] - mn) / (mx - mn) * 2.0 - 1.0

    # ── Build node list ───────────────────────────────────────────────────────
    nodes = [
        {
            'id':       g,
            'x':        float(xy[i, 0]),
            'y':        float(xy[i, 1]),
            'cluster':  query_info[g]['cluster'],
            'class':    query_info[g]['class'],
            'is_novel': query_info[g]['is_novel'],
            'family':   query_info[g]['family'],
            'genus':    query_info[g]['genus'],
        }
        for i, g in enumerate(q_genomes)
    ]

    # ── Build edge list ───────────────────────────────────────────────────────
    q_coo = q_sub.tocoo()
    seen  = set()
    edges = []
    for r, c, w in zip(q_coo.row, q_coo.col, q_coo.data):
        if r >= c:
            continue
        if w < min_edge_weight:
            continue
        key = (r, c)
        if key in seen:
            continue
        seen.add(key)
        edges.append({
            'source': q_genomes[r],
            'target': q_genomes[c],
            'weight': round(float(w), 2),
        })

    print(f"[vc3_layout] Nodes: {len(nodes)}, Edges: {len(edges)}", flush=True)

    with open(out_json, 'w', encoding='utf-8') as fout:
        json.dump({"nodes": nodes, "edges": edges}, fout, separators=(',', ':'))
    print(f"[vc3_layout] Saved → {out_json}", flush=True)


if __name__ == '__main__':
    if len(sys.argv) < 4:
        sys.exit("Usage: vcontact3_network_layout.py <h5> <csv> <out_json>")
    main(sys.argv[1], sys.argv[2], sys.argv[3])
