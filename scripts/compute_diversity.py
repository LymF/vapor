#!/usr/bin/env python3
"""
compute_diversity.py — Alpha + beta + combined diversity from CoverM abundance tables.

Inputs (via Snakemake):
  snakemake.input.viral_tables  — list of viral_abundance.tsv (one per sample)
  snakemake.input.prok_tables   — list of prok_abundance.tsv  (one per sample)

Outputs:
  diversity/alpha_diversity.tsv      — Shannon, Simpson, Chao1, Richness per sample×domain
  diversity/beta_pcoord_viral.tsv    — PCoA coords (Bray-Curtis, viral)
  diversity/beta_pcoord_prok.tsv     — PCoA coords (Bray-Curtis, prokaryotic)
  diversity/beta_pcoord_combined.tsv — PCoA coords (viral+prok merged features)
  diversity/procrustes_coords.tsv    — Procrustes alignment of viral vs prok ordinations
  diversity/diversity_done.txt       — sentinel
"""

import csv
import math
import os
import sys
from collections import defaultdict

import numpy as np

# scipy imports — conditionally used
try:
    from scipy.spatial import procrustes as scipy_procrustes
    from scipy.linalg import eigh
    SCIPY_OK = True
except ImportError:
    SCIPY_OK = False
    print("[compute_diversity] WARNING: scipy not available — beta diversity skipped",
          file=sys.stderr)


# ──────────────────────────────────────────────────────────────────────────────
#  I/O helpers
# ──────────────────────────────────────────────────────────────────────────────

def load_abundance_table(path, abundance_col=None):
    """
    Read a CoverM abundance TSV.
    Returns (feature_ids, {feature: abundance}) where abundance is numeric.
    First column = feature name; remaining columns = sample abundances.
    If multiple sample columns exist, sums them (multi-sample BAM edge case).
    """
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return [], {}
    rows = []
    with open(path) as f:
        reader = csv.DictReader(f, delimiter="\t")
        headers = reader.fieldnames or []
        feature_col = headers[0] if headers else "Contig"
        value_cols = [h for h in headers[1:] if abundance_col is None
                      or abundance_col.lower() in h.lower()]
        if not value_cols and len(headers) > 1:
            value_cols = headers[1:]
        for row in reader:
            feat = row.get(feature_col, "").strip()
            if not feat or feat == "unmapped":
                continue
            vals = []
            for c in value_cols:
                try:
                    vals.append(float(row.get(c, 0) or 0))
                except ValueError:
                    pass
            rows.append((feat, sum(vals) if vals else 0.0))
    ids   = [r[0] for r in rows]
    abund = {r[0]: r[1] for r in rows}
    return ids, abund


def build_matrix(tables, samples, method_hint="rpkm"):
    """
    Build abundance matrix: rows = features, cols = samples.
    Returns (sorted_feature_ids, {sample: {feature: abundance}}).
    """
    sample_data = {}
    all_features = set()
    for path, sample in zip(tables, samples):
        _, abund = load_abundance_table(path, abundance_col=method_hint)
        sample_data[sample] = abund
        all_features.update(abund.keys())
    features = sorted(all_features)
    return features, sample_data


# ──────────────────────────────────────────────────────────────────────────────
#  Alpha diversity
# ──────────────────────────────────────────────────────────────────────────────

def shannon(abundances):
    total = sum(abundances)
    if total == 0:
        return 0.0
    h = 0.0
    for a in abundances:
        if a > 0:
            p = a / total
            h -= p * math.log(p)
    return round(h, 6)


def simpson(abundances):
    total = sum(abundances)
    if total <= 1:
        return 0.0
    d = sum(a * (a - 1) for a in abundances)
    return round(1.0 - d / (total * (total - 1)), 6)


def chao1(abundances):
    """Chao1 richness estimator. Requires integer counts; floors floats."""
    counts = [int(a) for a in abundances if a > 0]
    n_obs = len(counts)
    if n_obs == 0:
        return 0
    f1 = sum(1 for c in counts if c == 1)
    f2 = sum(1 for c in counts if c == 2)
    if f2 == 0:
        chao = n_obs + f1 * (f1 - 1) / 2
    else:
        chao = n_obs + (f1 ** 2) / (2 * f2)
    return round(chao, 1)


def compute_alpha(features, sample_data, domain):
    rows = []
    for sample, abund in sample_data.items():
        vals = [abund.get(f, 0.0) for f in features]
        richness = sum(1 for v in vals if v > 0)
        rows.append({
            "sample":   sample,
            "domain":   domain,
            "richness": richness,
            "shannon":  shannon(vals),
            "simpson":  simpson(vals),
            "chao1":    chao1(vals),
        })
    return rows


# ──────────────────────────────────────────────────────────────────────────────
#  Beta diversity — Bray-Curtis + PCoA
# ──────────────────────────────────────────────────────────────────────────────

def bray_curtis_matrix(features, sample_data):
    """Return (sample_list, distance_matrix) — Bray-Curtis dissimilarity."""
    samples = sorted(sample_data.keys())
    if len(samples) < 2:
        return samples, np.zeros((len(samples), len(samples)))
    mat = np.array([
        [sample_data[s].get(f, 0.0) for f in features]
        for s in samples
    ], dtype=float)
    n = len(samples)
    dist = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            num = np.sum(np.abs(mat[i] - mat[j]))
            den = np.sum(mat[i]) + np.sum(mat[j])
            d = num / den if den > 0 else 0.0
            dist[i, j] = dist[j, i] = d
    return samples, dist


def pcoa(dist_matrix):
    """
    PCoA (Classical MDS) from a symmetric distance matrix.
    Returns (eigenvalues, coordinates) — rows = samples, cols = PC axes.
    """
    n = dist_matrix.shape[0]
    if n < 2:
        return np.array([0.0]), np.zeros((n, 1))
    # Double centering
    d2 = dist_matrix ** 2
    row_mean = d2.mean(axis=1, keepdims=True)
    col_mean = d2.mean(axis=0, keepdims=True)
    grand_mean = d2.mean()
    g = -0.5 * (d2 - row_mean - col_mean + grand_mean)
    # Eigendecomposition (symmetric)
    eigenvalues, eigenvectors = eigh(g)
    # Sort descending
    idx = np.argsort(eigenvalues)[::-1]
    eigenvalues = eigenvalues[idx]
    eigenvectors = eigenvectors[:, idx]
    # Keep only positive eigenvalues
    pos = eigenvalues > 1e-10
    coords = eigenvectors[:, pos] * np.sqrt(np.maximum(eigenvalues[pos], 0))
    return eigenvalues[pos], coords


def pcoa_to_tsv(samples, dist_matrix, domain, n_axes=5):
    """Return list of dicts for TSV output (up to n_axes PCoA axes)."""
    if len(samples) < 2:
        return [{"sample": s, "domain": domain, "PC1": 0, "PC2": 0} for s in samples]
    evals, coords = pcoa(dist_matrix)
    rows = []
    n_out = min(n_axes, coords.shape[1])
    for i, s in enumerate(samples):
        row = {"sample": s, "domain": domain}
        for ax in range(n_out):
            row[f"PC{ax+1}"] = round(float(coords[i, ax]), 6)
        # Variance explained
        total = max(evals.sum(), 1e-9)
        for ax in range(n_out):
            row[f"PC{ax+1}_var"] = round(float(evals[ax]) / total * 100, 2)
        rows.append(row)
    return rows


# ──────────────────────────────────────────────────────────────────────────────
#  Procrustes — co-variation between viral and prokaryotic ordinations
# ──────────────────────────────────────────────────────────────────────────────

def run_procrustes(viral_pcoa_rows, prok_pcoa_rows):
    """
    Align viral and prokaryotic PCoA coordinates via Procrustes.
    Only uses samples present in both ordinations.
    Returns list of dicts with aligned coordinates + disparity score.
    """
    if not SCIPY_OK:
        return [], None
    v_dict = {r["sample"]: r for r in viral_pcoa_rows}
    p_dict = {r["sample"]: r for r in prok_pcoa_rows}
    common = sorted(set(v_dict) & set(p_dict))
    if len(common) < 3:
        print(f"[compute_diversity] Procrustes skipped: only {len(common)} common samples",
              file=sys.stderr)
        return [], None
    pc_cols = [k for k in v_dict[common[0]] if k.startswith("PC") and not k.endswith("_var")]
    pc_cols = sorted(pc_cols, key=lambda x: int(x[2:]))[:2]
    if len(pc_cols) < 2:
        return [], None
    mat_v = np.array([[v_dict[s].get(c, 0.0) for c in pc_cols] for s in common])
    mat_p = np.array([[p_dict[s].get(c, 0.0) for c in pc_cols] for s in common])
    try:
        _, mat_v_t, disparity = scipy_procrustes(mat_p, mat_v)
    except Exception as e:
        print(f"[compute_diversity] Procrustes error: {e}", file=sys.stderr)
        return [], None
    rows = []
    for i, s in enumerate(common):
        rows.append({
            "sample":     s,
            "viral_PC1":  round(float(mat_v_t[i, 0]), 6),
            "viral_PC2":  round(float(mat_v_t[i, 1]), 6),
            "prok_PC1":   round(float(mat_p[i, 0]), 6),
            "prok_PC2":   round(float(mat_p[i, 1]), 6),
            "disparity":  round(float(disparity), 6),
        })
    print(f"[compute_diversity] Procrustes disparity: {disparity:.4f}", file=sys.stderr)
    return rows, float(disparity)


# ──────────────────────────────────────────────────────────────────────────────
#  Write TSV
# ──────────────────────────────────────────────────────────────────────────────

def write_tsv(path, rows, fieldnames=None):
    if not rows:
        open(path, "w").close()
        return
    if fieldnames is None:
        fieldnames = list(rows[0].keys())
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


# ──────────────────────────────────────────────────────────────────────────────
#  MAIN
# ──────────────────────────────────────────────────────────────────────────────

def main():
    viral_tables = list(snakemake.input.viral_tables)
    prok_tables  = list(snakemake.input.prok_tables)
    samples      = list(snakemake.params.samples)
    method       = snakemake.params.method
    out_alpha    = snakemake.output.alpha
    out_pcoa_v   = snakemake.output.pcoa_v
    out_pcoa_p   = snakemake.output.pcoa_p
    out_pcoa_c   = snakemake.output.pcoa_c
    out_procrust = snakemake.output.procrust
    out_done     = snakemake.output.done
    log_path     = snakemake.log[0]

    os.makedirs(os.path.dirname(out_alpha), exist_ok=True)
    log_f = open(log_path, "w")

    def log(msg):
        print(msg, file=log_f)
        print(msg, file=sys.stderr)

    log(f"[compute_diversity] Samples: {samples}")
    log(f"[compute_diversity] Abundance method: {method}")

    # ── Load abundance matrices ───────────────────────────────────────
    vfeats, v_data = build_matrix(viral_tables, samples, method)
    pfeats, p_data = build_matrix(prok_tables,  samples, method)
    log(f"[compute_diversity] Viral features: {len(vfeats)}")
    log(f"[compute_diversity] Prokaryotic features: {len(pfeats)}")

    # ── Alpha diversity ───────────────────────────────────────────────
    alpha_rows = []
    alpha_rows += compute_alpha(vfeats, v_data, "viral")
    alpha_rows += compute_alpha(pfeats, p_data, "prokaryotic")

    # Combined: merge feature sets
    all_feats = sorted(set(vfeats) | set(pfeats))
    comb_data = {}
    for s in samples:
        comb_data[s] = {}
        comb_data[s].update(v_data.get(s, {}))
        comb_data[s].update(p_data.get(s, {}))
    alpha_rows += compute_alpha(all_feats, comb_data, "combined")
    write_tsv(out_alpha, alpha_rows,
              ["sample", "domain", "richness", "shannon", "simpson", "chao1"])
    log(f"[compute_diversity] Alpha: {len(alpha_rows)} rows")

    if not SCIPY_OK:
        # Write empty beta outputs so Snakemake targets exist
        for p in [out_pcoa_v, out_pcoa_p, out_pcoa_c, out_procrust]:
            open(p, "w").close()
        open(out_done, "w").close()
        log_f.close()
        return

    # ── Beta diversity — viral ────────────────────────────────────────
    v_samples, v_dist = bray_curtis_matrix(vfeats, v_data)
    v_pcoa_rows = pcoa_to_tsv(v_samples, v_dist, "viral")
    write_tsv(out_pcoa_v, v_pcoa_rows)
    log(f"[compute_diversity] Viral PCoA: {len(v_pcoa_rows)} samples")

    # ── Beta diversity — prokaryotic ──────────────────────────────────
    p_samples, p_dist = bray_curtis_matrix(pfeats, p_data)
    p_pcoa_rows = pcoa_to_tsv(p_samples, p_dist, "prokaryotic")
    write_tsv(out_pcoa_p, p_pcoa_rows)
    log(f"[compute_diversity] Prok PCoA: {len(p_pcoa_rows)} samples")

    # ── Beta diversity — combined ─────────────────────────────────────
    c_samples, c_dist = bray_curtis_matrix(all_feats, comb_data)
    c_pcoa_rows = pcoa_to_tsv(c_samples, c_dist, "combined")
    write_tsv(out_pcoa_c, c_pcoa_rows)
    log(f"[compute_diversity] Combined PCoA: {len(c_pcoa_rows)} samples")

    # ── Procrustes ────────────────────────────────────────────────────
    proc_rows, disparity = run_procrustes(v_pcoa_rows, p_pcoa_rows)
    write_tsv(out_procrust, proc_rows)
    if disparity is not None:
        log(f"[compute_diversity] Procrustes disparity: {disparity:.4f} "
            f"(0=perfect co-variation, 1=no correlation)")

    open(out_done, "w").close()
    log_f.close()
    print(f"[compute_diversity] Done → {os.path.dirname(out_alpha)}", file=sys.stderr)


if __name__ == "__main__":
    main()
