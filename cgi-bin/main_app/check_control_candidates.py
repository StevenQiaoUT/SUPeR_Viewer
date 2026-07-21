#!/usr/bin/env python3
"""
Steven Qiao | BAR Lab | University of Toronto

Evaluate candidate "control" expression values for relative (log₂ ratio)
coloring in the SUPeR Viewer.

For datasets without an obvious biological control condition, we need a
single scalar per gene per cell type to serve as the denominator in:

    log₂(expression / control)

This script computes several candidate control metrics for every dataset,
both globally and per cell type, so we can compare and pick the most
sensible approach.

Candidates
----------
A) Global mean         – mean of ALL values (incl. zeros)
B) Global median       – median of ALL values (almost always 0 for sparse data)
C) Non-zero mean       – mean of values > 0 only
D) Non-zero median     – median of values > 0 only
E) Mean-of-CT-means    – average across per-cell-type pseudobulk means
                         (i.e. average over the values that actually color
                         the SVG — what the viewer already uses)
F) Percentiles of non-zero values (25th, 50th, 75th, 90th)

Per-cell-type breakdown reports how each cell type compares to each
candidate, which helps visualise what the relative map would look like.

Usage
-----
    python3 check_control_candidates.py all
    python3 check_control_candidates.py arabidopsis_nat
    python3 check_control_candidates.py /path/to/file.h5ad CellType
"""

from __future__ import annotations

import sys
import time
import warnings

import anndata as ad
import numpy as np
import scipy.sparse as sp

warnings.filterwarnings("ignore", category=FutureWarning)

# ── Dataset presets ──────────────────────────────────────────────────────
PRESETS = {
    "arabidopsis_nat": ("/mnt/home/sqiao/h5ad_files/arabidopsis_nat.h5ad", "label_major"),
    "rice":            ("/mnt/home/sqiao/h5ad_files/RiceOW.h5ad",          "CellAnnotation"),
    "at_root_rs":      ("/mnt/home/sqiao/h5ad_files/at_root_shahan.h5ad",  "Celltype"),
    "at_seed_martin":  ("/mnt/home/sqiao/h5ad_files/at_seed_martin.h5ad",  "level_2_annotation_timed"),
    "at_flower_lee":   ("/mnt/home/sqiao/h5ad_files/at_flower_lee.h5ad",   "CellType"),
    "at_silique_lee":  ("/mnt/home/sqiao/h5ad_files/at_silique_lee.h5ad",  "CellType"),
    "at_stem_lee":     ("/mnt/home/sqiao/h5ad_files/at_stem_lee.h5ad",     "CellType"),
    "at_shoot_zhang":  ("/mnt/home/sqiao/h5ad_files/at_shoot_zhang.h5ad",  "celltype_after"),
}


# ── Helpers ──────────────────────────────────────────────────────────────

def _dense_col(X, j):
    """Extract column j as a dense 1-D array."""
    col = X[:, j]
    if sp.issparse(col):
        return col.toarray().ravel()
    return np.asarray(col).ravel()


def _ensure_csc(X):
    """Convert to CSC for fast column slicing."""
    if sp.issparse(X) and not isinstance(X, sp.csc_matrix):
        return X.tocsc()
    return X


def _fmt(v, width=10):
    """Format a float for the table, or '—' if None/nan."""
    if v is None or (isinstance(v, float) and not np.isfinite(v)):
        return "—".center(width)
    return f"{v:.4f}".rjust(width)


# ── Core analysis ────────────────────────────────────────────────────────

def analyse(path: str, ct_col: str, ds_key: str = ""):
    t0 = time.perf_counter()
    adata = ad.read_h5ad(path)

    # logcounts fallback (same logic as check_sparsity / mini_app)
    if "logcounts" in adata.layers:
        x = adata.X
        if x is None or (sp.issparse(x) and x.nnz == 0) or (
            not sp.issparse(x) and np.asarray(x).sum() == 0
        ):
            print("  (X empty — using logcounts layer)")
            adata.X = adata.layers["logcounts"]

    X = _ensure_csc(adata.X)
    n_cells, n_genes = X.shape
    load_s = time.perf_counter() - t0

    print(f"\n{'=' * 88}")
    print(f"  {ds_key + '  ' if ds_key else ''}{path}")
    print(f"  {n_cells:,} cells × {n_genes:,} genes   ct_col: {ct_col}   loaded in {load_s:.1f}s")
    print(f"{'=' * 88}")

    # ── Global metrics (computed over ALL elements) ──────────────────────
    if sp.issparse(X):
        stored = X.data.copy()
        nnz = len(stored)
        total = n_cells * n_genes
        n_zero = total - nnz
        global_sum = float(stored.sum())
        global_mean = global_sum / total
        # For global median: if >50% are zeros, median = 0
        frac_nonzero = nnz / total
        if frac_nonzero <= 0.5:
            global_median = 0.0
        else:
            global_median = float(np.median(np.append(stored, np.zeros(n_zero))))
        nz_vals = stored[stored > 0]
    else:
        flat = np.asarray(X).ravel()
        total = len(flat)
        n_zero = int(np.sum(flat == 0))
        nnz = total - n_zero
        global_mean = float(flat.mean())
        global_median = float(np.median(flat))
        nz_vals = flat[flat > 0]

    sparsity = n_zero / total * 100

    nz_mean   = float(nz_vals.mean()) if len(nz_vals) else 0.0
    nz_median = float(np.median(nz_vals)) if len(nz_vals) else 0.0
    nz_p25    = float(np.percentile(nz_vals, 25)) if len(nz_vals) else 0.0
    nz_p75    = float(np.percentile(nz_vals, 75)) if len(nz_vals) else 0.0
    nz_p90    = float(np.percentile(nz_vals, 90)) if len(nz_vals) else 0.0
    nz_min    = float(nz_vals.min()) if len(nz_vals) else 0.0
    nz_max    = float(nz_vals.max()) if len(nz_vals) else 0.0

    print(f"\n  Sparsity: {sparsity:.2f}% zeros   "
          f"({nnz:,} non-zero / {total:,} total)")
    print(f"  Non-zero range: [{nz_min:.4g}, {nz_max:.4g}]")

    print(f"\n  ┌─────────────────────────────────────────────────────┐")
    print(f"  │  GLOBAL CONTROL CANDIDATES (across all genes)       │")
    print(f"  ├─────────────────────────────────────────────────────┤")
    print(f"  │  A) Global mean (incl. zeros)     {global_mean:>14.6f}   │")
    print(f"  │  B) Global median (incl. zeros)   {global_median:>14.6f}   │")
    print(f"  │  C) Non-zero mean                 {nz_mean:>14.6f}   │")
    print(f"  │  D) Non-zero median               {nz_median:>14.6f}   │")
    print(f"  │  ─── Non-zero percentiles ───                       │")
    print(f"  │     P25                           {nz_p25:>14.6f}   │")
    print(f"  │     P50 (= D above)               {nz_median:>14.6f}   │")
    print(f"  │     P75                           {nz_p75:>14.6f}   │")
    print(f"  │     P90                           {nz_p90:>14.6f}   │")
    print(f"  └─────────────────────────────────────────────────────┘")

    # ── Per-cell-type pseudobulk means ───────────────────────────────────
    if ct_col not in adata.obs.columns:
        print(f"\n  WARNING: column '{ct_col}' not found.")
        print(f"  Available: {list(adata.obs.columns)}")
        return

    ct_labels = adata.obs[ct_col].astype(str).values
    unique_cts = sorted(set(ct_labels))
    n_ct = len(unique_cts)

    # Pseudobulk mean per cell type (average over cells, then report per
    # cell type — same thing the SVG coloring uses)
    ct_means = {}        # cell type -> mean expression (across all genes)
    ct_nz_medians = {}   # cell type -> median of non-zero values
    ct_nz_means = {}     # cell type -> mean of non-zero values
    ct_n_cells = {}      # cell type -> number of cells

    for ct in unique_cts:
        mask = ct_labels == ct
        sub = X[mask]
        n = int(mask.sum())
        ct_n_cells[ct] = n

        if sp.issparse(sub):
            ct_total = sub.shape[0] * sub.shape[1]
            ct_sum = float(sub.data.sum()) if len(sub.data) else 0.0
            ct_means[ct] = ct_sum / ct_total if ct_total else 0.0
            nzv = sub.data[sub.data > 0]
        else:
            flat = np.asarray(sub).ravel()
            ct_means[ct] = float(flat.mean())
            nzv = flat[flat > 0]

        ct_nz_means[ct] = float(nzv.mean()) if len(nzv) else 0.0
        ct_nz_medians[ct] = float(np.median(nzv)) if len(nzv) else 0.0

    # Candidate E: mean of per-cell-type means
    mean_of_ct_means = float(np.mean(list(ct_means.values())))
    # Also: median of per-cell-type means
    median_of_ct_means = float(np.median(list(ct_means.values())))

    mean_of_ct_nz_medians = float(np.mean(list(ct_nz_medians.values())))

    print(f"\n  ┌─────────────────────────────────────────────────────┐")
    print(f"  │  AGGREGATED PER-CELL-TYPE CANDIDATES                │")
    print(f"  ├─────────────────────────────────────────────────────┤")
    print(f"  │  E) Mean of CT means              {mean_of_ct_means:>14.6f}   │")
    print(f"  │  F) Median of CT means             {median_of_ct_means:>13.6f}   │")
    print(f"  │  G) Mean of CT non-zero medians   {mean_of_ct_nz_medians:>14.6f}   │")
    print(f"  └─────────────────────────────────────────────────────┘")

    # ── Per-cell-type table ──────────────────────────────────────────────
    max_name = max(len(ct) for ct in unique_cts)
    max_name = max(max_name, 9)  # at least "Cell Type"

    print(f"\n  Per-cell-type breakdown:")
    hdr = (f"  {'Cell Type':<{max_name}}  {'Cells':>7}"
           f"  {'CT Mean':>10}  {'NZ Mean':>10}  {'NZ Median':>10}"
           f"  {'log₂(Mean/C)':>12}  {'log₂(Mean/D)':>12}")
    print(hdr)
    print(f"  {'-' * (max_name + 72)}")

    # For the log₂ ratio columns, use candidate C (global NZ mean) and
    # D (global NZ median) as the denominator to preview what relative
    # coloring would look like
    for ct in unique_cts:
        m = ct_means[ct]
        nzm = ct_nz_means[ct]
        nzmed = ct_nz_medians[ct]

        if m > 0 and nz_mean > 0:
            ratio_c = np.log2(m / nz_mean)
        else:
            ratio_c = None
        if m > 0 and nz_median > 0:
            ratio_d = np.log2(m / nz_median)
        else:
            ratio_d = None

        rc = _fmt(ratio_c, 12)
        rd = _fmt(ratio_d, 12)

        print(f"  {ct:<{max_name}}  {ct_n_cells[ct]:>7,}"
              f"  {m:>10.4f}  {nzm:>10.4f}  {nzmed:>10.4f}"
              f"  {rc}  {rd}")

    print(f"  {'-' * (max_name + 72)}")

    # ── Spread summary ───────────────────────────────────────────────────
    ct_mean_arr = np.array(list(ct_means.values()))
    ct_nzm_arr  = np.array(list(ct_nz_medians.values()))

    def spread_stats(arr, label):
        mn, mx = arr.min(), arr.max()
        fold = mx / mn if mn > 0 else float("inf")
        cv = float(np.std(arr) / np.mean(arr) * 100) if np.mean(arr) > 0 else 0
        print(f"  {label}:  min={mn:.4f}  max={mx:.4f}  "
              f"fold-range={fold:.1f}×  CV={cv:.1f}%")

    print(f"\n  Spread across cell types:")
    spread_stats(ct_mean_arr, "CT means        ")
    spread_stats(ct_nzm_arr,  "CT NZ medians   ")

    # ── Gene-level sampling ──────────────────────────────────────────────
    # Pick 5 random genes and show per-gene control candidates to give a
    # sense of how much the "best" control varies gene-to-gene
    rng = np.random.default_rng(42)
    sample_n = min(5, n_genes)
    sample_idx = rng.choice(n_genes, size=sample_n, replace=False)
    sample_idx.sort()

    print(f"\n  Gene-level sample (random {sample_n} genes):")
    ghdr = (f"  {'Gene':<20}  {'Gene Mean':>10}  {'Gene NZ Mean':>12}"
            f"  {'Gene NZ Med':>11}  {'% Expressing':>12}")
    print(ghdr)
    print(f"  {'-' * 72}")

    for gi in sample_idx:
        col = _dense_col(X, gi)
        gene_id = adata.var_names[gi]
        gene_mean = float(col.mean())
        nzv = col[col > 0]
        gene_nz_mean = float(nzv.mean()) if len(nzv) else 0.0
        gene_nz_med = float(np.median(nzv)) if len(nzv) else 0.0
        pct_expr = len(nzv) / len(col) * 100

        print(f"  {gene_id:<20}  {gene_mean:>10.4f}  {gene_nz_mean:>12.4f}"
              f"  {gene_nz_med:>11.4f}  {pct_expr:>11.1f}%")

    print()


# ── CLI ──────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage:")
        print(f"  {sys.argv[0]} <preset_key|all>")
        print(f"  {sys.argv[0]} <file.h5ad> <ct_column>")
        print(f"\nPresets: {', '.join(PRESETS.keys())}")
        sys.exit(1)

    arg = sys.argv[1]

    if arg == "all":
        for key, (path, col) in PRESETS.items():
            try:
                analyse(path, col, ds_key=key)
            except Exception as e:
                print(f"\n  SKIPPED {key}: {e}")
    elif arg in PRESETS:
        path, col = PRESETS[arg]
        analyse(path, col, ds_key=arg)
    elif len(sys.argv) >= 3:
        analyse(sys.argv[1], sys.argv[2])
    else:
        print(f"Unknown preset '{arg}'. Available: {', '.join(PRESETS.keys())}, all")
        sys.exit(1)
