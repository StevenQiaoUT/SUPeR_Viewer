#!/usr/bin/env python3
"""Report the percentage of zero, one, and other values in an H5AD
expression matrix, overall and broken down by cell type.
Also verifies all values fall within [0, 1]."""

import sys
import anndata as ad
import scipy.sparse as sp
import numpy as np

# ── Dataset presets: ds_key → (h5ad path, cell type column) ──────────────
PRESETS = {
    "arabidopsis_nat":  ("/mnt/home/sqiao/h5ad_files/arabidopsis_nat.h5ad",  "label_major"),
    "rice":             ("/mnt/home/sqiao/h5ad_files/RiceOW.h5ad",           "CellAnnotation"),
    "at_root_rs":       ("/mnt/home/sqiao/h5ad_files/at_root_shahan.h5ad",   "Celltype"),
    "at_seed_martin":   ("/mnt/home/sqiao/h5ad_files/at_seed_martin.h5ad",   "level_2_annotation_timed"),
    "at_flower_lee":    ("/mnt/home/sqiao/h5ad_files/at_flower_lee.h5ad",    "CellType"),
    "at_silique_lee":   ("/mnt/home/sqiao/h5ad_files/at_silique_lee.h5ad",   "CellType"),
    "at_stem_lee":      ("/mnt/home/sqiao/h5ad_files/at_stem_lee.h5ad",      "CellType"),
    "at_shoot_zhang":   ("/mnt/home/sqiao/h5ad_files/at_shoot_zhang.h5ad",   "celltype_after"),
}

def get_data_array(X):
    """Return the raw values as a 1-D array (sparse-aware)."""
    if sp.issparse(X):
        return X.data
    return np.asarray(X).ravel()

def range_check(X, label=""):
    """Check whether all values fall in [0, 1]. Print PASS/FAIL."""
    vals = get_data_array(X)
    if len(vals) == 0:
        print(f"  {label}RANGE CHECK:  no data")
        return True
    vmin, vmax = float(vals.min()), float(vals.max())
    n_neg   = int(np.sum(vals < 0))
    n_above = int(np.sum(vals > 1))
    ok = n_neg == 0 and n_above == 0
    if ok:
        print(f"  {label}RANGE CHECK:  PASS  all values in [0, 1]   min={vmin:.6g}  max={vmax:.6g}")
    else:
        parts = []
        if n_neg:   parts.append(f"{n_neg:,} below 0")
        if n_above: parts.append(f"{n_above:,} above 1")
        print(f"  {label}RANGE CHECK:  FAIL  min={vmin:.6g}  max={vmax:.6g}  ({', '.join(parts)})")
    return ok

def count_categories(X):
    """Return (total, n_zero, n_one, n_other) for a matrix."""
    total = X.shape[0] * X.shape[1]
    if sp.issparse(X):
        data = X.data
        n_nonzero = len(data)
        n_one = int(np.sum(data == 1))
        n_zero = total - n_nonzero
    else:
        flat = np.asarray(X).ravel()
        n_zero = int(np.sum(flat == 0))
        n_one = int(np.sum(flat == 1))
    n_other = total - n_zero - n_one
    return total, n_zero, n_one, n_other

def report(path, ct_col):
    adata = ad.read_h5ad(path)

    # Use logcounts if X is empty
    if 'logcounts' in adata.layers:
        x = adata.X
        if x is None or (sp.issparse(x) and x.nnz == 0) or (not sp.issparse(x) and x.sum() == 0):
            print("  (X empty — using logcounts layer)")
            adata.X = adata.layers['logcounts']

    total, n_zero, n_one, n_other = count_categories(adata.X)
    pz  = n_zero  / total * 100
    p1  = n_one   / total * 100
    po  = n_other / total * 100

    print(f"\n{'=' * 78}")
    print(f"  {path}")
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes    column: {ct_col}")
    print(f"  Overall:  {pz:.2f}% zeros  |  {p1:.2f}% ones  |  {po:.2f}% other")
    print(f"{'=' * 78}")

    # ── Global range check ────────────────────────────────────────────────
    global_ok = range_check(adata.X, label="Overall ")

    if ct_col not in adata.obs.columns:
        print(f"  WARNING: column '{ct_col}' not found. Available: {list(adata.obs.columns)}")
        return

    ct_labels = adata.obs[ct_col].astype(str)
    unique_cts = sorted(ct_labels.unique())

    # ── Per-cell-type range check (only print failures) ───────────────────
    any_fail = False
    for ct in unique_cts:
        mask = (ct_labels == ct).values
        sub = adata.X[mask]
        vals = get_data_array(sub)
        if len(vals) == 0:
            continue
        vmin, vmax = float(vals.min()), float(vals.max())
        n_neg   = int(np.sum(vals < 0))
        n_above = int(np.sum(vals > 1))
        if n_neg or n_above:
            any_fail = True
            parts = []
            if n_neg:   parts.append(f"{n_neg:,} below 0")
            if n_above: parts.append(f"{n_above:,} above 1")
            print(f"  FAIL  {ct:<30}  min={vmin:.6g}  max={vmax:.6g}  ({', '.join(parts)})")
    if not any_fail:
        print(f"  Per-cell-type range: all PASS")

    # ── Per-cell-type distribution table ──────────────────────────────────
    rows = []
    for ct in unique_cts:
        mask = (ct_labels == ct).values
        sub = adata.X[mask]
        t, nz, n1, no = count_categories(sub)
        n_cells = int(mask.sum())
        rows.append((ct, n_cells, nz/t*100, n1/t*100, no/t*100))

    max_name = max(len(r[0]) for r in rows)
    header = f"  {'Cell Type':<{max_name}}  {'Cells':>7}  {'% Zeros':>9}  {'% Ones':>9}  {'% Other':>9}"
    print(header)
    print(f"  {'-' * (max_name + 42)}")

    other_vals = [r[4] for r in rows]
    mean_o = np.mean(other_vals)
    std_o  = np.std(other_vals)

    for ct, n_cells, pz, p1, po in rows:
        flag = " ***" if abs(po - mean_o) > 2 * std_o else ""
        print(f"  {ct:<{max_name}}  {n_cells:>7,}  {pz:>8.2f}%  {p1:>8.2f}%  {po:>8.2f}%{flag}")

    print(f"  {'-' * (max_name + 42)}")
    mean_z = np.mean([r[2] for r in rows])
    mean_1 = np.mean([r[3] for r in rows])
    print(f"  {'Mean':<{max_name}}  {'':>7}  {mean_z:>8.2f}%  {mean_1:>8.2f}%  {mean_o:>8.2f}%")
    print(f"  {'Std':<{max_name}}  {'':>7}  {'':>9}  {'':>9}  {std_o:>8.2f}%")
    cv = (std_o / mean_o * 100) if mean_o else 0
    print(f"  CV for '% Other' (std/mean):  {cv:.1f}%   {'(fairly even)' if cv < 30 else '(uneven)'}")
    print(f"  *** = >2 std from the mean\n")


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage:")
        print(f"  {sys.argv[0]} <preset_key|all>          # use preset datasets")
        print(f"  {sys.argv[0]} <file.h5ad> <ct_column>   # custom file")
        print(f"\nPresets: {', '.join(PRESETS.keys())}")
        sys.exit(1)

    arg = sys.argv[1]

    if arg == "all":
        for key, (path, col) in PRESETS.items():
            try:
                report(path, col)
            except Exception as e:
                print(f"\n  SKIPPED {key}: {e}")
    elif arg in PRESETS:
        path, col = PRESETS[arg]
        report(path, col)
    elif len(sys.argv) >= 3:
        report(sys.argv[1], sys.argv[2])
    else:
        print(f"Unknown preset '{arg}'. Available: {', '.join(PRESETS.keys())}, all")
        sys.exit(1)