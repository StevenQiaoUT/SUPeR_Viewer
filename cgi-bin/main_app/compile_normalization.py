#!/usr/bin/env python3
"""For each dataset: identify what X contains, count true zeros,
and report data shape for discussion with Nick."""

import sys
import anndata as ad
import scipy.sparse as sp
import numpy as np

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

def classify_dataset(key, path, ct_col):
    adata = ad.read_h5ad(path)
    layers = list(adata.layers.keys()) if adata.layers else []

    # ── Determine active matrix ──────────────────────────────────────────
    using = "X"
    if 'logcounts' in layers:
        x = adata.X
        if x is None or (sp.issparse(x) and x.nnz == 0):
            using = "logcounts"
            adata.X = adata.layers['logcounts']

    X = adata.X
    vals = X.data if sp.issparse(X) else np.asarray(X).ravel()
    total_entries = X.shape[0] * X.shape[1]
    n_nonzero = len(vals) if sp.issparse(X) else int(np.count_nonzero(vals))
    n_matrix_zeros = total_entries - n_nonzero

    vmin = float(vals.min()) if len(vals) > 0 else 0.0
    vmax = float(vals.max()) if len(vals) > 0 else 0.0
    is_integer = np.allclose(vals[:10000], np.round(vals[:10000]))  # sample
    has_negatives = vmin < 0
    pct_ones = np.sum(vals == 1) / len(vals) * 100 if len(vals) else 0

    # ── Match X[0].sum() to obs columns ──────────────────────────────────
    row0_sum = float(X[0].sum())
    matches = {}
    for col in ['nCount_RNA', 'nCount_SCT']:
        if col in adata.obs.columns:
            ref = float(adata.obs[col].iloc[0])
            matches[col] = (ref, abs(row0_sum - ref) < 0.5)

    # ── Classify normalization ───────────────────────────────────────────
    if has_negatives:
        norm_type = "SCALED (z-scored / SCT residuals)"
        true_zeros_in_X = "UNRELIABLE — 0 means average, not absence"
        # Check if raw counts exist in a layer
        raw_layer = None
        for lname in ['counts', 'raw', 'spliced']:
            if lname in layers:
                raw_layer = lname
                break
        if raw_layer:
            raw_X = adata.layers[raw_layer]
            raw_total = raw_X.shape[0] * raw_X.shape[1]
            if sp.issparse(raw_X):
                raw_nz = raw_X.nnz
            else:
                raw_nz = int(np.count_nonzero(raw_X))
            raw_zeros = raw_total - raw_nz
            true_zeros_in_X = f"{raw_zeros:,} true zeros found in '{raw_layer}' layer ({raw_zeros/raw_total*100:.2f}%)"
        else:
            true_zeros_in_X += f" — no raw layer available (layers: {layers})"
    elif is_integer and pct_ones > 1:
        norm_type = "RAW COUNTS (UMI)"
        true_zeros_in_X = f"{n_matrix_zeros:,} ({n_matrix_zeros/total_entries*100:.2f}%) — all zeros are true zeros"
    elif vmax <= 1 and pct_ones < 0.1:
        norm_type = "SCTRANSFORM corrected (continuous, [0,1])"
        true_zeros_in_X = f"{n_matrix_zeros:,} ({n_matrix_zeros/total_entries*100:.2f}%) — zeros = true zeros (corrected counts)"
    elif vmin >= 0 and vmax > 1 and not is_integer:
        if pct_ones > 1:
            norm_type = "RAW or LIGHTLY NORMALIZED COUNTS"
        else:
            norm_type = "LOG1P-TRANSFORMED (log1p(0)=0, zeros preserved)"
        true_zeros_in_X = f"{n_matrix_zeros:,} ({n_matrix_zeros/total_entries*100:.2f}%) — all zeros are true zeros"
    else:
        norm_type = "UNKNOWN"
        true_zeros_in_X = f"{n_matrix_zeros:,} ({n_matrix_zeros/total_entries*100:.2f}%) — interpretation unclear"

    # ── Print report ─────────────────────────────────────────────────────
    print(f"\n{'=' * 78}")
    print(f"  {key}")
    print(f"  {path}")
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"{'=' * 78}")
    print(f"  Using:          {using}")
    print(f"  Layers:         {layers if layers else '(none)'}")
    print(f"  Value range:    min={vmin:.6g}   max={vmax:.6g}")
    print(f"  Integer values: {'yes' if is_integer else 'no'}    % exact 1s (of nonzero): {pct_ones:.2f}%")
    for col, (ref, ok) in matches.items():
        mark = "MATCH" if ok else "no match"
        print(f"  X[0].sum()={row0_sum:.1f}  vs  {col}[0]={ref:.1f}  → {mark}")
    print(f"  ──────────────────────────────────────────────────────────────")
    print(f"  NORMALIZATION:  {norm_type}")
    print(f"  TRUE ZEROS:     {true_zeros_in_X}")
    print(f"  Has negatives:  {'YES — zeros are NOT true zeros' if has_negatives else 'no'}")

    return {
        'key': key, 'cells': adata.n_obs, 'genes': adata.n_vars,
        'norm': norm_type, 'vmin': vmin, 'vmax': vmax,
        'matrix_zeros': n_matrix_zeros, 'total': total_entries,
    }


if __name__ == "__main__":
    if len(sys.argv) < 2:
        print(f"Usage:")
        print(f"  {sys.argv[0]} <preset_key|all>")
        print(f"  {sys.argv[0]} <file.h5ad> <ct_column>")
        print(f"\nPresets: {', '.join(PRESETS.keys())}")
        sys.exit(1)

    arg = sys.argv[1]
    results = []

    if arg == "all":
        for key, (path, col) in PRESETS.items():
            try:
                r = classify_dataset(key, path, col)
                results.append(r)
            except Exception as e:
                print(f"\n  SKIPPED {key}: {e}")
    elif arg in PRESETS:
        path, col = PRESETS[arg]
        results.append(classify_dataset(arg, path, col))
    elif len(sys.argv) >= 3:
        results.append(classify_dataset(sys.argv[1], sys.argv[1], sys.argv[2]))
    else:
        print(f"Unknown preset '{arg}'.")
        sys.exit(1)

    # ── Summary table ────────────────────────────────────────────────────
    if len(results) > 1:
        print(f"\n\n{'=' * 90}")
        print(f"  SUMMARY FOR NICK")
        print(f"{'=' * 90}")
        max_k = max(len(r['key']) for r in results)
        print(f"  {'Dataset':<{max_k}}  {'Cells':>8}  {'Genes':>7}  {'min':>10}  {'max':>10}  {'% Zero':>8}  Normalization")
        print(f"  {'-' * (max_k + 70)}")
        for r in results:
            pz = r['matrix_zeros'] / r['total'] * 100
            print(f"  {r['key']:<{max_k}}  {r['cells']:>8,}  {r['genes']:>7,}  {r['vmin']:>10.4f}  {r['vmax']:>10.4f}  {pz:>7.2f}%  {r['norm']}")
