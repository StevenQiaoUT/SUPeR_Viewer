#!/mnt/home/sqiao/venv/bin/python3
"""
generate_pseudobulk_dumps.py
=============================
Reads H5AD files for all SUPeR Viewer datasets, computes per-gene
per-cell-type average and standard deviation, and writes MySQL dump files
compatible with the BAR eFP structure.

Output files (one per dataset):
    arabidopsis_NIE_pseudobulk_dump.sql     <- Illouz-Eliaz drought/leaf
    rice_OW_pseudobulk_dump.sql             <- Robertson stress
    arabidopsis_root_rs_pseudobulk_dump.sql <- Shahan root
    arabidopsis_seed_martin_pseudobulk_dump.sql
    arabidopsis_flower_lee_pseudobulk_dump.sql
    arabidopsis_silique_lee_pseudobulk_dump.sql
    arabidopsis_stem_lee_pseudobulk_dump.sql
    arabidopsis_seed_0d_lee_pseudobulk_dump.sql
    arabidopsis_rosette_21d_lee_pseudobulk_dump.sql
    arabidopsis_shoot_zhang_pseudobulk_dump.sql
    arabidopsis_rosette_30d_lee_pseudobulk_dump.sql
    arabidopsis_seedling_3d_lee_pseudobulk_dump.sql
    arabidopsis_seedling_6d_lee_pseudobulk_dump.sql
    arabidopsis_seedling_12d_lee_pseudobulk_dump.sql

Usage:
    Option 1 - datasets already registered in dataset_config.yaml:
        python generate_pseudobulk_dumps.py                      # all datasets
        python generate_pseudobulk_dumps.py --ds rice            # one dataset
        python generate_pseudobulk_dumps.py --ds at_root rice    # several datasets
        python generate_pseudobulk_dumps.py --ds at_flower --outdir /path/to/output/

    Option 2 - any H5AD file, straight from the command line (no YAML editing):
        python generate_pseudobulk_dumps.py --h5ad lee_2025_at_root.h5ad
        python generate_pseudobulk_dumps.py --h5ad robertson_2026_rice.h5ad \\
            --cell-type-col CellAnnotation --cell-type-col2 Condition \\
            --db-name rice_OW_pseudobulk --outdir ~/dumps

        Defaults for Option 2:
            --cell-type-col  CellType  (the column rds_to_h5ad.r always writes)
            --db-name        <h5ad file name>_pseudobulk
            --cell-type-col2 / --gene-id-col: not used
            --use-logcounts: off

    Quick test on a single gene (works with both options):
        python generate_pseudobulk_dumps.py --ds rice -n 1 --no-std

    Available --ds keys:
        arabidopsis, rice, at_root, at_seed, at_flower, at_silique, at_stem,
        at_seed_0d_lee, at_rosette_21d_lee, at_shoot_zhang, at_rosette_30d_lee,
        at_seedling_3d_lee, at_seedling_6d_lee, at_seedling_12d_lee

Dataset-specific notes
----------------------
arabidopsis:
    Keys are the label_majorXcondition obs values (e.g. "W0_Guard").
    Synthetic "Phloem average" rows are appended for each condition
    (average of Phloem companion + Phloem Parenchyma means).

rice:
    Keys are "{Condition}_{CellAnnotation}" — Condition first so that
    splitting on the first underscore always gives (cond, ct), matching
    the lookup pattern in mini_app.py.

at_root:
    The raw obs gives {Genotype}_{Celltype} combinations that do NOT
    directly match SVG element IDs.  post_process_root() replicates
    the h5ad_to_id expansion and synthetic xylem/phloem averaging from
    mini_app.py so that data_bot_id values match SVG element IDs exactly.
    Final data_bot_id values are things like col0_top_endodermis,
    shr2_xpp_circle, col0_top_xylem (synthetic), etc. — 54 total.

    On top of that, post_process_root() ALSO emits a second family of
    data_bot_id values: the raw obs {Genotype}_{Celltype} combination
    unchanged except the Genotype text is swapped for its short SVG-style
    code (e.g. "Col-0 (shr2)_Root endodermis" -> "shr2_Root endodermis").
    This preserves every raw obs Celltype category, including ones with
    no SVG shape at all (e.g. "col0_G1/G0 phase", "col0_Root hair"), and
    matches the "{genotype_code}_{Celltype}" key generate_umap_dumps.py
    writes for this dataset's UMAP cell_type column — so the UMAP
    cell-type dropdown can look up an average/std for any selected
    category, not just the ones drawn on the eFP diagram.

at_seed:
    Keys are the level_2_annotation_timed obs values directly
    (e.g. "3DAP_Endosperm"), matching SVG element IDs.
    Uses adata.X directly (not logcounts layer).

at_flower / at_silique / at_stem:
    Keys are CellType obs values, matching SVG element IDs directly.
    X matrix may be empty; falls back to logcounts layer automatically.

at_seed_0d_lee / at_rosette_21d_lee / at_shoot_zhang / at_rosette_30d_lee /
at_seedling_3d_lee / at_seedling_6d_lee / at_seedling_12d_lee:
    Newly wired-in datasets, added from dataset_config.yaml with no code
    changes needed here beyond registering their ds_keys below — the
    per-dataset h5ad path / cell_type_col / use_logcounts / db_name all
    come from the shared config. Each is single-column (no condition or
    genotype), so cell_type_col2 is None and there's no combined-key or
    post_process step: unique_cts is just the raw CellType (or
    celltype_after, for Zhang) values, plus Mean_CTRL. No dedicated
    relative-view baseline existed for these before; Mean_CTRL (see
    below) now provides one automatically, same as every other dataset.

The script uses chunked computation (CHUNK=500 genes at a time) to
avoid out-of-memory errors on large datasets.

Mean_CTRL
---------
Every dataset now also gets a synthetic "Mean_CTRL" data_bot_id: the
per-gene average (and std) across ALL cells in the dataset, regardless
of cell type/condition. This is computed the same way as any other
per-cell-type column (it's just a mask that selects every cell), so it
is written into `sample_data` as one more row per gene, alongside the
existing per-tissue rows. It gives the SUPeR Viewer Gene Expression
endpoint a single dataset-wide baseline value per gene (an ~N-gene x 1
"mean expression table") without needing a second SQL table. Pass
--no-mean-ctrl to skip it.

Pass --no-std to skip std computation (sets data_signal_std=0) for
a quick test run.  Pass --gzip to compress output as .sql.gz.
"""

import os
import sys
import gzip
import argparse
import warnings
from datetime import datetime

warnings.filterwarnings("ignore", category=FutureWarning)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dataset_config"))
from dataset_config import (
    pseudobulk_columns, resolve_ds_key, pseudobulk_ready_keys,
    ROOT_GENO_PREFIX, ROOT_H5AD_TO_ID,
)

# ──────────────────────────────────────────────────────────────────────────────
# Dataset config now lives in dataset_config.yaml (see dataset_config.py),
# shared with mini_app.py and generate_umap_dumps.py. --ds still accepts
# the original short keys below (arabidopsis, at_root, ...); they're
# resolved to the canonical ds_key via legacy_key_aliases in the YAML.
# ──────────────────────────────────────────────────────────────────────────────

_DS_CHOICES = ["arabidopsis", "rice", "at_root", "at_seed",
               "at_flower", "at_silique", "at_stem"]

# Newly wired-in datasets (see dataset_config.yaml) — registered under
# their canonical ds_key directly, since none of them has a legacy short
# alias the way the original 7 above do.
_DS_CHOICES += [
    "at_seed_0d_lee", "at_rosette_21d_lee", "at_shoot_zhang",
    "at_rosette_30d_lee", "at_seedling_3d_lee", "at_seedling_6d_lee",
    "at_seedling_12d_lee",
]

# data_bot_id used for the synthetic "average across all cells" row that is
# now added to every dataset (see Mean_CTRL note in the module docstring).
MEAN_CTRL_BOT_ID = "Mean_CTRL"

# For each condition prefix, the two source cell types to average and the
# synthetic bot_id to emit.  Only applied for the arabidopsis leaf dataset.
PHLOEM_AVERAGE_RULES = [
    (f"{cond}_Phloem companion", f"{cond}_Phloem Parenchyma", f"{cond}_Phloem average")
    for cond in ("W0", "D0", "R15", "W15")
]

# ──────────────────────────────────────────────────────────────────────────────
# Root dataset post-processing
# ──────────────────────────────────────────────────────────────────────────────

# _ROOT_H5AD_TO_ID / _ROOT_GENO_PREFIX now come from dataset_config.yaml
# (root_h5ad_to_svg_ids / root_genotype_prefix) via the imports above —
# this used to be a local copy that generate_umap_dumps.py and mini_app.py
# each hardcoded separately, with a comment reminding you to keep them in
# sync by hand.
_ROOT_H5AD_TO_ID = ROOT_H5AD_TO_ID
_ROOT_GENO_PREFIX = ROOT_GENO_PREFIX


def post_process_root(gene_ids, raw_cell_types, avg_matrix, std_matrix):
    """
    Convert raw compute_pseudobulk() output for the Shahan root dataset into
    the final data_bot_id-keyed arrays, made of TWO families:

    Family 1 — SVG-element-ID-keyed columns (unchanged from before):
    exactly replicates the h5ad_to_id expansion and synthetic xylem/phloem
    logic from mini_app.py load_all (Genotype branch). e.g.
    col0_top_endodermis, shr2_xpp_circle, col0_top_xylem (synthetic).

    Family 2 — raw genotype-code + celltype columns (new): every raw obs
    {Genotype}_{Celltype} combination, carried through unchanged except the
    raw Genotype text is swapped for its short SVG-style code via
    _ROOT_GENO_PREFIX (e.g. "Col-0 (shr2)_Root endodermis" ->
    "shr2_Root endodermis"). This preserves every raw obs Celltype
    category — including ones with no SVG shape at all, like
    "col0_G1/G0 phase" or "col0_Root hair" — so the UMAP cell-type
    dropdown (keyed on "{genotype_code}_{Celltype}", see
    generate_umap_dumps.py) can look up an average/std for whatever a
    user selects, not just the categories drawn on the eFP diagram.

    Steps
    -----
    0. If a MEAN_CTRL_BOT_ID ("Mean_CTRL") column is present in
       raw_cell_types, copy it through unchanged — it's an all-cell
       dataset-wide stat, not part of either family's genotype/celltype
       mapping.
    1. For each raw key in _ROOT_H5AD_TO_ID that exists in raw_cell_types,
       copy its column to every target SVG ID column (family 1). Raw keys
       not in _ROOT_H5AD_TO_ID are skipped here (no SVG representation) —
       they still appear in family 2 below.
    2. For each genotype prefix (shr2, scr4, col0), compute synthetic
       family-1 columns:
         {pfx}_top_xylem  = {pfx}_side_xylem  = ({pfx}_proto_circle + {pfx}_meta_circle) / 2
         {pfx}_top_phloem = {pfx}_side_phloem = ({pfx}_xpp_circle   + {pfx}_ppp_circle)  / 2
       These overwrite any values set in step 1 for the phloem IDs, matching
       the mini_app.py behaviour exactly.
    3. For every raw "{Genotype}_{Celltype}" key (excluding Mean_CTRL), add
       a family-2 "{genotype_code}_{Celltype}" column carrying the raw
       average/std through unchanged. A Genotype value not found in
       _ROOT_GENO_PREFIX falls back to its raw text (with a warning)
       rather than being silently dropped.

    For std of a synthetic average: std_avg = sqrt((std_a² + std_b²) / 4),
    which is the std of the mean of two independent groups (equal-weight).

    Returns
    -------
    data_bot_ids : list[str]  — family 1 (54 SVG element IDs) followed by
                   family 2 (one per raw genotype/celltype combination)
    avg_matrix   : np.ndarray, shape (n_genes, len(data_bot_ids))
    std_matrix   : np.ndarray, shape (n_genes, len(data_bot_ids))
    """
    import numpy as np

    raw_idx = {ct: i for i, ct in enumerate(raw_cell_types)}
    n_genes = len(gene_ids)

    # Ordered dict so the final cell_types list is deterministic
    svg_col_avg = {}   # svg_id -> np.ndarray shape (n_genes,)
    svg_col_std = {}

    # ── Step 0: carry Mean_CTRL through untouched ────────────────────────────
    # It's a dataset-wide (all-cell) stat, not a per-genotype/celltype one, so
    # it isn't part of _ROOT_H5AD_TO_ID and needs no SVG-ID expansion.
    if MEAN_CTRL_BOT_ID in raw_idx:
        col_idx = raw_idx[MEAN_CTRL_BOT_ID]
        svg_col_avg[MEAN_CTRL_BOT_ID] = avg_matrix[:, col_idx].copy()
        svg_col_std[MEAN_CTRL_BOT_ID] = std_matrix[:, col_idx].copy()

    # ── Step 1: h5ad_to_id expansion ─────────────────────────────────────────
    for h5ad_key, svg_ids in _ROOT_H5AD_TO_ID.items():
        if h5ad_key not in raw_idx:
            print(
                f"  [root] WARNING: obs key '{h5ad_key}' not found in H5AD — "
                f"skipping {svg_ids}",
                file=sys.stderr,
            )
            continue
        col_idx = raw_idx[h5ad_key]
        avg_col = avg_matrix[:, col_idx].copy()
        std_col = std_matrix[:, col_idx].copy()
        for svg_id in svg_ids:
            svg_col_avg[svg_id] = avg_col
            svg_col_std[svg_id] = std_col

    # ── Step 2: synthetic xylem and phloem (overwrite) ────────────────────────
    for pfx in ('shr2', 'scr4', 'col0'):
        proto_key = f'{pfx}_proto_circle'
        meta_key  = f'{pfx}_meta_circle'
        xpp_key   = f'{pfx}_xpp_circle'
        ppp_key   = f'{pfx}_ppp_circle'

        # Synthetic xylem (new IDs — not overwriting anything meaningful)
        if proto_key in svg_col_avg and meta_key in svg_col_avg:
            xylem_avg = (svg_col_avg[proto_key] + svg_col_avg[meta_key]) / 2.0
            xylem_std = np.sqrt(
                (svg_col_std[proto_key] ** 2 + svg_col_std[meta_key] ** 2) / 4.0
            )
            svg_col_avg[f'{pfx}_top_xylem']  = xylem_avg
            svg_col_std[f'{pfx}_top_xylem']  = xylem_std
            svg_col_avg[f'{pfx}_side_xylem'] = xylem_avg
            svg_col_std[f'{pfx}_side_xylem'] = xylem_std
        else:
            print(
                f"  [root] WARNING: cannot build xylem synthetic for '{pfx}' "
                f"— missing {proto_key!r} or {meta_key!r}",
                file=sys.stderr,
            )

        # Synthetic phloem (overwrites the h5ad_to_id Phloem expansion above)
        if xpp_key in svg_col_avg and ppp_key in svg_col_avg:
            phloem_avg = (svg_col_avg[xpp_key] + svg_col_avg[ppp_key]) / 2.0
            phloem_std = np.sqrt(
                (svg_col_std[xpp_key] ** 2 + svg_col_std[ppp_key] ** 2) / 4.0
            )
            svg_col_avg[f'{pfx}_top_phloem']  = phloem_avg
            svg_col_std[f'{pfx}_top_phloem']  = phloem_std
            svg_col_avg[f'{pfx}_side_phloem'] = phloem_avg
            svg_col_std[f'{pfx}_side_phloem'] = phloem_std
        else:
            print(
                f"  [root] WARNING: cannot build phloem synthetic for '{pfx}' "
                f"— missing {xpp_key!r} or {ppp_key!r}",
                file=sys.stderr,
            )

    print(
        f"  [root] Expanded {len(raw_cell_types)} raw obs groups "
        f"→ {len(svg_col_avg)} SVG element IDs",
        file=sys.stderr,
    )

    # ── Step 3: raw genotype-code + celltype columns (family 2) ──────────────
    # Keep every raw obs Celltype value, per genotype, under its own
    # "{genotype_code}_{Celltype}" id — including categories with no SVG
    # shape at all — so the UMAP cell-type dropdown can look up stats for
    # any of them, matching generate_umap_dumps.py's cell_type key exactly.
    raw_col_avg = {}
    raw_col_std = {}
    _genotypes_by_len = sorted(_ROOT_GENO_PREFIX, key=len, reverse=True)
    unmatched_keys = []

    for raw_key in raw_cell_types:
        if raw_key == MEAN_CTRL_BOT_ID:
            continue
        matched_geno = next(
            (g for g in _genotypes_by_len if raw_key.startswith(g + '_')),
            None,
        )
        if matched_geno is None:
            unmatched_keys.append(raw_key)
            continue
        celltype = raw_key[len(matched_geno) + 1:]
        geno_code = _ROOT_GENO_PREFIX[matched_geno]
        new_key = f"{geno_code}_{celltype}"
        col_idx = raw_idx[raw_key]
        raw_col_avg[new_key] = avg_matrix[:, col_idx].copy()
        raw_col_std[new_key] = std_matrix[:, col_idx].copy()

    if unmatched_keys:
        print(
            f"  [root] WARNING: {len(unmatched_keys)} raw obs key(s) did not "
            f"match a known genotype prefix in _ROOT_GENO_PREFIX, skipped "
            f"from family 2: {unmatched_keys[:5]}"
            f"{'...' if len(unmatched_keys) > 5 else ''}",
            file=sys.stderr,
        )

    print(
        f"  [root] Added {len(raw_col_avg)} raw genotype-code + celltype "
        f"columns (family 2), e.g. 'col0_Root hair', 'shr2_G1/G0 phase'",
        file=sys.stderr,
    )

    # ── Combine both families ─────────────────────────────────────────────
    data_bot_ids = list(svg_col_avg.keys()) + list(raw_col_avg.keys())
    all_avg = {**svg_col_avg, **raw_col_avg}
    all_std = {**svg_col_std, **raw_col_std}

    new_avg = np.column_stack([all_avg[bid] for bid in data_bot_ids])
    new_std = np.column_stack([all_std[bid] for bid in data_bot_ids])

    return data_bot_ids, new_avg, new_std


# ──────────────────────────────────────────────────────────────────────────────
# Dump header / footer  (matches embryo_efp_feb_6_2025_dump.sql style)
# ──────────────────────────────────────────────────────────────────────────────

DUMP_HEADER = """\
-- MySQL dump 10.13  Distrib 9.4.0, for Linux (x86_64)
--
-- Host: localhost    Database: {db_name}
-- ------------------------------------------------------
-- Server version\t9.4.0

/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!50503 SET NAMES utf8mb4 */;
/*!40103 SET @OLD_TIME_ZONE=@@TIME_ZONE */;
/*!40103 SET TIME_ZONE='+00:00' */;
/*!40014 SET @OLD_UNIQUE_CHECKS=@@UNIQUE_CHECKS, UNIQUE_CHECKS=0 */;
/*!40014 SET @OLD_FOREIGN_KEY_CHECKS=@@FOREIGN_KEY_CHECKS, FOREIGN_KEY_CHECKS=0 */;
/*!40101 SET @OLD_SQL_MODE=@@SQL_MODE, SQL_MODE='NO_AUTO_VALUE_ON_ZERO' */;
/*!40111 SET @OLD_SQL_NOTES=@@SQL_NOTES, SQL_NOTES=0 */;

--
-- Current Database: `{db_name}`
--

CREATE DATABASE /*!32312 IF NOT EXISTS*/ `{db_name}` /*!40100 DEFAULT CHARACTER SET latin1 */ /*!80016 DEFAULT ENCRYPTION='N' */;

USE `{db_name}`;

--
-- Table structure for table `sample_data`
--

DROP TABLE IF EXISTS `sample_data`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `sample_data` (
  `data_probeset_id` varchar(32) NOT NULL,
  `data_signal`      float       DEFAULT '0',
  `data_signal_std`  float       DEFAULT '0',
  `data_bot_id`      varchar(64) NOT NULL,
  UNIQUE KEY `uq_probeset_bot` (`data_probeset_id`,`data_bot_id`),
  KEY `data_probeset_id` (`data_probeset_id`,`data_bot_id`,`data_signal`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

--
-- Dumping data for table `sample_data`
--

LOCK TABLES `sample_data` WRITE;
/*!40000 ALTER TABLE `sample_data` DISABLE KEYS */;
"""

DUMP_FOOTER = """\
/*!40000 ALTER TABLE `sample_data` ENABLE KEYS */;
UNLOCK TABLES;
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;
/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
-- Dump completed on {timestamp}
"""


# ──────────────────────────────────────────────────────────────────────────────
# Core computation
# ──────────────────────────────────────────────────────────────────────────────

def compute_pseudobulk(h5ad_path, cell_type_col, cell_type_col2=None,
                       gene_id_col=None, use_logcounts=False, compute_std=True,
                       add_mean_ctrl=True):
    """
    Load H5AD and compute per-(gene, cell_type) average and std.

    For at_root, cell_type_col2='Genotype' produces raw combined keys
    "{Genotype}_{Celltype}".  These are NOT the final data_bot_id values —
    call post_process_root() afterwards (Mean_CTRL, if present, is carried
    through separately — see post_process_root()).

    Parameters
    ----------
    h5ad_path      : path to .h5ad file
    cell_type_col  : primary obs column
    cell_type_col2 : secondary obs column; when set, keys are "{col2}_{col1}"
    gene_id_col    : var column to use for gene IDs; None -> var_names
    use_logcounts  : fall back to logcounts layer if adata.X is empty
    compute_std    : if False, std_matrix is all zeros (fast test mode)
    add_mean_ctrl  : if True (default), append a synthetic MEAN_CTRL_BOT_ID
                     "cell type" whose mask selects every cell in the
                     dataset, so its avg/std columns are the per-gene mean
                     and std across all cells regardless of cell type or
                     condition.

    Returns
    -------
    gene_ids   : list[str]
    cell_types : list[str]  — raw combined labels (may need post-processing);
                 includes MEAN_CTRL_BOT_ID as the last entry when
                 add_mean_ctrl is True
    avg_matrix : np.ndarray, shape (n_genes, n_cell_types)
    std_matrix : np.ndarray, shape (n_genes, n_cell_types)
    """
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp
    import pandas as pd

    print(f"  Loading H5AD: {h5ad_path}", file=sys.stderr)
    adata = ad.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes", file=sys.stderr)

    # ── Optionally fall back to logcounts layer (Lee datasets) ───────────────
    if use_logcounts and 'logcounts' in adata.layers:
        x = adata.X
        x_empty = (
            x is None
            or (sp.issparse(x) and x.nnz == 0)
            or (not sp.issparse(x) and float(x.sum()) == 0.0)
        )
        if x_empty:
            print("  X matrix empty — using logcounts layer", file=sys.stderr)
            adata.X = adata.layers['logcounts']

    # ── Resolve gene IDs ──────────────────────────────────────────────────────
    if gene_id_col and gene_id_col in adata.var.columns:
        gene_ids = [
            str(t).strip() if (isinstance(t, str) and t.strip()) else var
            for var, t in zip(adata.var_names, adata.var[gene_id_col])
        ]
        print(f"  Using '{gene_id_col}' column for gene IDs.", file=sys.stderr)
    else:
        gene_ids = list(adata.var_names)
        print(f"  Using var_names as gene IDs.", file=sys.stderr)

    # ── Resolve cell-type labels ──────────────────────────────────────────────
    if cell_type_col not in adata.obs.columns:
        raise ValueError(
            f"Column '{cell_type_col}' not found in obs. "
            f"Available: {list(adata.obs.columns)}"
        )

    col1 = adata.obs[cell_type_col].astype(str).to_numpy()

    if cell_type_col2 and cell_type_col2 in adata.obs.columns:
        col2 = adata.obs[cell_type_col2].astype(str).to_numpy()
        # "{col2}_{col1}" ordering:
        #   rice  -> "{Condition}_{CellAnnotation}" (split on first _ gives cond, ct)
        #   root  -> "{Genotype}_{Celltype}"        (matches _ROOT_H5AD_TO_ID keys)
        combined = np.array([f"{c}_{t}" for c, t in zip(col2, col1)])
        cell_type_labels = combined
        print(
            f"  Combined cell-type key: '{cell_type_col2}' + '{cell_type_col}'",
            file=sys.stderr,
        )
    else:
        cell_type_labels = col1

    unique_cts = list(pd.unique(cell_type_labels))

    # ── Mean_CTRL: one more "cell type" whose mask covers every cell ─────────
    # Computed by the exact same masked-mean/std loop below, so it needs no
    # special-cased math — the mask just happens to select everything.
    if add_mean_ctrl:
        unique_cts.append(MEAN_CTRL_BOT_ID)

    n_genes = adata.n_vars
    n_cts = len(unique_cts)
    print(f"  {n_cts} unique cell types / conditions found.", file=sys.stderr)

    # ── Densify X in gene chunks to avoid OOM ────────────────────────────────
    avg_matrix = np.zeros((n_genes, n_cts), dtype=np.float32)
    std_matrix = np.zeros((n_genes, n_cts), dtype=np.float32)
    masks = [
        (np.ones(len(cell_type_labels), dtype=bool) if ct == MEAN_CTRL_BOT_ID
         and add_mean_ctrl else cell_type_labels == ct)
        for ct in unique_cts
    ]

    CHUNK = 500
    for g_start in range(0, n_genes, CHUNK):
        g_end = min(g_start + CHUNK, n_genes)
        chunk = adata.X[:, g_start:g_end]
        if sp.issparse(chunk):
            chunk = chunk.toarray()
        chunk = chunk.astype(np.float32)

        for j, mask in enumerate(masks):
            sub = chunk[mask]
            avg_matrix[g_start:g_end, j] = sub.mean(axis=0)
            if compute_std:
                std_matrix[g_start:g_end, j] = sub.std(axis=0)

        if (g_start // CHUNK) % 10 == 0:
            print(
                f"  Processed genes {g_start}–{g_end} / {n_genes}",
                file=sys.stderr,
            )

    return gene_ids, unique_cts, avg_matrix, std_matrix


# ──────────────────────────────────────────────────────────────────────────────
# Phloem average synthetic rows  (arabidopsis leaf only)
# ──────────────────────────────────────────────────────────────────────────────

def build_phloem_averages(cell_types, avg_matrix, std_matrix):
    """
    For each rule in PHLOEM_AVERAGE_RULES, if both source cell types exist,
    compute element-wise mean of avg and propagated std, and return a list of
    (bot_id, avg_col, std_col) tuples to append to the dump.

    std of the average of two groups: sqrt((std_a² + std_b²) / 4)
    """
    import numpy as np

    ct_index = {ct: i for i, ct in enumerate(cell_types)}
    extras = []

    for ct_a, ct_b, out_id in PHLOEM_AVERAGE_RULES:
        if ct_a not in ct_index or ct_b not in ct_index:
            print(
                f"  Skipping phloem average for '{out_id}': "
                f"missing '{ct_a}' or '{ct_b}'",
                file=sys.stderr,
            )
            continue

        i_a = ct_index[ct_a]
        i_b = ct_index[ct_b]
        avg_col = (avg_matrix[:, i_a] + avg_matrix[:, i_b]) / 2.0
        std_col = np.sqrt(
            (std_matrix[:, i_a] ** 2 + std_matrix[:, i_b] ** 2) / 4.0
        )
        extras.append((out_id, avg_col, std_col))
        print(f"  Added synthetic cell type: '{out_id}'", file=sys.stderr)

    return extras


# ──────────────────────────────────────────────────────────────────────────────
# SQL dump writer
# ──────────────────────────────────────────────────────────────────────────────

ROWS_PER_INSERT = 500


def escape_sql(s):
    """Minimal SQL string escaping for latin1 context."""
    return str(s).replace("\\", "\\\\").replace("'", "\\'")


def write_dump(out_path, db_name, gene_ids, cell_types, avg_matrix, std_matrix,
               phloem_extras=None):
    """
    Write a mysqldump-compatible .sql (or .sql.gz) file.

    phloem_extras : list of (bot_id, avg_col, std_col) synthetic rows;
                    only used for the arabidopsis leaf dataset.
    """
    n_genes = len(gene_ids)
    n_cts = len(cell_types)
    n_extras = len(phloem_extras) if phloem_extras else 0
    total = n_genes * (n_cts + n_extras)

    print(f"  Writing {total:,} rows to {out_path} ...", file=sys.stderr)

    open_fn = gzip.open if out_path.endswith('.gz') else open
    mode = 'wt' if out_path.endswith('.gz') else 'w'

    with open_fn(out_path, mode, encoding='latin-1', errors='replace') as fh:
        fh.write(DUMP_HEADER.format(db_name=db_name))

        row_buf = []

        def flush_buf():
            if not row_buf:
                return
            fh.write("INSERT INTO `sample_data` VALUES ")
            fh.write(",".join(row_buf))
            fh.write(";\n")
            row_buf.clear()

        count = 0

        for g_idx, gene_id in enumerate(gene_ids):
            safe_gene = escape_sql(gene_id)
            for ct_idx, ct in enumerate(cell_types):
                avg = float(avg_matrix[g_idx, ct_idx])
                std = float(std_matrix[g_idx, ct_idx])
                safe_ct = escape_sql(ct)
                row_buf.append(
                    f"('{safe_gene}',{avg:.6g},{std:.6g},'{safe_ct}')"
                )
                count += 1
                if len(row_buf) >= ROWS_PER_INSERT:
                    flush_buf()

            if g_idx % 1000 == 0 and g_idx > 0:
                print(
                    f"    {g_idx:,} / {n_genes:,} genes written",
                    file=sys.stderr,
                )

        flush_buf()

        # ── Synthetic phloem average rows (arabidopsis leaf only) ─────────────
        if phloem_extras:
            print(
                f"  Writing {n_extras} synthetic phloem average cell type(s)...",
                file=sys.stderr,
            )
            for bot_id, avg_col, std_col in phloem_extras:
                safe_ct = escape_sql(bot_id)
                for g_idx, gene_id in enumerate(gene_ids):
                    safe_gene = escape_sql(gene_id)
                    avg = float(avg_col[g_idx])
                    std = float(std_col[g_idx])
                    row_buf.append(
                        f"('{safe_gene}',{avg:.6g},{std:.6g},'{safe_ct}')"
                    )
                    count += 1
                    if len(row_buf) >= ROWS_PER_INSERT:
                        flush_buf()
                flush_buf()

        ts = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
        fh.write(DUMP_FOOTER.format(timestamp=ts))

    print(f"  Done — {count:,} rows written.", file=sys.stderr)


# ──────────────────────────────────────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────────────────────────────────────

def build_custom_config(args):
    """
    Build a one-off dataset config from command-line arguments, so a new
    H5AD can be dumped without adding it to dataset_config.yaml first.
    Returns a dict with the same keys as pseudobulk_columns() returns.
    """
    h5ad_path = os.path.expanduser(args.h5ad)

    # Default database name: the H5AD file name + "_pseudobulk"
    # e.g. lee_2025_at_root.h5ad -> lee_2025_at_root_pseudobulk
    if args.db_name:
        db_name = args.db_name
    else:
        file_name = os.path.basename(h5ad_path)
        file_stem = os.path.splitext(file_name)[0]
        db_name = file_stem + "_pseudobulk"

    config = {
        "ds_key": "custom",
        "h5ad": h5ad_path,
        "cell_type_col": args.cell_type_col,
        "cell_type_col2": args.cell_type_col2,
        "gene_id_col": args.gene_id_col,
        "use_logcounts": args.use_logcounts,
        "post_process": None,
        "db_name": db_name,
    }
    return config


def main():
    parser = argparse.ArgumentParser(
        description="Generate MySQL pseudobulk dump files from H5AD datasets."
    )

    # ── Option 1: registered datasets from dataset_config.yaml ───────────────
    parser.add_argument(
        "--ds",
        nargs="+",
        choices=_DS_CHOICES + ["all"],
        default=["all"],
        help="One or more registered datasets to process (default: all)",
    )

    # ── Option 2: any H5AD file, no YAML editing needed ──────────────────────
    parser.add_argument(
        "--h5ad",
        default=None,
        help="Path to an H5AD file to process directly (ignores --ds)",
    )
    parser.add_argument(
        "--cell-type-col",
        default="CellType",
        help="With --h5ad: obs column holding the cell types (default: CellType)",
    )
    parser.add_argument(
        "--cell-type-col2",
        default=None,
        help="With --h5ad: optional second obs column (condition/genotype); "
             "data_bot_id becomes '{col2}_{col1}'",
    )
    parser.add_argument(
        "--gene-id-col",
        default=None,
        help="With --h5ad: optional var column holding gene IDs (e.g. TAIR_ID); "
             "default uses var_names",
    )
    parser.add_argument(
        "--use-logcounts",
        action="store_true",
        help="With --h5ad: use the 'logcounts' layer if X is empty",
    )
    parser.add_argument(
        "--db-name",
        default=None,
        help="With --h5ad: MySQL database name (default: <file name>_pseudobulk)",
    )
    parser.add_argument(
        "--outdir",
        default=".",
        help="Output directory for .sql dump files (default: current dir)",
    )
    parser.add_argument(
        "--no-std",
        action="store_true",
        help="Skip std computation (sets data_signal_std=0); much faster for testing",
    )
    parser.add_argument(
        "--gzip",
        action="store_true",
        help="Compress output files with gzip (.sql.gz)",
    )
    parser.add_argument(
        "--no-mean-ctrl",
        action="store_true",
        help="Skip the synthetic Mean_CTRL (all-cell average) row per gene",
    )
    parser.add_argument(
        "-n", "--ngenes",
        type=int,
        default=None,
        help="Limit to the first N genes (e.g. -n 1 for a quick test)",
    )
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    compute_std = not args.no_std
    add_mean_ctrl = not args.no_mean_ctrl

    # ── Collect the config for every dataset we need to process ──────────────
    configs = []
    if args.h5ad:
        configs.append(build_custom_config(args))
    elif "all" in args.ds:
        for ds_key in _DS_CHOICES:
            configs.append(pseudobulk_columns(ds_key))
    else:
        for ds_key in args.ds:
            configs.append(pseudobulk_columns(ds_key))

    for cfg in configs:
        ds_key         = cfg["ds_key"]
        db_name        = cfg["db_name"]
        h5ad_path      = cfg["h5ad"]
        cell_type_col  = cfg["cell_type_col"]
        cell_type_col2 = cfg["cell_type_col2"]
        gene_id_col    = cfg["gene_id_col"]
        use_logcounts  = cfg["use_logcounts"]
        post_process   = cfg["post_process"]

        print(f"\n{'=' * 60}", file=sys.stderr)
        print(f"Processing dataset: {ds_key} -> {db_name}", file=sys.stderr)
        print(f"{'=' * 60}", file=sys.stderr)

        if not os.path.exists(h5ad_path):
            print(f"  ERROR: H5AD not found at {h5ad_path}", file=sys.stderr)
            continue

        gene_ids, cell_types, avg_matrix, std_matrix = compute_pseudobulk(
            h5ad_path=h5ad_path,
            cell_type_col=cell_type_col,
            cell_type_col2=cell_type_col2,
            gene_id_col=gene_id_col,
            use_logcounts=use_logcounts,
            compute_std=compute_std,
            add_mean_ctrl=add_mean_ctrl,
        )

        # ── Optional gene-count limit (for quick testing) ───────────────────
        if args.ngenes is not None:
            n = min(args.ngenes, len(gene_ids))
            gene_ids   = gene_ids[:n]
            avg_matrix = avg_matrix[:n]
            std_matrix = std_matrix[:n]
            print(f"  --ngenes {args.ngenes}: truncated to {n} gene(s)", file=sys.stderr)

        # ── Root: expand raw obs keys to SVG element IDs ──────────────────────
        if post_process == "root":
            cell_types, avg_matrix, std_matrix = post_process_root(
                gene_ids, cell_types, avg_matrix, std_matrix
            )

        # ── Arabidopsis leaf: append synthetic phloem average rows ────────────
        phloem_extras = None
        if resolve_ds_key(ds_key) == "arabidopsis_nat":
            phloem_extras = build_phloem_averages(cell_types, avg_matrix, std_matrix)

        ext = ".sql.gz" if args.gzip else ".sql"
        out_file = os.path.join(args.outdir, f"{db_name}_dump{ext}")

        write_dump(
            out_file, db_name, gene_ids, cell_types, avg_matrix, std_matrix,
            phloem_extras=phloem_extras,
        )
        print(f"  Saved: {out_file}", file=sys.stderr)

    print("\nAll done.", file=sys.stderr)


if __name__ == "__main__":
    main()