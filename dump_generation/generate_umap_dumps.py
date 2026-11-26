#!/usr/bin/env python3
"""
Steven Qiao | BAR Lab | University of Toronto

Generate MySQL SQL dumps for UMAP coordinate data from H5AD files.

Two-table design:
    umap_coords     — sampled cell coordinates, stored ONCE per dataset,
                      shared across all genes (no redundancy)
    umap_expression — one row per gene with a JSON object of expression
                      values mapping cell_id -> non-zero expression.

Positional matching guarantee:
    umap_coords row with cell_id=X  <->  expression["X"] (if present)
    Zero expression values are entirely dropped to save space.

Key design decisions:
    - X, Y, and cell_type do NOT change per gene — stored once in umap_coords
    - Expression is a sparse JSON object per gene: {"cell_id_1": 12.5, "cell_id_2": 5.0, ...}
    - Sampling is done ONCE per dataset — same cells reused for every gene
    - Expression JSON is built in-memory per gene and flushed immediately
    - Inner loop vectorized with numpy fancy indexing (no Python for loops)

Cell-type label convention
---------------------------
For datasets where the pseudobulk SVG viewer keys expression on a
condition/genotype PLUS a cell type (rather than cell type alone), the
`cell_type` column written here is the same combined key, so a hover on
the SVG can find matching points in the UMAP via a plain string match.

    rice ("CellAnnotation" + "Condition"):
        combined as "{Condition}_{CellAnnotation}", using the raw obs
        values directly — these already match generate_pseudobulk_dumps.py's
        data_bot_id cell-type portion (both scripts read the same obs
        columns with no relabeling), including the '.'-separated spelling
        already present in the source h5ad (e.g. "Mild.Drought", "Phloem.SE").

    at_root_rs ("Celltype" + "Genotype"):
        combined as "{genotype_code}_{Celltype}", where genotype_code is
        the SHORT SVG-style code (shr2 / scr4 / col0) — NOT the raw
        Genotype obs value ("Col-0 (shr2)", "Ler (scr4)", "Col-0"). The
        SVG element ids and generate_pseudobulk_dumps.py's post-processed
        data_bot_id values (e.g. "shr2_top_endodermis") use these short
        codes, so raw-value concatenation would silently never match
        anything. _ROOT_GENO_PREFIX below mirrors the identical table in
        generate_pseudobulk_dumps.py.

Usage:
    Option 1 - datasets already registered in dataset_config.yaml:
        python3 generate_umap_dumps.py                                  # all datasets
        python3 generate_umap_dumps.py --datasets rice arabidopsis_nat
        python3 generate_umap_dumps.py --max-cells 5000 --out-dir /tmp/umap_dumps

    Option 2 - any H5AD file, straight from the command line (no YAML editing):
        python3 generate_umap_dumps.py --h5ad lee_2025_at_root.h5ad
        python3 generate_umap_dumps.py --h5ad robertson_2026_rice.h5ad \\
            --umap-col CellAnnotation --umap-col2 Condition \\
            --db-name rice_OW_umap --max-cells 15000

        Defaults for Option 2:
            --umap-col  CellType  (the column rds_to_h5ad.r always writes)
            --db-name   <h5ad file name>_umap
            --umap-col2 / --gene-id-col: not used
        The H5AD must contain obsm['X_umap'] (rds_to_h5ad.r writes this).

    Quick test on a couple of genes (works with both options):
        python3 generate_umap_dumps.py --datasets rice --genes LOC_Os01g01010 LOC_Os01g01019

    --outdir is accepted as another spelling of --out-dir.

Output:
    One .sql file per dataset in --out-dir (default: ./umap_dumps/)

DO NOT commit the output .sql files to git.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import warnings
import time
from datetime import datetime
from pathlib import Path

warnings.filterwarnings("ignore", category=FutureWarning)
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dataset_config"))
from dataset_config import umap_columns, umap_ready_keys, ROOT_GENO_PREFIX

# ---------------------------------------------------------------------------
# Dataset registry — now sourced from dataset_config.yaml (see
# dataset_config.py), shared with mini_app.py and
# generate_pseudobulk_dumps.py, instead of a local copy here.
# ---------------------------------------------------------------------------

DATASETS = {ds_key: umap_columns(ds_key) for ds_key in umap_ready_keys()}

# Genotype prefix in SVG element ids for each raw obs Genotype value.
# Now sourced from dataset_config.yaml's root_genotype_prefix — previously
# a local copy of the identical table in generate_pseudobulk_dumps.py, kept
# in sync by hand via a comment in both files.
_ROOT_GENO_PREFIX = ROOT_GENO_PREFIX

# ---------------------------------------------------------------------------
# SQL templates — mysqldump-compatible format
# ---------------------------------------------------------------------------

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

"""

DUMP_FOOTER = """\
/*!40103 SET TIME_ZONE=@OLD_TIME_ZONE */;
/*!40101 SET SQL_MODE=@OLD_SQL_MODE */;
/*!40014 SET FOREIGN_KEY_CHECKS=@OLD_FOREIGN_KEY_CHECKS */;
/*!40014 SET UNIQUE_CHECKS=@OLD_UNIQUE_CHECKS */;
/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
/*!40111 SET SQL_NOTES=@OLD_SQL_NOTES */;
"""

TABLE_DDL_COORDS = """\
--
-- Table structure for table `umap_coords`
--

DROP TABLE IF EXISTS `umap_coords`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `umap_coords` (
  `cell_id`   INT          NOT NULL,
  `umap_1`    FLOAT        NOT NULL,
  `umap_2`    FLOAT        NOT NULL,
  `cell_type` VARCHAR(128) NOT NULL,
  PRIMARY KEY (`cell_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

"""

TABLE_DDL_EXPRESSION = """\
--
-- Table structure for table `umap_expression`
--

DROP TABLE IF EXISTS `umap_expression`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `umap_expression` (
  `gene_id`    VARCHAR(32) NOT NULL,
  `expression` JSON        NOT NULL,
  PRIMARY KEY (`gene_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

"""


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _escape(s: str) -> str:
    return s.replace("\\", "\\\\").replace("'", "\\'")


def stratified_sample(labels, max_cells: int, rng):
    """Sample cell indices once, stratified by cell type."""
    import numpy as np
    n_cells = len(labels)
    if n_cells <= max_cells:
        return np.arange(n_cells)

    unique_cts, counts = np.unique(labels, return_counts=True)
    n_types = len(unique_cts)
    sampled = []
    remaining = max_cells

    for i, (ct, count) in enumerate(zip(unique_cts, counts)):
        quota = min(remaining // (n_types - i), count)
        ct_indices = np.where(labels == ct)[0]
        sampled.append(rng.choice(ct_indices, size=quota, replace=False))
        remaining -= quota

    return np.sort(np.concatenate(sampled))


def write_coords(fh, sampled_idx, umap, labels, chunk_size: int) -> int:
    """Write umap_coords INSERT statements. Returns number of rows written."""
    fh.write("--\n-- Dumping data for table `umap_coords`\n--\n\n")
    fh.write("LOCK TABLES `umap_coords` WRITE;\n")
    fh.write("/*!40000 ALTER TABLE `umap_coords` DISABLE KEYS */;\n")

    rows = []
    for cell_id, idx in enumerate(sampled_idx):
        u1, u2 = float(umap[idx, 0]), float(umap[idx, 1])
        if not (u1 == u1 and u2 == u2 and abs(u1) < 1e10 and abs(u2) < 1e10):
            continue
        rows.append(f"({cell_id},{u1:.6f},{u2:.6f},'{_escape(labels[idx])}')")

    header = "INSERT INTO `umap_coords` VALUES "
    for start in range(0, len(rows), chunk_size):
        chunk = rows[start: start + chunk_size]
        fh.write(header)
        fh.write(",".join(chunk))
        fh.write(";\n")

    fh.write("/*!40000 ALTER TABLE `umap_coords` ENABLE KEYS */;\n")
    fh.write("UNLOCK TABLES;\n\n")

    return len(rows)


def build_expr_tuple(gene_id: str, sampled_idx, X_col) -> str:
    """
    Build a single INSERT value tuple for one gene.

    Filters out zeros and stores only non-zero expression values in a
    compact JSON object: {"cell_id": value, "cell_id": value, ...}

    Position i corresponds to cell_id in `umap_coords`.
    """
    import numpy as np

    escaped_gene = _escape(gene_id)

    # 1. Slice expression values for sampled cells
    expr_values = X_col[sampled_idx].astype(float)
    expr_values = np.where(np.isfinite(expr_values), expr_values, 0.0)

    # 2. Vectorized mask to select only non-zero entries
    nonzero_mask = expr_values != 0.0

    # 3. Get the corresponding virtual cell_ids (0 to len(sampled_idx)-1)
    cell_ids = np.where(nonzero_mask)[0]
    filtered_values = expr_values[nonzero_mask]

    # 4. Map them to a dictionary { "cell_id": rounded_val }
    # (JSON keys must be strings, so we convert cell_id to str)
    expr_dict = {
        str(cell_id): round(float(val), 6)
        for cell_id, val in zip(cell_ids, filtered_values)
    }

    # 5. Build compact JSON string — no whitespaces
    json_str = json.dumps(expr_dict, separators=(",", ":"))
    json_str_escaped = json_str.replace("\\", "\\\\").replace("'", "\\'")

    return f"('{escaped_gene}','{json_str_escaped}')"


# ---------------------------------------------------------------------------
# Main dump writer — streaming, one gene at a time
# ---------------------------------------------------------------------------

def write_dump(
        ds_key: str,
        cfg: dict,
        out_dir: Path,
        chunk_size: int,
        gene_ids: list[str] | None,
        max_cells: int,
        seed: int,
) -> None:
    import anndata as ad
    import numpy as np
    import scipy.sparse as sp

    db_name = cfg["db_name"]
    out_path = out_dir / f"{db_name}.sql"

    print(f"\n=== {ds_key} -> {db_name} ===", file=sys.stderr)

    print(f"  Loading {cfg['h5ad']} ...", file=sys.stderr)
    t0 = time.perf_counter()
    adata = ad.read_h5ad(cfg["h5ad"])
    print(
        f"  Loaded: {adata.n_obs:,} cells x {adata.n_vars:,} genes "
        f"({time.perf_counter() - t0:.1f}s)",
        file=sys.stderr,
    )

    if sp.issparse(adata.X):
        if not isinstance(adata.X, sp.csc_matrix):
            print("  [performance] Converting expression matrix to CSC format for rapid column slicing...",
                  file=sys.stderr)
            t_conv = time.perf_counter()
            adata.X = adata.X.tocsc()
            print(f"  Converted in {time.perf_counter() - t_conv:.1f}s", file=sys.stderr)

    if "X_umap" not in adata.obsm:
        print("  [error] No X_umap -- skipping", file=sys.stderr)
        return

    umap = np.asarray(adata.obsm["X_umap"])[:, :2].astype(float)

    # ── Cell-type label(s) ────────────────────────────────────────────────
    umap_col = cfg["umap_col"]
    umap_col2 = cfg.get("umap_col2")

    if umap_col not in adata.obs.columns:
        fallback = adata.obs.columns[0]
        print(f"  [warn] umap_col '{umap_col}' not found, using '{fallback}'", file=sys.stderr)
        umap_col = fallback

    col1 = adata.obs[umap_col].astype(str).str.strip().to_numpy()

    if umap_col2 and umap_col2 in adata.obs.columns:
        col2_raw = adata.obs[umap_col2].astype(str).str.strip().to_numpy()

        if ds_key == "at_root_rs":
            # Map the raw Genotype obs value to its short SVG-style code
            # (shr2 / scr4 / col0) before combining — the SVG element ids
            # and the pseudobulk data_bot_id values use these short codes,
            # not the raw obs text. Unrecognized genotype values fall back
            # to the raw text so nothing silently disappears; a warning is
            # printed so the mapping can be extended if new genotypes show up.
            unmapped = sorted(set(col2_raw) - set(_ROOT_GENO_PREFIX))
            if unmapped:
                print(
                    f"  [warn] Genotype value(s) not in _ROOT_GENO_PREFIX, "
                    f"using raw text: {unmapped}",
                    file=sys.stderr,
                )
            col2 = np.array([_ROOT_GENO_PREFIX.get(g, g) for g in col2_raw])
        else:
            col2 = col2_raw

        # "{col2}_{col1}" ordering — mirrors generate_pseudobulk_dumps.py's
        # combined data_bot_id key exactly (condition/genotype-code first,
        # then cell type), so this UMAP dump's `cell_type` column agrees
        # with the pseudobulk viewer's cell-type keys.
        labels = np.array([f"{c}_{t}" for c, t in zip(col2, col1)])
        print(f"  Combined cell-type key: '{umap_col2}' + '{umap_col}'", file=sys.stderr)
    elif umap_col2:
        print(f"  [warn] umap_col2 '{umap_col2}' not found — using '{umap_col}' alone", file=sys.stderr)
        labels = col1
    else:
        labels = col1

    # Sample cell indices ONCE for the whole dataset
    rng = np.random.default_rng(seed)
    sampled_idx = stratified_sample(labels, max_cells, rng)
    print(f"  Sampled {len(sampled_idx):,} / {adata.n_obs:,} cells", file=sys.stderr)

    # Resolve gene list — use gene_id_col if configured, otherwise var_names
    gene_id_col = cfg.get("gene_id_col")
    var_index = adata.var_names

    if gene_id_col and gene_id_col in adata.var.columns:
        # Build lookup: TAIR_ID (or equivalent) -> positional index in var
        id_series = adata.var[gene_id_col].astype(str).str.strip()
        # gene_id_to_varidx: maps external gene ID -> integer position
        gene_id_to_varidx = {}
        for pos, ext_id in enumerate(id_series):
            if ext_id and ext_id.upper() not in gene_id_to_varidx:
                gene_id_to_varidx[ext_id.upper()] = pos
        all_gene_ids = list(dict.fromkeys(id_series))  # deduplicated, ordered
        print(f"  Using '{gene_id_col}' column for gene IDs.", file=sys.stderr)
    else:
        gene_id_to_varidx = {v.upper(): i for i, v in enumerate(var_index)}
        all_gene_ids = list(var_index)

    target_genes = all_gene_ids if gene_ids is None else gene_ids
    n_genes = len(target_genes)

    print(f"  Writing to {out_path} ...", file=sys.stderr)

    with open(out_path, "w", encoding="utf-8", buffering=8 * 1024 * 1024) as fh:
        # --- Header ---
        fh.write(DUMP_HEADER.format(db_name=db_name))

        # --- umap_coords DDL + data ---
        fh.write(TABLE_DDL_COORDS)
        n_coords = write_coords(fh, sampled_idx, umap, labels, chunk_size)
        print(f"  Coords written: {n_coords:,} rows", file=sys.stderr)

        # --- umap_expression DDL ---
        fh.write(TABLE_DDL_EXPRESSION)
        fh.write("--\n-- Dumping data for table `umap_expression`\n--\n\n")
        fh.write("LOCK TABLES `umap_expression` WRITE;\n")
        fh.write("/*!40000 ALTER TABLE `umap_expression` DISABLE KEYS */;\n")

        # --- Stream expression: one row per gene, chunked INSERTs ---
        t_gene_start = time.perf_counter()
        genes_written = 0
        insert_header = "INSERT INTO `umap_expression` VALUES "
        value_buf: list[str] = []

        for g_num, gene_id in enumerate(target_genes, 1):
            varidx = gene_id_to_varidx.get(gene_id.upper())
            if varidx is None:
                print(f"  [warn] Gene '{gene_id}' not found -- skipping", file=sys.stderr)
                continue

            X_col = adata.X[:, varidx]
            if sp.issparse(X_col):
                X_col = X_col.toarray().ravel()
            else:
                X_col = np.asarray(X_col).ravel()

            value_buf.append(build_expr_tuple(gene_id, sampled_idx, X_col))
            genes_written += 1

            # Flush chunk to disk
            if len(value_buf) >= chunk_size:
                fh.write(insert_header)
                fh.write(",".join(value_buf))
                fh.write(";\n")
                value_buf.clear()

            if g_num % 500 == 0 or g_num == n_genes:
                elapsed = time.perf_counter() - t_gene_start
                rate = g_num / elapsed
                eta = (n_genes - g_num) / rate if rate > 0 else 0
                print(
                    f"  [{g_num:>6}/{n_genes}] {gene_id} | "
                    f"{rate:.0f} genes/s | ETA {eta / 60:.1f} min",
                    file=sys.stderr,
                )

        # Flush remaining
        if value_buf:
            fh.write(insert_header)
            fh.write(",".join(value_buf))
            fh.write(";\n")
            value_buf.clear()

        fh.write("/*!40000 ALTER TABLE `umap_expression` ENABLE KEYS */;\n")
        fh.write("UNLOCK TABLES;\n")

        # --- Footer ---
        fh.write(DUMP_FOOTER)
        fh.write(f"-- Dump completed on {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")

    size_mb = out_path.stat().st_size / (1024 * 1024)
    print(
        f"  Done: {n_coords:,} coord rows, {genes_written} genes, {size_mb:.1f} MB",
        file=sys.stderr,
    )


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_custom_config(args) -> dict:
    """
    Build a one-off dataset config from command-line arguments, so a new
    H5AD can be dumped without adding it to dataset_config.yaml first.
    Returns a dict with the same keys as umap_columns() returns.
    """
    h5ad_path = os.path.expanduser(args.h5ad)

    # Default database name: the H5AD file name + "_umap"
    # e.g. lee_2025_at_root.h5ad -> lee_2025_at_root_umap
    if args.db_name:
        db_name = args.db_name
    else:
        file_name = os.path.basename(h5ad_path)
        file_stem = os.path.splitext(file_name)[0]
        db_name = file_stem + "_umap"

    config = {
        "ds_key": "custom",
        "h5ad": h5ad_path,
        "umap_col": args.umap_col,
        "umap_col2": args.umap_col2,
        "gene_id_col": args.gene_id_col,
        "db_name": db_name,
    }
    return config


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Generate two-table UMAP SQL dumps from H5AD files for BAR API",
    )

    # ── Option 1: registered datasets from dataset_config.yaml ───────────────
    parser.add_argument(
        "--datasets", nargs="+", choices=list(DATASETS.keys()),
        default=list(DATASETS.keys()),
        help="Which registered datasets to process (default: all)",
    )

    # ── Option 2: any H5AD file, no YAML editing needed ──────────────────────
    parser.add_argument(
        "--h5ad", default=None,
        help="Path to an H5AD file to process directly (ignores --datasets)",
    )
    parser.add_argument(
        "--umap-col", default="CellType",
        help="With --h5ad: obs column used to label UMAP points (default: CellType)",
    )
    parser.add_argument(
        "--umap-col2", default=None,
        help="With --h5ad: optional second obs column (condition/genotype); "
             "cell_type becomes '{col2}_{col1}'",
    )
    parser.add_argument(
        "--gene-id-col", default=None,
        help="With --h5ad: optional var column holding gene IDs (e.g. TAIR_ID); "
             "default uses var_names",
    )
    parser.add_argument(
        "--db-name", default=None,
        help="With --h5ad: MySQL database name (default: <file name>_umap)",
    )
    parser.add_argument(
        "--genes", nargs="+", default=None, metavar="GENE_ID",
        help="Specific gene IDs to dump (default: ALL genes)",
    )
    parser.add_argument(
        "--max-cells", type=int, default=10000, metavar="N",
        help="Max sampled cells per dataset (default: 10000)",
    )
    parser.add_argument(
        "--chunk-size", type=int, default=1000, metavar="N",
        help="Rows per INSERT statement (default: 1000)",
    )
    parser.add_argument(
        "--out-dir", "--outdir", dest="out_dir", default="umap_dumps", metavar="DIR",
        help="Output directory (default: ./umap_dumps/)",
    )
    parser.add_argument(
        "--seed", type=int, default=42,
        help="Random seed for reproducible sampling (default: 42)",
    )
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.genes is None:
        print(
            "WARNING: No --genes specified -- dumping ALL genes.\n"
            "Use --genes GENE_ID ... to test with a subset first.\n",
            file=sys.stderr,
        )

    # ── Collect the config for every dataset we need to process ──────────────
    configs = []
    if args.h5ad:
        configs.append(build_custom_config(args))
    else:
        for ds_key in args.datasets:
            configs.append(DATASETS[ds_key])

    for cfg in configs:
        ds_key = cfg["ds_key"]
        if not os.path.exists(cfg["h5ad"]):
            print(f"[skip] {ds_key}: H5AD not found at {cfg['h5ad']}", file=sys.stderr)
            continue
        write_dump(ds_key, cfg, out_dir, args.chunk_size, args.genes, args.max_cells, args.seed)

    print("\nAll done.", file=sys.stderr)
    print(
        "\nNext steps:\n"
        "  1. Copy .sql files to config/databases/\n"
        "  2. Add load lines to config/init.sh\n"
        "  3. Add binds to config/BAR_API.cfg (host: BAR_mysqldb)\n"
        "  4. Add db names to DATABASE_SPECIES in gene_id_utils.py\n"
        "  5. Update umap_expression.py to query JSON expression object keys\n"
        "  6. DO NOT git add the .sql files",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()

