"""
dataset_config.py
==================

Loads dataset_config.yaml and exposes everything mini_app.py,
generate_pseudobulk_dumps.py, and generate_umap_dumps.py need from it, so
a new dataset (or a column-name fix) only ever gets added in one place.

Usage
-----
    from dataset_config import DATASETS, resolve_columns

    # mini_app.py: replaces its old three-times-duplicated
    # ds_key/col_key if/elif dispatch chain
    col1, col2, svg_path = resolve_columns("at_seed_martin",
                                            col_key="level_3_annotation_full_timed")

    from dataset_config import resolve_ds_key, pseudobulk_columns

    # generate_pseudobulk_dumps.py: resolves legacy short --ds keys
    # (arabidopsis, at_root, ...) to the canonical ds_key, then reads
    # the same columns mini_app.py and generate_umap_dumps.py use
    ds_key = resolve_ds_key("at_root")   # -> "at_root_rs"
    cols = pseudobulk_columns(ds_key)
    # cols == {"h5ad": ..., "cell_type_col": "Celltype",
    #          "cell_type_col2": "Genotype", "gene_id_col": None,
    #          "use_logcounts": False, "post_process": "root",
    #          "db_name": "arabidopsis_root_rs_pseudobulk"}

    from dataset_config import umap_columns, ROOT_GENO_PREFIX, ROOT_H5AD_TO_ID

    # generate_umap_dumps.py
    cols = umap_columns("at_root_rs")
    # cols == {"h5ad": ..., "umap_col": "Celltype", "umap_col2": "Genotype",
    #          "gene_id_col": None, "db_name": "arabidopsis_root_shahan_umap"}
"""
import os

import yaml

_CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "dataset_config.yaml")


def _load(path=_CONFIG_PATH):
    with open(path, "r") as f:
        raw = yaml.safe_load(f)
    return raw


_raw = _load()
DATASETS = _raw.get("datasets", {})
LEGACY_KEY_ALIASES = _raw.get("legacy_key_aliases", {})
ROOT_GENO_PREFIX = _raw.get("root_genotype_prefix", {})
ROOT_H5AD_TO_ID = _raw.get("root_h5ad_to_svg_ids", {})


def resolve_ds_key(ds_key):
    """
    Map a legacy short key (arabidopsis, at_root, at_seed, ...) to the
    canonical ds_key used by mini_app.py / generate_umap_dumps.py. Returns
    ds_key unchanged if it's already canonical or unrecognized (callers
    still get a clear KeyError from DATASETS[ds_key] in that case).
    """
    return LEGACY_KEY_ALIASES.get(ds_key, ds_key)


def resolve_columns(ds_key, col_key=None):
    """
    Return (cell_type_col, cell_type_col_2, svg_path) for a given dataset
    and (optional) requested column key. Used by mini_app.py in place of
    its old hardcoded if/elif ds_key/col_key dispatch chain.

    - col_key=None (or equal to the dataset's own cell_type_col) returns
      the dataset's default columns.
    - col_key matching a key under `alt_columns` in the YAML returns that
      alternate's columns instead (rice's CAxCondition, at_root_rs's
      TypexGenotype, at_seed_martin's level_3 view).
    - Any other col_key is passed straight through as cell_type_col with
      no secondary column (single-column datasets like arabidopsis_nat).
    """
    if ds_key not in DATASETS:
        raise KeyError(f"Unknown dataset key: {ds_key!r}. "
                        f"Known keys: {sorted(DATASETS)}")

    cfg = DATASETS[ds_key]
    svg_path = cfg["svg"]

    alt = cfg.get("alt_columns", {})
    if col_key and col_key in alt:
        alt_cfg = alt[col_key]
        col1 = alt_cfg.get("cell_type_col", cfg.get("cell_type_col"))
        col2 = alt_cfg.get("cell_type_col_2")
        svg_key = alt_cfg.get("svg")
        if svg_key:
            svg_path = cfg.get(svg_key, svg_path)
        return col1, col2, svg_path

    if col_key is None or col_key == cfg.get("cell_type_col"):
        return cfg.get("cell_type_col"), cfg.get("cell_type_col_2"), svg_path

    return col_key, None, svg_path


def umap_col_for(ds_key):
    return DATASETS[ds_key].get("umap_col")


def pseudobulk_columns(ds_key):
    """
    Return the dict generate_pseudobulk_dumps.py needs for one dataset,
    accepting either a canonical ds_key or one of its legacy short-key
    aliases (arabidopsis, at_root, at_seed, at_flower, at_silique, at_stem).
    """
    key = resolve_ds_key(ds_key)
    if key not in DATASETS:
        raise KeyError(f"Unknown dataset key: {ds_key!r}. "
                        f"Known keys: {sorted(DATASETS)} "
                        f"(legacy aliases: {sorted(LEGACY_KEY_ALIASES)})")
    cfg = DATASETS[key]
    return {
        "ds_key": key,
        "h5ad": cfg["h5ad"],
        "cell_type_col": cfg.get("cell_type_col"),
        "cell_type_col2": cfg.get("cell_type_col_2"),
        "gene_id_col": cfg.get("gene_id_col"),
        "use_logcounts": cfg.get("use_logcounts", False),
        "post_process": cfg.get("post_process"),
        "db_name": cfg.get("db_name_pseudobulk"),
    }


def umap_columns(ds_key):
    """Return the dict generate_umap_dumps.py needs for one dataset."""
    key = resolve_ds_key(ds_key)
    if key not in DATASETS:
        raise KeyError(f"Unknown dataset key: {ds_key!r}. "
                        f"Known keys: {sorted(DATASETS)}")
    cfg = DATASETS[key]
    return {
        "ds_key": key,
        "h5ad": cfg["h5ad"],
        "umap_col": cfg.get("umap_col"),
        "umap_col2": cfg.get("umap_col2"),
        "gene_id_col": cfg.get("gene_id_col"),
        "db_name": cfg.get("db_name_umap"),
    }


def pseudobulk_ready_keys():
    """ds_keys that have a db_name_pseudobulk configured."""
    return [k for k, v in DATASETS.items() if v.get("db_name_pseudobulk")]


def umap_ready_keys():
    """ds_keys that have a db_name_umap configured."""
    return [k for k, v in DATASETS.items() if v.get("db_name_umap")]


def reload():
    """Re-read dataset_config.yaml from disk (useful after editing it
    without restarting a long-running process, e.g. in a notebook)."""
    global _raw, DATASETS, LEGACY_KEY_ALIASES, ROOT_GENO_PREFIX, ROOT_H5AD_TO_ID
    _raw = _load()
    DATASETS = _raw.get("datasets", {})
    LEGACY_KEY_ALIASES = _raw.get("legacy_key_aliases", {})
    ROOT_GENO_PREFIX = _raw.get("root_genotype_prefix", {})
    ROOT_H5AD_TO_ID = _raw.get("root_h5ad_to_svg_ids", {})
    return DATASETS
