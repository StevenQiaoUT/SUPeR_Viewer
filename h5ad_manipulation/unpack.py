#!/usr/bin/env python3
"""
unpack.py — h5ad metadata sniffer
==================================

Every new dataset export names its cell-type column differently
("CellType", "label_major", "label_majorXcondition", "Celltype",
"celltype_after", ...), so a fully-automated `unpack.py <file> --gene X`
CLI just trades one failure mode ("wrong hardcoded column") for another
("guessed the wrong column silently"). This script does NOT try to guess
and run — it's a read-only report: point it at an .h5ad and it prints
every obs column, flags the ones that *look* like cell-type/tissue labels
using a small plant-biology vocabulary, and shows sample values so you can
pick the right one by eye in a few seconds instead of trial-and-error runs
against the full pipeline.

Usage
-----
    python unpack.py path/to/file.h5ad
    python unpack.py path/to/file.h5ad --top 8
    python unpack.py path/to/file.h5ad --vocab suspensor cambial
    python unpack.py path/to/file.h5ad --yaml-stub my_dataset_key

Nothing here writes to the h5ad or touches X/layers beyond checking their
names — this is inspection only.
"""
import argparse
import sys

# ── Vocabulary ──────────────────────────────────────────────────────────────
# Deliberately a flat, case-insensitive substring list, not an ontology.
# Extend this as new datasets introduce new tissue/cell-type terms —
# the point is to catch *most* obs columns worth looking at, not all of them.
CELL_TYPE_VOCAB = [
    "mesophyll", "epidermis", "epidermal", "guard cell", "guard_cell",
    "trichome", "phloem", "xylem", "cortex", "endodermis", "pericycle",
    "cambium", "cambial", "parenchyma", "meristem", "meristematic", "sieve",
    "companion cell", "stele", "columella", "root cap", "root_cap",
    "procambium", "vascular", "stomat", "collenchyma", "sclerenchyma",
    "tracheid", "protoxylem", "metaxylem", "seed coat", "endosperm",
    "embryo", "cotyledon", "hypocotyl", "shoot apical", "suspensor",
    "unlabeled", "unlabelled", "unknown", "quiescent",
]

# Columns whose *names* strongly suggest cell-type content, even before
# checking values — bumps a column up the ranking without requiring a
# vocab hit in the values (useful for coded/abbreviated labels).
NAME_HINTS = ["celltype", "cell_type", "cell type", "annotation", "label", "type"]

# Columns whose names suggest they're never a cell-type column, no matter
# what their values look like — QC/technical metadata gets deprioritized
# even if a stray vocab word shows up in a free-text field.
NAME_DEPRIORITIZE = [
    "n_genes", "n_counts", "pct_counts", "total_counts", "doublet",
    "barcode", "sample", "batch", "condition", "genotype", "replicate",
    "timepoint", "index", "leiden", "louvain", "cluster",
]


def score_column(values, vocab):
    """Return (# unique values containing a vocab hit, matched vocab terms)."""
    hits = set()
    matched_terms = set()
    for v in values:
        vl = str(v).lower()
        for term in vocab:
            if term in vl:
                hits.add(v)
                matched_terms.add(term)
    return len(hits), sorted(matched_terms)


def name_bonus(col_name):
    cl = col_name.lower()
    if any(bad in cl for bad in NAME_DEPRIORITIZE):
        return -1
    if any(good in cl for good in NAME_HINTS):
        return 1
    return 0


def main():
    parser = argparse.ArgumentParser(
        description="Sniff an h5ad's obs columns for likely cell-type labels."
    )
    parser.add_argument("h5ad_path")
    parser.add_argument("--top", type=int, default=6,
                         help="Show this many top candidate columns in detail (default 6)")
    parser.add_argument("--vocab", nargs="*", default=[],
                         help="Extra vocabulary terms to check for, on top of the built-in list")
    parser.add_argument("--max-samples", type=int, default=12,
                         help="Max sample values to print per column (default 12)")
    parser.add_argument("--yaml-stub", metavar="DATASET_KEY", default=None,
                         help="Print a dataset_config.yaml stanza for the "
                              "top-ranked column under this dataset key")
    args = parser.parse_args()

    try:
        import anndata as ad
    except ImportError:
        print("This script needs anndata (`pip install anndata`).", file=sys.stderr)
        sys.exit(1)

    vocab = CELL_TYPE_VOCAB + [v.lower() for v in args.vocab]

    print(f"Reading {args.h5ad_path} ...", file=sys.stderr)
    # backed="r" avoids loading the (often huge) X matrix into memory just
    # to look at obs — this is a metadata report, not an analysis step.
    adata = ad.read_h5ad(args.h5ad_path, backed="r")
    obs = adata.obs

    print(f"\n{adata.n_obs:,} cells x {adata.n_vars:,} genes")
    print(f"obsm keys : {list(adata.obsm.keys())}")
    print(f"layers    : {list(adata.layers.keys())}")
    print(f"uns keys  : {list(adata.uns.keys())}")

    print(f"\n{'=' * 78}\nobs column report ({len(obs.columns)} columns)\n{'=' * 78}")

    scored = []
    for col in obs.columns:
        series = obs[col]
        dtype = series.dtype
        is_textlike = dtype == object or str(dtype).startswith("category")
        n_unique = series.nunique(dropna=True)

        hit_count, matched_terms = (0, [])
        if is_textlike:
            uniques = series.dropna().astype(str).unique().tolist()
            hit_count, matched_terms = score_column(uniques, vocab)
        else:
            uniques = []

        bonus = name_bonus(col)
        # Cardinality sanity check: cell-type columns are rarely 2-3
        # categories (that's usually condition/genotype) or thousands
        # (that's usually a barcode or a numeric ID cast to string).
        cardinality_ok = is_textlike and 3 <= n_unique <= 80

        rank_score = (hit_count * 2) + bonus + (1 if cardinality_ok else 0)
        scored.append({
            "col": col, "dtype": dtype, "n_unique": n_unique,
            "is_textlike": is_textlike, "hit_count": hit_count,
            "matched_terms": matched_terms, "uniques": uniques,
            "rank_score": rank_score,
        })

    scored.sort(key=lambda r: r["rank_score"], reverse=True)

    print(f"\nTop {min(args.top, len(scored))} candidates, ranked by vocabulary "
          f"match + column-name hints + cardinality:\n")
    for r in scored[:args.top]:
        marker = "*" if r["hit_count"] > 0 else " "
        print(f"{marker} {r['col']!r:38s} dtype={str(r['dtype']):10s} "
              f"n_unique={r['n_unique']:5d}  score={r['rank_score']}")
        if r["hit_count"]:
            shown = ", ".join(r["matched_terms"][:6])
            more = "..." if len(r["matched_terms"]) > 6 else ""
            print(f"      vocab hits: {r['hit_count']} values matched ({shown}{more})")
        if r["is_textlike"]:
            sample = sorted(r["uniques"])[:args.max_samples]
            trailing = "" if len(r["uniques"]) <= args.max_samples else ", ..."
            print(f"      sample values: {sample}{trailing}")
        print()

    if len(scored) > args.top:
        print(f"{'-' * 78}\nOther columns (lower-ranked)\n{'-' * 78}")
        for r in scored[args.top:]:
            print(f"  {r['col']!r:38s} dtype={str(r['dtype']):10s} n_unique={r['n_unique']}")

    if args.yaml_stub:
        best = scored[0] if scored else None
        key = args.yaml_stub
        col_name = best["col"] if best else "REPLACE_ME"
        print(f"\n{'=' * 78}\nSuggested dataset_config.yaml stanza\n{'=' * 78}")
        print(f"""
{key}:
  label: "REPLACE_ME (e.g. Species Tissue (Author et al., Year))"
  h5ad: "{args.h5ad_path}"
  svg: "REPLACE_ME path/to/template.svg"
  cell_type_col: "{col_name}"   # <- verify this is right before committing
  umap_col: "{col_name}"
  opacity: 1.0
""".strip())
        print("\nDouble-check the sample values above actually read as cell "
              "types before using this column name — this is a starting "
              "guess, not a guarantee.")


if __name__ == "__main__":
    main()
