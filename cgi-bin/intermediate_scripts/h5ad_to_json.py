#!/mnt/home/sqiao/venv/bin/python3

"""
Convert H5AD to JSON with Average Expression by Cell Type

Computes average gene expression for each cell type and exports to JSON format.
Can be run as a command-line script or as a CGI script on a web server.

Usage (command-line):
    python h5ad_to_json.py <input_h5ad> <output_json> <cell_type_column> [gene1 gene2 ...]

Usage (CGI / browser):
    http://localhost:8080/cgi-bin/h5ad_to_json.py?file=/Users/stevenqiao/Desktop/at.h5ad&col=label_majorXcondition&gene=AT3G05727

Arguments (command-line):
    input_h5ad        : Path to input .h5ad file
    output_json       : Path to output JSON file
    cell_type_column  : Column name in adata.obs containing cell type labels
    gene1 gene2 ...   : (Optional) List of genes to include. If not provided, uses all genes.

Query Parameters (CGI):
    file              : Filename of the .h5ad file
    col               : Column name in adata.obs containing cell type labels
    gene              : Gene ID to include (repeat for multiple genes)

Examples (command-line):
    # Single gene
    python h5ad_to_json.py at.h5ad output.json label_majorXcondition AT3G05727

    # Multiple genes
    python h5ad_to_json.py at.h5ad output.json label_majorXcondition AT3G05727 AT1G01010 AT5G12345

    # All genes (no gene list specified)
    python h5ad_to_json.py at.h5ad output.json label_v2

    # Different paths
    python h5ad_to_json.py /path/to/data.h5ad /path/to/output.json cell_type
"""

import sys
import os
import json


def h5ad_to_json_avg(input_h5ad, output_json, cell_type_column="cell_type", gene_list=None):
    """
    Convert h5ad file to JSON with averaged expression per cell type.

    Parameters
    ----------
    input_h5ad : str
        Path to input .h5ad file
    output_json : str
        Path to output JSON file
    cell_type_column : str
        Column name in adata.obs containing cell type labels
    gene_list : list of str or None
        List of genes to include. If None, uses all genes.

    Returns
    -------
    list of dict
        The output records (also written to output_json).
    """
    import anndata as ad
    import numpy as np

    print(f"Loading data from: {input_h5ad}", file=sys.stderr)
    try:
        adata = ad.read_h5ad(input_h5ad, backed='r')
    except Exception as e:
        print(f"Error reading h5ad file: {e}", file=sys.stderr)
        sys.exit(1)

    print(f"  Total cells: {adata.n_obs:,}", file=sys.stderr)
    print(f"  Total genes: {adata.n_vars:,}", file=sys.stderr)

    # Validate cell type column
    if cell_type_column not in adata.obs.columns:
        print(f"\nError: Column '{cell_type_column}' not found in dataset", file=sys.stderr)
        print(f"Available columns: {list(adata.obs.columns)}", file=sys.stderr)
        sys.exit(1)

    # Validate and filter gene list
    if gene_list is None:
        gene_list = list(adata.var_names)
        print(f"  Using all {len(gene_list):,} genes", file=sys.stderr)
    else:
        missing = [g for g in gene_list if g not in adata.var_names]
        if missing:
            print(f"\nWarning: Genes not found in dataset: {missing}", file=sys.stderr)
            print(f"Available genes (first 10): {list(adata.var_names[:10])}", file=sys.stderr)
        gene_list = [g for g in gene_list if g in adata.var_names]
        if not gene_list:
            print("\nError: No valid genes specified", file=sys.stderr)
            sys.exit(1)
        print(f"  Using {len(gene_list)} gene(s): {gene_list}", file=sys.stderr)

    # KEY FIX: slice only the needed gene columns BEFORE densifying.
    # This avoids loading the full 144K x 2000 matrix into memory —
    # for a single gene query it's ~2000x less data.
    print("\nSlicing gene columns...", file=sys.stderr)
    gene_indices = [adata.var_names.get_loc(g) for g in gene_list]
    X_slice = adata.X[:, gene_indices]
    if not isinstance(X_slice, __import__('numpy').ndarray):
        X_slice = X_slice.toarray()

    cell_types = adata.obs[cell_type_column].values
    unique_cell_types = adata.obs[cell_type_column].unique()
    print(f"  Cell types found: {len(unique_cell_types)}", file=sys.stderr)

    # Compute average expression per cell type
    print("Computing average expression per cell type...", file=sys.stderr)
    output = []
    for ct in unique_cell_types:
        mask = cell_types == ct
        avg = X_slice[mask].mean(axis=0)
        for i, gene in enumerate(gene_list):
            output.append({
                "cell_type": str(ct),
                "gene": gene,
                "average_expression": float(avg[i])
            })

    print(f"  Computed {len(output):,} records", file=sys.stderr)

    # Write JSON
    try:
        with open(output_json, "w") as f:
            json.dump(output, f, indent=4)
        print(f"\n✓ Saved to: {output_json}", file=sys.stderr)
    except Exception as e:
        print(f"\nError writing JSON: {e}", file=sys.stderr)
        sys.exit(1)

    return output


def main():
    # Detect CGI mode: Apache sets REQUEST_METHOD in the environment
    is_cgi = "REQUEST_METHOD" in os.environ

    if is_cgi:
        # --------------------------------------------------------------- #
        # CGI mode                                                          #
        # CRITICAL: flush headers IMMEDIATELY so Apache does not time out  #
        # before the heavy H5AD loading is done.                           #
        # --------------------------------------------------------------- #
        sys.stdout.write("Content-Type: application/json\r\n\r\n")
        sys.stdout.flush()

        # Parse query string without the deprecated `cgi` module
        from urllib.parse import parse_qs
        query_string = os.environ.get("QUERY_STRING", "")
        params = parse_qs(query_string)

        input_h5ad       = params.get("file", [None])[0]
        cell_type_column = params.get("col",  [None])[0]
        gene_list        = params.get("gene", None)  # list or None

        if not input_h5ad or not cell_type_column:
            json.dump(
                {"error": "Missing required query parameters: file, col"},
                sys.stdout
            )
            sys.stdout.flush()
            return

        # Run computation and stream JSON directly to stdout
        import anndata as ad
        import numpy as np

        print(f"Loading data from: {input_h5ad}", file=sys.stderr)
        try:
            adata = ad.read_h5ad(input_h5ad, backed='r')
        except Exception as e:
            json.dump({"error": f"Could not read H5AD file: {e}"}, sys.stdout)
            sys.stdout.flush()
            return

        print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes", file=sys.stderr)

        if cell_type_column not in adata.obs.columns:
            json.dump(
                {"error": f"Column '{cell_type_column}' not found",
                 "available_columns": list(adata.obs.columns)},
                sys.stdout
            )
            sys.stdout.flush()
            return

        if gene_list:
            missing = [g for g in gene_list if g not in adata.var_names]
            if missing:
                print(f"Warning: genes not found: {missing}", file=sys.stderr)
            gene_list = [g for g in gene_list if g in adata.var_names]
        else:
            gene_list = list(adata.var_names)

        if not gene_list:
            json.dump({"error": "No valid genes specified"}, sys.stdout)
            sys.stdout.flush()
            return

        print(f"  Genes: {gene_list}", file=sys.stderr)

        # Slice only needed gene columns before densifying
        gene_indices = [adata.var_names.get_loc(g) for g in gene_list]
        X_slice = adata.X[:, gene_indices]
        if not isinstance(X_slice, np.ndarray):
            X_slice = X_slice.toarray()

        cell_types = adata.obs[cell_type_column].values
        unique_cell_types = adata.obs[cell_type_column].unique()

        print("Computing averages...", file=sys.stderr)
        output = []
        for ct in unique_cell_types:
            mask = cell_types == ct
            avg = X_slice[mask].mean(axis=0)
            for i, gene in enumerate(gene_list):
                output.append({
                    "cell_type": str(ct),
                    "gene": gene,
                    "average_expression": float(avg[i])
                })

        json.dump(output, sys.stdout, indent=4)
        sys.stdout.flush()
        print(f"\n✓ Returned {len(output):,} records", file=sys.stderr)

    else:
        # --------------------------------------------------------------- #
        # Command-line mode                                                 #
        # --------------------------------------------------------------- #
        if len(sys.argv) < 4:
            print("Usage: python h5ad_to_json.py <input_h5ad> <output_json> <cell_type_column> [gene1 gene2 ...]")
            sys.exit(1)

        input_h5ad       = sys.argv[1]
        output_json      = sys.argv[2]
        cell_type_column = sys.argv[3]
        gene_list        = sys.argv[4:] if len(sys.argv) > 4 else None

        if not os.path.exists(input_h5ad):
            print(f"Error: Input file '{input_h5ad}' does not exist")
            sys.exit(1)
        if not input_h5ad.endswith('.h5ad'):
            print(f"Warning: '{input_h5ad}' does not have .h5ad extension")

        print("=" * 60)
        print("H5AD to JSON Converter — Average Expression")
        print("=" * 60)
        print(f"Input file:       {input_h5ad}")
        print(f"Output file:      {output_json}")
        print(f"Cell type column: {cell_type_column}")
        print(f"Gene list:        {', '.join(gene_list) if gene_list else 'All genes'}")
        print("=" * 60)
        print()

        h5ad_to_json_avg(
            input_h5ad=input_h5ad,
            output_json=output_json,
            cell_type_column=cell_type_column,
            gene_list=gene_list
        )
        print("\nDone!")


if __name__ == "__main__":
    main()

