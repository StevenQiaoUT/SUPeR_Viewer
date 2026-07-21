#!/mnt/home/sqiao/venv/bin/python3
"""
eFP SVG Generator — CGI Script
================================
Reads an H5AD file, computes average gene expression per cell type, and
returns a colored SVG directly — no intermediate JSON file needed.

Usage (CGI / browser):
    http://142.150.215.219/~sqiao/cgi-bin/efp_svg.py?file=/path/to/data.h5ad&col=label_majorXcondition&gene=AT3G05727&svg=/path/to/new_grid.svg

Query Parameters:
    file   : Absolute path to the .h5ad file on the server
    col    : Column name in adata.obs containing cell type labels
    gene   : Gene ID (repeat for multiple genes; first gene is used for coloring)
    svg    : Absolute path to the template SVG file on the server
    opacity: (Optional) Fill opacity 0.0–1.0 (default 1.0)

Usage (command-line):
    python efp_svg.py <h5ad_file> <svg_template> <cell_type_col> <gene> [output.svg] [--opacity 0.9]
"""

import sys
import os
import json
from xml.etree import ElementTree as ET


# ──────────────────────────────────────────────────────────────────────────────
# Color utilities  (from color_guard_cells_cli.py)
# ──────────────────────────────────────────────────────────────────────────────

def expression_to_color(value, max_val, scheme='yellow_to_red'):
    if scheme == 'yellow_to_red':
        minColor = {'red': 255, 'green': 255, 'blue': 0}
    else:  # white_to_red
        minColor = {'red': 255, 'green': 255, 'blue': 255}
    maxColor = {'red': 255, 'green': 0, 'blue': 0}

    ratio = value / max_val if max_val > 0 else 0
    if ratio > 1:
        ratio = 1
    if not (0 <= ratio <= 1):
        return "#ffffff"

    red   = minColor['red']   + round((maxColor['red']   - minColor['red'])   * ratio)
    green = minColor['green'] + round((maxColor['green'] - minColor['green']) * ratio)
    blue  = minColor['blue']  + round((maxColor['blue']  - minColor['blue'])  * ratio)
    return f"#{red:02x}{green:02x}{blue:02x}"


def compute_global_shift(expression_dict):
    all_values = list(expression_dict.values())
    min_val = min(all_values) if all_values else 0
    return -min_val if min_val < 0 else 0


# ──────────────────────────────────────────────────────────────────────────────
# Legend + title helpers  (from color_guard_cells_cli.py, unchanged)
# ──────────────────────────────────────────────────────────────────────────────

def add_legend(root, gene_name, colour_ceiling, shift, n_boxes=10):
    SVG_NS = "http://www.w3.org/2000/svg"
    box_w        = 24
    box_h        = 20
    tick_gap     = 6
    font_size    = 10
    header_fs    = 11
    left_pad     = 10
    right_pad    = 20
    tick_label_w = 40
    legend_width = left_pad + box_w + tick_gap + tick_label_w + right_pad

    width_attr = root.get("width", "")
    try:
        old_w  = float("".join(c for c in width_attr if c.isdigit() or c == "."))
        w_unit = "".join(c for c in width_attr if c.isalpha())
        root.set("width", f"{old_w + legend_width}{w_unit}")
    except ValueError:
        pass

    vb = root.get("viewBox", "")
    if vb:
        try:
            parts    = vb.split()
            parts[2] = str(float(parts[2]) + legend_width)
            root.set("viewBox", " ".join(parts))
        except (IndexError, ValueError):
            pass

    existing = [c for c in list(root) if c.get("id") not in ("expression_legend",)]
    for child in existing:
        root.remove(child)
    shift_g = ET.Element(f"{{{SVG_NS}}}g")
    shift_g.set("transform", f"translate({legend_width}, 0)")
    for child in existing:
        shift_g.append(child)
    root.append(shift_g)

    x0 = left_pad
    height_attr = root.get("height", "")
    try:
        canvas_h = float("".join(c for c in height_attr if c.isdigit() or c == "."))
    except ValueError:
        try:
            canvas_h = float(root.get("viewBox", "0 0 800 600").split()[3])
        except (IndexError, ValueError):
            canvas_h = 600

    header_lines = 3
    header_h     = header_lines * (header_fs + 4)
    legend_h     = header_h + n_boxes * box_h + font_size + 6
    y0           = canvas_h - legend_h - 10

    g = ET.Element(f"{{{SVG_NS}}}g")
    g.set("id", "expression_legend")

    def txt(x, y, text, fs=font_size, anchor="start", weight="normal", style=""):
        el = ET.Element(f"{{{SVG_NS}}}text")
        el.set("x", str(x)); el.set("y", str(y))
        el.set("font-family", "Arial, sans-serif")
        el.set("font-size", str(fs))
        el.set("font-weight", weight)
        el.set("text-anchor", anchor)
        el.set("dominant-baseline", "auto")
        el.set("fill", "#222222")
        if style:
            el.set("font-style", style)
        el.text = text
        return el

    hy = y0 + header_fs
    g.append(txt(x0, hy, gene_name, fs=header_fs, weight="bold", style="italic"))
    hy += header_fs + 4
    g.append(txt(x0, hy, "Single Cell Max", fs=font_size))
    hy += font_size + 4
    g.append(txt(x0, hy, "Linear", fs=font_size))

    boxes_y0   = y0 + header_h
    raw_top    = colour_ceiling - shift
    raw_bottom = -shift

    for i in range(n_boxes):
        ratio_top    = 1.0 - i       / n_boxes
        ratio_bottom = 1.0 - (i + 1) / n_boxes
        ratio_mid    = (ratio_top + ratio_bottom) / 2
        color        = expression_to_color(ratio_mid * colour_ceiling, colour_ceiling, "yellow_to_red")

        bx = x0
        by = boxes_y0 + i * box_h
        rect = ET.Element(f"{{{SVG_NS}}}rect")
        rect.set("x", str(bx)); rect.set("y", str(by))
        rect.set("width", str(box_w)); rect.set("height", str(box_h))
        rect.set("fill", color); rect.set("stroke", "#888888")
        rect.set("stroke-width", "0.5")
        g.append(rect)

        raw_tick  = raw_bottom + ratio_top * (raw_top - raw_bottom)
        ty        = by + 4
        g.append(txt(bx + box_w + tick_gap, ty, f"{raw_tick:.2f}", fs=font_size - 1))

    bottom_y = boxes_y0 + n_boxes * box_h + 4
    g.append(txt(x0 + box_w + tick_gap, bottom_y, f"{raw_bottom:.2f}", fs=font_size - 1))
    root.append(g)


def add_gene_title(root, gene_name, font_size=18, padding_bottom=10):
    SVG_NS   = 'http://www.w3.org/2000/svg'
    offset_y = font_size + padding_bottom

    existing_title = root.find(f'{{{SVG_NS}}}title')
    if existing_title is not None:
        root.remove(existing_title)
    title_elem      = ET.Element(f'{{{SVG_NS}}}title')
    title_elem.text = gene_name
    root.insert(0, title_elem)

    height_attr = root.get('height', '')
    if height_attr:
        try:
            new_h = float(''.join(c for c in height_attr if c.isdigit() or c == '.')) + offset_y
            unit  = ''.join(c for c in height_attr if c.isalpha())
            root.set('height', f"{new_h}{unit}")
        except ValueError:
            pass

    vb = root.get('viewBox', '')
    if vb:
        try:
            parts    = vb.split()
            parts[3] = str(float(parts[3]) + offset_y)
            root.set('viewBox', ' '.join(parts))
        except (IndexError, ValueError):
            pass

    children = [c for c in list(root) if c is not title_elem]
    for child in children:
        root.remove(child)
    wrapper = ET.Element(f'{{{SVG_NS}}}g')
    wrapper.set('transform', f'translate(0, {offset_y})')
    for child in children:
        wrapper.append(child)
    root.append(wrapper)

    width_attr = root.get('width', '')
    cx = None
    if width_attr:
        try:
            cx = float(''.join(c for c in width_attr if c.isdigit() or c == '.')) / 2
        except ValueError:
            pass
    if cx is None:
        try:
            cx = float(vb.split()[2]) / 2
        except (IndexError, ValueError):
            cx = 300

    text_elem = ET.Element(f'{{{SVG_NS}}}text')
    text_elem.set('x', str(cx)); text_elem.set('y', str(font_size))
    text_elem.set('text-anchor', 'middle')
    text_elem.set('dominant-baseline', 'auto')
    text_elem.set('font-family', 'Arial, sans-serif')
    text_elem.set('font-size', str(font_size))
    text_elem.set('font-weight', 'bold')
    text_elem.set('fill', '#222222')
    text_elem.text = gene_name
    root.insert(1, text_elem)


# ──────────────────────────────────────────────────────────────────────────────
# Core coloring logic  (from color_guard_cells_cli.py, adapted to work in-memory)
# ──────────────────────────────────────────────────────────────────────────────

def color_svg(svg_file, expression_dict, gene_name=None, opacity=1.0):
    """
    Parse svg_file, apply expression colors, return the modified ET tree.
    No files are written — the caller serializes to wherever it needs.
    """
    tree = ET.parse(svg_file)
    root = tree.getroot()

    conditions = ['D0', 'W0', 'R15', 'W15']
    SVG_NS     = 'http://www.w3.org/2000/svg'

    shift        = compute_global_shift(expression_dict)
    shifted_dict = {k: v + shift for k, v in expression_dict.items()}

    visualized_keys = (
        [f"{c}_Guard"             for c in conditions] +
        [f"{c}_Mesophyll"         for c in conditions] +
        [f"{c}_Epidermal"         for c in conditions] +
        [f"{c}_Trichome"          for c in conditions] +
        [f"{c}_Vascular"          for c in conditions] +
        [f"{c}_Phloem Parenchyma" for c in conditions] +
        [f"{c}_Phloem companion"  for c in conditions]
    )
    vis_vals       = [shifted_dict[k] for k in visualized_keys if k in shifted_dict]
    colour_ceiling = max(vis_vals) if vis_vals else 1.0

    print(f"  Shift: {shift:.4f}  Colour ceiling: {colour_ceiling:.4f}", file=sys.stderr)

    def get_color(key):
        val = shifted_dict.get(key)
        return expression_to_color(val, colour_ceiling, 'yellow_to_red') if val is not None else None

    def raw_val(key):
        """Return the original (unshifted) average expression value, formatted for display."""
        val = expression_dict.get(key)
        return f"{val:.4f}" if val is not None else "N/A"

    def add_tooltip(elem, text):
        """Insert a <title> child as the first child so the browser shows it on hover."""
        title_el = ET.Element(f'{{{SVG_NS}}}title')
        title_el.text = text
        elem.insert(0, title_el)

    # Guard cells
    for condition in conditions:
        key   = f"{condition}_Guard"
        color = get_color(key)
        if not color: continue
        tooltip = f"{key}\nAvg expression: {raw_val(key)}"
        for elem in root.iter():
            if elem.get('id') == f"{condition}_guard":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    add_tooltip(path, tooltip)

    # Mesophyll
    for condition in conditions:
        key   = f"{condition}_Mesophyll"
        color = get_color(key)
        if not color: continue
        tooltip = f"{key}\nAvg expression: {raw_val(key)}"
        for elem in root.iter():
            if elem.get('id') == f"mesophyll_{condition}":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    add_tooltip(path, tooltip)

    # Epidermal
    for condition in conditions:
        key   = f"{condition}_Epidermal"
        color = get_color(key)
        if not color: continue
        tooltip = f"{key}\nAvg expression: {raw_val(key)}"
        for elem in root.iter():
            if elem.get('id') == f"{condition}_epidermal":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    add_tooltip(path, tooltip)

    # Trichome
    for condition in conditions:
        key   = f"{condition}_Trichome"
        color = get_color(key)
        if not color: continue
        tooltip = f"{key}\nAvg expression: {raw_val(key)}"
        for elem in root.iter():
            if elem.get('id') == f"{condition}_trichome":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    add_tooltip(path, tooltip)

    # Vascular circles
    for condition in conditions:
        key   = f"{condition}_Vascular"
        color = get_color(key)
        if not color: continue
        tooltip = f"{key}\nAvg expression: {raw_val(key)}"
        for elem in root.iter():
            if elem.get('id') == f"vascular_{condition}":
                for circle in elem.iter(f'{{{SVG_NS}}}circle'):
                    circle.set('fill', color)
                    circle.set('fill-opacity', str(opacity))
                    add_tooltip(circle, tooltip)

    # Phloem Parenchyma circles
    for condition in conditions:
        key   = f"{condition}_Phloem Parenchyma"
        color = get_color(key)
        if not color: continue
        tooltip = f"{key}\nAvg expression: {raw_val(key)}"
        for elem in root.iter():
            if elem.get('id') == f"phloem_parenchyma_{condition}" and elem.tag == f'{{{SVG_NS}}}circle':
                elem.set('fill', color)
                elem.set('fill-opacity', str(opacity))
                add_tooltip(elem, tooltip)

    # Phloem Companion circles
    for condition in conditions:
        key   = f"{condition}_Phloem companion"
        color = get_color(key)
        if not color: continue
        tooltip = f"{key}\nAvg expression: {raw_val(key)}"
        for elem in root.iter():
            if elem.get('id') == f"phloem_companion_{condition}" and elem.tag == f'{{{SVG_NS}}}circle':
                elem.set('fill', color)
                elem.set('fill-opacity', str(opacity))
                add_tooltip(elem, tooltip)

    # Phloem paths — average of Phloem Parenchyma + Phloem companion
    for condition in conditions:
        pp_key = f"{condition}_Phloem Parenchyma"
        pc_key = f"{condition}_Phloem companion"
        if pp_key not in expression_dict or pc_key not in expression_dict:
            continue
        pp_raw      = expression_dict[pp_key]
        pc_raw      = expression_dict[pc_key]
        avg_raw     = (pp_raw + pc_raw) / 2
        avg_shifted = avg_raw + shift
        color       = expression_to_color(avg_shifted, colour_ceiling, 'yellow_to_red')
        tooltip     = (f"{condition}_Phloem (averaged)\n"
                       f"Phloem Parenchyma: {pp_raw:.4f}\n"
                       f"Phloem Companion:  {pc_raw:.4f}\n"
                       f"Avg expression:    {avg_raw:.4f}")
        for elem in root.iter():
            if elem.get('id') == f"{condition}_phloem":
                for path in elem.iter(f'{{{SVG_NS}}}path'):
                    path.set('fill', color)
                    path.set('fill-opacity', str(opacity))
                    add_tooltip(path, tooltip)

    # Title + legend
    if gene_name:
        add_gene_title(root, gene_name)
    add_legend(root, gene_name or "Expression", colour_ceiling, shift)

    return tree


# ──────────────────────────────────────────────────────────────────────────────
# H5AD → expression dict  (from h5ad_to_json.py, returns dict directly)
# ──────────────────────────────────────────────────────────────────────────────

def compute_expression(h5ad_path, cell_type_column, gene_list):
    """
    Load h5ad, slice requested genes, compute per-cell-type averages.
    Returns a dict: { cell_type_label -> average_expression } for the first gene,
    plus the gene name used.
    Raises ValueError on bad inputs.
    """
    import anndata as ad
    import numpy as np

    print(f"Loading {h5ad_path}", file=sys.stderr)
    adata = ad.read_h5ad(h5ad_path, backed='r')
    print(f"  {adata.n_obs:,} cells × {adata.n_vars:,} genes", file=sys.stderr)

    if cell_type_column not in adata.obs.columns:
        raise ValueError(f"Column '{cell_type_column}' not found. "
                         f"Available: {list(adata.obs.columns)}")

    # Validate genes
    missing = [g for g in gene_list if g not in adata.var_names]
    if missing:
        print(f"Warning: genes not in dataset: {missing}", file=sys.stderr)
    gene_list = [g for g in gene_list if g in adata.var_names]
    if not gene_list:
        raise ValueError("No valid genes found in dataset.")

    print(f"  Genes: {gene_list}", file=sys.stderr)

    # Slice only needed columns before densifying (huge performance win)
    gene_indices = [adata.var_names.get_loc(g) for g in gene_list]
    X_slice      = adata.X[:, gene_indices]
    if not isinstance(X_slice, np.ndarray):
        X_slice = X_slice.toarray()

    cell_types        = adata.obs[cell_type_column].values
    unique_cell_types = adata.obs[cell_type_column].unique()

    print("Computing averages...", file=sys.stderr)

    # Build expression dict: cell_type -> avg expression (first gene only for SVG coloring)
    # If multiple genes were requested, we use only the first for the SVG.
    # The label format expected by the SVG coloring is e.g. "D0_Guard".
    expression_dict = {}
    for ct in unique_cell_types:
        mask = cell_types == ct
        avg  = X_slice[mask].mean(axis=0)
        for i, gene in enumerate(gene_list):
            # Store every gene, but we only use gene_list[0] for coloring
            expression_dict[str(ct)] = float(avg[0])   # first gene

    print(f"  {len(expression_dict)} cell types computed", file=sys.stderr)
    return expression_dict, gene_list[0]


# ──────────────────────────────────────────────────────────────────────────────
# Main entry point
# ──────────────────────────────────────────────────────────────────────────────

def main():
    is_cgi = "REQUEST_METHOD" in os.environ

    if is_cgi:
        # ── CGI mode ─────────────────────────────────────────────────────────
        # Flush SVG content-type header IMMEDIATELY to prevent Apache timeout
        sys.stdout.write("Content-Type: image/svg+xml\r\n\r\n")
        sys.stdout.flush()

        from urllib.parse import parse_qs
        params = parse_qs(os.environ.get("QUERY_STRING", ""))

        h5ad_path        = params.get("file",    [None])[0]
        cell_type_column = params.get("col",     [None])[0]
        gene_list        = params.get("gene",    None)
        svg_template     = params.get("svg",     [None])[0]
        opacity          = float(params.get("opacity", ["1.0"])[0])

        def svg_error(msg):
            """Return a minimal SVG with an error message."""
            sys.stdout.write(
                f'<svg xmlns="http://www.w3.org/2000/svg" width="600" height="60">'
                f'<text x="10" y="35" font-family="Arial" font-size="14" fill="red">'
                f'Error: {msg}</text></svg>'
            )
            sys.stdout.flush()

        if not h5ad_path or not cell_type_column or not gene_list or not svg_template:
            svg_error("Missing required params: file, col, gene, svg")
            return

        if not os.path.exists(h5ad_path):
            svg_error(f"H5AD file not found: {h5ad_path}")
            return
        if not os.path.exists(svg_template):
            svg_error(f"SVG template not found: {svg_template}")
            return

        try:
            expression_dict, gene_name = compute_expression(
                h5ad_path, cell_type_column, gene_list
            )
            tree = color_svg(svg_template, expression_dict,
                             gene_name=gene_name, opacity=opacity)

            ET.register_namespace('', 'http://www.w3.org/2000/svg')
            import io
            buf = io.BytesIO()
            tree.write(buf, encoding='utf-8', xml_declaration=True)
            sys.stdout.write(buf.getvalue().decode('utf-8'))
            sys.stdout.flush()

        except Exception as e:
            import traceback
            print(traceback.format_exc(), file=sys.stderr)
            svg_error(str(e))

    else:
        # ── Command-line mode ─────────────────────────────────────────────────
        import argparse
        parser = argparse.ArgumentParser(
            description="Compute expression from H5AD and color an SVG template."
        )
        parser.add_argument("h5ad_file",   help="Path to .h5ad file")
        parser.add_argument("svg_template", help="Path to input SVG template")
        parser.add_argument("col",          help="Cell type column in adata.obs")
        parser.add_argument("gene",         nargs='+', help="Gene ID(s) — first is used for coloring")
        parser.add_argument("output",       nargs='?', default="output_colored.svg",
                            help="Output SVG path (default: output_colored.svg)")
        parser.add_argument("--opacity",    type=float, default=1.0,
                            help="Fill opacity 0.0–1.0 (default 1.0)")
        args = parser.parse_args()

        if not os.path.exists(args.h5ad_file):
            print(f"Error: H5AD file not found: {args.h5ad_file}")
            sys.exit(1)
        if not os.path.exists(args.svg_template):
            print(f"Error: SVG template not found: {args.svg_template}")
            sys.exit(1)

        print("=" * 60)
        print("eFP SVG Generator")
        print("=" * 60)
        print(f"H5AD file:    {args.h5ad_file}")
        print(f"SVG template: {args.svg_template}")
        print(f"Column:       {args.col}")
        print(f"Gene(s):      {args.gene}")
        print(f"Output:       {args.output}")
        print(f"Opacity:      {args.opacity}")
        print("=" * 60)

        expression_dict, gene_name = compute_expression(
            args.h5ad_file, args.col, args.gene
        )
        tree = color_svg(args.svg_template, expression_dict,
                         gene_name=gene_name, opacity=args.opacity)

        ET.register_namespace('', 'http://www.w3.org/2000/svg')
        tree.write(args.output, encoding='utf-8', xml_declaration=True)
        print(f"\n✓ Written to: {args.output}")


if __name__ == "__main__":
    main()
