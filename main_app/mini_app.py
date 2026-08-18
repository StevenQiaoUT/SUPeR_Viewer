#!/mnt/home/sqiao/venv/bin/python3
"""
eFP Viewer — CGI Script
========================
Loads an H5AD file once, then produces a single HTML page containing:
  1. A colored eFP SVG (tissue expression browser)
  2. An interactive Plotly UMAP with cell-type dropdown

URL Parameters (minimal — dataset and column are chosen via in-page dropdowns):
    gene : Gene ID (e.g. AT3G05727)
    ds   : Dataset key — "arabidopsis_nat" or "rice" (default: arabidopsis_nat)
    col  : obs column key (default: label_majorXcondition)

All other parameters (file path, SVG template, umapcol, opacity)
are resolved server-side from the dataset key and are never exposed in the URL.

"""
import math
import sys
import os
import io
import json
import warnings
from xml.etree import ElementTree as ET
import re

warnings.filterwarnings("ignore", category=FutureWarning)

# ──────────────────────────────────────────────────────────────────────────────
# Dataset registry — all server-side config lives here, never in the URL
# ──────────────────────────────────────────────────────────────────────────────

DATASETS = {
    "arabidopsis_nat": {
        "label": "Arabidopsis Drought (Illouz-Eliaz et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/arabidopsis_nat.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_natanella.svg",
        "umap_col": "label_majorXcondition",
        "opacity": 1.0,
    },
    "rice": {
        "label": "Rice Drought-Salinity (Robertson et al., 2026)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/RiceOW.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/rice_template_ow.svg",
        "umap_col": "CellAnnotation",
        "opacity": 1.0,
    },
    "at_root_rs": {
        "label": "Arabidopsis Root (Shahan et al., 2022)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_root_shahan.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_root_template_shahan.svg",
        "umap_col": "Celltype",
        "opacity": 1.0,
    },
    "at_seed_martin": {
        "label": "Arabidopsis Seed (Martin et al., 2026)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_seed_martin.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_seed_template_martin.svg",
        "svg_l3": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_seed_l3_template_martin.svg",
        "umap_col": "level_2_annotation_timed",
        "opacity": 1.0,
    },
    "at_flower_lee": {
        "label": "Arabidopsis Flower (Lee et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_flower_lee.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_flower_template_lee.svg",
        "umap_col": "CellType",
        "opacity": 1.0,
    },
    "at_silique_lee": {
        "label": "Arabidopsis Silique (Lee et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_silique_lee.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_silique_template_lee.svg",
        "umap_col": "CellType",
        "opacity": 1.0,
    },
    "at_stem_lee": {
        "label": "Arabidopsis Stem (Lee et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_stem_lee.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_stem_template_lee.svg",
        "umap_col": "CellType",
        "opacity": 1.0,
    },
    "at_seed_0d_lee": {
        "label": "Arabidopsis Seed — 0 DAP (Lee et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_seed_0d_lee.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_seed_0d_lee.svg",
        "umap_col": "CellType",
        "opacity": 1.0,
    },
    "at_rosette_21d_lee": {
        "label": "Arabidopsis Rosette — 21 DAP (Lee et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_rosette_21d_lee.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_rosette_21d_lee.svg",
        "umap_col": "CellType",
        "opacity": 1.0,
    },
    "at_shoot_zhang": {
        "label": "Arabidopsis Shoot (Zhang et al., 2021)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_shoot_zhang.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_shoot_zhang.svg",
        "umap_col": "celltype_after",
        "opacity": 1.0,
    },
    "at_rosette_30d_lee": {
        "label": "Arabidopsis Rosette — 30 DAP (Lee et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_rosette_30d_lee.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_rosette_30d_lee.svg",
        "umap_col": "CellType",
        "opacity": 1.0,
    },
    "at_seedling_3d_lee": {
        "label": "Arabidopsis Seedling — 3 DAP (Lee et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_seedling_3d_lee.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_seedling_3d_lee.svg",
        "umap_col": "CellType",
        "opacity": 1.0,
    },
    "at_seedling_6d_lee": {
        "label": "Arabidopsis Seedling — 6 DAP (Lee et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_seedling_6d_lee.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_seedling_6d_lee.svg",
        "umap_col": "CellType",
        "opacity": 1.0,
    },
    "at_seedling_12d_lee": {
        "label": "Arabidopsis Seedling — 12 DAP (Lee et al., 2025)",
        "h5ad": "/mnt/home/sqiao/h5ad_files/at_seedling_12d_lee.h5ad",
        "svg": "/mnt/home/sqiao/public_html/cgi-bin/svg_templates/at_seedling_12d_lee.svg",
        "umap_col": "CellType",
        "opacity": 1.0,
    },
}

# Column registry — per-dataset map of obs column key → display label.
DATASET_COLUMNS = {
    "arabidopsis_nat": {
        "label_majorXcondition": "label_majorXcondition",
    },
    "rice": {
        "CAxCondition": "CAxCondition",
    },
    "at_root_rs": {
        "TypexGenotype": "TypexGenotype",
    },
    "at_seed_martin": {
        "level_2_annotation_timed": "TypexTimepoint",
        "level_3_annotation_full_timed": "Level 3",
    },
    "at_flower_lee": {
        "CellType": "cell type",
    },
    "at_silique_lee": {
        "CellType": "cell type",
    },
    "at_stem_lee": {
        "CellType": "cell type",
    },
    "at_seed_0d_lee": {
        "CellType": "cell type",
    },
    "at_rosette_21d_lee": {
        "CellType": "cell type",
    },
    "at_shoot_zhang": {
        "celltype_after": "cell type",
    },
    "at_rosette_30d_lee": {
        "CellType": "cell type",
    },
    "at_seedling_3d_lee": {
        "CellType": "cell type",
    },
    "at_seedling_6d_lee": {
        "CellType": "cell type",
    },
    "at_seedling_12d_lee": {
        "CellType": "cell type",
    },
}

COLUMNS = {k: k for cols in DATASET_COLUMNS.values() for k in cols}

# Default gene shown in the Gene ID box when each dataset is selected from
# the dropdown. Keyed by DATASETS key (i.e. the same key used by ds-select).
DATASET_DEFAULT_GENES = {
    "arabidopsis_nat": "SZF1",  # Illouz-Eliaz — SZF1/AT3G55980
    "rice": "Os10g0168500",  # Rice (Wilkins/Roberston)
    "at_root_rs": "AT1G79580",  # Shahan
    "at_seed_martin": "AT2G42840",  # Martin
    "at_flower_lee": "AT3G15510",  # Lee — flower
    "at_silique_lee": "AT3G24140",  # Lee — silique
    "at_stem_lee": "AT5G26000",  # Lee — stem
    "at_seed_0d_lee": "AT1G01010",  # Lee — seed 0 DAP — TODO: replace with a real marker gene
    "at_rosette_21d_lee": "AT1G01010",  # Lee — rosette 21 DAP — TODO: replace with a real marker gene
    "at_shoot_zhang": "AT1G65480",  # Zhang — shoot
    "at_rosette_30d_lee": "AT1G01010",  # Lee — rosette 30 DAP — TODO: replace with a real marker gene
    "at_seedling_3d_lee": "AT1G01010",  # Lee — seedling 3 DAP — TODO: replace with a real marker gene
    "at_seedling_6d_lee": "AT1G01010",  # Lee — seedling 6 DAP — TODO: replace with a real marker gene
    "at_seedling_12d_lee": "AT1G01010",  # Lee — seedling 12 DAP — TODO: replace with a real marker gene
}

# Curated marker genes per dataset — shown as autocomplete suggestions
# when the user types in the Gene ID field.
DATASET_MARKER_GENES = {
    "arabidopsis_nat": [
        {"gene": "AT1G17840", "label": "AT1G17840 — epidermal"},
        {"gene": "AT1G22430", "label": "AT1G22430 — mesophylls"},
        {"gene": "AT1G74730", "label": "AT1G74730 — W0 condition"},
        {"gene": "AT1G35780", "label": "AT1G35780 — phloem cell"},
        {"gene": "AT3G04810", "label": "AT3G04810 — guard cell"},
    ],
    "rice": [
        {"gene": "Os01g0291500", "label": "Os01g0291500 — fibre"},
        {"gene": "Os01g0610400", "label": "Os01g0610400 — mestome sheath"},
        {"gene": "Os01g0711400", "label": "Os01g0711400 — mesophyll"},
        {"gene": "Os01g0878900", "label": "Os01g0878900 — epidermis"},
        {"gene": "Os02g0627100", "label": "Os02g0627100 — xylem parenchyma"},
    ],
    "at_root_rs": [
        {"gene": "AT1G20823", "label": "AT1G20823 — phloem"},
        {"gene": "AT1G20010", "label": "AT1G20010 — cortex"},
        {"gene": "AT1G26810", "label": "AT1G26810 — metaxylem"},
        {"gene": "AT1G77690", "label": "AT1G77690 — procambium"},
    ],
    "at_seed_martin": [
        {"gene": "AT1G03860", "label": "AT1G03860 — endosperm"},
        {"gene": "AT1G15930", "label": "AT1G15930 — embryo"},
        {"gene": "AT1G12080", "label": "AT1G12080 — 7DAP MCE and PEN"},
        {"gene": "AT1G14350", "label": "AT1G14350 — seed coat"},
        {"gene": "AT1G15500", "label": "AT1G15500 — 5DAP embryo"},
    ],
    "at_flower_lee": [
        {"gene": "AT1G04270", "label": "AT1G04270 — epidermal"},
        {"gene": "AT1G17020", "label": "AT1G17020 — anther and pollen"},
        {"gene": "AT1G19950", "label": "AT1G19950 — male meiocyte"},
        {"gene": "AT1G15740", "label": "AT1G15740 — ovule"},
    ],
    "at_silique_lee": [
        {"gene": "AT1G51035", "label": "AT1G51035 — mature silique"},
        {"gene": "AT4G38420", "label": "AT4G38420 — young silique"},
        {"gene": "AT2G14690", "label": "AT2G14690 — seed"},
    ],
    "at_stem_lee": [
        {"gene": "AT5G41315", "label": "AT5G41315 — epidermis"},
        {"gene": "AT2G37260", "label": "AT2G37260 — guard"},
        {"gene": "AT1G07640", "label": "AT1G07640 — trichome"},
    ],
    "at_seed_0d_lee": [
        # TODO: replace with real marker genes for this dataset's cell types
        # (Epidermal, Guard, Meristematic, Seed_coat, Stele)
    ],
    "at_rosette_21d_lee": [
        # TODO: replace with real marker genes for this dataset's cell types
        # (Epidermal, Guard, Stele are drawn in the SVG; Meristematic and
        # Unannotated exist in obs but have no corresponding SVG shape)
    ],
    "at_shoot_zhang": [
        {"gene": "AT1G11840", "label": "AT1G11840 — shoot apical meristem"},
        {"gene": "AT1G01090", "label": "AT1G01090 — mesophyll"},
        {"gene": "AT5G53210", "label": "AT5G53210 — meristemoid and pavement"},
        {"gene": "AT2G25810", "label": "AT2G25810 — guard cells"},
        {"gene": "AT1G12880", "label": "AT1G12880 — companion cells"},
    ],
    "at_rosette_30d_lee": [
        # TODO: replace with real marker genes for this dataset's cell types
        # (Epidermal, Guard, Meristematic, Mesophyll, Stele; Unannotated
        # exists in obs but has no corresponding SVG shape)
    ],
    "at_seedling_3d_lee": [
        # TODO: replace with real marker genes for this dataset's cell types
        # (Vascular, Mesophyll, Epidermis, Phloem, Dividing, Trichoblast,
        # Phloem parenchyma; Unknown exists in obs but has no corresponding
        # SVG shape)
    ],
    "at_seedling_6d_lee": [
        # TODO: replace with real marker genes for this dataset's cell types
        # (Epidermal, Stele, Mesophyll, Meristematic, Guard all drawn in the SVG)
    ],
    "at_seedling_12d_lee": [
        # TODO: replace with real marker genes for this dataset's cell types
        # (Epidermal, Mesophyll, Stele, Guard drawn in the SVG; Unannotated
        # exists in obs but has no corresponding SVG shape)
    ],
}

DEFAULT_DATASET = "arabidopsis_nat"
DEFAULT_GENE = "SZF1"
DEFAULT_COLUMN = "label_majorXcondition"

_SEED_SKIP_IDS = {
    "3DAP unlabeled", "5DAP unlabeled", "7DAP unlabeled",
    "3DAP L1", "3DAP L2L", "3DAP L2R",
    "5DAP L1", "5DAP L2L", "5DAP L2R",
    "7DAP L1", "7DAP L2L", "7DAP L2R",
    "3DAP-unlabeled", "5DAP-unlabeled", "7DAP-unlabeled",
    "3DAP", "5DAP", "7DAP",
    "top_bar", "layer1",
}

_VASCULAR_CIRCLE_ID = {
    "D0": "D0_vascular_circle",
    "R15": "R15_vascular_circle",
    "W0": "W0_vascular_circle",
    "W15": "W15_vascular_circle",
}

_ROOT_GENOTYPE_PREFIXES = ('shr2_', 'scr4_', 'col0_')
_ROOT_POSITION_PREFIXES = ('top_', 'side_', 'tip_')
_ROOT_SUFFIX_LABELS = {
    'endodermis': 'Root endodermis',
    'columella': 'Collumella root cap',
    'lateral_root_cap': 'Lateral root cap',
    'phloem': 'Phloem',
    'procambium': 'Root procambium',
    'cortex': 'Root cortex',
    'xpp_circle': 'Xylem pole pericycle',
    'proto_circle': 'Protoxylem',
    'ppp_circle': 'Phloem pole pericycle',
    'meta_circle': 'Metaxylem',
    'xylem': 'Xylem (combined)',
}

# Aliases for the "well-watered" / control condition across rice obs
# spellings. The SVG panel header text is always literally "Well Watered"
# (group id "Well.Watered"), but the underlying obs `Condition` value can be
# spelled differently depending on how the dataset was exported (e.g. "WW",
# "Control", "CK", "Mock"). Direct normalised-string matching against the
# SVG label fails whenever the obs value doesn't textually resemble "well
# watered", which is exactly what was happening here — every other panel
# (Mild Drought, Moderate Drought, Mild/Moderate Salinity) coincidentally
# matched the obs spelling, but the well-watered control did not, so its
# `cond_val` resolved to None and the whole panel was skipped (left
# uncoloured/outline-only).
_RICE_CONTROL_ALIASES = {
    'well watered', 'wellwatered', 'control', 'ck', 'ww', 'mock', 'untreated',
}


def _resolve_rice_control(norm_to_cond):
    """
    Resolve the obs Condition value that represents the well-watered
    control. Tries the literal 'well watered' spelling first (handles the
    common case directly), then falls back to scanning all obs condition
    values for one whose normalised form matches a known control alias.
    Returns None if nothing matches.
    """
    direct = norm_to_cond.get('well watered')
    if direct is not None:
        return direct
    for norm_val, actual in norm_to_cond.items():
        if norm_val in _RICE_CONTROL_ALIASES:
            return actual
    return None


def root_celltype_label(svg_id):
    """
    Strip the genotype and position prefixes off a Shahan-root SVG id to
    recover a human-readable cell-type label usable for highlighting
    (e.g. 'shr2_top_endodermis' -> 'Root endodermis').
    """
    rest = svg_id
    for pre in _ROOT_GENOTYPE_PREFIXES:
        if rest.startswith(pre):
            rest = rest[len(pre):]
            break
    for pre in _ROOT_POSITION_PREFIXES:
        if rest.startswith(pre):
            rest = rest[len(pre):]
            break
    return _ROOT_SUFFIX_LABELS.get(rest, rest.replace('_', ' ').title())


def seed_celltype_label(svg_id):
    """Strip the DAP timepoint prefix off a Martin-seed SVG id."""
    for tp in ('3DAP_', '5DAP_', '7DAP_'):
        if svg_id.startswith(tp):
            return svg_id[len(tp):]
    return svg_id


# ──────────────────────────────────────────────────────────────────────────────
# Color utilities
# ──────────────────────────────────────────────────────────────────────────────

def expression_to_color(value, max_val, scheme='yellow_to_red', ratio=None):
    if scheme == 'yellow_to_red':
        minColor = {'red': 255, 'green': 255, 'blue': 0}
        maxColor = {'red': 255, 'green': 0, 'blue': 0}
    elif scheme == 'yellow_to_blue':
        minColor = {'red': 255, 'green': 255, 'blue': 0}
        maxColor = {'red': 0, 'green': 0, 'blue': 255}
    else:  # blue_to_yellow (kept for back-compat)
        minColor = {'red': 0, 'green': 0, 'blue': 255}
        maxColor = {'red': 255, 'green': 255, 'blue': 0}

    if ratio is None:
        ratio = min(value / max_val, 1.0) if max_val > 0 else 0

    if not (0 <= ratio <= 1):
        return "#ffffff"

    red = minColor['red'] + round((maxColor['red'] - minColor['red']) * ratio)
    green = minColor['green'] + round((maxColor['green'] - minColor['green']) * ratio)
    blue = minColor['blue'] + round((maxColor['blue'] - minColor['blue']) * ratio)
    return f"#{red:02x}{green:02x}{blue:02x}"


def expression_to_color_relative(value, control_val, log2_max):
    """
    Diverging color scale centered on the control condition.
      log2 ratio > 0  →  yellow → red   (upregulated vs control)
      log2 ratio = 0  →  yellow          (same as control)
      log2 ratio < 0  →  yellow → blue  (downregulated vs control)
    norm is clamped to [0, 1] using log2_max as the symmetric ceiling.

    Special cases:
      value is None or control_val is None → white (no data)
      control_val == 0 and value == 0      → yellow (neither expressed)
      control_val == 0 and value > 0       → full red (expressed only in sample)
      value == 0 and control_val > 0       → full blue (expressed only in control)
    """
    if value is None or control_val is None:
        return "#ffffff"
    if control_val == 0 and value == 0:
        # Neither condition has expression — show yellow (no change / none)
        return expression_to_color(0.0, 1.0, scheme='yellow_to_red', ratio=0.0)
    if control_val == 0 and value > 0:
        # Gene expressed in sample but absent in control → maximum upregulation
        return expression_to_color(1.0, 1.0, scheme='yellow_to_red', ratio=1.0)
    if value == 0:
        # Gene expressed in control but absent in sample → maximum downregulation
        return expression_to_color(1.0, 1.0, scheme='yellow_to_blue', ratio=1.0)
    try:
        log2_ratio = value - control_val   # already log-normalized: difference = fold change
    except (ValueError, ZeroDivisionError):
        return "#ffffff"

    clamped = max(-log2_max, min(log2_max, log2_ratio))
    norm = abs(clamped) / log2_max  # 0 = no change, 1 = max deviation

    if clamped >= 0:
        return expression_to_color(norm, 1.0, scheme='yellow_to_red', ratio=norm)
    else:
        return expression_to_color(norm, 1.0, scheme='yellow_to_blue', ratio=norm)


# ──────────────────────────────────────────────────────────────────────────────
# SVG legend + title
# ──────────────────────────────────────────────────────────────────────────────

def add_legend(root, gene_name, colour_ceiling, n_boxes=10):
    SVG_NS = "http://www.w3.org/2000/svg"
    box_w = 24
    box_h = 20
    tick_gap = 6
    font_size = 10
    header_fs = 11
    left_pad = 10
    tick_label_w = 40
    right_pad = 20
    legend_width = left_pad + box_w + tick_gap + tick_label_w + right_pad

    for attr in ["width", "viewBox"]:
        val = root.get(attr, "")
        if not val:
            continue
        try:
            if attr == "width":
                num = float("".join(c for c in val if c.isdigit() or c == "."))
                unit = "".join(c for c in val if c.isalpha())
                root.set(attr, f"{num + legend_width}{unit}")
            else:
                parts = val.split()
                parts[2] = str(float(parts[2]) + legend_width)
                root.set(attr, " ".join(parts))
        except (ValueError, IndexError):
            pass

    existing = [c for c in list(root) if c.get("id") != "expression_legend"]
    for child in existing:
        root.remove(child)

    shift_g = ET.Element(f"{{{SVG_NS}}}g")
    shift_g.set("transform", f"translate({legend_width}, 0)")
    for child in existing:
        shift_g.append(child)
    root.append(shift_g)

    x0 = left_pad
    try:
        h_attr = root.get("height", "")
        canvas_h = float("".join(c for c in h_attr if c.isdigit() or c == "."))
    except ValueError:
        canvas_h = 600

    header_h = 3 * (header_fs + 4)
    legend_h = header_h + n_boxes * box_h + font_size + 6
    y0 = canvas_h - legend_h - 10

    g = ET.Element(f"{{{SVG_NS}}}g")
    g.set("id", "expression_legend")

    def txt(x, y, text, fs=font_size, anchor="start", weight="normal", style=""):
        el = ET.Element(f"{{{SVG_NS}}}text")
        el.set("x", str(x))
        el.set("y", str(y))
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
    g.append(txt(x0, hy, "Single Cell Max"))
    hy += font_size + 4
    g.append(txt(x0, hy, "Linear"))

    boxes_y0 = y0 + header_h

    for i in range(n_boxes):
        ratio_mid = 1.0 - (i + 0.5) / n_boxes
        color = expression_to_color(ratio_mid * colour_ceiling, colour_ceiling, "yellow_to_red")
        by = boxes_y0 + i * box_h
        rect = ET.Element(f"{{{SVG_NS}}}rect")
        rect.set("x", str(x0))
        rect.set("y", str(by))
        rect.set("width", str(box_w))
        rect.set("height", str(box_h))
        rect.set("fill", color)
        rect.set("stroke", "#888888")
        rect.set("stroke-width", "0.5")
        g.append(rect)

        ratio_top = 1.0 - i / n_boxes
        tick_val = ratio_top * colour_ceiling
        g.append(txt(x0 + box_w + tick_gap, by + 4, f"{tick_val:.2f}", fs=font_size - 1))

    g.append(txt(x0 + box_w + tick_gap, boxes_y0 + n_boxes * box_h + 4,
                 f"0.00", fs=font_size - 1))
    root.append(g)


def add_relative_legend(root, gene_name, log2_max, n_boxes=10):
    """Diverging legend for relative (log2 ratio) mode: red=up, yellow=no change, blue=down."""
    SVG_NS = "http://www.w3.org/2000/svg"
    box_w = 24
    box_h = 20
    tick_gap = 6
    font_size = 10
    header_fs = 11
    left_pad = 10
    tick_label_w = 50
    right_pad = 20
    legend_width = left_pad + box_w + tick_gap + tick_label_w + right_pad

    for attr in ["width", "viewBox"]:
        val = root.get(attr, "")
        if not val:
            continue
        try:
            if attr == "width":
                num = float("".join(c for c in val if c.isdigit() or c == "."))
                unit = "".join(c for c in val if c.isalpha())
                root.set(attr, f"{num + legend_width}{unit}")
            else:
                parts = val.split()
                parts[2] = str(float(parts[2]) + legend_width)
                root.set(attr, " ".join(parts))
        except (ValueError, IndexError):
            pass

    existing = [c for c in list(root) if c.get("id") != "expression_legend"]
    for child in existing:
        root.remove(child)

    shift_g = ET.Element(f"{{{SVG_NS}}}g")
    shift_g.set("transform", f"translate({legend_width}, 0)")
    for child in existing:
        shift_g.append(child)
    root.append(shift_g)

    x0 = left_pad
    try:
        canvas_h = float("".join(c for c in root.get("height", "") if c.isdigit() or c == "."))
    except ValueError:
        canvas_h = 600

    # Auto-scale n_boxes so the legend uses at most 60% of the canvas height.
    # Each box is box_h px; header is ~3*(header_fs+4) px; 2*n_boxes total boxes.
    if n_boxes is None:
        header_h_est = 3 * (header_fs + 4)
        max_legend_h = canvas_h * 0.6
        n_boxes = max(3, min(10, int((max_legend_h - header_h_est) // (box_h * 2))))

    total_boxes = n_boxes * 2  # positive half + negative half
    header_h = 3 * (header_fs + 4)
    legend_h = header_h + total_boxes * box_h + font_size + 6
    y0 = canvas_h - legend_h - 10

    g = ET.Element(f"{{{SVG_NS}}}g")
    g.set("id", "expression_legend")

    def txt(x, y, text, fs=font_size, anchor="start", weight="normal", style=""):
        el = ET.Element(f"{{{SVG_NS}}}text")
        el.set("x", str(x));
        el.set("y", str(y))
        el.set("font-family", "Arial, sans-serif");
        el.set("font-size", str(fs))
        el.set("font-weight", weight);
        el.set("text-anchor", anchor)
        el.set("dominant-baseline", "auto");
        el.set("fill", "#222222")
        if style:
            el.set("font-style", style)
        el.text = text
        return el

    hy = y0 + header_fs
    g.append(txt(x0, hy, gene_name, fs=header_fs, weight="bold", style="italic"))
    hy += header_fs + 4
    g.append(txt(x0, hy, "vs WW Control" if True else "vs W0 Control"))
    hy += font_size + 4
    g.append(txt(x0, hy, "Log\u2082 Ratio"))

    boxes_y0 = y0 + header_h

    # Draw top-to-bottom: i=0 is most positive (red), i=total_boxes-1 is most negative (blue)
    for i in range(total_boxes):
        if i < n_boxes:
            norm = 1.0 - i / n_boxes  # 1.0 → 0.0  (red → yellow)
            color = expression_to_color(norm, 1.0, scheme='yellow_to_red', ratio=norm)
        else:
            norm = (i - n_boxes) / n_boxes  # 0.0 → 1.0  (yellow → blue)
            color = expression_to_color(norm, 1.0, scheme='yellow_to_blue', ratio=norm)

        by = boxes_y0 + i * box_h
        rect = ET.Element(f"{{{SVG_NS}}}rect")
        rect.set("x", str(x0));
        rect.set("y", str(by))
        rect.set("width", str(box_w));
        rect.set("height", str(box_h))
        rect.set("fill", color);
        rect.set("stroke", "#888888");
        rect.set("stroke-width", "0.5")
        g.append(rect)

        if i == 0:
            g.append(txt(x0 + box_w + tick_gap, by + 4, f"+{log2_max:.2f}", fs=font_size - 1))
        elif i == n_boxes:
            g.append(txt(x0 + box_w + tick_gap, by + 4, "0.00", fs=font_size - 1))
        elif i == total_boxes - 1:
            g.append(txt(x0 + box_w + tick_gap, by + box_h, f"-{log2_max:.2f}", fs=font_size - 1))

    root.append(g)


def add_relative_legend_labeled(root, gene_name, log2_max, control_label, n_boxes=None):
    """
    Wrapper around add_relative_legend that patches the legend subtitle
    to show the correct control label (e.g. 'vs W0 Control', 'vs WW Control',
    'vs Col-0 Control') before delegating to the shared drawing routine.
    We rebuild the legend inline here so the label is correct.
    """
    SVG_NS = "http://www.w3.org/2000/svg"
    box_w = 24
    box_h = 20
    tick_gap = 6
    font_size = 10
    header_fs = 11
    left_pad = 10
    tick_label_w = 50
    right_pad = 20
    legend_width = left_pad + box_w + tick_gap + tick_label_w + right_pad

    for attr in ["width", "viewBox"]:
        val = root.get(attr, "")
        if not val:
            continue
        try:
            if attr == "width":
                num = float("".join(c for c in val if c.isdigit() or c == "."))
                unit = "".join(c for c in val if c.isalpha())
                root.set(attr, f"{num + legend_width}{unit}")
            else:
                parts = val.split()
                parts[2] = str(float(parts[2]) + legend_width)
                root.set(attr, " ".join(parts))
        except (ValueError, IndexError):
            pass

    existing = [c for c in list(root) if c.get("id") != "expression_legend"]
    for child in existing:
        root.remove(child)

    shift_g = ET.Element(f"{{{SVG_NS}}}g")
    shift_g.set("transform", f"translate({legend_width}, 0)")
    for child in existing:
        shift_g.append(child)
    root.append(shift_g)

    x0 = left_pad
    try:
        canvas_h = float("".join(c for c in root.get("height", "") if c.isdigit() or c == "."))
    except ValueError:
        canvas_h = 600

    # Auto-scale n_boxes so the legend uses at most 60% of the canvas height.
    # Each box is box_h px; header is ~3*(header_fs+4) px; 2*n_boxes total boxes.
    if n_boxes is None:
        header_h_est = 3 * (header_fs + 4)
        max_legend_h = canvas_h * 0.6
        n_boxes = max(3, min(10, int((max_legend_h - header_h_est) // (box_h * 2))))

    total_boxes = n_boxes * 2
    header_h = 3 * (header_fs + 4)
    legend_h = header_h + total_boxes * box_h + font_size + 6
    y0 = canvas_h - legend_h - 10

    g = ET.Element(f"{{{SVG_NS}}}g")
    g.set("id", "expression_legend")

    def txt(x, y, text, fs=font_size, anchor="start", weight="normal", style=""):
        el = ET.Element(f"{{{SVG_NS}}}text")
        el.set("x", str(x));
        el.set("y", str(y))
        el.set("font-family", "Arial, sans-serif");
        el.set("font-size", str(fs))
        el.set("font-weight", weight);
        el.set("text-anchor", anchor)
        el.set("dominant-baseline", "auto");
        el.set("fill", "#222222")
        if style:
            el.set("font-style", style)
        el.text = text
        return el

    hy = y0 + header_fs
    g.append(txt(x0, hy, gene_name, fs=header_fs, weight="bold", style="italic"))
    hy += header_fs + 4
    g.append(txt(x0, hy, f"vs {control_label} Control"))
    hy += font_size + 4
    g.append(txt(x0, hy, "Log\u2082 Ratio"))

    boxes_y0 = y0 + header_h

    for i in range(total_boxes):
        if i < n_boxes:
            norm = 1.0 - i / n_boxes
            color = expression_to_color(norm, 1.0, scheme='yellow_to_red', ratio=norm)
        else:
            norm = (i - n_boxes) / n_boxes
            color = expression_to_color(norm, 1.0, scheme='yellow_to_blue', ratio=norm)

        by = boxes_y0 + i * box_h
        rect = ET.Element(f"{{{SVG_NS}}}rect")
        rect.set("x", str(x0));
        rect.set("y", str(by))
        rect.set("width", str(box_w));
        rect.set("height", str(box_h))
        rect.set("fill", color);
        rect.set("stroke", "#888888");
        rect.set("stroke-width", "0.5")
        g.append(rect)

        if i == 0:
            g.append(txt(x0 + box_w + tick_gap, by + 4, f"+{log2_max:.2f}", fs=font_size - 1))
        elif i == n_boxes:
            g.append(txt(x0 + box_w + tick_gap, by + 4, "0.00", fs=font_size - 1))
        elif i == total_boxes - 1:
            g.append(txt(x0 + box_w + tick_gap, by + box_h, f"-{log2_max:.2f}", fs=font_size - 1))

    root.append(g)


def add_gene_title(root, gene_name, font_size=18, padding_bottom=10):
    SVG_NS = 'http://www.w3.org/2000/svg'
    offset_y = font_size + padding_bottom

    existing = root.find(f'{{{SVG_NS}}}title')
    if existing is not None:
        root.remove(existing)

    t = ET.Element(f'{{{SVG_NS}}}title')
    t.text = gene_name
    root.insert(0, t)

    for attr in ['height', 'viewBox']:
        val = root.get(attr, '')
        if not val:
            continue
        try:
            if attr == 'height':
                num = float(''.join(c for c in val if c.isdigit() or c == '.')) + offset_y
                unit = ''.join(c for c in val if c.isalpha())
                root.set(attr, f"{num}{unit}")
            else:
                parts = val.split()
                parts[3] = str(float(parts[3]) + offset_y)
                root.set(attr, ' '.join(parts))
        except (ValueError, IndexError):
            pass

    children = [c for c in list(root) if c is not t]
    for c in children:
        root.remove(c)

    wrap = ET.Element(f'{{{SVG_NS}}}g')
    wrap.set('transform', f'translate(0, {offset_y})')
    for c in children:
        wrap.append(c)
    root.append(wrap)

    try:
        cx = float(''.join(c for c in root.get('width', '') if c.isdigit() or c == '.')) / 2
    except ValueError:
        cx = 300

    tx = ET.Element(f'{{{SVG_NS}}}text')
    tx.set('x', str(cx))
    tx.set('y', str(font_size))
    tx.set('text-anchor', 'middle')
    tx.set('font-family', 'Arial, sans-serif')
    tx.set('font-size', str(font_size))
    tx.set('font-weight', 'bold')
    tx.set('fill', '#222222')
    tx.text = gene_name
    root.insert(1, tx)


# ──────────────────────────────────────────────────────────────────────────────
# SVG coloring
# ──────────────────────────────────────────────────────────────────────────────

def add_median_circle(root, gene_name, median_ctrl, radius=40):
    """
    Appends a reference circle to the SVG for the median-as-control relative view.
    The circle is always yellow — it IS the control (log2 ratio = 0).
    Expands the canvas downward to make room so the circle never overlaps content.
    Placed centred horizontally in the new bottom strip.
    """
    SVG_NS = "http://www.w3.org/2000/svg"
    font_size = 11
    n_label_lines = 2
    pad = 16
    strip_h = radius * 2 + pad * 2 + (font_size + 4) * n_label_lines

    try:
        w = float(''.join(c for c in root.get('width', '') if c.isdigit() or c == '.'))
    except ValueError:
        w = 600
    try:
        h = float(''.join(c for c in root.get('height', '') if c.isdigit() or c == '.'))
    except ValueError:
        h = 600

    # Expand canvas height to make room below existing content
    new_h = h + strip_h
    for attr in ['height', 'viewBox']:
        val = root.get(attr, '')
        if not val:
            continue
        try:
            if attr == 'height':
                unit = ''.join(c for c in val if c.isalpha())
                root.set(attr, f"{new_h}{unit}")
            else:
                parts = val.split()
                parts[3] = str(float(parts[3]) + strip_h)
                root.set(attr, ' '.join(parts))
        except (ValueError, IndexError):
            pass

    # Centre the circle horizontally in the full (legend-expanded) canvas width
    cx = w / 2
    cy = h + pad + radius

    color = expression_to_color(0.0, 1.0, 'yellow_to_red', ratio=0.0)

    g = ET.Element(f'{{{SVG_NS}}}g')
    g.set('id', 'median_summary_circle')

    circle = ET.Element(f'{{{SVG_NS}}}circle')
    circle.set('cx', str(cx))
    circle.set('cy', str(cy))
    circle.set('r', str(radius))
    circle.set('fill', color)
    circle.set('stroke', '#555555')
    circle.set('stroke-width', '1.5')
    circle.set('data-tooltip',
               f"Median control\nMedian of cell-type means: {median_ctrl:.4f}\n"
               f"(log\u2082 ratio = 0 reference)")
    circle.set('class', 'hoverable')
    g.append(circle)

    def lbl(text, dy_offset):
        t = ET.Element(f'{{{SVG_NS}}}text')
        t.set('x', str(cx))
        t.set('y', str(cy + radius + dy_offset))
        t.set('text-anchor', 'middle')
        t.set('font-family', 'Arial, sans-serif')
        t.set('font-size', str(font_size))
        t.set('fill', '#222222')
        t.text = text
        return t

    g.append(lbl("Median (control)", font_size + 4))
    g.append(lbl(f"{median_ctrl:.4f}", (font_size + 4) * 2))

    root.append(g)


def color_svg(svg_file, expression_dict, ds_key, gene_name=None, opacity=1.0, control=None):
    tree = ET.parse(svg_file)
    root = tree.getroot()
    SVG_NS = 'http://www.w3.org/2000/svg'
    conditions = ['D0', 'W0', 'R15', 'W15']

    # Tracks every distinct "cell type" label encountered while colouring,
    # used to populate the highlight dropdown for this dataset/view.
    cell_types_seen = set()

    if ds_key == 'arabidopsis_nat':
        vis_keys = (
                [f"{c}_Guard" for c in conditions] +
                [f"{c}_Mesophyll" for c in conditions] +
                [f"{c}_Epidermal" for c in conditions] +
                [f"{c}_Trichome" for c in conditions] +
                [f"{c}_Vascular" for c in conditions] +
                [f"{c}_Phloem Parenchyma" for c in conditions] +
                [f"{c}_Phloem companion" for c in conditions]
        )
        vis_vals = [expression_dict[k] for k in vis_keys if k in expression_dict]
    elif ds_key == 'at_seed_martin':
        svg_group_ids = {
            el.get('id')
            for el in root.iter()
            if el.get('id') and el.get('id') not in _SEED_SKIP_IDS
               and el.tag == f'{{{SVG_NS}}}g'
        }
        vis_vals = [expression_dict[k] for k in svg_group_ids if k in expression_dict]
    elif ds_key in ('at_flower_lee', 'at_silique_lee', 'at_stem_lee', 'at_seed_0d_lee', 'at_rosette_21d_lee', 'at_shoot_zhang', 'at_rosette_30d_lee', 'at_seedling_3d_lee', 'at_seedling_6d_lee', 'at_seedling_12d_lee'):
        # Build vis_vals using exactly the same element-matching logic as the
        # Lee coloring loop: only elements whose id is in expression_dict,
        # restricted to shape tags (<path/circle/rect/ellipse/polygon>) or
        # <g> groups. This prevents structural SVG IDs (layer groups, defs,
        # etc.) from accidentally matching expression_dict keys and skewing
        # the median. (at_rosette_21d_lee draws its cell-type shapes as
        # <polygon>, unlike the path-based templates used by the other Lee
        # datasets, hence 'polygon' being included here.)
        _lee_shape_tags = {
            f'{{{SVG_NS}}}path', f'{{{SVG_NS}}}circle',
            f'{{{SVG_NS}}}rect', f'{{{SVG_NS}}}ellipse',
            f'{{{SVG_NS}}}polygon', f'{{{SVG_NS}}}g',
        }
        _lee_vis_keys = {
            el.get('id') for el in root.iter()
            if el.get('id') and el.tag in _lee_shape_tags
               and el.get('id') in expression_dict
        }
        vis_vals = [expression_dict[k] for k in _lee_vis_keys]
        if not vis_vals:
            vis_vals = list(expression_dict.values())
        print(f"  Lee vis_vals: {len(vis_vals)} SVG-matched cell types: {sorted(_lee_vis_keys)}", file=sys.stderr)
        # For median_ctrl: exclude zero-expression cell types so the reference
        # is centred on the expressing population only, not dragged down by
        # cell types where the gene is absent.
        _lee_expressing_vals = [v for v in vis_vals if v > 0]
        if _lee_expressing_vals:
            vis_vals = _lee_expressing_vals
            print(f"  Lee median pool: {len(vis_vals)} expressing cell types (zeros excluded)", file=sys.stderr)
    else:
        svg_ids = {el.get('id') for el in root.iter() if el.get('id')}
        vis_vals = [v for k, v in expression_dict.items() if k in svg_ids]
        if not vis_vals:
            vis_vals = list(expression_dict.values())

    colour_ceiling = max(vis_vals) if vis_vals else 1.0

    # ── Compute log2_max from expression_dict directly ───────────────────────
    log2_max = 1.0  # fallback

    def _compute_log2_max(pairs):
        """pairs: list of (sample_val, control_val). Returns max abs fold-change.
        Data is already log-normalized, so fold-change = difference, not log2(ratio)."""
        diffs = []
        for v, ctrl in pairs:
            if v is None or ctrl is None:
                continue
            diffs.append(abs(v - ctrl))
        return max(diffs) if diffs else 1.0

    if control == "WW" and ds_key == "rice":
        # expression_dict keys are "{cond}_{ct}" built from obs columns.
        # The obs Condition value for "well watered" may be "Well-Watered",
        # "WW", "Control", "CK", "Mock", etc. — find it via alias-aware
        # matching rather than assuming it textually resembles the SVG
        # panel label "Well.Watered".
        def _norm_rice(s):
            return ' '.join(s.lower().replace('.', ' ').replace('-', ' ').split())

        _rice_conds = set()
        for k in expression_dict:
            parts = k.split('_', 1)
            if len(parts) == 2:
                _rice_conds.add(parts[0])

        _rice_norm_to_cond = {_norm_rice(c): c for c in _rice_conds}
        _ww_cond_val = _resolve_rice_control(_rice_norm_to_cond)
        print(f"  Rice WW cond value resolved to: {_ww_cond_val!r}", file=sys.stderr)
        print(f"  Rice all cond values: {sorted(_rice_conds)}", file=sys.stderr)

        pairs = []
        for key, v in expression_dict.items():
            parts = key.split('_', 1)
            if len(parts) != 2:
                continue
            cond, ct = parts
            if cond != _ww_cond_val:
                ctrl = expression_dict.get(f'{_ww_cond_val}_{ct}') if _ww_cond_val else None
                if ctrl is not None:
                    pairs.append((v, ctrl))
        log2_max = _compute_log2_max(pairs)
        print(f"  Rice relative log2_max={log2_max:.4f} from {len(pairs)} pairs", file=sys.stderr)

    elif control == "Col-0" and ds_key == "at_root_rs":
        pairs = []
        for svg_id, v in expression_dict.items():
            for mut_prefix in ("shr2", "scr4"):
                if svg_id.startswith(mut_prefix + "_"):
                    suffix = svg_id[len(mut_prefix):]
                    ctrl = expression_dict.get(f"col0{suffix}")
                    pairs.append((v, ctrl))
        log2_max = _compute_log2_max(pairs)

    elif control == "W0" and ds_key == "arabidopsis_nat":
        _NAT_CELL_TYPES = [
            "Guard", "Mesophyll", "Epidermal", "Trichome",
            "Vascular", "Phloem Parenchyma", "Phloem companion",
        ]
        pairs = []
        for ct in _NAT_CELL_TYPES:
            ctrl = expression_dict.get(f"W0_{ct}")
            for c in conditions:
                v = expression_dict.get(f"{c}_{ct}")
                pairs.append((v, ctrl))
        log2_max = _compute_log2_max(pairs)

    # ── Median-as-control (seed / Lee datasets) ───────────────────────────────
    # median_ctrl = median of vis_vals — the same SVG-visible per-cell-type means
    # used for colour_ceiling. This excludes non-SVG keys (e.g. L1 annotation
    # entries mixed into the seed expression_dict) so the reference is exactly
    # the population being displayed.
    median_ctrl = None
    if control == "median":
        import statistics as _stats
        _vis = [v for v in vis_vals if v is not None and v >= 0]
        median_ctrl = float(_stats.median(_vis)) if _vis else 0.0
        pairs = [(v, median_ctrl) for v in _vis]
        log2_max = _compute_log2_max(pairs)
        print(
            f"  Median control (of {len(_vis)} SVG-visible cell-type means): {median_ctrl:.4f}  log2_max={log2_max:.4f}",
            file=sys.stderr)

    def get_color(key):
        v = expression_dict.get(key)
        return expression_to_color(v, colour_ceiling, 'yellow_to_red') if v is not None else None

    def raw_val(key):
        v = expression_dict.get(key)
        return f"{v:.4f}" if v is not None else "N/A"

    def add_tooltip(elem, text):
        for child in list(elem):
            if child.tag == f'{{{SVG_NS}}}title':
                elem.remove(child)
        elem.set('data-tooltip', text)
        elem.set('class', (elem.get('class', '') + ' hoverable').strip())

    def tag_celltype(elem, ct_label):
        """Tag an element with its highlight-dropdown cell-type label."""
        if ct_label:
            elem.set('data-celltype', ct_label)
            cell_types_seen.add(ct_label)

    def set_fill(elem, color, opacity_val):
        """Set fill, handling inline style="fill:none" and bare fill="none" attribute."""
        style = elem.get('style', '')
        if 'fill' in style:
            style = re.sub(r'fill\s*:[^;]+', f'fill:{color}', style)
            elem.set('style', style)
        else:
            elem.set('fill', color)
        # Also override a bare fill="none" attribute that would win via cascade
        if elem.get('fill', '').lower() == 'none':
            elem.set('fill', color)
        elem.set('fill-opacity', str(opacity_val))

    _NAT_SHAPE_TAGS = (
        f'{{{SVG_NS}}}path', f'{{{SVG_NS}}}circle',
        f'{{{SVG_NS}}}rect', f'{{{SVG_NS}}}ellipse',
    )

    def color_paths_in_element(el, color, tip, ct_label=None, umap_ct_label=None):
        """Color all drawable shape children of el (path/circle/rect/ellipse).

        Natanella's SVG renders several cell types (Vascular, Phloem
        Parenchyma, Phloem companion — and the guard/mesophyll/epidermal/
        trichome cells too) as literal <circle> elements rather than <path>
        groups. Only matching <path> here meant those circles never got
        set_fill/add_tooltip/tag_celltype applied, so they had no fill,
        no data-tooltip, and no data-celltype — i.e. no hover behavior at
        all, even though the color-selection logic upstream matched fine.

        umap_ct_label carries the *condition-qualified* label (e.g.
        "D0_Mesophyll") used by natanella's UMAP (umap_col =
        label_majorXcondition combines condition + cell type). It's kept
        separate from data-celltype/ct_label — which stays condition-agnostic
        ("Mesophyll") so the existing SVG-side Highlight dropdown still
        highlights that cell type across all four condition panels at once.
        """
        for tag in _NAT_SHAPE_TAGS:
            for p in el.iter(tag):
                set_fill(p, color, opacity)
                add_tooltip(p, tip)
                tag_celltype(p, ct_label)
                if umap_ct_label:
                    p.set('data-umap-celltype', umap_ct_label)

    if ds_key == 'rice':
        # ── Runtime-derived SVG ID → expression_dict key matching ────────────
        # expression_dict keys are "{Condition}_{CellAnnotation}" built from
        # obs columns. We don't hardcode the obs value spellings; instead we
        # derive the mapping at runtime so any spelling works.
        #
        # Step 1: extract all unique (cond, ct) pairs from expression_dict.
        # Keys are "{cond}_{ct}"; cond has no underscore (WW, MS, S, MD, D)
        # and ct may have spaces (e.g. "Large Parenchyma").
        _all_conds = set()
        _all_cts = set()
        for k in expression_dict:
            parts = k.split('_', 1)
            if len(parts) == 2:
                _all_conds.add(parts[0])
                _all_cts.add(parts[1])

        print(f"  Rice conds in expr_dict: {sorted(_all_conds)}", file=sys.stderr)
        print(f"  Rice cts   in expr_dict: {sorted(_all_cts)}", file=sys.stderr)

        def _normalise(s):
            """Lower-case, replace dots/hyphens with spaces, collapse whitespace."""
            return ' '.join(s.lower().replace('.', ' ').replace('-', ' ').split())

        # Step 2: normalised-string → actual obs value lookup.
        _norm_to_cond = {_normalise(c): c for c in _all_conds}
        _norm_to_ct = {_normalise(c): c for c in _all_cts}

        def _strip_digit_suffix(svg_id):
            """'Phloem.CC-7-0' -> 'Phloem.CC',  'Epidermis-39' -> 'Epidermis'"""
            return re.sub(r'(-\d+)+$', '', svg_id)

        # Resolve the actual obs Condition value that means "well watered".
        # Try matching the normalised form of the SVG group name
        # "Well.Watered" directly against the obs values first, then fall
        # back to a set of common control aliases. This is only a first
        # guess — if the obs spelling isn't covered by either of those, the
        # elimination logic below (run once panel/group labels are known)
        # will resolve it definitively regardless of spelling.
        _ww_cond_val = _resolve_rice_control(_norm_to_cond)
        print(f"  Rice coloring: _ww_cond_val (alias guess)={_ww_cond_val!r}", file=sys.stderr)

        def _resolve_panel_conds(label_norms):
            """
            label_norms: list of normalised label strings (or None) for
            every panel/group found in the SVG, in iteration order.

            Pass 1: direct match of each label against obs Condition values
            — this is what already resolves Mild/Moderate Drought and
            Mild/Moderate Salinity correctly, since their SVG text and obs
            spelling agree.

            Pass 2 (elimination): whichever obs Condition value was *not*
            claimed by any directly-matched label must be the control —
            the well-watered panel is drawn twice in the SVG but there's
            only one obs value left unclaimed once the other four panels
            are matched. This works no matter how the control happens to
            be spelled in obs ("WW", "Control", "CK", "Mock", "Baseline",
            ...) — no alias list needed.

            Returns (resolved, elim): resolved is a list aligned with
            label_norms giving the directly-matched cond_val (or None);
            elim is the single leftover obs value if exactly one remains,
            else None.
            """
            resolved = []
            matched_conds = set()
            for label_norm in label_norms:
                cond_val = _norm_to_cond.get(label_norm) if label_norm else None
                if cond_val:
                    matched_conds.add(cond_val)
                resolved.append(cond_val)
            leftover = _all_conds - matched_conds
            elim = next(iter(leftover)) if len(leftover) == 1 else None
            print(f"  Rice: directly-matched conds: {sorted(matched_conds)}", file=sys.stderr)
            print(f"  Rice: leftover (unmatched) conds: {sorted(leftover)}", file=sys.stderr)
            if elim:
                print(f"  Rice: control resolved via elimination -> {elim!r}", file=sys.stderr)
            return resolved, elim

        def _colour_rice_ct_group(ct_g, cond_val):
            """Colour every <path> in one cell-type group, given the
            already-resolved obs Condition value for the panel it's in."""
            ct_prefix = _strip_digit_suffix(ct_g.get('id', ''))
            ct_val = _norm_to_ct.get(_normalise(ct_prefix))
            if ct_val is None:
                return

            # Record this cell type for the highlight dropdown regardless
            # of whether the current gene happens to colour it.
            cell_types_seen.add(ct_val)

            key = f'{cond_val}_{ct_val}'
            v = expression_dict.get(key)
            if v is None:
                return

            if control == "WW":
                ctrl_key = f'{_ww_cond_val}_{ct_val}' if _ww_cond_val else None
                ctrl_val = expression_dict.get(ctrl_key) if ctrl_key else None
                color = expression_to_color_relative(v, ctrl_val, log2_max)
                if ctrl_val is not None:
                    try:
                        log2_r = v - ctrl_val
                        tip = (
                            f"{key}\nAvg expression: {v:.4f}\n"
                            f"WW expression: {ctrl_val:.4f}\n"
                            f"Log fold-change: {log2_r:.4f}"
                        )
                    except (ValueError, ZeroDivisionError):
                        tip = f"{key}\nAvg expression: {v:.4f}"
                else:
                    tip = f"{key}\nAvg expression: {v:.4f}\nWW expression: N/A"
            else:
                color = get_color(key)
                tip = f"{key}\nAvg expression: {raw_val(key)}"

            if not color or color == "#ffffff":
                return
            for p in ct_g.iter(f'{{{SVG_NS}}}path'):
                set_fill(p, color, opacity)
                add_tooltip(p, tip)
                tag_celltype(p, ct_val)

        samples_g = root.find(f'{{{SVG_NS}}}g[@id="samples"]')

        if samples_g is not None:
            # ── Current eFP-standard template ────────────────────────────────
            # All cell-type groups for *every* condition sit flat as direct
            # children of <g id="samples"> (no per-condition wrapper group),
            # repeated once per panel — one panel per condition, laid out in
            # a row x column grid (e.g. top row = drought panels, bottom row
            # = salinity panels, with a "Well Watered" control panel re-drawn
            # in both rows). We figure out which condition each panel belongs
            # to from where it sits on the canvas, matched against the
            # condition names drawn as text in <g id="labels">.
            try:
                _vb = [float(x) for x in root.get('viewBox', '').split()]
                _canvas_w, _canvas_h = _vb[2], _vb[3]
            except (ValueError, IndexError):
                _canvas_w = float(re.sub(r'[^0-9.]', '', root.get('width', '') or '') or 1000)
                _canvas_h = float(re.sub(r'[^0-9.]', '', root.get('height', '') or '') or 1000)

            _m_re = re.compile(r'^\s*[Mm]\s*([-\d.eE]+)[,\s]+([-\d.eE]+)')

            def _first_point(d_str):
                """Read the (x, y) of a path's first (absolute) moveto."""
                m = _m_re.match(d_str or '')
                if not m:
                    return None
                try:
                    return float(m.group(1)), float(m.group(2))
                except ValueError:
                    return None

            def _quadrant(x, y):
                col = 0 if x < _canvas_w / 3 else (1 if x < 2 * _canvas_w / 3 else 2)
                row = 'top' if y < _canvas_h / 2 else 'bottom'
                return row, col

            # (row, col) -> condition label, read straight off the aria-label
            # text drawn in the SVG, so panel↔condition assignment follows
            # the template even if the layout, labels, or the wrapper
            # group's id (seen as both "labels" and "label" across export
            # versions) change later. We deliberately don't look this up by
            # a specific group id — just scan the whole tree for any path
            # carrying an aria-label.
            _quad_to_label = {}
            for el in root.iter(f'{{{SVG_NS}}}path'):
                label = (el.get('aria-label') or '').strip()
                if not label:
                    continue
                pt = _first_point(el.get('d', ''))
                if pt is not None:
                    _quad_to_label[_quadrant(*pt)] = label
            print(f"  Rice quad->label: {_quad_to_label}", file=sys.stderr)

            # Chunk the flat list of cell-type groups under <g id="samples">
            # into panels: a new panel starts every time the first cell-type
            # category ("Epidermis") shows up again.
            panels = []
            current = []
            for ct_g in samples_g:
                if not hasattr(ct_g, 'get') or ct_g.tag != f'{{{SVG_NS}}}g':
                    continue
                if _strip_digit_suffix(ct_g.get('id', '')) == 'Epidermis' and current:
                    panels.append(current)
                    current = []
                current.append(ct_g)
            if current:
                panels.append(current)

            # Pass 1: figure out each panel's centroid/label without
            # colouring anything yet, so we can run elimination across all
            # panels before committing to any single panel's cond_val.
            _panel_labels = []  # normalised label per panel, aligned with `panels`
            for panel in panels:
                pts = []
                for g in panel:
                    for p in g.iter(f'{{{SVG_NS}}}path'):
                        pt = _first_point(p.get('d', ''))
                        if pt is not None:
                            pts.append(pt)
                if not pts:
                    _panel_labels.append(None)
                    continue
                cx = sum(p[0] for p in pts) / len(pts)
                cy = sum(p[1] for p in pts) / len(pts)
                label = _quad_to_label.get(_quadrant(cx, cy))
                _panel_labels.append(_normalise(label) if label else None)

            _resolved, _elim = _resolve_panel_conds(_panel_labels)
            if _elim and _ww_cond_val is None:
                _ww_cond_val = _elim

            # Pass 2: colour each panel using its directly-matched cond_val,
            # falling back to the alias guess or the elimination result for
            # whichever panel(s) couldn't be matched directly (this is
            # almost always the well-watered control, drawn twice).
            for panel, label_norm, cond_val in zip(panels, _panel_labels, _resolved):
                if cond_val is None:
                    if label_norm in (_RICE_CONTROL_ALIASES | {'well watered'}) and _ww_cond_val:
                        cond_val = _ww_cond_val
                    elif _elim is not None:
                        cond_val = _elim
                if cond_val is None:
                    print(f"  Rice: unresolved panel condition (label_norm={label_norm!r})", file=sys.stderr)
                    continue
                for ct_g in panel:
                    _colour_rice_ct_group(ct_g, cond_val)

        else:
            # ── Old template fallback ────────────────────────────────────────
            #   layer1 > <g id="Well.Watered-1"> > <g id="Epidermis-39"> > <path>
            layer1 = root.find(f'{{{SVG_NS}}}g[@id="layer1"]')
            cond_groups = [
                cond_g for cond_g in (layer1 if layer1 is not None else root)
                if hasattr(cond_g, 'get') and cond_g.get('id', '')
            ]

            _group_labels = [
                _normalise(_strip_digit_suffix(cond_g.get('id', '')))
                for cond_g in cond_groups
            ]
            _resolved, _elim = _resolve_panel_conds(_group_labels)
            if _elim and _ww_cond_val is None:
                _ww_cond_val = _elim

            for cond_g, label_norm, cond_val in zip(cond_groups, _group_labels, _resolved):
                if cond_val is None:
                    if label_norm in (_RICE_CONTROL_ALIASES | {'well watered'}) and _ww_cond_val:
                        cond_val = _ww_cond_val
                    elif _elim is not None:
                        cond_val = _elim
                if cond_val is None:
                    continue  # top_bar, left_line, layer1, etc.

                for ct_g in cond_g:
                    if not hasattr(ct_g, 'get') or ct_g.tag != f'{{{SVG_NS}}}g':
                        continue
                    _colour_rice_ct_group(ct_g, cond_val)


    elif ds_key == 'arabidopsis_nat':
        id_to_el = {el.get('id'): el for el in root.iter() if el.get('id')}

        _CT_CONFIG = [
            ("Guard", lambda c: f"{c}_guard_cell"),
            ("Mesophyll", lambda c: f"{c}_mesophyll"),
            ("Epidermal", lambda c: f"{c}_epidermal_cell"),
            ("Trichome", lambda c: f"{c}_trichome"),
            ("Vascular", lambda c: _VASCULAR_CIRCLE_ID[c]),
            ("Phloem Parenchyma", lambda c: f"{c}_phloem_parenchyma_circle"),
            ("Phloem companion", lambda c: f"{c}_phloem_companion_circle"),
        ]

        _NAT_CELL_TYPES = [
            "Guard", "Mesophyll", "Epidermal", "Trichome",
            "Vascular", "Phloem Parenchyma", "Phloem companion",
        ]

        for ct, svg_id_fn in _CT_CONFIG:
            cell_types_seen.add(ct)
            w0_val = expression_dict.get(f"W0_{ct}")
            for c in conditions:
                key = f"{c}_{ct}"
                v = expression_dict.get(key)
                if v is None:
                    continue

                if control == "W0":
                    color = expression_to_color_relative(v, w0_val, log2_max)
                    try:
                        log2_r = v - w0_val
                        tip = f"{key}\nAvg expression: {v:.4f}\nW0 expression: {w0_val:.4f}\nLog fold-change: {log2_r:.4f}"
                    except (ValueError, ZeroDivisionError):
                        tip = f"{key}\nAvg expression: {v:.4f}"
                else:
                    color = get_color(key)
                    tip = f"{key}\nAvg expression: {raw_val(key)}"

                if not color or color == "#ffffff":
                    continue
                el = id_to_el.get(svg_id_fn(c))
                if el is not None:
                    color_paths_in_element(el, color, tip, ct, umap_ct_label=key)

        # Phloem (averaged)
        cell_types_seen.add("Phloem (avg)")
        for c in conditions:
            pp_key = f"{c}_Phloem Parenchyma"
            pc_key = f"{c}_Phloem companion"
            if pp_key not in expression_dict or pc_key not in expression_dict:
                continue
            pp_raw = expression_dict[pp_key]
            pc_raw = expression_dict[pc_key]
            avg_raw = (pp_raw + pc_raw) / 2

            if control == "W0":
                w0_pp = expression_dict.get("W0_Phloem Parenchyma")
                w0_pc = expression_dict.get("W0_Phloem companion")
                if w0_pp and w0_pc:
                    ctrl_avg = (w0_pp + w0_pc) / 2
                    color = expression_to_color_relative(avg_raw, ctrl_avg, log2_max)
                    try:
                        log2_r = avg_raw - ctrl_avg
                        tip = (
                            f"{c}_Phloem (averaged)\n"
                            f"Phloem Parenchyma: {pp_raw:.4f}\n"
                            f"Phloem Companion: {pc_raw:.4f}\n"
                            f"Avg expression: {avg_raw:.4f}\n"
                            f"Log fold-change: {log2_r:.4f}"
                        )
                    except (ValueError, ZeroDivisionError):
                        tip = f"{c}_Phloem (averaged)\nAvg expression: {avg_raw:.4f}"
                else:
                    continue
            else:
                color = expression_to_color(avg_raw, colour_ceiling, 'yellow_to_red')
                tip = (
                    f"{c}_Phloem (averaged)\n"
                    f"Phloem Parenchyma: {pp_raw:.4f}\n"
                    f"Phloem Companion: {pc_raw:.4f}\n"
                    f"Avg expression: {avg_raw:.4f}"
                )

            if not color or color == "#ffffff":
                continue
            el = id_to_el.get(f"{c}_phloem")
            if el is not None:
                color_paths_in_element(el, color, tip, "Phloem (avg)")

    elif ds_key == 'at_root_rs':
        for el in root.iter():
            el_id = el.get('id')
            if el_id not in expression_dict:
                continue
            v = expression_dict[el_id]
            ct_label = root_celltype_label(el_id)
            cell_types_seen.add(ct_label)

            if control == "Col-0":
                ctrl_val = None
                for mut_prefix in ("shr2", "scr4"):
                    if el_id.startswith(mut_prefix + "_"):
                        suffix = el_id[len(mut_prefix):]
                        ctrl_val = expression_dict.get(f"col0{suffix}")
                        break
                else:
                    ctrl_val = v
                color = expression_to_color_relative(v, ctrl_val, log2_max)
                if ctrl_val is not None:
                    try:
                        log2_r = v - ctrl_val
                        tip = f"{el_id}\nAvg expression: {v:.4f}\ncol0 expression: {ctrl_val:.4f}\nLog fold-change: {log2_r:.4f}"
                    except (ValueError, ZeroDivisionError):
                        tip = f"{el_id}\nAvg expression: {v:.4f}"
                else:
                    tip = f"{el_id}\nAvg expression: {v:.4f}\ncol0 expression: N/A"
            else:
                color = get_color(el_id)
                tip = f"{el_id}\nAvg expression: {raw_val(el_id)}"

            if not color or color == "#ffffff":
                continue

            if el.tag in (f'{{{SVG_NS}}}path', f'{{{SVG_NS}}}circle',
                          f'{{{SVG_NS}}}rect', f'{{{SVG_NS}}}ellipse'):
                set_fill(el, color, opacity)
                add_tooltip(el, tip)
                tag_celltype(el, ct_label)

            elif el.tag == f'{{{SVG_NS}}}g':
                for child in el.iter():
                    if child is el:
                        continue
                    if child.tag in (f'{{{SVG_NS}}}path', f'{{{SVG_NS}}}circle',
                                     f'{{{SVG_NS}}}rect', f'{{{SVG_NS}}}ellipse'):
                        set_fill(child, color, opacity)
                        add_tooltip(child, tip)
                        tag_celltype(child, ct_label)

    elif ds_key == 'at_seed_martin':
        for el in root.iter():
            el_id = el.get('id', '')
            if not el_id or el_id in _SEED_SKIP_IDS:
                continue
            if el.tag != f'{{{SVG_NS}}}g':
                continue
            if el_id not in expression_dict:
                continue

            ct_label = seed_celltype_label(el_id)
            cell_types_seen.add(ct_label)

            v = expression_dict[el_id]
            if control == 'median' and median_ctrl is not None:
                color = expression_to_color_relative(v, median_ctrl, log2_max)
                log2_r = v - median_ctrl
                tip = (f"{el_id}\nAvg expression: {v:.4f}\n"
                       f"Median control: {median_ctrl:.4f}\n"
                       f"Log fold-change: {log2_r:.4f}")
            else:
                color = get_color(el_id)
                tip = f"{el_id}\nAvg expression: {raw_val(el_id)}"

            if not color:
                continue

            for child in el.iter():
                if child is el:
                    continue
                if child.tag in (f'{{{SVG_NS}}}path', f'{{{SVG_NS}}}circle',
                                 f'{{{SVG_NS}}}rect', f'{{{SVG_NS}}}ellipse'):
                    set_fill(child, color, opacity)
                    add_tooltip(child, tip)
                    tag_celltype(child, ct_label)

    elif ds_key in ('at_flower_lee', 'at_silique_lee', 'at_stem_lee', 'at_seed_0d_lee', 'at_rosette_21d_lee', 'at_shoot_zhang', 'at_rosette_30d_lee', 'at_seedling_3d_lee', 'at_seedling_6d_lee', 'at_seedling_12d_lee'):
        _lee_paint_tags = (f'{{{SVG_NS}}}path', f'{{{SVG_NS}}}circle',
                           f'{{{SVG_NS}}}rect', f'{{{SVG_NS}}}ellipse',
                           f'{{{SVG_NS}}}polygon')
        for el in root.iter():
            el_id = el.get('id', '')
            if not el_id or el_id not in expression_dict:
                continue

            ct_label = el_id
            cell_types_seen.add(ct_label)

            v = expression_dict[el_id]
            if control == 'median' and median_ctrl is not None:
                color = expression_to_color_relative(v, median_ctrl, log2_max)
                log2_r = v - median_ctrl
                tip = (f"{el_id}\nAvg expression: {v:.4f}\n"
                       f"Median control: {median_ctrl:.4f}\n"
                       f"Log fold-change: {log2_r:.4f}")
            else:
                color = get_color(el_id)
                tip = f"{el_id}\nAvg expression: {raw_val(el_id)}"

            if not color:
                continue

            if el.tag in _lee_paint_tags:
                set_fill(el, color, opacity)
                add_tooltip(el, tip)
                tag_celltype(el, ct_label)

            elif el.tag == f'{{{SVG_NS}}}g':
                # Clear group-level fill="none" so it doesn't cascade and override children
                if el.get('fill', '').lower() == 'none':
                    el.attrib.pop('fill')
                for child in el.iter():
                    if child is el:
                        continue
                    if child.tag in _lee_paint_tags:
                        set_fill(child, color, opacity)
                        add_tooltip(child, tip)
                        tag_celltype(child, ct_label)

    if gene_name:
        add_gene_title(root, gene_name)

    # ── Legend — use per-control label ──────────────────────────────────────
    _CONTROL_LABELS = {
        'W0': 'W0',
        'WW': 'WW',
        'Col-0': 'Col-0',
    }
    if control in _CONTROL_LABELS:
        add_relative_legend_labeled(root, gene_name or "Expression", log2_max,
                                    control_label=_CONTROL_LABELS[control])
    elif control == 'median':
        # Seed / Lee: diverging scale centred on median of cell-type means.
        # No reference circle here — 'median' isn't a real experimental
        # control (unlike W0/WW/Col-0), so there's nothing to anchor a
        # circle to.
        add_relative_legend_labeled(root, gene_name or "Expression", log2_max,
                                    control_label="Median")
    else:
        add_legend(root, gene_name or "Expression", colour_ceiling)

    return tree, sorted(cell_types_seen)


# ──────────────────────────────────────────────────────────────────────────────
# H5AD loader
# ──────────────────────────────────────────────────────────────────────────────

def load_all(h5ad_path, gene_list, cell_type_column_1, cell_type_column_2=None, umap_col='label_majorXcondition'):
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    print(f"Loading full H5AD from {h5ad_path}", file=sys.stderr)
    adata = ad.read_h5ad(h5ad_path)
    print(f"  {adata.n_obs:,} cells x {adata.n_vars:,} genes", file=sys.stderr)

    # ── Use logcounts layer as X if X is empty (e.g. Lee datasets) ───────
    if 'logcounts' in adata.layers:
        import scipy.sparse as sp_check
        x = adata.X
        if x is None or (sp_check.issparse(x) and x.nnz == 0) or (not sp_check.issparse(x) and x.sum() == 0):
            print("  X matrix empty — using logcounts layer", file=sys.stderr)
            adata.X = adata.layers['logcounts']

    if cell_type_column_1 not in adata.obs.columns:
        raise ValueError(
            f"Column '{cell_type_column_1}' not found. Available: {list(adata.obs.columns)}"
        )

    tair_col = None
    for candidate in ("TAIR_ID", "tair_id", "gene_id"):
        if candidate in adata.var.columns:
            tair_col = candidate
            break

    sym_to_tair = {}
    tair_to_sym = {}
    if tair_col:
        for sym, tair in zip(adata.var_names, adata.var[tair_col]):
            if isinstance(tair, str) and tair:
                sym_to_tair[sym.upper()] = tair.upper()
                tair_to_sym[tair.upper()] = sym

    var_names_upper = {v.upper(): v for v in adata.var_names}

    def resolve(query):
        q = query.strip().upper()
        if q in var_names_upper:
            return var_names_upper[q]
        if q in tair_to_sym:
            return tair_to_sym[q]
        return None

    resolved = [r for g in gene_list if (r := resolve(g))]
    if not resolved:
        raise ValueError("No valid genes found in dataset.")

    gene = resolved[0]
    tair_id = sym_to_tair.get(gene.upper(), "")
    display_name = f"{gene}/{tair_id}" if tair_id and tair_id.upper() != gene.upper() else gene

    gene_idx = adata.var_names.get_loc(gene)

    if cell_type_column_2 == 'level_1_annotation_timed':
        X_col = adata.X[:, gene_idx]
        if sp.issparse(X_col):
            X_col = X_col.toarray().ravel()
        else:
            X_col = np.asarray(X_col).ravel()
        X_col = X_col.astype(float)

        cell_types_l2 = adata.obs[cell_type_column_1].astype(str).to_numpy()
        expression_dict = {
            str(ct): float(np.mean(X_col[cell_types_l2 == ct]))
            for ct in pd.unique(cell_types_l2)
            if (cell_types_l2 == ct).sum() > 0
        }

        if cell_type_column_2 in adata.obs.columns:
            cell_types_l1 = adata.obs[cell_type_column_2].astype(str).to_numpy()
            for ct in pd.unique(cell_types_l1):
                mask = cell_types_l1 == ct
                if mask.sum() > 0:
                    expression_dict[str(ct)] = float(np.mean(X_col[mask]))

        print(f"  Gene: {display_name}  ({len(expression_dict)} groups, L1+L2)", file=sys.stderr)

    elif cell_type_column_2 == 'Condition':
        X_col = adata.X[:, gene_idx]
        if sp.issparse(X_col):
            X_col = X_col.toarray().ravel()
        else:
            X_col = np.asarray(X_col).ravel()
        X_col = X_col.astype(float)

        cell_types = adata.obs[cell_type_column_1].astype(str).to_numpy()
        condition = adata.obs[cell_type_column_2].astype(str).to_numpy()

        # Keys: "{Condition}_{CellAnnotation}" — Condition is the first part
        # so that splitting on the first underscore always gives (cond, ct).
        expression_dict = {
            f"{cond}_{ct}": float(np.mean(X_col[(cell_types == ct) & (condition == cond)]))
            for ct in pd.unique(cell_types)
            for cond in pd.unique(condition)
            if ((cell_types == ct) & (condition == cond)).sum() > 0
        }
        print(f"  Gene: {display_name}  ({len(expression_dict)} groups, CAxCondition)", file=sys.stderr)
        print(f"  Sample keys: {list(expression_dict.keys())[:6]}", file=sys.stderr)

    elif cell_type_column_2 == 'Genotype':
        X_col = adata.X[:, gene_idx]
        if sp.issparse(X_col):
            X_col = X_col.toarray().ravel()
        else:
            X_col = np.asarray(X_col).ravel()
        X_col = X_col.astype(float)

        expression_dict = {}

        h5ad_to_id = {
            'Col-0 (shr2)_Root endodermis': ['shr2_top_endodermis', 'shr2_side_endodermis', 'shr2_tip_endodermis'],
            'Col-0 (shr2)_Collumella root cap': ['shr2_tip_columella'],
            'Col-0 (shr2)_Lateral root cap': ['shr2_tip_lateral_root_cap'],
            'Col-0 (shr2)_Phloem': ['shr2_top_phloem', 'shr2_side_phloem'],
            'Col-0 (shr2)_Root procambium': ['shr2_top_procambium', 'shr2_side_procambium'],
            'Col-0 (shr2)_Root cortex': ['shr2_tip_cortex', 'shr2_top_cortex', 'shr2_side_cortex'],
            'Col-0 (shr2)_Xylem pole pericycle': ['shr2_xpp_circle'],
            'Col-0 (shr2)_Protoxylem': ['shr2_proto_circle'],
            'Col-0 (shr2)_Phloem pole pericycle': ['shr2_ppp_circle'],
            'Col-0 (shr2)_Metaxylem': ['shr2_meta_circle'],
            'Ler (scr4)_Root endodermis': ['scr4_top_endodermis', 'scr4_side_endodermis', 'scr4_tip_endodermis'],
            'Ler (scr4)_Collumella root cap': ['scr4_tip_columella'],
            'Ler (scr4)_Lateral root cap': ['scr4_tip_lateral_root_cap'],
            'Ler (scr4)_Phloem': ['scr4_top_phloem', 'scr4_side_phloem'],
            'Ler (scr4)_Root procambium': ['scr4_top_procambium', 'scr4_side_procambium'],
            'Ler (scr4)_Root cortex': ['scr4_tip_cortex', 'scr4_top_cortex', 'scr4_side_cortex'],
            'Ler (scr4)_Xylem pole pericycle': ['scr4_xpp_circle'],
            'Ler (scr4)_Protoxylem': ['scr4_proto_circle'],
            'Ler (scr4)_Phloem pole pericycle': ['scr4_ppp_circle'],
            'Ler (scr4)_Metaxylem': ['scr4_meta_circle'],
            'Col-0_Root endodermis': ['col0_top_endodermis', 'col0_side_endodermis', 'col0_tip_endodermis'],
            'Col-0_Collumella root cap': ['col0_tip_columella'],
            'Col-0_Lateral root cap': ['col0_tip_lateral_root_cap'],
            'Col-0_Phloem': ['col0_top_phloem', 'col0_side_phloem'],
            'Col-0_Root procambium': ['col0_top_procambium', 'col0_side_procambium'],
            'Col-0_Root cortex': ['col0_tip_cortex', 'col0_top_cortex', 'col0_side_cortex'],
            'Col-0_Xylem pole pericycle': ['col0_xpp_circle'],
            'Col-0_Protoxylem': ['col0_proto_circle'],
            'Col-0_Phloem pole pericycle': ['col0_ppp_circle'],
            'Col-0_Metaxylem': ['col0_meta_circle'],
        }

        cell_types = adata.obs[cell_type_column_1].astype(str).to_numpy()
        genotype = adata.obs[cell_type_column_2].astype(str).to_numpy()
        intermediate_dict = {
            f"{geno}_{ct}": float(np.mean(X_col[(cell_types == ct) & (genotype == geno)]))
            for ct in pd.unique(cell_types)
            for geno in pd.unique(genotype)
            if ((cell_types == ct) & (genotype == geno)).sum() > 0 and f"{geno}_{ct}" in h5ad_to_id.keys()
        }

        for key in intermediate_dict.keys():
            for pov in h5ad_to_id[key]:
                expression_dict[pov] = intermediate_dict[key]

        for prefix in ('shr2', 'scr4', 'col0'):
            if f'{prefix}_proto_circle' not in expression_dict or f'{prefix}_meta_circle' not in expression_dict:
                continue
            expression_dict[f'{prefix}_top_xylem'] = (
                    (expression_dict[f'{prefix}_proto_circle'] + expression_dict[f'{prefix}_meta_circle']) / 2
            )
            expression_dict[f'{prefix}_side_xylem'] = expression_dict[f'{prefix}_top_xylem']
            expression_dict[f'{prefix}_top_phloem'] = (
                    (expression_dict[f'{prefix}_xpp_circle'] + expression_dict[f'{prefix}_ppp_circle']) / 2
            )
            expression_dict[f'{prefix}_side_phloem'] = expression_dict[f'{prefix}_top_phloem']

    else:
        X_col = adata.X[:, gene_idx]
        if sp.issparse(X_col):
            X_col = X_col.toarray().ravel()
        else:
            X_col = np.asarray(X_col).ravel()
        X_col = X_col.astype(float)

        cell_types = adata.obs[cell_type_column_1].astype(str).to_numpy()
        unique_cell_types = pd.unique(cell_types)

        expression_dict = {
            str(ct): float(np.mean(X_col[cell_types == ct]))
            for ct in unique_cell_types
            if (cell_types == ct).sum() > 0
        }

    print(f"  Gene: {display_name}  ({len(expression_dict)} cell-type groups)", file=sys.stderr)

    # ── UMAP ─────────────────────────────────────────────────────────────────
    umap_data = None

    if 'X_umap' in adata.obsm:
        umap_col_actual = umap_col if umap_col in adata.obs.columns else cell_type_column_1
        umap = np.asarray(adata.obsm['X_umap'])[:, :2].astype(float)
        labels = adata.obs[umap_col_actual].astype(str).to_numpy()

        X_col_full = adata.X[:, gene_idx]
        if sp.issparse(X_col_full):
            X_col_full = X_col_full.toarray().ravel()

        umap_df = pd.DataFrame({
            'UMAP_1': umap[:, 0],
            'UMAP_2': umap[:, 1],
            'cell_type': labels,
            'gene_expression': X_col_full.astype(float),
        })
        umap_data = umap_df.replace([np.inf, -np.inf], np.nan).dropna(
            subset=['UMAP_1', 'UMAP_2', 'gene_expression']
        )
        print(f"  UMAP rows loaded: {len(umap_data):,}", file=sys.stderr)
    else:
        print("  Warning: X_umap not found — UMAP panel will be skipped", file=sys.stderr)

    print(f"  Final UMAP dataframe size: {len(umap_data) if umap_data is not None else 0:,}", file=sys.stderr)
    return expression_dict, display_name, umap_data


# ──────────────────────────────────────────────────────────────────────────────
# UMAP Plotly figure → HTML div string
# ──────────────────────────────────────────────────────────────────────────────

def build_umap_div(umap_df, gene_name):
    import plotly.graph_objects as go
    import numpy as np

    if umap_df is None or len(umap_df) == 0:
        return "<div style='padding:20px;color:#b00'>No UMAP points available.</div>"

    x_all = umap_df['UMAP_1'].to_numpy(dtype=float)
    y_all = umap_df['UMAP_2'].to_numpy(dtype=float)
    expr_all = umap_df['gene_expression'].to_numpy(dtype=float)
    labels_all = umap_df['cell_type'].astype(str).to_numpy()

    cell_types = sorted(umap_df['cell_type'].astype(str).unique())

    fig = go.Figure()

    fig.add_trace(go.Scattergl(
        x=x_all, y=y_all,
        mode='markers',
        marker=dict(
            size=6,
            color=expr_all,
            colorscale=[[0.0, '#ffff00'], [1.0, '#ff0000']],
            showscale=True,
            colorbar=dict(title=f"{gene_name}<br>Expression"),
            opacity=0.6,
            line=dict(width=0.5, color='rgba(0,0,0,0.35)')
        ),
        # selected/unselected styling lets JS spotlight one cell type via
        # selectedpoints without needing a separate trace per cell type.
        selected=dict(marker=dict(opacity=0.95)),
        unselected=dict(marker=dict(opacity=0.08)),
        customdata=np.column_stack([labels_all, expr_all]),
        hovertemplate=(
            "<b>%{customdata[0]}</b><br>"
            "UMAP1: %{x:.2f}<br>"
            "UMAP2: %{y:.2f}<br>"
            "Expression: %{customdata[1]:.4f}<extra></extra>"
        ),
        name='All Cells',
        visible=True
    ))

    for ct in cell_types:
        ct_mask = (labels_all == ct)
        fig.add_trace(go.Scattergl(
            x=x_all[ct_mask], y=y_all[ct_mask],
            mode='markers',
            marker=dict(
                size=9,
                color='#1a6fcc',
                line=dict(width=0.5, color='#0a3d6b'),
                opacity=1.0
            ),
            customdata=np.column_stack([labels_all[ct_mask], expr_all[ct_mask]]),
            hovertemplate=(
                "<b>%{customdata[0]}</b><br>"
                "UMAP1: %{x:.2f}<br>"
                "UMAP2: %{y:.2f}<br>"
                "Expression: %{customdata[1]:.4f}<extra></extra>"
            ),
            name=ct,
            visible=False
        ))

    pad = 0.5
    x_min = float(np.min(x_all)) - pad
    x_max = float(np.max(x_all)) + pad
    y_min = float(np.min(y_all)) - pad
    y_max = float(np.max(y_all)) + pad

    buttons = [{
        "label": "All Cells",
        "method": "update",
        "args": [
            {"visible": [True] + [False] * len(cell_types)},
            {
                "title": {"text": f"UMAP - All Cells - {gene_name}"},
                "xaxis": {"range": [x_min, x_max], "title": "UMAP_1"},
                "yaxis": {"range": [y_min, y_max], "title": "UMAP_2"},
            }
        ]
    }]

    for i, ct in enumerate(cell_types):
        visible = [True] + [False] * len(cell_types)
        visible[i + 1] = True

        buttons.append({
            "label": ct,
            "method": "update",
            "args": [
                {"visible": visible},
                {
                    "title": {"text": f"UMAP - {ct} highlighted - {gene_name}"},
                    "xaxis": {"range": [x_min, x_max], "title": "UMAP_1"},
                    "yaxis": {"range": [y_min, y_max], "title": "UMAP_2"},
                }
            ]
        })

    fig.update_layout(
        updatemenus=[dict(
            buttons=buttons,
            direction="down",
            pad={"r": 10, "t": 10},
            showactive=True,
            x=0.01, xanchor="left",
            y=0.99, yanchor="top",
            bgcolor="white",
            bordercolor="gray",
            borderwidth=1
        )],
        title={"text": f"UMAP - All Cells - {gene_name}"},
        xaxis=dict(title="UMAP_1", range=[x_min, x_max]),
        yaxis=dict(title="UMAP_2", range=[y_min, y_max]),
        template="plotly_white",
        showlegend=False,
        height=650,
        margin=dict(t=60, b=40, l=40, r=40),
        plot_bgcolor = "#f2f2f2",
        paper_bgcolor = "#ffffff",
    )

    return fig.to_html(full_html=False, include_plotlyjs='cdn', div_id='umap-plot')


def build_violin_div(umap_df, gene_name, cell_types_order=None):
    """
    Plotly violin plot showing per-cell expression distributions by cell type.
    Complements the pseudobulk eFP SVG by revealing distributional shape
    (bimodality, outliers, zero-inflation, spread).

    Handles zero-inflated scRNA-seq data by:
      - Using bandwidth > 0 so the KDE doesn't collapse to a line
      - Showing the % of expressing cells in hover annotations
      - Providing interactive buttons to toggle between all cells,
        expressing-only view, and log-scale

    Parameters
    ----------
    umap_df : pd.DataFrame
        Must have columns 'cell_type' and 'gene_expression'.
    gene_name : str
        Gene display name for the title.
    cell_types_order : list[str] | None
        If given, show only these cell types in this order (matches the
        SVG highlight dropdown).  Otherwise sorted by descending median.

    Returns an HTML div string (no full page, no plotly.js include —
    relies on the UMAP section having already loaded plotly via CDN).
    """
    import plotly.graph_objects as go
    import numpy as np

    if umap_df is None or len(umap_df) == 0:
        return ""

    ct_col = umap_df['cell_type'].astype(str)
    expr_col = umap_df['gene_expression'].to_numpy(dtype=float)

    # ── Determine cell-type ordering ─────────────────────────────────────
    if cell_types_order:
        available = set(ct_col.unique())
        cell_types = [ct for ct in cell_types_order if ct in available]
        if not cell_types:
            cell_types = sorted(available)
    else:
        cell_types = sorted(ct_col.unique())

    # Sort by median expression descending so the most-expressed cell
    # types appear on the left — easier to scan visually.
    ct_medians = {}
    for ct in cell_types:
        vals = expr_col[ct_col == ct]
        ct_medians[ct] = float(np.median(vals)) if len(vals) else 0.0
    cell_types = sorted(cell_types, key=lambda c: ct_medians[c], reverse=True)

    if not cell_types:
        return ""

    # ── Compute stats for annotations ────────────────────────────────────
    ct_stats = {}
    for ct in cell_types:
        vals = expr_col[ct_col == ct]
        n = len(vals)
        n_expr = int(np.sum(vals > 0))
        pct = (n_expr / n * 100) if n else 0.0
        ct_stats[ct] = {
            'n': n, 'n_expr': n_expr, 'pct': pct,
            'mean': float(np.mean(vals)) if n else 0.0,
            'median': float(np.median(vals)) if n else 0.0,
            'max': float(np.max(vals)) if n else 0.0,
        }

    # ── Compute a sensible global bandwidth ──────────────────────────────
    # When data is zero-inflated, the default Plotly bandwidth collapses
    # to near-zero, rendering a flat line instead of a violin shape.
    # We set a minimum bandwidth based on the global range of non-zero
    # expression values so the KDE always produces a visible curve.
    all_nonzero = expr_col[expr_col > 0]
    if len(all_nonzero) > 1:
        global_range = float(np.ptp(all_nonzero))
        global_std = float(np.std(all_nonzero))
        # Silverman-like minimum: ensure visible width even for sparse data
        min_bw = max(global_range * 0.05, global_std * 0.3, 0.01)
    else:
        min_bw = 0.1

    # ── Build figure ─────────────────────────────────────────────────────
    fig = go.Figure()

    colors = [
        '#636EFA', '#EF553B', '#00CC96', '#AB63FA',
        '#FFA15A', '#19D3F3', '#FF6692', '#B6E880',
        '#FF97FF', '#FECB52',
    ]

    for i, ct in enumerate(cell_types):
        mask = (ct_col == ct)
        vals = expr_col[mask]
        s = ct_stats[ct]
        color = colors[i % len(colors)]

        # Per-cell-type bandwidth: use the larger of the per-group
        # Silverman estimate and the global minimum so the violin is
        # always visible.
        nz = vals[vals > 0]
        if len(nz) > 1:
            ct_bw = max(float(np.std(nz)) * 0.5, min_bw)
        else:
            ct_bw = min_bw

        fig.add_trace(go.Violin(
            y=vals,
            name=ct,
            box_visible=True,
            meanline_visible=True,
            points='outliers',
            pointpos=0,
            jitter=0.4,
            bandwidth=ct_bw,
            hoveron='violins+points',
            hoverinfo='y+name',
            scalemode='width',
            spanmode='soft',
            marker=dict(size=3, color=color, opacity=0.5),
            line=dict(color=color),
            fillcolor=color,
            opacity=0.7,
            customdata=[[s['n'], s['n_expr'], s['pct'],
                         s['mean'], s['median'], s['max']]]*len(vals),
            hovertemplate=(
                '<b>%{fullData.name}</b><br>'
                'Expression: %{y:.3f}<br>'
                'Cells: %{customdata[0]:,} '
                '(%{customdata[2]:.1f}% expressing)<br>'
                'Mean: %{customdata[3]:.3f} · '
                'Median: %{customdata[4]:.3f}<br>'
                '<extra></extra>'
            ),
        ))

    # ── Annotation bar: % expressing per cell type ───────────────────────
    annotations = []
    for i, ct in enumerate(cell_types):
        s = ct_stats[ct]
        annotations.append(dict(
            x=i, y=s['max'] + (s['max'] * 0.08 if s['max'] > 0 else 0.05),
            text=f"<b>{s['pct']:.0f}%</b>",
            showarrow=False,
            font=dict(size=10, color='#555'),
            yanchor='bottom',
        ))

    # ── Interactive buttons ──────────────────────────────────────────────
    fig.update_layout(
        title=dict(
            text=(f"Per-Cell Expression Distribution — {gene_name}"
                  f"<br><sup style='color:#888'>Percentages show fraction "
                  f"of expressing cells (> 0) per cell type</sup>"),
            font=dict(size=15),
        ),
        yaxis_title="Expression",
        xaxis_title="Cell Type",
        template="plotly_white",
        showlegend=False,
        height=520,
        margin=dict(t=90, b=140, l=60, r=40),
        xaxis=dict(tickangle=-45),
        violinmode='group',
        annotations=annotations,
        # Toggle buttons for log scale and point display
        updatemenus=[
            # Points toggle
            dict(
                type='buttons',
                direction='left',
                x=1.0, xanchor='right',
                y=1.15, yanchor='top',
                showactive=True,
                buttons=[
                    dict(label='Outliers',
                         method='restyle',
                         args=[{'points': 'outliers',
                                'jitter': 0.4, 'pointpos': 0}]),
                    dict(label='All Points',
                         method='restyle',
                         args=[{'points': 'all',
                                'jitter': 0.4, 'pointpos': 0}]),
                    dict(label='No Points',
                         method='restyle',
                         args=[{'points': False}]),
                ],
                font=dict(size=11),
                pad=dict(r=5, t=0),
            ),
            # Log-scale toggle
            dict(
                type='buttons',
                direction='left',
                x=0.0, xanchor='left',
                y=1.15, yanchor='top',
                showactive=True,
                buttons=[
                    dict(label='Linear',
                         method='relayout',
                         args=[{'yaxis.type': 'linear'}]),
                    dict(label='Log',
                         method='relayout',
                         args=[{'yaxis.type': 'log'}]),
                ],
                font=dict(size=11),
                pad=dict(r=5, t=0),
            ),
        ],
    )

    # include_plotlyjs=False because the UMAP div above already loads it
    # via CDN.  If UMAP is absent, we fall back to 'cdn' so the script
    # tag is still emitted (see build_html integration below).
    return fig.to_html(full_html=False, include_plotlyjs=False)


TOOLTIP_JS = """
<style>
  #efp-tooltip {
    position: fixed;
    background: rgba(30,30,30,0.92);
    color: #fff;
    padding: 7px 11px;
    border-radius: 5px;
    font-size: 13px;
    font-family: monospace;
    pointer-events: none;
    white-space: pre;
    z-index: 9999;
    display: none;
    box-shadow: 0 2px 8px rgba(0,0,0,0.35);
    line-height: 1.5;
  }
  .hoverable {
    cursor: pointer;
  }
  .group-hover-active {
    transition: stroke-width 0.1s, stroke 0.1s;
  }
  .celltype-dim {
    opacity: 0.12;
    transition: opacity 0.15s;
  }
  .celltype-highlight {
    stroke: #111111 !important;
    stroke-width: 2.5px !important;
    filter: drop-shadow(0 0 4px rgba(0,0,0,0.55));
    transition: opacity 0.15s;
  }
</style>
<div id="efp-tooltip"></div>
<script>
  const tooltip = document.getElementById('efp-tooltip');

  document.addEventListener('mousemove', function(e) {
    const el = e.target;
    const tip = el.getAttribute('data-tooltip');
    if (tip) {
      tooltip.textContent = tip;
      tooltip.style.display = 'block';
      const pad = 14;
      let x = e.clientX + pad;
      let y = e.clientY + pad;
      if (x + tooltip.offsetWidth > window.innerWidth)  x = e.clientX - tooltip.offsetWidth - pad;
      if (y + tooltip.offsetHeight > window.innerHeight) y = e.clientY - tooltip.offsetHeight - pad;
      tooltip.style.left = x + 'px';
      tooltip.style.top  = y + 'px';
    } else {
      tooltip.style.display = 'none';
    }
  });

  document.addEventListener('mouseleave', function() {
    tooltip.style.display = 'none';
  }, true);

  function attachHoverableListeners() {
    // Listen only on <path>/<circle> that were actually coloured (tagged
    // 'hoverable' by add_tooltip during SVG colouring) — this excludes
    // purely structural elements like the "outlines" and "label" groups,
    // which are never in expression_dict and so never get tooltipped.
    document.querySelectorAll('#svg-string path.hoverable, #svg-string circle.hoverable, #svg-string polygon.hoverable').forEach(function(el) {
      if (el.dataset.hoverBound) return;
      el.dataset.hoverBound = '1';

      el.addEventListener('mouseenter', function() {
        // Immediate parent <g> — this is the cell-type group
        var group = el.parentElement;
        if (!group || group.tagName.toLowerCase() !== 'g') return;

        // Highlight all direct child paths/circles in this <g>
        group.querySelectorAll(':scope > path, :scope > circle, :scope > polygon').forEach(function(sib) {
          if (!sib.dataset.origStroke) {
            sib.dataset.origStroke = sib.getAttribute('stroke-width') || '1';
            sib.dataset.origStrokeColor = sib.getAttribute('stroke') || '';
          }
          sib.setAttribute('stroke-width', parseFloat(sib.dataset.origStroke) + 2);
          sib.setAttribute('stroke', '#333');
          sib.classList.add('group-hover-active');
        });

        // Mirror the hover onto the UMAP below (debounced — see
        // scheduleUmapHighlight — so fast mouse movement across many
        // SVG shapes doesn't fire a Plotly restyle per pixel).
        // data-umap-celltype (condition-qualified, e.g. "D0_Mesophyll")
        // takes priority over the plain data-celltype: natanella's UMAP
        // clusters are labeled condition+celltype together, so the bare
        // cell-type name alone ("Mesophyll") would never match a cluster.
        scheduleUmapHighlight(el.getAttribute('data-umap-celltype') || el.getAttribute('data-celltype'));
      });

      el.addEventListener('mouseleave', function() {
        document.querySelectorAll('.group-hover-active').forEach(function(sib) {
          sib.setAttribute('stroke-width', sib.dataset.origStroke || '1');
          if (sib.dataset.origStrokeColor) {
            sib.setAttribute('stroke', sib.dataset.origStrokeColor);
          } else {
            sib.removeAttribute('stroke');
          }
          sib.classList.remove('group-hover-active');
        });

        scheduleUmapHighlight(null);
      });
    });
  }
  attachHoverableListeners();

  /* ── SVG-hover → UMAP highlight ─────────────────────────────────────
     Debounced so rapidly sweeping the mouse across many SVG shapes
     coalesces into a single Plotly.restyle call instead of one per
     mouseenter/mouseleave, which is what actually causes the UI to
     lag when someone hovers quickly. */
  let __umapCellIndex = null;   // { cellType: [pointIndices...] }
  let __umapHoverTimer = null;
  const UMAP_HOVER_DEBOUNCE_MS = 45;

  function buildUmapCellIndex() {
    const gd = document.getElementById('umap-plot');
    if (!gd || !gd.data || !gd.data[0] || !gd.data[0].customdata) return;
    const cd = gd.data[0].customdata;
    const idx = {};
    for (let i = 0; i < cd.length; i++) {
      const ct = cd[i][0];
      (idx[ct] || (idx[ct] = [])).push(i);
    }
    __umapCellIndex = idx;
  }

  function scheduleUmapHighlight(cellType) {
    const gd = document.getElementById('umap-plot');
    if (!gd) return;
    if (__umapHoverTimer) clearTimeout(__umapHoverTimer);
    __umapHoverTimer = setTimeout(function() {
      __umapHoverTimer = null;
      if (!__umapCellIndex) buildUmapCellIndex();
      if (!__umapCellIndex) return;
      if (!cellType) {
        Plotly.restyle(gd, { selectedpoints: [null] }, [0]);
        return;
      }
      const points = __umapCellIndex[cellType] || [];
      Plotly.restyle(gd, { selectedpoints: [points] }, [0]);
    }, UMAP_HOVER_DEBOUNCE_MS);
  }
</script>
"""


# ──────────────────────────────────────────────────────────────────────────────
# HTML page builder
# ──────────────────────────────────────────────────────────────────────────────

def build_html(svg_string, umap_div, gene_name, active_ds_key, active_col,
               relative_svg_string=None, cell_types=None, violin_div=None):
    umap_section = ""
    if umap_div:
        umap_section = f"""
            <div class="section-title">Single-Cell UMAP - {gene_name}</div>
            <div id="umap-container">
                {umap_div}
            </div>
        """

    violin_section = ""
    if violin_div:
        # If there's no UMAP (and thus no plotly CDN tag), we need to load
        # plotly ourselves for the violin to render.
        plotly_cdn = ""
        if not umap_div:
            plotly_cdn = '<script src="https://cdn.plot.ly/plotly-latest.min.js"></script>'
        violin_section = f"""
            {plotly_cdn}
            <div class="section-title">Expression Distribution — {gene_name}</div>
            <div id="violin-container" style="background:white; border:1px solid #ddd;
                 border-radius:6px; padding:8px; overflow:hidden;">
                {violin_div}
            </div>
        """

    ds_options = ""
    for key, meta in DATASETS.items():
        sel = ' selected' if key == active_ds_key else ''
        ds_options += f'            <option value="{key}"{sel}>{meta["label"]}</option>\n'

    js_col_map = json.dumps(DATASET_COLUMNS)
    js_default_genes = json.dumps(DATASET_DEFAULT_GENES)
    js_marker_genes = json.dumps(DATASET_MARKER_GENES)

    active_ds_cols = DATASET_COLUMNS.get(active_ds_key, {active_col: active_col})
    col_options = ""
    for key, label in active_ds_cols.items():
        sel = ' selected' if key == active_col else ''
        col_options += f'            <option value="{key}"{sel}>{label}</option>\n'

    active_col_json = json.dumps(active_col)

    svg_absolute_js = json.dumps(svg_string)
    svg_relative_js = json.dumps(relative_svg_string or "")

    absrel_btn = (
        '<button id="absrel" onclick="toggleView(this.textContent)">absolute</button>'
        if relative_svg_string
        else ''
    )

    highlight_dropdown = ""
    if cell_types:
        ct_options = "\n".join(
            f'                <option value="{ct}">{ct}</option>' for ct in cell_types
        )
        highlight_dropdown = f"""
            <div class="ctrl-group" id="highlight-group">
                <label for="highlight-select">Highlight:</label>
                <select id="highlight-select" onchange="currentHighlight=this.value; applyHighlight(this.value);">
                    <option value="">None</option>
{ct_options}
                </select>
            </div>
"""

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>eFP Viewer - {gene_name}</title>
    <style>
        * {{ box-sizing: border-box; margin: 0; padding: 0; }}
        body {{
            font-family: Arial, sans-serif;
            background: #f5f5f5;
            color: #222;
        }}
        header {{
            background: #2c5f2e;
            color: white;
            padding: 16px 32px;
            display: flex;
            align-items: center;
            gap: 24px;
            flex-wrap: wrap;
        }}
        .header-title {{
            font-size: 20px;
            font-weight: bold;
            letter-spacing: 0.5px;
            white-space: nowrap;
        }}
        .header-title span {{
            font-weight: normal;
            font-style: italic;
            margin-left: 12px;
            font-size: 16px;
            opacity: 0.85;
        }}
        .header-controls {{
            display: flex;
            align-items: center;
            gap: 16px;
            margin-left: auto;
            flex-wrap: wrap;
        }}
        .ctrl-group {{
            display: flex;
            align-items: center;
            gap: 6px;
        }}
        .ctrl-group label {{
            font-size: 13px;
            opacity: 0.9;
            white-space: nowrap;
        }}
        .ctrl-group select {{
            padding: 6px 10px;
            font-size: 14px;
            border: none;
            border-radius: 4px;
            background: white;
            color: #222;
            outline: none;
            cursor: pointer;
            min-width: 130px;
        }}
        .ctrl-group select:focus {{
            box-shadow: 0 0 0 2px #a8d5a2;
        }}
        #gene-input {{
            padding: 6px 10px;
            font-size: 14px;
            border: none;
            border-radius: 4px;
            width: 160px;
            outline: none;
            font-family: monospace;
            letter-spacing: 0.5px;
        }}
        #gene-input:focus {{
            box-shadow: 0 0 0 2px #a8d5a2;
        }}
        #gene-input.invalid {{
            box-shadow: 0 0 0 2px #e05c5c;
            background: #fff5f5;
        }}
        /* ── Autocomplete dropdown ─────────────────────────────────── */
        .gene-input-wrapper {{
            position: relative;
        }}
        #gene-suggestions {{
            display: none;
            position: absolute;
            top: 100%;
            left: 0;
            right: 0;
            background: #fff;
            border: 1px solid #ccc;
            border-top: none;
            border-radius: 0 0 6px 6px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.15);
            z-index: 1000;
            max-height: 240px;
            overflow-y: auto;
            font-family: monospace;
            font-size: 13px;
            min-width: 280px;
        }}
        #gene-suggestions .sg-item {{
            padding: 7px 10px;
            cursor: pointer;
            border-bottom: 1px solid #f0f0f0;
            white-space: nowrap;
        }}
        #gene-suggestions .sg-item:last-child {{
            border-bottom: none;
        }}
        #gene-suggestions .sg-item:hover,
        #gene-suggestions .sg-item.sg-active {{
            background: #e8f5e9;
        }}
        #gene-suggestions .sg-item .sg-gene {{
            font-weight: bold;
            color: #2e7d32;
        }}
        #gene-suggestions .sg-item .sg-label {{
            color: #777;
            margin-left: 4px;
        }}
        #gene-suggestions .sg-header {{
            padding: 5px 10px;
            font-size: 11px;
            color: #999;
            background: #fafafa;
            border-bottom: 1px solid #eee;
            font-family: sans-serif;
        }}
        .search-btn {{
            padding: 6px 14px;
            background: #4a8f4c;
            color: white;
            border: none;
            border-radius: 4px;
            font-size: 14px;
            cursor: pointer;
            white-space: nowrap;
            transition: background 0.15s;
        }}
        .search-btn:hover {{ background: #3a7a3c; }}
        .container {{
            max-width: 1600px;
            margin: 0 auto;
            padding: 24px 32px;
        }}
        .section-title {{
            font-size: 15px;
            font-weight: bold;
            color: #444;
            margin: 28px 0 10px 0;
            padding-bottom: 6px;
            border-bottom: 2px solid #ddd;
        }}
        #svg-toolbar {{
            display: flex;
            align-items: center;
            gap: 16px;
            margin-bottom: 10px;
            flex-wrap: wrap;
        }}
        #svg-container {{
            background: white;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 16px;
            min-height: 500px;
            overflow-x: auto;
        }}
        #svg-container svg {{
            max-height: 100%;
            max-width: 100%;
            display: block;
            height: auto;
            margin: 0 auto;
        }}
        #umap-container {{
            background: white;
            border: 1px solid #ddd;
            border-radius: 6px;
            padding: 8px;
            overflow: hidden;
        }}
        #absrel {{
            margin-bottom: 0;
            padding: 4px 12px;
            font-size: 13px;
            cursor: pointer;
            border: 1px solid #aaa;
            border-radius: 4px;
            background: #f0f0f0;
        }}
        #absrel:hover {{ background: #e0e0e0; }}
        #highlight-group select {{
            padding: 4px 10px;
            font-size: 13px;
            border: 1px solid #aaa;
            border-radius: 4px;
            background: white;
            color: #222;
            cursor: pointer;
            min-width: 160px;
        }}
        #highlight-group label {{
            font-size: 13px;
            color: #444;
        }}

        /* Hide the popup by default */
        .popup-overlay {{
          display: none;
          position: fixed;
          top: 0;
          left: 0;
          width: 100%;
          height: 100%;
          background-color: rgba(0, 0, 0, 0.5); /* Semi-transparent background */
          justify-content: center;
          align-items: center;
          z-index: 1000;
        }}

        /* Show the popup */
        .popup-overlay.active {{
          display: flex;
        }}

        /* Popup Box Styling */
        .popup-content {{
          background: white;
          padding: 30px;
          border-radius: 8px;
          position: relative;
          width: 90%;
          max-width: 400px;
          text-align: center;
        }}

        /* Close Button Styling */
        .close-btn {{
          position: absolute;
          top: 10px;
          right: 15px;
          font-size: 24px;
          background: none;
          border: none;
          cursor: pointer;
        }}

    </style>
</head>

<body>
    <header>
        <div class="header-title">
            SUPeR Viewer
            <span>{gene_name}</span>
        </div>

        <form class="header-controls" method="GET" onsubmit="return validateGene()">

            <div class="ctrl-group">
                <label for="ds-select">Dataset:</label>
                <select id="ds-select" name="ds" onchange="updateColumnDropdown(this.value); updateDefaultGene(this.value);">
{ds_options}                </select>
            </div>

            <div class="ctrl-group">
                <label for="col-select">Column:</label>
                <select id="col-select" name="col">
{col_options}                </select>
            </div>

            <div class="ctrl-group">
                <label for="gene-input">Gene ID:</label>
                <div class="gene-input-wrapper">
                    <input
                        type="text"
                        id="gene-input"
                        name="gene"
                        value="{gene_name.split('/')[0]}"
                        placeholder="e.g. AT3G05727"
                        spellcheck="false"
                        autocomplete="off"
                        autocorrect="off"
                        autocapitalize="characters"
                    >
                    <div id="gene-suggestions"></div>
                </div>
            </div>

            <button type="submit" class="search-btn">&#x1F50D; Search</button>
        </form>
    </header>

    <div class="container">
        <div class="section-title">Tissue eFP - {gene_name}</div>

        <div id="svg-toolbar">
            {absrel_btn}
            {highlight_dropdown}
        </div>

        <div id="svg-container">
            <!-- Trigger Button to Open Popup -->
            <button id="openPopupBtn">Open Popup</button>

            <!-- Popup Overlay and Box -->
            <div id="popupOverlay" class="popup-overlay">
              <div class="popup-content">
                <button id="closePopupBtn" class="close-btn">&times;</button>
                <h2>Hello There!</h2>
                <p>This is a custom webpage popup window.</p>
              </div>
            </div>

            <div id="svg-string">
                {svg_string}
            </div>
        </div>
        {umap_section}
        {violin_section}
    </div>

    {TOOLTIP_JS}

    <script>
      const DATASET_COLUMNS = {js_col_map};
      const DATASET_DEFAULT_GENES = {js_default_genes};
      const DATASET_MARKER_GENES = {js_marker_genes};
      const ACTIVE_COL = {active_col_json};
      const SVG_ABSOLUTE = {svg_absolute_js};
      const SVG_RELATIVE = {svg_relative_js};
      let currentHighlight = "";

      /* ── Gene autocomplete / marker suggestions ──────────────── */
      (function() {{
        const input = document.getElementById('gene-input');
        const dropdown = document.getElementById('gene-suggestions');
        const dsSelect = document.getElementById('ds-select');
        let activeIdx = -1;

        function getMarkers() {{
          return DATASET_MARKER_GENES[dsSelect.value] || [];
        }}

        function renderSuggestions(items, showAll) {{
          if (!items.length) {{ dropdown.style.display = 'none'; return; }}
          let html = '<div class="sg-header">'
                   + (showAll ? 'Marker genes for this dataset'
                              : 'Matching markers')
                   + '</div>';
          items.forEach(function(m, i) {{
            html += '<div class="sg-item" data-gene="' + m.gene + '" data-idx="' + i + '">'
                  + '<span class="sg-gene">' + m.gene + '</span>'
                  + '<span class="sg-label">' + m.label.replace(m.gene + ' — ', '— ') + '</span>'
                  + '</div>';
          }});
          dropdown.innerHTML = html;
          dropdown.style.display = 'block';
          activeIdx = -1;
        }}

        function pick(gene) {{
          input.value = gene;
          input.classList.remove('invalid');
          dropdown.style.display = 'none';
        }}

        function setActive(idx, items) {{
          const nodes = dropdown.querySelectorAll('.sg-item');
          nodes.forEach(function(n) {{ n.classList.remove('sg-active'); }});
          if (idx >= 0 && idx < nodes.length) {{
            nodes[idx].classList.add('sg-active');
            nodes[idx].scrollIntoView({{ block: 'nearest' }});
          }}
          activeIdx = idx;
        }}

        input.addEventListener('focus', function() {{
          const q = input.value.trim().toUpperCase();
          const markers = getMarkers();
          if (!q) {{ renderSuggestions(markers, true); return; }}
          const filtered = markers.filter(function(m) {{
            return m.gene.toUpperCase().indexOf(q) !== -1
                || m.label.toUpperCase().indexOf(q) !== -1;
          }});
          renderSuggestions(filtered.length ? filtered : markers, !filtered.length);
        }});

        input.addEventListener('input', function() {{
          const q = input.value.trim().toUpperCase();
          const markers = getMarkers();
          if (!q) {{ renderSuggestions(markers, true); return; }}
          const filtered = markers.filter(function(m) {{
            return m.gene.toUpperCase().indexOf(q) !== -1
                || m.label.toUpperCase().indexOf(q) !== -1;
          }});
          if (filtered.length) {{
            renderSuggestions(filtered, false);
          }} else {{
            dropdown.style.display = 'none';
          }}
        }});

        input.addEventListener('keydown', function(e) {{
          const items = dropdown.querySelectorAll('.sg-item');
          if (!items.length || dropdown.style.display === 'none') return;
          if (e.key === 'ArrowDown') {{
            e.preventDefault();
            setActive(Math.min(activeIdx + 1, items.length - 1), items);
          }} else if (e.key === 'ArrowUp') {{
            e.preventDefault();
            setActive(Math.max(activeIdx - 1, 0), items);
          }} else if (e.key === 'Enter' && activeIdx >= 0) {{
            e.preventDefault();
            pick(items[activeIdx].getAttribute('data-gene'));
          }} else if (e.key === 'Escape') {{
            dropdown.style.display = 'none';
          }}
        }});

        dropdown.addEventListener('mousedown', function(e) {{
          // mousedown (not click) so it fires before blur
          const item = e.target.closest('.sg-item');
          if (item) {{
            e.preventDefault();
            pick(item.getAttribute('data-gene'));
          }}
        }});

        input.addEventListener('blur', function() {{
          // Small delay so click on dropdown can fire first
          setTimeout(function() {{ dropdown.style.display = 'none'; }}, 150);
        }});
      }})();

      function updateColumnDropdown(dsKey) {{
        const sel = document.getElementById('col-select');
        const cols = DATASET_COLUMNS[dsKey] || {{}};
        const prev = sel.value;
        sel.innerHTML = '';
        Object.entries(cols).forEach(function([value, label]) {{
          const opt = document.createElement('option');
          opt.value = value;
          opt.textContent = label;
          if (value === prev || value === ACTIVE_COL) opt.selected = true;
          sel.appendChild(opt);
        }});
      }}

      // Pre-fills the Gene ID box with a sensible default gene for whichever
      // dataset was just selected. Only fires on user-initiated dataset
      // changes (onchange), so it never overwrites the gene on initial load.
      function updateDefaultGene(dsKey) {{
        const geneInput = document.getElementById('gene-input');
        const defaultGene = DATASET_DEFAULT_GENES[dsKey];
        if (defaultGene) {{
          geneInput.value = defaultGene;
          geneInput.classList.remove('invalid');
        }}
      }}

      function toggleView(curr) {{
        const container = document.getElementById('svg-string');
        const btn = document.getElementById('absrel');
        if (curr === 'absolute') {{
          container.innerHTML = SVG_RELATIVE;
          btn.textContent = 'relative';
        }} else {{
          container.innerHTML = SVG_ABSOLUTE;
          btn.textContent = 'absolute';
        }}
        attachHoverableListeners();
        applyHighlight(currentHighlight);
      }}

      // Highlights every SVG element tagged with the chosen cell type
      // (data-celltype attribute) and dims everything else that was coloured.
      function applyHighlight(ct) {{
        document.querySelectorAll('#svg-string [data-celltype]').forEach(function(el) {{
          el.classList.remove('celltype-highlight', 'celltype-dim');
          if (ct) {{
            if (el.getAttribute('data-celltype') === ct) {{
              el.classList.add('celltype-highlight');
            }} else {{
              el.classList.add('celltype-dim');
            }}
          }}
        }});

        // Also fade purely structural top-level groups (e.g. "outlines",
        // "label") into the background whenever a highlight is active.
        // These never carry data-celltype themselves (they aren't coloured
        // cell-type data), so they're identified generically as: a direct
        // child of the SVG root, with an id, that contains no
        // data-celltype-tagged element anywhere inside it. This leaves the
        // "samples" wrapper (and flat cell-type groups on datasets with no
        // such wrapper) alone, since those DO contain data-celltype elements
        // and are already handled per cell-type above.
        var svgRoot = document.querySelector('#svg-string > svg');
        if (svgRoot) {{
          Array.from(svgRoot.children).forEach(function(el) {{
            if (!el.id) return;
            el.classList.remove('celltype-dim');
            if (ct && !el.querySelector('[data-celltype]')) {{
              el.classList.add('celltype-dim');
            }}
          }});
        }}
      }}

      const openBtn = document.getElementById('openPopupBtn');
      const closeBtn = document.getElementById('closePopupBtn');
      const overlay = document.getElementById('popupOverlay');

      // Show popup on button click
      openBtn.addEventListener('click', () => {{
        overlay.classList.add('active');
      }});

      // Hide popup on close button click
      closeBtn.addEventListener('click', () => {{
        overlay.classList.remove('active');
      }});

      // Hide popup when clicking outside the box
      overlay.addEventListener('click', (e) => {{
        if (e.target === overlay) {{
          overlay.classList.remove('active');
        }}
      }});

    </script>
</body>
</html>"""


# ──────────────────────────────────────────────────────────────────────────────
# Main
# ──────────────────────────────────────────────────────────────────────────────

def main():
    is_cgi = "REQUEST_METHOD" in os.environ

    if is_cgi:
        sys.stdout.write("Content-Type: text/html\r\n\r\n")
        sys.stdout.flush()

        from urllib.parse import parse_qs
        params = parse_qs(os.environ.get("QUERY_STRING", ""))

        ds_key = params.get("ds", [DEFAULT_DATASET])[0]
        col_key = params.get("col", [DEFAULT_COLUMN])[0]
        gene_list = params.get("gene", [DEFAULT_GENE])

        if ds_key not in DATASETS:
            ds_key = DEFAULT_DATASET

        allowed_cols = DATASET_COLUMNS.get(ds_key, {})
        if col_key not in allowed_cols:
            col_key = next(iter(allowed_cols), DEFAULT_COLUMN)

        cfg = DATASETS[ds_key]
        h5ad_path = cfg["h5ad"]
        svg_template = cfg["svg"]
        umap_col = cfg["umap_col"]
        opacity = cfg["opacity"]
        relative_svg_string = None

        def html_error(msg):
            sys.stdout.write(
                f"<html><body><pre style='color:red'>Error: {msg}</pre></body></html>"
            )
            sys.stdout.flush()

        try:
            if not os.path.exists(h5ad_path):
                html_error(f"H5AD file not found: {h5ad_path}")
                return
            if not os.path.exists(svg_template):
                html_error(f"SVG template not found: {svg_template}")
                return

            if ds_key == 'rice' and col_key == 'CAxCondition':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellAnnotation', 'Condition', umap_col
                )
            elif ds_key == 'at_root_rs' and col_key == 'TypexGenotype':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'Celltype', 'Genotype', umap_col
                )
            elif ds_key == 'at_seed_martin':
                if col_key == 'level_3_annotation_full_timed':
                    svg_template = cfg.get("svg_l3", svg_template)
                    expression_dict, gene_name, umap_df = load_all(
                        h5ad_path, gene_list,
                        'level_3_annotation_full_timed', None,
                        umap_col
                    )
                else:
                    expression_dict, gene_name, umap_df = load_all(
                        h5ad_path, gene_list,
                        'level_2_annotation_timed', 'level_1_annotation_timed',
                        umap_col
                    )
            elif ds_key == 'at_flower_lee':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellType', None, umap_col
                )
            elif ds_key == 'at_silique_lee':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellType', None, umap_col
                )
            elif ds_key == 'at_stem_lee':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellType', None, umap_col
                )
            elif ds_key == 'at_seed_0d_lee':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellType', None, umap_col
                )
            elif ds_key == 'at_rosette_21d_lee':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellType', None, umap_col
                )
            elif ds_key == 'at_shoot_zhang':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'celltype_after', None, umap_col
                )
            elif ds_key == 'at_rosette_30d_lee':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellType', None, umap_col
                )
            elif ds_key == 'at_seedling_3d_lee':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellType', None, umap_col
                )
            elif ds_key == 'at_seedling_6d_lee':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellType', None, umap_col
                )
            elif ds_key == 'at_seedling_12d_lee':
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, 'CellType', None, umap_col
                )
            else:
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, gene_list, col_key, None, umap_col
                )

            # ── Absolute SVG ─────────────────────────────────────────────────
            ET.register_namespace('', 'http://www.w3.org/2000/svg')
            svg_tree, cell_types = color_svg(
                svg_template,
                expression_dict,
                ds_key,
                gene_name=gene_name,
                opacity=opacity
            )
            buf = io.BytesIO()
            svg_tree.write(buf, encoding='utf-8', xml_declaration=False)
            svg_string = buf.getvalue().decode('utf-8')

            # ── Relative SVG ─────────────────────────────────────────────────
            _RELATIVE_CONTROL = {
                'arabidopsis_nat': 'W0',
                'rice': 'WW',
                'at_root_rs': 'Col-0',
                'at_seed_martin': 'median',
                'at_flower_lee': 'median',
                'at_silique_lee': 'median',
                'at_stem_lee': 'median',
                'at_seed_0d_lee': 'median',
                'at_rosette_21d_lee': 'median',
                'at_shoot_zhang': 'median',
                'at_rosette_30d_lee': 'median',
                'at_seedling_3d_lee': 'median',
                'at_seedling_6d_lee': 'median',
                'at_seedling_12d_lee': 'median',
            }
            rel_control = _RELATIVE_CONTROL.get(ds_key)
            if rel_control:
                ET.register_namespace('', 'http://www.w3.org/2000/svg')
                svg_tree_rel, _ = color_svg(
                    svg_template,
                    expression_dict,
                    ds_key,
                    gene_name,
                    opacity,
                    control=rel_control,
                )
                buf_rel = io.BytesIO()
                svg_tree_rel.write(buf_rel, encoding='utf-8', xml_declaration=False)
                relative_svg_string = buf_rel.getvalue().decode('utf-8')

            umap_div = build_umap_div(umap_df, gene_name) if umap_df is not None else None
            violin_div = build_violin_div(umap_df, gene_name,
                                          cell_types_order=cell_types) if umap_df is not None else None
            html = build_html(svg_string, umap_div, gene_name,
                              active_ds_key=ds_key, active_col=col_key,
                              relative_svg_string=relative_svg_string,
                              cell_types=cell_types,
                              violin_div=violin_div)

            sys.stdout.write(html)
            sys.stdout.flush()

        except Exception as e:
            import traceback
            print(traceback.format_exc(), file=sys.stderr)
            html_error(str(e))

    else:
        import argparse

        parser = argparse.ArgumentParser(description="eFP SVG + UMAP HTML viewer")
        parser.add_argument(
            "ds_key", choices=list(DATASETS.keys()),
            help="Dataset key: " + ", ".join(DATASETS.keys())
        )
        parser.add_argument("gene", nargs='+', help="Gene ID(s); first valid one is used")
        parser.add_argument("output", nargs='?', default="efp_viewer.html")
        parser.add_argument(
            "--col",
            choices=list(COLUMNS.keys()),
            default=DEFAULT_COLUMN,
            help="obs column to colour by (default: %(default)s)"
        )
        args = parser.parse_args()

        cfg = DATASETS[args.ds_key]
        h5ad_path = cfg["h5ad"]

        allowed_cols = DATASET_COLUMNS.get(args.ds_key, {})
        col_key = args.col if args.col in allowed_cols else next(iter(allowed_cols), DEFAULT_COLUMN)

        svg_template = cfg["svg"]
        umap_col = cfg["umap_col"]
        opacity = cfg["opacity"]
        relative_svg_string = None

        for f in [h5ad_path, svg_template]:
            if not os.path.exists(f):
                print(f"Error: file not found: {f}")
                sys.exit(1)

        if args.ds_key == 'rice' and col_key == 'CAxCondition':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellAnnotation', 'Condition', umap_col
            )
        elif args.ds_key == 'at_root_rs' and col_key == 'TypexGenotype':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'Celltype', 'Genotype', umap_col
            )
        elif args.ds_key == 'at_seed_martin':
            if col_key == 'level_3_annotation_full_timed':
                svg_template = cfg.get("svg_l3", svg_template)
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, args.gene,
                    'level_3_annotation_full_timed', None,
                    umap_col
                )
            else:
                expression_dict, gene_name, umap_df = load_all(
                    h5ad_path, args.gene,
                    'level_2_annotation_timed', 'level_1_annotation_timed',
                    umap_col
                )
        elif args.ds_key == 'at_flower_lee':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellType', None, umap_col
            )
        elif args.ds_key == 'at_silique_lee':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellType', None, umap_col
            )
        elif args.ds_key == 'at_stem_lee':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellType', None, umap_col
            )
        elif args.ds_key == 'at_seed_0d_lee':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellType', None, umap_col
            )
        elif args.ds_key == 'at_rosette_21d_lee':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellType', None, umap_col
            )
        elif args.ds_key == 'at_shoot_zhang':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'celltype_after', None, umap_col
            )
        elif args.ds_key == 'at_rosette_30d_lee':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellType', None, umap_col
            )
        elif args.ds_key == 'at_seedling_3d_lee':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellType', None, umap_col
            )
        elif args.ds_key == 'at_seedling_6d_lee':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellType', None, umap_col
            )
        elif args.ds_key == 'at_seedling_12d_lee':
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, 'CellType', None, umap_col
            )
        else:
            expression_dict, gene_name, umap_df = load_all(
                h5ad_path, args.gene, col_key, None, umap_col
            )

        # ── Absolute SVG ─────────────────────────────────────────────────────
        ET.register_namespace('', 'http://www.w3.org/2000/svg')
        svg_tree, cell_types = color_svg(
            svg_template,
            expression_dict,
            args.ds_key,
            gene_name=gene_name,
            opacity=opacity
        )
        buf = io.BytesIO()
        svg_tree.write(buf, encoding='utf-8', xml_declaration=False)
        svg_string = buf.getvalue().decode('utf-8')

        # ── Relative SVG ─────────────────────────────────────────────────────
        _RELATIVE_CONTROL = {
            'arabidopsis_nat': 'W0',
            'rice': 'WW',
            'at_root_rs': 'Col-0',
            'at_seed_martin': 'median',
            'at_flower_lee': 'median',
            'at_silique_lee': 'median',
            'at_stem_lee': 'median',
            'at_seed_0d_lee': 'median',
            'at_rosette_21d_lee': 'median',
            'at_shoot_zhang': 'median',
            'at_rosette_30d_lee': 'median',
            'at_seedling_3d_lee': 'median',
            'at_seedling_6d_lee': 'median',
            'at_seedling_12d_lee': 'median',
        }
        rel_control = _RELATIVE_CONTROL.get(args.ds_key)
        if rel_control:
            ET.register_namespace('', 'http://www.w3.org/2000/svg')
            svg_tree_rel, _ = color_svg(
                svg_template,
                expression_dict,
                args.ds_key,
                gene_name=gene_name,
                opacity=opacity,
                control=rel_control,
            )
            buf_rel = io.BytesIO()
            svg_tree_rel.write(buf_rel, encoding='utf-8', xml_declaration=False)
            relative_svg_string = buf_rel.getvalue().decode('utf-8')

        umap_div = build_umap_div(umap_df, gene_name) if umap_df is not None else None
        violin_div = build_violin_div(umap_df, gene_name, cell_types_order=cell_types) if umap_df is not None else None
        html = build_html(svg_string, umap_div, gene_name,
                          active_ds_key=args.ds_key, active_col=col_key,
                          relative_svg_string=relative_svg_string,
                          cell_types=cell_types,
                          violin_div=violin_div)

        with open(args.output, 'w', encoding='utf-8') as f:
            f.write(html)

        print(f"Written to: {args.output}")


if __name__ == "__main__":
    main()