# SUPeR Viewer

## Table of Contents

- [Software and Package Requirements](#software-and-package-requirements)
- [Background Information](#background-information)
- [SVG Illustration](#svg-illustration)
- [RDS/H5AD Data Manipulation](#rdsh5ad-data-manipulation)
  - [RDS to H5AD Conversion](#rds-to-h5ad-conversion)
  - [MEX to H5AD Conversion](#mex-to-h5ad-conversion)
  - [Reducing H5AD Size](#reducing-h5ad-size)
  - [How Can I Convert from Raw Counts and Z-Score to Log-Normalization?](#how-can-i-convert-from-raw-counts-and-z-score-to-log-normalization)
  - [MySQL Dump Generation](#mysql-dump-generation)

---

## Software and Package Requirements

- **R**: version 4.5.1
  - Packages: Seurat, Matrix, anndata
- **Python**: version 3.12.1
  - Packages: json, xml, anndata, pathlib
- **Inkscape**: version 1.4.2
  - Apply Transform extension for transform removal (instructions linked in original doc)

---

## Background Information

Currently, two techniques are commonly practiced for gene expression analysis: bulk RNA-seq and single cell RNA-seq (scRNA-seq).

For a given gene, conventional bulk RNA-seq returns one average value across an entire tissue sample containing millions of cells. Although it preserves information about the tissue of origin, only taking one average value masks cellular heterogeneity and potentially important variations amidst different cell types. For example, a plant tissue might contain guard cells, mesophyll cells, epidermal cells, and bundle sheath cells, each with a distinct expression identity; however, bulk RNA-seq collapses all these signals into one averaged measurement.

Despite the loss of anatomical context, scRNA-seq yields information at a higher resolution by establishing gene expression profiles for individual cells. It maintains cellular heterogeneity while allowing accurate identification of the original cell types. As a result, researchers can generate large-scale datasets containing thousands to millions of cells and gene expression features, including specific gene expression profiles, cell type identities, and the distribution of cell types within a tissue.

<img width="468" height="180" alt="image" src="https://github.com/user-attachments/assets/1eeeea7a-4705-416d-b776-1b0b5668f84d" />

In plant scRNA-seq, there are four major steps:

- Tissue dissociation and isolation
- Single cell capture and barcoding
- Amplification and sequencing
- Data Analysis

The plant cell walls are digested using cellulase and pectinase enzyme to create protoplasts, which are known to be "naked cells" that is only held together by its cell membrane. Flow cytometry (FANS) is then applied to filter out debris and sort/enrich DAPI-stained nuclei.

Afterwards, the isolated protoplasts are loaded into 10x Genomics chromium platforms to partition individual cells into tiny droplets. This process appends a unique molecular barcode while converting RNA to cDNA so that researchers can trace where the RNA molecule came from, thereby making scRNA-seq possible.

Further down the line, Illumina reads the actual tagged nucleotide sequences to generate the raw sequence counts. Online processing platforms such as Cell Ranger take Illumina's raw sequencing output as input to produce a UMI raw count matrix as output, which is what we see within some of our H5AD files before log-normalization. Lastly, some quality control processes are conducted before eventually using dimensionality reduction methods (tSNE, PCA, and UMAP) along with marker genes to group cells and identify specific cell types and plant development.

This project is described as a **pseudo-bulk viewer** because conventional bulk RNA-seq viewers take the average expression value across all cell types within a tissue. In contrast, this tool computes a separate average for each individual cell type using scRNA-seq data, so each cell type has its own color based on its expression level intensity, hence the attachment of the "pseudo-" prefix.

<img width="241" height="118" alt="image" src="https://github.com/user-attachments/assets/2a83f332-1af6-4971-81f3-7ba16806c2ca" />

<img width="208" height="118" alt="image" src="https://github.com/user-attachments/assets/8a6a8311-9d3b-4c63-9917-60bf3848807b" />

---

## SVG Illustration

Use an SVG program with a graphical user interface to create an SVG from scratch. Inkscape is highly recommended, and the following steps describe this program.

<img width="349" height="186" alt="image" src="https://github.com/user-attachments/assets/8d31a1aa-6237-4c38-976c-4ce6170b3d04" />

Open Inkscape. If you have a reference image from which to work, open the image or paste it on the canvas.

### Document Properties Setup

<img width="184" height="264" alt="image" src="https://github.com/user-attachments/assets/9473307b-d6cc-4cf7-b78d-c5240ab16d1c" />

First, open **File → Document Properties** (shortcut: `CMD/CTRL + SHIFT + D`). Make the following changes:

- **Format Unit**: change to px
- **Display Units**: change to px
- **Scale**: change to 1.000000
- **Orientation**: horizontal/landscape (preferably, but depending on the original layout)
- **Important**: expand the Viewbox dropdown to ensure that the view box origin coordinate is 0 for box X and Y
- **Extra tip**: as you organize your tissues and cell type groups, there will be a lot of extra space left on the edges; whenever this is the case, open the window again and click "Resize to content" to remove the extra space on the outer edges.

### Drawing Paths

In Inkscape, you draw images with paths connected by nodes. Draw your image with the *Draw Bezier curves and straight lines* tool (shortcut: `Shift + F6`).

<img width="104" height="155" alt="image" src="https://github.com/user-attachments/assets/430908cc-6963-4363-90de-8b8e5b9d7701" />

Drawing with Bezier curves is the ideal method to maximize detail while minimizing the number of coordinates needed to represent the path. Begin with your first node by clicking anywhere in the canvas and continue adding nodes by subsequent clicks elsewhere in the canvas. Finish your path by either closing the path (looping back to your first node) or pressing `ESC`.

If your path is closed, the path fill colour (selected through the bottom colour strip toolbar) will be enclosed within the loop. Otherwise, the fill will be enclosed in the loop that would be created when the first and last nodes are connected by a straight line.

You can select and resize your paths or reference image with the *Select and transform object* tool (shortcut: `F1` or `S` on macOS) and edit the nodes within a path with the *Edit nodes by path* tool (shortcut: `F2` or `N` on macOS).

### Finding Reference Images

A decent sample layer image is necessary for making your SVG tissues. Some publications will provide a decent PNG/JPG of all their cell types in their publication; if so, you can directly adapt their image for tracing.

 <img width="468" height="336" alt="image" src="https://github.com/user-attachments/assets/7f6ae683-435a-495d-97c0-8186cd57acb8" />

A publication with good Arabidopsis root image (top left corner)

Unfortunately, a lot of publications do not provide PNG/JPG of the cell types within their H5AD, so you would have to search the web for accurate and well-illustrated ones. Make the best attempt to find images of **TISSUES** so that it tells the structure and how multiple cell types spatially arrange relative to one another. If you cannot find entire tissues or certain cell types, you can draw the individual cell type shapes and put the shape of the entire plant as a visual reference.

<img width="468" height="499" alt="image" src="https://github.com/user-attachments/assets/022e288d-4cc5-4c5d-8a94-6cf9638122b7" />

An attempt to search for good anatomical image to produce SVG template

<img width="468" height="306" alt="image" src="https://github.com/user-attachments/assets/ff5f0767-d0dc-48dc-a54f-c559cb1623fa" />

Illustration of full plant (left) and individual cell types (right) due to a lack of quality tissue images

Copying and pasting are more convenient than uploading or letting Inkscape open a file from your directory. Once pasted onto the canvas layer, you can directly trace on top of it.

<img width="467" height="263" alt="image" src="https://github.com/user-attachments/assets/1a32295c-3b1b-4214-b75e-de6cbd30169b" />

<img width="280" height="160" alt="image" src="https://github.com/user-attachments/assets/2d3a5360-b7bd-42a6-a342-35342ca44894" />

The Bezier tool has many shortcuts. For example, holding while clicking will bend the path parabolically when you drag, but holding Shift while doing this will adjust the axis to which the parabola relatively bends.

### Grouping Tissues

You will want to create paths in your SVG that will eventually be colour-filled to represent gene expression levels. To accomplish this, group all cells that are within the same tissue using **Object → Group** (shortcut: `CTRL/CMD + G`). You can select multiple by pressing `CTRL/SHIFT` and clicking on the paths.

<img width="472" height="307" alt="image" src="https://github.com/user-attachments/assets/e2c57862-786d-450a-b224-8cdff61aa9a6" />

Ensure you have exactly one group for every tissue type your dataset includes, including when the tissue is only represented with a single path (i.e., make sure to group tissues with only one path too).

### Naming Groups via XML Editor

Open the XML editor with **Edit → XML Editor** (shortcut: `CTRL/CMD + SHIFT + X`).

<img width="472" height="266" alt="image" src="https://github.com/user-attachments/assets/2bee060a-aa9c-46f5-80fe-c91e2e54937c" />

XML panel on the bottom right section; renaming the <g> id attribute

If you are unsure if the path is closed properly, try filling the group with a colour from the bottom palette bar to make sure it fills the tissue representation as expected. However, make sure to change the FILL back to none afterwards. Key steps:

1. Double click on the group folder in Layers and Objects until all paths are highlighted in blue. There should be a large, dotted rectangle that includes all the paths in the group.(left side) <img width="558" height="308" alt="image" src="https://github.com/user-attachments/assets/a7f88f45-8711-49ad-a6e7-f4a88de4cbd1" />

2. Select a high-contrast color as the FILL on the bottom palette bar (e.g., Navy). <img width="544" height="287" alt="image" src="https://github.com/user-attachments/assets/7103d8e9-dbef-4a38-b361-5a88d132e267" />

3. Click the <img width="10" height="10" alt="image" src="https://github.com/user-attachments/assets/98708f42-6286-445c-a1b8-faf307072399" />
"remove fill" icon (bottom left corner) to remove fill.
<img width="580" height="326" alt="image" src="https://github.com/user-attachments/assets/e0f2803c-083b-4e6c-b45a-1f2945292c79" />

<img width="461" height="300" alt="image" src="https://github.com/user-attachments/assets/9a294600-1524-4f7a-9e6e-e2f72a639d21" />

The H5AD cell type names identical to the SVG `<g>` id attributes
Epidermis, Mesophyll, Phloem parenchyma, and etc.


From here, you can label your group with an ID that will then be matched to the contents of your XML file or the cell type names of the unpacked H5AD file.

Ensure your id labels begin with a letter of the alphabet and do not contain any special characters (except underscore `_` to replace spaces). The distinction between space (" ") vs underscore ("_") depends on the cell type names you see while unpacking the H5AD obs metadata. For example, the H5AD file might add an underscore as "phloem_parenchyma" or directly write "phloem parenchyma".

<img width="468" height="91" alt="image" src="https://github.com/user-attachments/assets/985d58d4-8bdc-4a99-b389-e2555f08cf7a" />

<img width="108" height="128" alt="image" src="https://github.com/user-attachments/assets/97945e37-8e50-4249-b142-f033c25a47f1" />

In this example, the H5AD uses a space for phloem parenchyma, thus "Phloem parenchyma" should be the ID name in the SVG to ensure seamless coloring. **In short, ALWAYS unpack the H5AD files and match SVG IDs with H5AD cell types to prevent further headaches down the line.**

**Important:** make sure that your SVG `<g>` tags do not have spaces in their id attributes because spaces might cause further issues down the line with HTML JavaScript selectors. If you insert an underscore where the H5AD has a space, keep in mind you'll need a regex check so the paths are colored correctly.

### Removing Transforms

Follow the <a href="https://github.com/Klowner/inkscape-applytransforms">linked instructions </a> to download Klowner's Inkscape *Apply Transform* tool at GitHub (you will need to restart Inkscape the first time you download this to use this tool). Deselect any selections you currently have and apply the tool through **Extensions → Modify Path → Apply Transform**. This removes transformations in every group by having their children inherit their transformations. This step is necessary because ePlant can only parse SVG tissue groups without transforms.

- **Note:** It is very possible to have H5AD-to-SVG cell type mismatches (e.g., "Dividing" and "Unknown" categories not represented visually). You should always try your best to include all cell types from the H5AD. It is acceptable to omit a cell type if you really struggle to find a good cell type image, or if the author explicitly asks not to include it. One solution is simply drawing a circle with the cell type caption underneath.

<img width="105" height="101" alt="image" src="https://github.com/user-attachments/assets/c8f18dc0-54b8-4776-bd6f-84a776efbe03" />

- **Note:** The *Apply Transform* tool requires all objects and shapes in the SVG to be paths, so you might get a warning that all elements need to be paths before transformation. To convert everything to a path, select all elements and go **Path → Object to Path** (shortcut: `CTRL/CMD + SHIFT + C`).

<img width="361" height="198" alt="image" src="https://github.com/user-attachments/assets/698ec336-acd1-4f8b-9c35-1fabc683be35" />

“Apply Transform” with no selections will apply transforms directly to the paths of the entire SVG. With selections, it will apply transforms only to the paths within the selection.

<img width="298" height="203" alt="image" src="https://github.com/user-attachments/assets/70fa001f-2052-4669-bf94-05c0981dc9ea" />
<img width="159" height="213" alt="image" src="https://github.com/user-attachments/assets/8ef47a37-ab5c-4607-911d-373e4c70267e" />

"Apply Transform" with no selections will apply transforms directly to the paths of the entire SVG. With selections, it will apply transforms only to the paths within the selection.

### Canonical ePlant SVG Structure

<img width="586" height="422" alt="image" src="https://github.com/user-attachments/assets/efb68656-ea18-4b88-b53a-edbbc50a885e" />

To produce an SVG that complies with the ePlant standard, complete the following:

- Organize the SVG into three major groups with IDs: **label**, **samples**, and **outlines**
  - **Label**: any text (captions, information, or subtitles) describing plant/tissue structure. Not data-fillable — purely informational.
  - **Samples**: cell type paths/shapes that you are coloring (expression values available). Direct children of "samples" are the cell-type groups (e.g., Vascular, Mesophyll) — these get colored by expression values from the H5AD file.
  - **Outlines**: non-fillable paths that are still informative in the SVG (e.g., arrows, guide lines, structural paths without expression data, or tissues/cells available in the illustration but without expression values in the H5AD file).
- The data-fillable cell type groups within "samples" (the direct children of samples) must not contain any further nested groups — they only contain `<path>` elements as direct children. If nesting happens during export, ungroup those paths before saving.
- Non-filled objects are never grouped with data-fillable groups. Within the raw script, element order equals paint order: move "outlines" content that should sit underneath (e.g., background guides) so it precedes "samples," ensuring it never covers a fillable path.

Example canonical structure (from `at_seedling_3d_opt.svg`, Lee et al. 2025 seedling dataset):

```
<svg id="svg1">          (3 top-level groups)
 ├─ <g id="outlines">    (73 <path> — non-fillable)
 ├─ <g id="samples">     (6 cell-type groups — fillable)
 │   ├─ Phloem parenchyma (1 path)
 │   ├─ Mesophyll (12 paths)
 │   ├─ Vascular (27 paths)
 │   ├─ Phloem (20 paths)
 │   ├─ Trichoblast (8 paths)
 │   └─ Epidermal (20 paths)
 └─ <g id="label">       (10 <path> — captions / text)
```

### Example SVG Snippet

```xml
<?xml version="1.0" encoding="UTF-8"?>
<svg id="svg1" width="417.31" height="320.45" version="1.1" viewBox="0 0 417.31 320.45"
xmlns="http://www.w3.org/2000/svg">
 <g id="outlines" fill="none" stroke="#000000">
  <path id="path2" d="..."/>
  <!-- additional outline paths -->
 </g>
 <g id="samples">
  <g id="Young_silique">
   <path id="path1" d="..." fill="none" stroke="#000000"/>
  </g>
  <g id="Mature_silique">
   <path id="path3" d="..." fill="none" stroke="#000000"/>
  </g>
  <g id="Seed_(silique)" fill="none" stroke="#000000">
   <path id="path9-7-80" d="..." stroke-width="1.0718"/>
   <!-- additional seed paths -->
  </g>
 </g>
</svg>
```

### Saving as Optimized SVG

Save your file normally under a different name, then save again under another different name as an **Optimized SVG** (under the "Save as type:" dropdown menu). You will be prompted to select parameters for your Optimized SVG. Recommended settings:

<img width="152" height="157" alt="image" src="https://github.com/user-attachments/assets/2e105c07-c50e-4173-b5d6-c23f110aa070" />
<img width="152" height="157" alt="image" src="https://github.com/user-attachments/assets/d77b70a8-62ac-4c26-90dd-0dd7aa957664" />
<img width="152" height="156" alt="image" src="https://github.com/user-attachments/assets/0faf724f-147b-442c-8972-fbaec1fad449" />

- Remove metadata: ✔
- Embed raster images: ✔
- Format output with line-breaks and indentation: ✔
- Strip the "xml:space" attribute from the root SVG element: ✔
- Convert CSS attributes to XML attributes: ✔
- Work around renderer bugs: ✔
- Preserve manually created IDs not ending with digits: ✔

**Note:** *Convert CSS attributes to XML attributes* matters because a leftover inline `style="fill:#..."` overrides the attribute the viewer sets, and the tissue then never colors.

Lastly, open the SVG product in the browser and view page source. Since it is optimized, the SVG should **NOT** contain information relevant to Inkscape or `transform="..."` (double check this with `CMD/CTRL + F`). The product should be clean and minimalistic.

<img width="468" height="262" alt="image" src="https://github.com/user-attachments/assets/f5dd5a6b-5b0c-4d46-a2d9-8ece0da2c3ee" />

As a reference, here is what it looks like WITHOUT the optimizations below, so you know what to avoid:

<img width="468" height="243" alt="image" src="https://github.com/user-attachments/assets/e16fa21c-8fd0-4e97-9d16-286db9fe96e3" />

Viola! You can now move on to the next step of SVG coloring.

---

## RDS/H5AD Data Manipulation

With the SUPeR Viewer, we work with H5AD datasets, and our CellxGene matrices all follow a log-normalized standard with the gene expression values.

For log-normalization, we standardize the per-cell count to 10000 and then apply `log1p`. The `target_sum=1e4` can be understood as:

```
expr(gene X, cell Y) = ln[1 + (count of gene X in Y × 10000) / (total gene counts in cell Y, i.e. library size)]
```

This ensures that the total gene count in every cell sums to 10000.

**Note:** in the pseudobulk viewer, the mean cell type values are calculated from log1p values from ALL cells, NOT just the 10000 subsampled cells.

### RDS to H5AD Conversion

Whether it is from a researcher or the GEO database, you should convert the RDS file into an H5AD file. The script `rds_to_h5ad.r` extracts two components: 1) the RAW COUNT CellxGene matrix and 2) the UMAP coordinates, to produce a minimized H5AD file.

#### Command Line Usage

```bash
Rscript rds_to_h5ad.R [rds dataset name] [output h5ad dataset name]
Rscript rds_to_h5ad.R lee_rosette_30d.rds lee_rosette_30d_min.h5ad
```

#### Script: `rds_to_h5ad.R`

```r
library(Seurat)
library(Matrix)
library(anndata)

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript rds_to_h5ad.R <input.rds> <output.h5ad> [assay] [reduction] [celltype_col]")
}

rds_path <- args[1]
h5ad_path <- args[2]
assay_name <- if (length(args) >= 3) args[3] else "RNA"
reduction_name <- if (length(args) >= 4) args[4] else "umap"
celltype_col <- if (length(args) >= 5) args[5] else "CellType"

message("Loading: ", rds_path)
obj <- readRDS(rds_path)

message("Extracting raw counts from assay '", assay_name, "'")
counts <- obj[[assay_name]]$counts
counts_t <- Matrix::t(counts)
counts_t <- as(counts_t, "CsparseMatrix")
gene_ids <- rownames(counts)
cell_ids <- colnames(counts)
rm(counts); gc()

lib_sizes <- Matrix::rowSums(counts_t)
counts_t@x <- counts_t@x / lib_sizes[counts_t@i + 1L] * 1e4
counts_t@x <- log1p(counts_t@x)

if (!celltype_col %in% colnames(obj@meta.data)) {
  stop("Column '", celltype_col, "' not found in obj@meta.data. Available: ",
       paste(colnames(obj@meta.data), collapse = ", "))
}

obs_df <- data.frame(
  CellType = as.character(obj@meta.data[[celltype_col]]), row.names = cell_ids,
  stringsAsFactors = FALSE
)
var_df <- data.frame(row.names = gene_ids)

message("Extracting reduction '", reduction_name, "'")
umap_mat <- Embeddings(obj, reduction = reduction_name)
umap_mat <- umap_mat[cell_ids, , drop = FALSE]
colnames(umap_mat) <- c("UMAP1", "UMAP2")
rm(obj); gc()

message("Building AnnData (", nrow(counts_t), " cells x ", ncol(counts_t), " genes)")
ad_py <- reticulate::import("anndata")
ad <- ad_py$AnnData(
  X = counts_t,
  obs = obs_df,
  var = var_df,
  obsm = list(X_umap = umap_mat)
)

message("Writing: ", h5ad_path)
ad$write_h5ad(h5ad_path)
message("Done.")
```

*(This R script conducts RDS → H5AD conversion with raw counts and UMAP coordinates.)*

#### Mini Python Normalization Script (when given a RAW COUNT H5AD directly, no RDS)

```python
import scanpy as sc

adata = sc.read_h5ad("at_seed_martin.h5ad")
adata1 = sc.read_h5ad("at_flower_lee.h5ad")
adata2 = sc.read_h5ad("at_stem_lee.h5ad")
adata3 = sc.read_h5ad("at_shoot_zhang.h5ad")

sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.normalize_total(adata1, target_sum=1e4)
sc.pp.normalize_total(adata2, target_sum=1e4)
sc.pp.normalize_total(adata3, target_sum=1e4)

sc.pp.log1p(adata)
sc.pp.log1p(adata1)
sc.pp.log1p(adata2)
sc.pp.log1p(adata3)

adata.write_h5ad("at_seed_martin_normalized.h5ad")
adata1.write_h5ad("at_flower_lee_normalized.h5ad")
adata2.write_h5ad("at_stem_lee_normalized.h5ad")
adata3.write_h5ad("at_shoot_zhang_normalized.h5ad")
```

### MEX to H5AD Conversion

*(Another file format — section not elaborated in the source document.)*

### Reducing H5AD Size

This is a relatively simple process requiring a short Python script with a few packages and poking around within each H5AD dataset.

Typically, what you need to construct the SUPeR Viewer are the `X_umap` from `obsm` and `CellType` from `obs`. The other columns and information can be deleted with multiple `del` statements or a loop that iterates through a list of attributes to be deleted. This step is necessary to minimize the space required to save the H5AD files.

**Note:** When completing this step, make sure to add `compression='gzip'` in the `.write_h5ad` method. This parameter is crucial because it most effectively minimizes H5AD file sizes.

Example workflow:

```python
import anndata as ad
import numpy as np
import scipy.sparse as sp

adata1 = ad.read_h5ad("stem.h5ad")

print(adata1)
print(adata1.obs['CellType'].unique())

adata1.raw = None
del adata1.obs['nCount_RNA']
del adata1.obs['nFeature_RNA']
del adata1.obs['orig.ident']
del adata1.obs['percent.mt']
del adata1.obs['percent.cp']
del adata1.obs['integrated_snn_res.0.5']
del adata1.uns['seurat_clusters']
del adata1.uns
del adata1.layers

adata1.write_h5ad("stem_slim.h5ad", compression='gzip')
```

To accelerate the process of cell type identification, `unpack.py` can be used to quickly find which metadata column most likely contains the conditions and cell types. It works from a flat, case-insensitive vocabulary list of cell-type terms (e.g., "mesophyll", "epidermis", "guard cell", "trichome", "phloem", "xylem", "cortex", "endodermis", "pericycle", "cambium", "parenchyma", "meristem", "sieve", "companion cell", "stele", "columella", "root cap", "procambium", "vascular", "stomat", "collenchyma", "sclerenchyma", "tracheid", "protoxylem", "metaxylem", "seed coat", "endosperm", "embryo", "cotyledon", "hypocotyl", "shoot apical", "suspensor", "unlabeled", "unlabelled", "unknown", "quiescent"), plus name hints ("celltype", "cell_type", "cell type", "annotation", "label", "type") that bump a column's ranking, and a deprioritization list for columns that are never cell-type columns regardless of values (e.g., "n_genes", "n_counts", "pct_counts", "total_counts", "doublet", "barcode", "sample", "batch", "condition", "genotype", "replicate", "timepoint", "index", "leiden", "louvain", "cluster").

Example command line usage:

```bash
python3 unpack.py lee_seedling_6d_min.h5ad
```

Example output:

```
Reading lee_seedling_6d_min.h5ad ...
41,314 cells x 22,376 genes
obsm keys : ['X_umap']
layers    : []
uns keys  : []

===========================================================
obs column report (1 columns)
===========================================================

Top 1 candidates, ranked by vocabulary match + column-name hints + cardinality:

* 'CellType'                          dtype=category  n_unique= 5   score=10
    vocab hits: 4 values matched (epidermal, meristem, meristematic, mesophyll, stele)
    sample values: ['Epidermal', 'Guard', 'Meristematic', 'Mesophyll', 'Stele']
```

<img width="328" height="204" alt="image" src="https://github.com/user-attachments/assets/4ae7350f-6220-4f58-9dbb-d27e218047fe" />

### How Can I Convert from Raw Counts and Z-Score to Log-Normalization?

As a great man of science once said, data is dirty in the real world.

While working with different datasets from many publications, the H5AD CellxGene matrix expression values might be manipulated in different ways, but the SUPeR viewer requires library-size log-normalization.

<img width="115" height="113" alt="image" src="https://github.com/user-attachments/assets/7493b63b-fb94-477f-895c-d5b86b163e1a" />
<img width="144" height="113" alt="image" src="https://github.com/user-attachments/assets/a1fb1de6-33fc-4784-b984-1fab51ba4106" />
<img width="177" height="113" alt="image" src="https://github.com/user-attachments/assets/58bd3764-b73e-4a64-a5fd-8155f193e84b" />

#### Command Line Usage

To generate a report of the statistical overview:

```bash
python3 aggregate_control_candidates.py
cd ~/h5ad_files
/mnt/home/sqiao/venv/bin/python3 aggregate_control_candidates.py
```

Patterns to get a quick first judgement:

- **Raw counts**: very big positive integers as max, and 0 as minimum
- **Log normalized**: decimals with not extremely large max values, and minimums close to 0
- **Z-Score**: negative minimums

You will likely find yourself in two situations:

- You have the RDS file close by, OR the H5AD CellxGene matrix is RAW COUNTS
- You only have the H5AD file, and it is Z-SCORE

If you have the RDS file or raw count H5AD, that is good news — you can use the R (`rds_to_h5ad.r`) or Python script above to directly work with the raw count layer to obtain the standard CellxGene matrix.

If you have done some H5AD analysis and realize it has been z-score transformed, you can look for the RDS file on GEO, scPlantDB, or PlantscRNAdb to get a hold of the raw counts. To search effectively, use the publication information such as BioProject ID, Author Name, or GEO accession numbers.

### MySQL Dump Generation

The SUPeR viewer requires two different scripts (`generate_umap_dumps.py` and `generate_pseudobulk_dumps.py`) to generate the expression values for SVG coloring and UMAP coordinates.

#### Command Line Usage

*Note: no arguments means dumps for all H5ADs and n=10000 UMAP subsample.*

```bash
python generate_pseudobulk_dumps.py --ds [dataset names/paths] --outdir [output directory path]
python generate_umap_dumps.py --datasets [dataset names/paths] --outdir [output directory path] --max-cells [# of cells to subsample (e.g. 5000 or 15000)]

/mnt/home/sqiao/venv/bin/python3 ~/public_html/cgi-bin/db_dump/generate_pseudobulk_dumps.py --ds at_root rice
```

The annotations and comments must conform to the current existing sample dumps on BAR API.

#### Pseudobulk Dump

The pseudobulk dump has 1 table containing 4 columns:

- **Data_probeset_id**: gene id (e.g., AT1G01010)
- **Data_signal**: mean of log-normalized values across all cells
- **Data_signal_std**: standard deviation of values from the population, NOT subsample
- **Data_bot_id**: cell type name corresponding with H5AD (e.g., D0_mesophyll, vascular, etc.)

**Note:** For `data_probeset_id`, gene models are recommended over gene symbols for a more standardized organization (e.g., prefer AT3G55980 over SZF1). There are mapping dictionaries in the BAR and many other public resources containing this conversion.

Pseudobulk analysis is a way of collapsing single-cell RNA-seq data back down to something that behaves like classic bulk RNA-seq, but broken out by cell type instead of by whole-tissue sample: for each gene, you group all the individual cells belonging to a given cell type/condition (e.g., R15_Vascular, W0_Epidermal) and compute a summary statistic (mean, standard deviation, etc.) rather than working with the noisy, sparse single-cell counts directly. This is done because single-cell data is extremely noisy and sparse at the individual-cell level (dropout, low counts per cell, technical variation), so aggregating to the cell-type level produces a much more statistically stable and interpretable signal — exactly what's needed for something like the SUPeR viewer, where one color is needed per SVG cell type.

```sql
DROP TABLE IF EXISTS `sample_data`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `sample_data` (
  `data_probeset_id` varchar(32) NOT NULL,
  `data_signal`       float        DEFAULT '0',
  `data_signal_std`   float        DEFAULT '0',
  `data_bot_id`       varchar(64) NOT NULL,
  UNIQUE KEY `uq_probeset_bot` (`data_probeset_id`,`data_bot_id`),
  KEY `data_probeset_id` (`data_probeset_id`,`data_bot_id`,`data_signal`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;
```

#### UMAP Dump

The UMAP dump has two separate tables:

```sql
DROP TABLE IF EXISTS `umap_coords`;
/*!40101 SET @saved_cs_client     = @@character_set_client */;
/*!50503 SET character_set_client = utf8mb4 */;
CREATE TABLE `umap_coords` (
  `cell_id`    INT          NOT NULL,
  `umap_1`     FLOAT        NOT NULL,
  `umap_2`     FLOAT        NOT NULL,
  `cell_type`  VARCHAR(128) NOT NULL,
  PRIMARY KEY (`cell_id`)
) ENGINE=InnoDB DEFAULT CHARSET=latin1;
/*!40101 SET character_set_client = @saved_cs_client */;

TABLE_DDL_EXPRESSION = """
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
```

**The UMAP coordinates table:**

- `cell_id`: a number (indexed 0–9999) designated to each of the 10000 proportionally sampled cells, to know the expression value for any corresponding gene.
  - Example: if the original H5AD has a ratio of 5:4:1 among cell types 1, 2, and 3, then the 10000 random sample would consist of 5000 type 1 cells, 4000 type 2 cells, and 1000 type 3 cells.
- `umap_1` and `umap_2`: coordinates of the UMAP points to know where to plot
- `cell_type`: cell type name corresponding with H5AD (e.g., D0_mesophyll, vascular, etc.)

**The UMAP expression table:**

- `gene_id`: the name of the gene (e.g., AT1G01010)
- `expression`: a JSON object of expression values for the given gene_id
  - Example: `{"43":1.152292,"44":1.546603,"46":1.392931,"76":0.911963}`
  - The key is the `cell_id` from the other table
  - The value is the expression level
  - Indices that are missing (47 to 75 in this case) are assumed to have the value 0

**Very Important:** regardless of the package or method used to unpack the H5AD, always convert the H5AD from COMPRESSED SPARSE ROW (CSR) to COMPRESSED SPARSE COLUMN (CSC). This step reduces the dump generation process per dataset from ~60 minutes to ~1 minute. This conversion should be incorporated within both `generate_umap_dumps.py` and `generate_pseudobulk_dumps.py`, since these scripts directly extract values from the H5AD.

```python
if sp.issparse(adata.X):
    if not isinstance(adata.X, sp.csc_matrix):
        print("  [performance] Converting expression matrix to CSC format for rapid column slicing...",
              file=sys.stderr)
        t_conv = time.perf_counter()
        adata.X = adata.X.tocsc()
        print(f"    Converted in {time.perf_counter() - t_conv:.1f}s", file=sys.stderr)
```

---

*Note: All embedded screenshots and figures from the original PDF have been omitted; this document preserves the text content only.*


## References and Publications

## License

MIT License

## Acknowledgements

I would like to acknowledge Prof. Nicholas J. Provart, Vincent Lau, Asher Pasha, and Reena Obmina for supervising and setting up the BAR API.

## Contact

Steven Qiao — University of Toronto
- NSERC USRA Award
- Lab: Provart Lab, Department of Cell & Systems Biology
- Email: steven.qiao@mail.utoronto.ca
