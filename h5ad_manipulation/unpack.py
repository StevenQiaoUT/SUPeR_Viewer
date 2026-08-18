import anndata as ad
import numpy as np
import scipy.sparse as sp

adata1 = ad.read_h5ad("control_full.h5ad")

print(adata1)
print(adata1.obs['CellType'].unique())


"""
del adata1.obs['nCount_RNA']
del adata1.obs['nFeature_RNA']
del adata1.obs['orig.ident']
del adata1.obs['percent.mt']
del adata1.obs['percent.cp']
del adata1.obs['integrated_snn_res.0.5']
del adata1.obs['seurat_clusters']
del adata1.obs['ident']
del adata1.uns
del adata1.layers
del adata1.obsm['PCA']
del adata1.obsm['UMAP']

adata1.write_h5ad("silique_slim.h5ad", compression='gzip')

adata1.raw = None

del adata1.obs['nCount_RNA']
del adata1.obs['nFeature_RNA']
del adata1.obs['orig.ident']
del adata1.obs['percent.mt']
del adata1.obs['percent.cp']
del adata1.obs['integrated_snn_res.0.5']
del adata1.obs['seurat_clusters']
del adata1.obs['ident']
del adata1.uns
del adata1.layers

adata1.write_h5ad("flower_slim.h5ad", compression='gzip')

"""
"""
adata = sc.read_h5ad("my_file.h5ad")

print(adata)
print(adata.X)
for cat in adata.obs['level_1_annotation_timed'].cat.categories:
    print(cat)
for cat in adata.obs['level_2_annotation_timed'].cat.categories:
    print(cat)
for cat in adata.obs['level_3_annotation_abbr'].cat.categories:
    print(cat)
for cat in adata.obs['level_3_annotation_full_timed'].cat.categories:
    print(cat)
print(adata.layers['logcounts'])

cols_to_keep = [
    'level_1_annotation_timed',
    'level_2_annotation_timed',
    'level_3_annotation_full_timed',
    'level_3_annotation_abbr'
]

adata.obs = adata.obs[cols_to_keep]

del adata.uns
del adata.layers

adata.write_h5ad("my_file_clean.h5ad", compression='gzip')
print(adata)

adata1.raw = None

del adata1.obs['Orig.ident']
del adata1.obs['nCount_RNA']
del adata1.obs['nFeature_RNA']
del adata1.obs['Percent.mt']
del adata1.obs['Seurat_clusters']
del adata1.obs['Dataset']
del adata1.obs['Tissue']
del adata1.obs['Organ']
del adata1.obs['Condition']
del adata1.obs['Libraries']
del adata1.obs['ACE']

adata1.write_h5ad('slim_root.h5ad', compression='gzip')


adata_1 = sc.read_h5ad("OsLeafStressIntegrated.h5ad")

print(adata_1)
print(adata_1.obs)
print(adata_1.var)
print(adata_1.obs['CellAnnotation'].unique())
print(adata_1.obs['Condition'].unique())
print(adata_1.X[:5, :5])        # what are the actual values?
print(adata_1.layers.keys())    # are there other layers (raw, counts, normalized)?
print(adata_1.uns.keys())       # any metadata about normalization?

adata_1.raw = None

adata_1.obs = adata_1.obs.drop(columns=[
    'Identity',
    'nCount_RNA',
    'nFeature_RNA',
    'percent.mt',
    'percent.pt',
    'TissueSystem',
    'Pseudotime'

del adata_1.obs['Identity']
del adata_1.obs['nCount_RNA']
del adata_1.obs['nFeature_RNA']
del adata_1.obs['percent.mt']
del adata_1.obs['percent.pt']
del adata_1.obs['TissueSystem']
del adata_1.obs['Pseudotime']
del adata_1.obsm['X_pca']
del adata_1.varm['PCs']
del adata_1.obsp['distances']
del adata_1.uns['neighbors']

adata_1.write_h5ad('slim.h5ad', compression='gzip')

adata_2 = sc.read_h5ad('brapa.h5ad')
print(adata_2)

adata_2.raw = None

del adata_2.obs['orig.ident']
del adata_2.obs['nCount_RNA']
del adata_2.obs['nFeature_RNA']
del adata_2.obs['nCount_SCT']
del adata_2.obs['nFeature_SCT']
del adata_2.obs['seurat_clusters']
del adata_2.var['features']
del adata_2.uns['neighbors']
del adata_2.obsm['X_pca']
del adata_2.varm['PCs']
del adata_2.obsp['distances']

adata_2.write_h5ad("brapa_slim.h5ad", compression='gzip')

adata_3 = sc.read_h5ad('maizleaf.h5ad')
adata_3.raw = None

print(adata_3)

del adata_3.obs['Orig.ident']
del adata_3.obs['nCount_RNA']
del adata_3.obs['nFeature_RNA']
del adata_3.obs['Percent.mt']
del adata_3.obs['Seurat_clusters']
del adata_3.obs['Dataset']
del adata_3.obs['ACE']

adata_3.write_h5ad('maizleaf_slim.h5ad', compression='gzip')
"""
