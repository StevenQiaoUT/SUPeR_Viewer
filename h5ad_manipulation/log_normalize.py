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
