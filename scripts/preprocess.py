import sys, os
import scanpy as sc
import tempfile

os.chdir("/well/ocallaghan/users/nem384/software/pyflow_seurat_parameter")

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

adata = sc.read_h5ad(snakemake.input[0])
k = int(snakemake.wildcards.k)
pc_use = int(snakemake.wildcards.pc)
resolution = float(snakemake.wildcards.resolution)

sc.pp.neighbors(adata, n_neighbors=k, n_pcs=pc_use, use_rep = "")
sc.tl.leiden(adata, resolution=resolution, key_added='leiden')
adata.obs['leiden'] = adata.obs['leiden'].astype(str)
adata.obs['pc.use'] = pc_use
adata = adata[:, 0].copy()

adata.write_h5ad(f"full_sample_preprocess/full_sample_k_{k}_resolution_{resolution}_PC_{pc_use}.h5ad")

sys.exit(0)
