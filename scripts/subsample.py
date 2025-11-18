import sys, os
import scanpy as sc
import tempfile
import numpy as np
import pandas as pd
import subprocess
from pathlib import Path
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri

log = open(snakemake.log[0], "w")
sys.stdout = log
sys.stderr = log

adata = sc.read_h5ad(snakemake.input[0])
k = int(snakemake.wildcards.k)
pc_use = int(snakemake.wildcards.pc)
resolution = float(snakemake.wildcards.resolution)
run_id = snakemake.wildcards.run_id
rate = float(snakemake.params.rate)
seed = None

ncells = adata.n_obs
ncells_subsample = int(round(ncells * float(rate)))
rng = np.random.RandomState(None if seed is None else int(seed))
selected_cells = rng.choice(adata.obs_names, size=ncells_subsample, replace=False)
adata = adata[:, 0].copy()
idx = adata.obs_names.get_indexer(selected_cells) 
subset_adata = adata[idx, :].copy()
original_ident_list = subset_adata.obs["leiden"].astype(str).tolist()

sc.pp.neighbors(subset_adata, n_neighbors=k, n_pcs=pc_use, use_rep = "X_scANVI")
sc.tl.leiden(subset_adata, resolution=resolution, key_added='recluster_leiden')
recluster_ident_list = subset_adata.obs["recluster_leiden"].astype(str).tolist()

res_df = pd.DataFrame({
    "pc": [pc_use],
    "resolution": [resolution],
    "k.param": [k],
    "original_ident": [original_ident_list],
    "recluster_ident": [recluster_ident_list],
    "round": [run_id]
})

outfile = Path(f"subsample/subsample_k_{k}_resolution_{resolution}_PC_{pc_use}_round_{run_id}.rds")
with (ro.default_converter + pandas2ri.converter).context():
    r_obj = ro.conversion.py2rpy(res_df)
    ro.r["saveRDS"](r_obj, str(outfile))
print(f"Saved RDS via rpy2: {outfile}")

sys.exit(0)
