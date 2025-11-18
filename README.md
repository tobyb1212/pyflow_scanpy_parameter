# Snakemake workflow for subsampling and repeat clustering

A snakemake pipeline to scatter and gather cluster identities by subsampling the cells and repeat for 
multiple times. This is useful for evaluating the cluster stability using different parameters. All credit for the development of the original code
goes to @crazyhottommy.

For now, three parameters are tested.

* `k.param` and `components` for `sc.pp.neighbors`
* `resolution` for `sc.tl.leiden` 

on `odyssey` cluster(SLURM):

```bash
ssh odyssey

## start a screen session
screen

git clone https://github.com/tobyb1212/pyflow_scanpy_parameter

conda install -c conda-forge mamba
mamba create -c conda-forge -c bioconda -n snakemake snakemake

source activate snakemake

# R3.5.1, make sure you load R after source activate conda environment
module load R

#hdf5, seurat needs a more recent hdf5 to be able to install
module load hdf5

R
>install.package("Seurat")
>install.package("schard")
q()

python
>pip install scanpy numpy pandas rpy2

```

Users are required to use an R-reticulate environment to enable python packages to be called from R. This envionment must contain the necessary packages and the environment location should be referenced in scripts/gather_fullsample.R.

```
R
> if (!requireNamespace("reticulate", quietly = TRUE)) {
  install.packages("reticulate")
}

> library(reticulate)

# create a new Python environment named "r-reticulate"
> reticulate::virtualenv_create("r-reticulate")

# install Python packages inside that environment
> reticulate::virtualenv_install("r-reticulate", 
  packages = c("scanpy", "numpy", "pandas", "rpy2")
)

# use this environment in R
> reticulate::use_virtualenv("r-reticulate", required = TRUE)

```

copy your .h5ad object into the `pyflow_scanpy_parameter` folder and 

open the `config.yaml` file to edit some configurations.



```bash
# dry run
snakemake -np 

# if on bioinfo1 or bioinfo2 (there are 64 cores avaiable on each node)
snakemake -j 40

# if submitting job to queue (for slurm only)

./pyflow-scBoot.sh
```

The `rds` output from this workflow can be directly used by the [scclusteval R package](https://github.com/crazyhottommy/scclusteval) for downstream analysis, with the following code.

```R
subsample_idents <- subsample_idents %>%
  mutate(
    original_ident = map(original_ident, ~ {
      names(.x) <- paste0("cell_", seq_along(.x))
      .x
    }),
    recluster_ident = map2(recluster_ident, original_ident, ~ {
      names(.x) <- names(.y)
      .x
    })
  )

subsample_idents_list <- subsample_idents %>% 
  group_by(pc, resolution, k.param) %>% 
  nest()

# Assign Stable Clusters

stable_clusters<- subsample_idents_list %>%
  mutate(stable_cluster = map(data, ~ AssignStableCluster(.x$original_ident,
                                                          .x$recluster_ident,
                                                          jaccard_cutoff = 0.6,
                                                          method = "jaccard_percent", 
                                                          percent_cutoff = 0.6)))

# Parameter Scatter Plot

ParameterSetScatterPlot(stable_clusters = stable_clusters,
                        fullsample_idents = fullsample_idents,
                        x_var = "k.param",
                        y_var = "number",
                        facet_rows = "resolution",
                        facet_cols = "pc") + theme_bw()
ggsave(plot = last_plot(), filename = "parameter_scatter.jpeg", dpi = 300, width = 10, height = 10)

# Calculate Percentage of Cells in Stable Clusters

ParameterSetScatterPlot(stable_clusters = stable_clusters,
                        fullsample_idents = fullsample_idents,
                        x_var = "k.param",
                        y_var = "percentage",
                        facet_rows = "resolution",
                        facet_cols = "pc") +
  ggtitle("percentage of cells in stable clusters") + theme_bw()
ggsave(plot = last_plot(), filename = "cells_in_cluster_scatter.jpeg", dpi = 300, width = 10, height = 10)

```
If there are any questions regarding the updated scanpy adaptation of this workflow, please contact the author of this adapted package. 

