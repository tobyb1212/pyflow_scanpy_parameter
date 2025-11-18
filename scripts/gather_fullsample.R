library(tidyverse)
library(Seurat)
Sys.setenv(RETICULATE_PYTHON = "")
library(reticulate)
library(schard)
library(anndata)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

rdss<- snakemake@input[["rds"]]

get_idents<- function(rds){
        x<- schard::h5ad2seurat(rds)
        k<- gsub("full_sample_k_([0-9]+)_resolution_([0-9\\.]+)_PC_([0-9]+).h5ad", "\\1", basename(rds))
        resolution<- gsub("full_sample_k_([0-9]+)_resolution_([0-9\\.]+)_PC_([0-9]+).h5ad", "\\2", basename(rds))
        pc.use<- gsub("full_sample_k_([0-9]+)_resolution_([0-9\\.]+)_PC_([0-9]+).h5ad", "\\3", basename(rds))
        df<- tibble::tibble(pc = pc.use, resolution = resolution, k_param = k, original_ident_full = list(x$leiden))
        return(df)
}

dat.list<- lapply(rdss, get_idents)

gather_idents<- do.call(bind_rows, dat.list)
saveRDS(gather_idents, file = snakemake@output[[1]])

