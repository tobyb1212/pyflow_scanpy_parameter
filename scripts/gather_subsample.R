.libPaths(c("/well/ocallaghan/users/nem384/R/4.3/skylake", .libPaths()))
library(tidyverse)

log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

rdss<- snakemake@input[["rds"]]

get_df <- function(rds) {
  df <- readRDS(rds)
  
  # Convert numeric/integer columns to character
  df$pc <- as.character(df$pc)
  df$resolution <- as.character(df$resolution)
  df$k.param <- as.character(df$k.param)
  df$round <- as.character(df$round)
  
  # Helper to parse the string columns into factor vectors
  parse <- function(s) {
    if (is.na(s) || grepl("__truncated__", s)) stop("string truncated or missing")
    s <- gsub("^\\[|\\]$", "", s)
    s <- gsub("['\"]", "", s)
    v <- trimws(strsplit(s, ",")[[1]])
    factor(v)
  }
  
  # Convert stringified list-columns into actual list-columns of factors
  df$original_ident  <- lapply(df$original_ident,  parse)
  df$recluster_ident <- lapply(df$recluster_ident, parse)
  
  # Return as tibble to match test_1 format
  tibble::as_tibble(df)
}

dat.list<- lapply(rdss, get_df)

gather_idents<- do.call(bind_rows, dat.list)
saveRDS(gather_idents, file = "gather_subsample.rds")

