library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(tictoc)

synLogin()
syn <- synDownloader(here("tempdl"))

# set directories & import files -----------------------------------------------
###############################################################################T

release <- "chembl_v25"
dir_release <- here(release)

rscores <- read_rds(file.path(dir_release, "pheno_data_rscores.rds"))

rs <- rscores$data[[2]]

# PFP similarity calculation ---------------------------------------------------
###############################################################################T

# Remove any compounds with less than 6 assays
rs_selected <- rs %>%
  as.data.table() %>%
  .[is.finite(rscore_tr), .(lspci_id, assay_id, rscore_tr)] %>%
  .[, if (.N >= 6) .SD, by = lspci_id]
# Set index on the assay and lspci id for efficient joining
setkey(rs_selected, assay_id, lspci_id)

pfp_correlation <- function(query_id, pfp_data) {
  tic()
  query_pfps <- pfp_data[lspci_id == query_id]
  res <- pfp_data[
    query_pfps,
    on = "assay_id",
    nomatch = NULL
  ][
    ,
    mask := abs(rscore_tr) >= 2.5 | abs(i.rscore_tr) >= 2.5
  ][
    ,
    if(sum(mask) >= 6) .(
      "correlation" = cor(rscore_tr, i.rscore_tr),
      "n_prior" = .N,
      "n" = sum(mask)
    ),
    by = lspci_id
  ]
  toc()
  res
}

x <- pfp_correlation(42, rs_selected)

