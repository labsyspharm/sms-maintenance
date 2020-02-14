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

tas <- syn("syn20982111") %>%
  read_csv() %>%
  select(lspci_id, entrez_gene_id, tas) %>%
  mutate_all(as.integer)

# TAS similarity calculation ---------------------------------------------------
###############################################################################T

# Remove any compounds with less than 6 assays
tas_selected <- tas %>%
  as.data.table() %>%
  .[, if (.N >= 6) .SD, by = lspci_id]
# Set index on the entrez_gene_id and lspci id for efficient joining
setkey(tas_selected, entrez_gene_id, lspci_id)

tas_weighted_jaccard <- function(query_id, tas_data) {
  tic()
  query_tas <- tas_data[lspci_id == query_id, .(entrez_gene_id, tas)]
  res <- tas_data[
    query_tas,
    on = "entrez_gene_id",
    nomatch = NULL
  ][
    ,
    mask := tas < 10 | i.tas < 10
  ][
    ,
    if (sum(mask) >= 6) .(
      "tas_similarity" = sum(pmin(tas[mask], i.tas[mask])) / sum(pmax(tas[mask], i.tas[mask])),
      "n" = sum(mask),
      "n_prior" = .N
    ),
    by = "lspci_id"
  ]
  toc()
  res
}

x <- tas_weighted_jaccard(1379, tas_selected)


