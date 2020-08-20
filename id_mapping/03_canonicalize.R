# Run on O2 cluster

library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(batchtools)
library(lspcheminf)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Loading files ----------------------------------------------------------------
###############################################################################T

inputs <- list(
  chembl_compounds = synPluck(syn_release, "raw_data", "chembl_compounds_raw.rds"),
  hmsl_compounds = synPluck(syn_release, "raw_data", "hmsl_compounds_raw.rds")
)

raw_compounds <- tribble(
  ~source, ~data,
  "chembl", inputs[["chembl_compounds"]] %>%
    syn() %>%
    read_rds() %>%
    distinct(id = chembl_id, inchi = standard_inchi),
  "hsml", inputs[["hmsl_compounds"]] %>%
    syn() %>%
    read_rds() %>%
    distinct(id = hms_id, inchi)
)

# Canonicalize compounds -------------------------------------------------------
###############################################################################T

all_inchis <- raw_compounds %>%
  unnest(data) %>%
  select(inchi) %>%
  distinct()

all_inchis_dfs <- all_inchis %>%
  chunk_df(150, seed = 1)

# Set up jobs


# Store to synapse -------------------------------------------------------------
###############################################################################T

cmpd_wrangling_activity <- Activity(
  name = "Map and wrangle compound IDs",
  used = c(
    "syn20692501",
    "syn20692443",
    "syn20692440",
    "syn20692514",
    "syn20830516",
    "syn21901782",
    "syn22080194",
    "syn21901780"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/04_compound_id_mapping.R"
)

syn_id_mapping <- Folder("id_mapping", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "canonical_table.rds"),
  file.path(dir_release, "lspci_id_compound_inchi_mapping.rds")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = cmpd_wrangling_activity)

