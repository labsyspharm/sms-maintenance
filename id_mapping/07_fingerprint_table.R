## this script makes a table of compoundfingerprints

library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(data.table)

synLogin()
syn <- synDownloader(here("tempdl"))

source(here("utils", "load_save.R"))

# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")


inputs <- list(
  all_fp = c("fingerprints", "all_compounds_fingerprints.csv.gz"),
  canonical_members = c("compounds_processed", "lspci_id_canonical_members.csv.gz")
) %>%
  pluck_inputs(syn_parent = syn_release)

input_data <- inputs %>%
  load_input_data(syn = syn)

# Prepare fingerprint database -------------------------------------------------
###############################################################################T

fingerprint_table <- input_data[["canonical_members"]][
  ,
  .(
    lspci_id,
    inchi_id
  )
] %>%
  merge(
    input_data[["all_fp"]],
    all = FALSE,
    by = "inchi_id"
  ) %>%
  setkey(lspci_id, fp_name)


fwrite(
  fingerprint_table,
  file.path(dir_release, "lspci_id_fingerprints.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

fp_table_activity <- Activity(
  name = "Create table of molecular fingerprints",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/07_fingerprint_table.R"
)

fp_folder <- synMkdir(syn_release, "compounds_processed")

c(
  file.path(dir_release, "lspci_id_fingerprints.csv.gz")
) %>%
  synStoreMany(parent = fp_folder, activity = fp_table_activity)
