library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

source(here("utils", "load_save.R"))

# Set directories, import files ------------------------------------------------
###############################################################################T

inputs <- list(
  references = c("reference_table")
) %>%
  pluck_inputs(syn_parent = syn_release)

# Make reference file ----------------------------------------------------------
###############################################################################T

lsp_references <- inputs[["references"]] %>%
  {synTableQuery(sprintf("SELECT * FROM %s", .))} %>%
  as.data.frame() %>%
  select(
    reference_id = ROW_ID, reference_type, reference_value, url
  ) %>%
  distinct() %>%
  setDT()

fwrite(
  lsp_references,
  file.path(dir_release, "references.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

hmsl_activity <- Activity(
  name = "Wrangle reference table.",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/02_reference_table.R"
)

syn_id_mapping <- synMkdir(syn_release, "references")

c(
  file.path(dir_release, "references.csv.gz")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = hmsl_activity)

