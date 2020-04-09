## this script makes a table of compoundfingerprints

library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Prepare fingerprint database -------------------------------------------------
###############################################################################T


all_fp <- syn("syn20692501") %>%
  read_rds()

canonical_compounds <- syn("syn20835543") %>%
  read_rds()

all_fp_flat <- all_fp %>%
  unnest(data)

canonical_fp <- canonical_compounds %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        select(lspci_id, chembl_id, hms_id) %>%
        gather("source", "id", hms_id, chembl_id) %>%
        # Prefer Chembl if present
        mutate(source = factor(source, levels = c("chembl_id", "hms_id"))) %>%
        group_by(lspci_id) %>%
        arrange(source, .by_group = TRUE) %>%
        slice(1) %>%
        ungroup() %>%
        left_join(all_fp_flat, by = "id") %>%
        select(lspci_id, fp_name, fp_type, fingerprint)
    )
  )

write_rds(
  canonical_fp,
  file.path(dir_release, "canonical_fingerprints.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

fp_table_activity <- Activity(
  name = "Create table of molecular fingerprints",
  used = c(
    "syn20692501",
    "syn20835543"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/07_fingerprint_table.R"
)

fp_folder <- Folder("fingerprints", syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "canonical_fingerprints.rds")
) %>%
  synStoreMany(parent = fp_folder, activity = fp_table_activity)

