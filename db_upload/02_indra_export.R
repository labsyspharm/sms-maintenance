library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

db_version <- "morgan_normal"

# read tables ------------------------------------------------------------------
###############################################################################T

name_map <- syn("syn22035396") %>%
  read_rds() %>%
  filter(fp_name == db_version) %>%
  chuck("data", 1)

tas_annotated <- syn("syn21664452") %>%
  read_rds()

id_map <- syn("syn20830516") %>%
  read_rds() %>%
  filter(fp_name == db_version) %>%
  chuck("data", 1) %>%
  transmute(
    lspci_id = eq_class,
    id_source = if_else(str_starts(.[["id"]], fixed("CHEMBL")), "chembl", "hmsl"),
    compound_id = id
  ) %>%
  arrange(lspci_id)

tas_table <- tas_annotated %>%
  filter(fp_name == db_version) %>%
  chuck("data", 1) %>%
  drop_na(lspci_id, entrez_gene_id, tas)

write_csv(
  tas_table,
  file.path(dir_release, "indra_tas.csv.gz")
)

write_csv(
  id_map,
  file.path(dir_release, "indra_id_map.csv.gz")
)

write_csv(
  name_map,
  file.path(dir_release, "indra_name_map.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

export_tas_indra <- Activity(
  name = "Export TAS data to INDRA",
  used = c(
    "syn22035396",
    "syn21664452",
    "syn20830516"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/db_upload/02_indra_export.R"
)

syn_export <- Folder("export", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "indra_tas.csv.gz"),
  file.path(dir_release, "indra_id_map.csv.gz"),
  file.path(dir_release, "indra_name_map.csv.gz")
) %>%
  synStoreMany(parentId = syn_export, activity = export_tas_indra)


