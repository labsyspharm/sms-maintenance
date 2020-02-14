library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

syn_tables <- "syn20981852"

canonical_table <- syn("syn20835543") %>%
  read_rds()

# Name mapping table -----------------------------------------------------------
###############################################################################T

gather_names <- function(d) {
  d %>%
    as.data.table() %>%
    .[
      ,
      .(
        lspci_id,
        name = map2(
          pref_name, alt_names,
          ~unique(c(if (!is.na(.x)) .x, if (!is.null(.y)) .y))
        )
      )
      ] %>%
    .[!map_lgl(name, is.null)] %>%
    .[, .(name = unlist(name)), by = lspci_id] %>%
    as_tibble()
}

name_mapping <- canonical_table %>%
  mutate(
    data = map(
      data,
      gather_names
    )
  )

activity <- Activity(
  name = "Wrangle name mapping table",
  used = c(
    "syn20835543"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/11_prepare_website_files.R"
)

pwalk(
  name_mapping,
  function(fp_name, data, ...) {
    write_csv(
      data,
      file.path(dir_release, paste0("lspci_id_name_mapping_", fp_name, ".csv.gz"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("lspci_id_name_mapping_", fp_name, ".csv.gz")),
      parent = syn_parent,
      name = "lspci_id_name_map.csv.gz"
    ) %>%
      synStore(activity = activity)
  }
)

c(
  file.path(dir_release, "all_compounds_similarity.csv.gz"),
  file.path(dir_release, "compound_masses.csv.gz")
) %>%
  synStoreMany(parentId = syn_id_mapping, activity = activity)
