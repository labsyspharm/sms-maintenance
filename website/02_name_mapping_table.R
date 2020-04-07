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

NAME_SOURCE_PREFERENCE <- c("pref_name", "alt_name", "chembl_id", "hms_id")

gather_names <- function(d) {
  d %>%
    as.data.table() %>%
    {
      .[, list(lspci_id, pref_name, alt_name = alt_names, chembl_id, hms_id)]
    } %>%
    melt(
      id.vars = "lspci_id", variable.name = "source", value.name = "name", variable.factor = FALSE
    ) %>%
    {
      .[
        , list(name = unlist(name)), by = list(lspci_id, source)
      ][
        , source := factor(source, levels = NAME_SOURCE_PREFERENCE)
      ][
        !is.na(name)
      ][
        , name := if_else(str_length(name) > 1000, str_sub(str_split_fixed(name, "\\s", n = 2)[, 1], end = 1000), name)
      ]
    } %>%
    setkey(lspci_id, source)
}

make_lspci_name_mapping <- function(d) {
  d %>%
    {
      .[, head(.SD, 1), by = lspci_id][, list(lspci_id, name)]
    } %>%
    setkey(lspci_id)
}

make_name_lspci_mapping <- function(d) {
  d %>%
    # We want higher priority sources to come up first
    setkey(source, lspci_id) %>%
    {
      .[, list(lspci_id, name)]
    }
}

name_mapping <- canonical_table %>%
  mutate(
    data = map(
      data,
      gather_names
    ),
    name_lspci_id_map = map(data, make_name_lspci_mapping),
    lspci_id_name_map = map(data, make_lspci_name_mapping)
  )

activity <- Activity(
  name = "Wrangle name mapping table",
  used = c(
    "syn20835543"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/02_name_mapping_table.R"
)

pwalk(
  name_mapping,
  function(fp_name, data, name_lspci_id_map, lspci_id_name_map, ...) {
    write_csv(
      data,
      file.path(dir_release, paste0("lspci_id_name_map.csv.gz"))
    )
    write_fst(
      name_lspci_id_map,
      file.path(dir_release, paste0("name_lspci_id_map.fst"))
    )
    write_fst(
      lspci_id_name_map,
      file.path(dir_release, paste0("lspci_id_name_map.fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    c(
      # file.path(dir_release, paste0("lspci_id_name_map.csv.gz")),
      file.path(dir_release, paste0("name_lspci_id_map.fst"))
      # file.path(dir_release, paste0("lspci_id_name_map.fst"))
    ) %>%
      synStoreMany(syn_parent, activity = activity)
  }
)
