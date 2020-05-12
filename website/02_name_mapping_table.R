library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(fst)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

syn_tables <- "syn20981852"

eq_classes <- syn("syn20830516") %>%
  read_rds()

all_names <- syn("syn22035396") %>%
  read_rds()

# Compound name mapping table --------------------------------------------------
###############################################################################T

gather_names <- function(names, ids, ...) {
  list(
    names,
    ids %>%
      transmute(
        lspci_id = eq_class,
        name = id,
        source = as.factor(paste0(source, "_id")),
        source_collapsed = as.factor("secondary")
      )
  ) %>%
    rbindlist(use.names = TRUE) %>%
    setkey(lspci_id, source_collapsed)
}

make_lspci_name_mapping <- function(d) {
  d %>%
    {
      .[, head(.SD, 1), by = lspci_id][, list(lspci_id, name)]
    } %>%
    unique()
}

name_mapping <- all_names %>%
  rename(names = data) %>%
  inner_join(
    rename(eq_classes, ids = data)
  ) %>%
  transmute(
    fp_name, fp_type,
    data = pmap(
      .,
      gather_names
    )
  ) %>%
  mutate(
    lspci_id_name_map = map(
      data,
      make_lspci_name_mapping
    )
  )

activity <- Activity(
  name = "Wrangle name mapping table",
  used = c(
    "syn20830516",
    "syn22035396"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/02_name_mapping_table.R"
)

pwalk(
  name_mapping,
  function(fp_name, data, lspci_id_name_map, ...) {
    write_csv(
      data,
      file.path(dir_release, paste0("all_names_lspci_id_map.csv.gz"))
    )
    write_fst(
      data,
      file.path(dir_release, paste0("all_names_lspci_id_map.fst"))
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
      file.path(dir_release, paste0("all_names_lspci_id_map.csv.gz")),
      file.path(dir_release, paste0("all_names_lspci_id_map.fst")),
      file.path(dir_release, paste0("lspci_id_name_map.fst"))
    ) %>%
      synStoreMany(syn_parent, activity = activity)
  }
)
