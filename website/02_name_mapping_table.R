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

all_names <- syn("syn22035396") %>%
  read_rds()

id_maps <- syn("syn20830516") %>%
  read_rds()

# Compound name mapping table --------------------------------------------------
###############################################################################T

make_lspci_name_mapping <- function(d) {
  d %>%
    {
      .[, head(.SD, 1), by = lspci_id][, list(lspci_id, name)]
    } %>%
    unique()
}

name_levels <- c("vendor", "chembl_pref", "hmsl_pref", "chembl_alt", "hmsl_alt", "chembl_id", "hmsl_id")

name_mapping <- all_names %>%
  inner_join(
    rename(id_maps, ids = data)
  ) %>%
  mutate(
    data = map2(
      data, ids,
      ~.x %>%
        # Remove duplicate names, only keep top source
        group_by(lspci_id, name) %>%
        slice(1) %>%
        ungroup() %>%
        mutate_at(vars(source), factor, levels = name_levels) %>%
        bind_rows(
          .y %>%
            transmute(
              lspci_id = eq_class,
              name = id,
              source = if_else(str_starts(id, fixed("H")), "hmsl_id", "chembl_id") %>%
                factor(
                  levels = name_levels
                ),
              source_collapsed = factor("secondary", levels = c("primary", "secondary"))
            )
        ) %>%
        setDT() %>%
        {
          .[
            # Prefer name with fewer spaces (no salts) and
            # shorter name if multiple sources with same priority are available
            order(lspci_id, source_collapsed, str_count(name, fixed(" ")), str_length(name))
          ][
            ,
            # For selectize.js, the values (lspci_ids) must be unique, duplicates
            # are silently discarded. Making a unique version by affixing -1, -2, ...
            lspci_id_unique := paste(lspci_id, 1:.N, sep = "-"),
            by = "lspci_id"
          ]
        } %>%
        setkey(lspci_id)
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
