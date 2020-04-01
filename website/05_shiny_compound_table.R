library(tidyverse)
library(here)
library(data.table)
library(fst)
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

compounds <- syn("syn20835543") %>%
  read_rds()

clinical_info <- syn("syn21064122") %>%
  read_rds()

# Compound table ---------------------------------------------------------------
###############################################################################T

compound_tables <- compounds %>%
  inner_join(
    clinical_info %>%
      rename(ci = data),
    by = c("fp_name", "fp_type")
  ) %>%
  mutate(
    data = map2(
      data, ci,
      ~.x %>%
        arrange(lspci_id) %>%
        left_join(
          .y %>%
            select(lspci_id, max_phase),
          by = "lspci_id"
        ) %>%
        distinct(
          lspci_id,
          chembl_id,
          hms_id,
          pref_name,
          alt_names = map_chr(
            alt_names,
            ~if (is.null(.x)) NA_character_ else paste(.x, collapse = "; ")
          ),
          max_phase
        ) %>%
        as.data.table()
    )
  )

activity <- Activity(
  name = "Convert compound table to fst",
  used = c(
    "syn20835543"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/05_shiny_compound_table.R"
)

pwalk(
  compound_tables,
  function(data, fp_name, ...) {
    setkey(data, lspci_id)
    write_fst(
      data,
      file.path(dir_release, paste0("shiny_compounds_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_compounds_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_compounds.fst"
    ) %>%
      synStore(activity = activity)
  }
)
