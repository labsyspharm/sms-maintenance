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

commercial_info <- syn("syn21049601") %>%
  read_rds()

# Commercially available compounds table ---------------------------------------
###############################################################################T

commercial_tables <- commercial_info %>%
  mutate(
    data = map(
      data,
      function(df) {
        df %>%
          distinct(lspci_id, vendor, vendor_id) %>%
          setDT(key = "lspci_id")
      }
    )
  )

activity <- Activity(
  name = "Convert commercial info table to fst",
  used = c(
    "syn21049601"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/07_shiny_commercial_availability.R"
)

pwalk(
  commercial_tables,
  function(data, fp_name, ...) {
    write_fst(
      data,
      file.path(dir_release, paste0("shiny_commercial_info_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_commercial_info_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_commercial_info.fst"
    ) %>%
      synStore(activity = activity)
  }
)
