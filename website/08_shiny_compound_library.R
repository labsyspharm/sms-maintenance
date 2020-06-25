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

optimal_library <- syn("syn22153734") %>%
  read_rds()

# Commercially available compounds table ---------------------------------------
###############################################################################T

optimal_library_tables <- optimal_library %>%
  mutate(
    data = map(
      data,
      function(df) {
        df %>%
          mutate(reason_included = factor(reason_included, levels = c("selectivity", "clinical"))) %>%
          select(gene_id, lspci_id, reason_included, rank) %>%
          setDT(key = c("gene_id", "rank"))
      }
    )
  )

activity <- Activity(
  name = "Convert optimal compound table to fst",
  used = c(
    "syn22153734"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/08_shiny_compound_library.R"
)

pwalk(
  optimal_library_tables,
  function(data, fp_name, ...) {
    write_fst(
      data,
      file.path(dir_release, paste0("shiny_optimal_compound_table_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_optimal_compound_table_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_optimal_compound_table.fst"
    ) %>%
      synStore(activity = activity)
  }
)
