library(tidyverse)
library(here)
library(fst)
library(data.table)
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

selectivity_classes <- syn("syn20836653") %>%
  read_rds()

# Selectivity table ------------------------------------------------------------
###############################################################################T

selectivity_table <- selectivity_classes %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        distinct(
          gene_id,
          lspci_id,
          chembl_id,
          hms_id,
          toolscore = tool_score,
          ontarget_IC50_Q1,
          n_measurement = ontarget_IC50_N,
          source = recode(
            selectivity_class,
            most_selective = "Most selective",
            semi_selective = "Semi-selective",
            poly_selective = "Poly-selective",
            unknown_selective = "Unknown",
            other_selective = "Other"
          ),
          symbol = gene_symbol
        ) %>%
        as.data.table()
    )
  )

activity <- Activity(
  name = "Wrangle selectivity table",
  used = c(
    "syn20836653"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/03_shiny_selectivity.R"
)

pwalk(
  selectivity_table,
  function(fp_name, data, ...) {
    setkey(data, lspci_id, gene_id)
    write_fst(
      data,
      file.path(dir_release, paste0("shiny_selectivity_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_selectivity_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_selectivity.fst"
    ) %>%
      synStore(activity = activity)
  }
)

