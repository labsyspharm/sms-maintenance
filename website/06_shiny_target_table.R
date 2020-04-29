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

target_dict <- syn("syn20693721") %>%
  read_csv()

# Target table -----------------------------------------------------------------
###############################################################################T

target_table <- target_dict %>%
  drop_na(entrez_gene_id) %>%
  # Replace nonsense symbol "-" with entrez id
  mutate(symbol = recode(symbol, `-` = paste0("uniprot_", uniprot_id))) %>%
  distinct(tax_id, gene_id = entrez_gene_id, symbol) %>%
  as.data.table()

setkey(target_table, gene_id)

activity <- Activity(
  name = "Wrangle target table",
  used = c(
    "syn20693721"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/06_shiny_target_table.R"
)

walk(
  synChildren(syn_tables) %>%
    names(),
  function(fp_name) {
    write_fst(
      target_table,
      file.path(dir_release, paste0("shiny_targets_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_targets_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_targets.fst"
    ) %>%
      synStore(activity = activity)
  }
)
