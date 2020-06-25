library(tidyverse)
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

optimal_library <- syn("syn22153734") %>%
  read_rds()

kinases <- syn("syn12617467") %>%
  read_csv(col_types = "iccllllic")

target_dict <- syn("syn20693721") %>%
  read_csv() %>%
  distinct(gene_id = entrez_gene_id, symbol = entrez_symbol)

# Generating gene sets ---------------------------------------------------------
###############################################################################T

full_liganded_genome <- optimal_library %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        distinct(gene_id) %>%
        inner_join(target_dict, by = "gene_id") %>%
        drop_na()
    )
  )

liganded_kinome <- full_liganded_genome %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        filter(gene_id %in% kinases[["gene_id"]])
    )
  )

gene_lists <- list(
  "Full_LigandedGenome" = full_liganded_genome,
  "Kinome" = liganded_kinome
) %>%
  bind_rows(.id = "gene_list")

activity <- Activity(
  name = "Generate liganded gene lists",
  used = c(
    "syn22153734",
    "syn12617467",
    "syn20693721"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/09_shiny_gene_lists.R"
)

pwalk(
  gene_lists,
  function(data, fp_name, gene_list, ...) {
    write_lines(
      data[["symbol"]] %>%
        unique(),
      file.path(dir_release, paste0(gene_list, ".txt"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0(gene_list, ".txt")),
      parent = syn_parent
    ) %>%
      synStore(activity = activity)
  }
)
