library(tidyverse)
library(data.table)
library(fst)
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

rscores <- syn("syn21590477") %>%
  read_rds()

tas <- syn("syn21664452") %>%
  read_rds()

# Phenotypic score and TAS tables ----------------------------------------------
###############################################################################T

# Remove any compounds with less than 2 assays
rscores_selected <- rscores %>%
  mutate(
    data = map(
      data,
      function(d) {
        d %>%
          mutate_at(vars(lspci_id, assay_id), as.integer) %>%
          as.data.table() %>%
          {
            .[
              is.finite(rscore_tr)
            ][
              , .(lspci_id, assay_id, rscore_tr)
            ][
              , if (.N >= 2) .SD, by = lspci_id
            ]
          } %>%
          # Set index on the assay and lspci id for efficient joining
          setkey(assay_id, lspci_id)
      }
    )
  )

activity <- Activity(
  name = "Wrangle phenotypic scores for website",
  used = c(
    "syn21590477"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/01_shiny_pfp_tas.R"
)

pwalk(
  rscores_selected,
  function(fp_name, data, ...) {
    write_fst(
      data,
      file.path(dir_release, paste0("shiny_phenotypic_rscore_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_phenotypic_rscore_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_phenotypic_rscore.fst"
    ) %>%
      synStore(activity = activity)
  }
)


# Remove any compounds with less than 2 assays
tas_selected <- tas %>%
  mutate(
    data = map(
      data,
      function(d) {
        d %>%
          mutate_at(vars(lspci_id, entrez_gene_id, tas), as.integer) %>%
          select(
            lspci_id, gene_id = entrez_gene_id, tas, source, measurement, unit, references
          ) %>%
          as.data.table() %>%
          {
            .[
              , if (.N >= 2) .SD, by = lspci_id
            ]
          } %>%
          # Set index on the entrez_gene_id and lspci id for efficient joining
          setkey(gene_id, lspci_id)
      }
    )
  )

activity <- Activity(
  name = "Wrangle TAS scores for website",
  used = c(
    "syn21664452"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/01_shiny_pfp_tas.R"
)

pwalk(
  tas_selected,
  function(fp_name, data, ...) {
    write_fst(
      data,
      file.path(dir_release, paste0("shiny_tas_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_tas_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_tas.fst"
    ) %>%
      synStore(activity = activity)
  }
)
