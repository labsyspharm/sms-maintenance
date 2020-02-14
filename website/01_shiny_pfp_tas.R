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

tas <- syn("syn20830939") %>%
  read_rds()

# Phenotypic score and TAS tables ----------------------------------------------
###############################################################################T

# Remove any compounds with less than 6 assays
rscores_selected <- rscores %>%
  mutate(
    data = map(
      data,
      function(d) {
        res <- d %>%
          as.data.table() %>%
          .[is.finite(rscore_tr), .(lspci_id, assay_id, rscore_tr)] %>%
          .[, if (.N >= 6) .SD, by = lspci_id]
        # Set index on the assay and lspci id for efficient joining
        setkey(res, assay_id, lspci_id)
        res
      }
    )
  )

activity <- Activity(
  name = "Wrangle phenotypic scores for website",
  used = c(
    "syn21590477"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/11_prepare_website_files.R"
)

pwalk(
  rscores_selected,
  function(fp_name, data, ...) {
    write_fst(
      data,
      file.path(dir_release, paste0("phenotypic_rscore_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("phenotypic_rscore_", fp_name, ".fst")),
      parent = syn_parent,
      name = "phenotypic_rscore.fst"
    ) %>%
      synStore(activity = activity)
  }
)


# Remove any compounds with less than 6 assays
tas_selected <- tas %>%
  mutate(
    data = map(
      data,
      function(d) {
        res <- d %>%
          as.data.table() %>%
          .[
            ,
            .(lspci_id, entrez_gene_id, tas)
          ] %>%
          .[
            ,
            if (.N >= 6) .SD,
            by = lspci_id
          ]
        # Set index on the entrez_gene_id and lspci id for efficient joining
        setkey(res, entrez_gene_id, lspci_id)
        res
      }
    )
  )

activity <- Activity(
  name = "Wrangle TAS scores for website",
  used = c(
    "syn20830939"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/11_prepare_website_files.R"
)

pwalk(
  tas_selected,
  function(fp_name, data, ...) {
    write_fst(
      data,
      file.path(dir_release, paste0("tas_score_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("tas_score_", fp_name, ".fst")),
      parent = syn_parent,
      name = "tas_score.fst"
    ) %>%
      synStore(activity = activity)
  }
)
