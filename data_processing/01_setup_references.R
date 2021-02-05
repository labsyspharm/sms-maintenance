library(tidyverse)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")


# Set directories, import files ------------------------------------------------
###############################################################################T

cols <- list(
  Column(name = "reference_type", columnType = "STRING", maximumSize = 20),
  Column(name = "reference_value", columnType = "STRING", maximumSize = 100),
  Column(name = "chembl_id_doc", columnType = "STRING", maximumSize = 20),
  Column(name = "url", columnType = "LARGETEXT")
)

schema <- Schema(name = "reference_table", columns = cols, parent = syn_release)

table <- Table(
  schema,
  data.frame(
    reference_type = character(),
    reference_value = character(),
    url = character()
  )
)

synStore(table)
