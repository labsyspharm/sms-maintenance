library(tidyverse)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)
library(furrr)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Checking for organic molecules -----------------------------------------------
###############################################################################T

canonical_compounds <- syn("syn20835543") %>%
  read_rds()

plan(multisession(workers = 6))
organic <- canonical_compounds %>%
  mutate(
    data = map(
      data,
      function(df) {
        df %>%
          drop_na(inchi) %>%
          chunk_df(12) %>%
          map(~set_names(.x[["inchi"]], .x[["lspci_id"]])) %>%
          future_map(is_organic, .progress = TRUE)
      }
    )
  )

organic_df <- organic %>%
  mutate(
    data = map(
      data,
      bind_rows
    ) %>%
      map(transmute, lspci_id = as.integer(compound), is_organic) %>%
      map(arrange, lspci_id)
  )

write_rds(
  organic_df,
  file.path(dir_release, "organic_compounds.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Find organic compounds",
  used = c(
    "syn20835543"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/06_organic_compounds.R"
)

syn_folder <- Folder("properties", syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "organic_compounds.rds")
) %>%
  synStoreMany(parent = syn_folder, activity = activity)

