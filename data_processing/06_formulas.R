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

options(future.globals.maxSize = 1024 * 1024 * 1000)
plan(multisession(workers = 3))
chemical_formulas <- canonical_compounds %>%
  mutate(
    data = future_map(
      data,
      function(df) {
        df %>%
          drop_na(inchi) %>%
          transmute(
            lspci_id,
            formula = inchi %>%
              str_split_fixed(fixed("/"), 3) %>%
              {.[, 2]},
            formula_vector = formula %>%
              str_match_all("([A-Z][a-z]?)([0-9]*)") %>%
              map(~set_names(replace_na(as.integer(.x[, 3]), 1L), .x[, 2]))
          )
      },
      .progress = TRUE
    )
  )

chemical_formulas_df <- chemical_formulas %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        mutate(
          is_organic = map_lgl(
            formula_vector,
            ~"C" %in% names(.x)
          )
        )
    )
  )

write_rds(
  chemical_formulas_df,
  file.path(dir_release, "chemical_formulas.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Compute compound formulas and find organic compounds",
  used = c(
    "syn20835543"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/06_formulas.R"
)

syn_folder <- Folder("properties", syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "chemical_formulas.rds")
) %>%
  synStoreMany(parent = syn_folder, activity = activity)

