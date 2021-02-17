library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(qs)

synLogin()
syn <- synDownloader(here("tempdl"))

source(here("utils", "load_save.R"))

# Set directories, import files ------------------------------------------------
###############################################################################T

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

inputs <- list(
  canonical_compounds = c("compounds_processed", "compound_dictionary.csv.gz")
)  %>%
  pluck_inputs(syn_parent = syn_release)

input_data <- inputs %>%
  load_input_data(syn = syn)

# Checking for organic molecules -----------------------------------------------
###############################################################################T

# options(future.globals.maxSize = 1024 * 1024 * 1000)

chemical_formulas_raw <- input_data[["canonical_compounds"]] %>%
  drop_na(canonical_inchi) %>%
  transmute(
    lspci_id,
    formula = canonical_inchi %>%
      str_split_fixed(fixed("/"), 3) %>%
      {.[, 2]},
    formula_vector = formula %>%
      str_match_all("([A-Z][a-z]?)([0-9]*)") %>%
      map(~set_names(replace_na(as.integer(.x[, 3]), 1L), .x[, 2]))
  )

chemical_formulas_df <- chemical_formulas_raw %>%
  mutate(
    is_organic = map_lgl(
      formula_vector,
      ~"C" %in% names(.x)
    )
  )

qsave(
  chemical_formulas_df,
  file.path(dir_release, "chemical_formulas.qs"),
  preset = "fast"
)

fwrite(
  chemical_formulas_df %>%
    select(-formula_vector),
  file.path(dir_release, "chemical_formulas.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Compute compound formulas and find organic compounds",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/06_formulas.R"
)

syn_folder <- synMkdir(syn_release, "properties")

c(
  file.path(dir_release, "chemical_formulas.qs"),
  file.path(dir_release, "chemical_formulas.csv.gz")
) %>%
  synStoreMany(parent = syn_folder, activity = activity)

