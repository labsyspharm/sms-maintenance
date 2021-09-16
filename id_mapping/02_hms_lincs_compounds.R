library(tidyverse)
library(httr)
library(furrr)
library(jsonlite)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)
library(data.table)
library(qs)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v29"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Retrieving list of LINCS compounds -------------------------------------------
###############################################################################T

# Attempt to get hms lincs ID map automagically from the HMS LINCS reagent tracker
# Set username and password in ~/.Renviron
# ECOMMONS_USERNAME=xxx
# ECOMMONS_PASSWORD=xxx

login_token <- POST(
  "https://reagenttracker.hms.harvard.edu/api/v0/login",
  body = list(
    username = Sys.getenv("ECOMMONS_USERNAME"),
    password = Sys.getenv("ECOMMONS_PASSWORD")
  ),
  encode = "json"
)

rt_response = GET(
  "https://reagenttracker.hms.harvard.edu/api/v0/search?q=",
  accept_json()
)

rt_df <- rt_response %>%
  content("parsed", type = "application/json", encoding = "UTF-8", simplifyVector = TRUE) %>%
  pluck("canonicals") %>%
  as_tibble() %>%
  filter(type == "small_molecule", name != "DEPRECATED") %>%
  transmute(
    hms_id = lincs_id,
    name,
    # Remove spurious empty lists
    alternate_names = map(alternate_names, ~if (length(.x) == 0) NULL else .x),
    smiles,
    inchi,
    inchi_key,
    chembl_id,
    n_batches = map_int(batches, length)
  )

# 2021-09-09
# SMILES for belvarafenib is wrong
setDT(rt_df)[
  name == "Belvarafenib",
  smiles := "CC1=C(C2=C(C=C1)C(=NC=C2)NC3=C(C(=CC=C3)Cl)F)NC(=O)C4=CSC5=C4N=CN=C5N"
]

# Convert smiles to inchi, where no inchi is provided
rt_df_inchi <- rt_df %>%
  mutate(
    inchi = map2_chr(
      inchi, smiles,
      # Only convert if inchi is unknown and smiles is known
      # otherwise use known inchi
      ~if (is.na(.x) && !is.na(.y))
        convert_compound_descriptor(compounds(.y, descriptor = "smiles"), target_descriptor = "inchi")[["compounds"]]
      else .x
    )
  ) %>%
  # Some inchis have newline characters at the end of line, remove them
  mutate(inchi = trimws(inchi))

qsave(
  rt_df_inchi,
  file.path(dir_release, "hmsl_compounds_raw.qs")
)
# rt_df_inchi <- read_rds(file.path(dir_release, "hmsl_compounds_raw.rds"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

fetch_hmsl_activity <- Activity(
  name = "Fetch HITS compound data and canonicalize",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/02_hms_lincs_compounds.R"
)

c(
  file.path(dir_release, "hmsl_compounds_raw.qs")
) %>%
  synStoreMany(parent = synMkdir(syn_release, "raw_data"), activity = fetch_hmsl_activity, forceVersion = FALSE)
