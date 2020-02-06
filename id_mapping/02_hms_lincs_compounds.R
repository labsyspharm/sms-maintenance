library(tidyverse)
library(httr)
library(furrr)
library(jsonlite)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
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

# Convert smiles to inchi, where no inchi is provided
rt_df_inchi <- rt_df %>%
  mutate(
    inchi = map2_chr(
      inchi, smiles,
      # Only convert if inchi is unknown and smiles is known
      # otherwise use known inchi
      ~if (is.na(.x) && !is.na(.y))
        convert_compound_identifier(.y, identifier = "smiles", target_identifier = "inchi")[["compounds"]]
      else .x
    )
  ) %>%
  # Some inchis have newline characters at the end of line, remove them
  mutate(inchi = trimws(inchi))

write_rds(
  rt_df_inchi,
  file.path(dir_release, "hmsl_compounds_raw.rds")
)
rt_df_inchi <- read_rds(file.path(dir_release, "hmsl_compounds_raw.rds"))

# Canonicalize LINCS compounds -------------------------------------------------
###############################################################################T

# Disregard the annotated Chembl ID in the HMS LINCS database because it may
# point to the salt compound and we always want the free base ID
# Only use annotated LINCS Chembl ID if we can't find it using the inchi_key

# We have to first generate the canonical tautomer for each compound
plan(multisession(workers = 4))
hms_lincs_compounds_canonical_inchis <- rt_df_inchi %>%
  drop_na(inchi) %>%
  select(hms_id, inchi) %>%
  chunk_df(12) %>%
  future_map(
    ~canonicalize_compound(set_names(.x$inchi, .x$hms_id)),
    .progress = TRUE
  )

hms_lincs_compounds_canonical_inchis_df <- hms_lincs_compounds_canonical_inchis %>%
  bind_rows() %>%
  distinct(hms_id = compound, inchi)

hms_lincs_compounds_canonical <- rt_df_inchi %>%
  select(hms_id, original_inchi = inchi) %>%
  left_join(
    hms_lincs_compounds_canonical_inchis_df,
    by = "hms_id"
  ) %>%
  mutate(
    # Whe canonicalization failed, use original inchi
    inchi = if_else(!is.na(inchi), inchi, original_inchi)
  )

write_csv(
  hms_lincs_compounds_canonical,
  file.path(dir_release, "hmsl_compounds_canonical.csv.gz")
)
hms_lincs_compounds_canonical <- read_csv(file.path(dir_release, "hmsl_compounds_canonical.csv.gz"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

fetch_hmsl_activity <- Activity(
  name = "Fetch HITS compound data and canonicalize",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/02_hms_lincs_compounds.R"
)

c(
  file.path(dir_release, "hmsl_compounds_raw.rds")
) %>%
  synStoreMany(parent = "syn21064123", activity = fetch_hmsl_activity)

c(
  file.path(dir_release, "hmsl_compounds_canonical.csv.gz")
) %>%
  synStoreMany(parent = "syn21093671", activity = fetch_hmsl_activity)
