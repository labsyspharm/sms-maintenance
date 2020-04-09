library(tidyverse)
library(httr)
library(furrr)
library(bit64)
library(jsonlite)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")


# Load compound data -----------------------------------------------------------
###############################################################################T

compound_sources <- tribble(
  ~source_name, ~synapse_id,
  "chembl", "syn20692439",
  "hmsl", "syn20692442"
)

compounds_canonical <- compound_sources %>%
  mutate(
    canonical = map(
      synapse_id,
      . %>%
        syn() %>%
        read_csv(col_types = "ccc")
    ) %>%
      map(
        rename_all,
        str_replace,
        "chembl_id|hms_id", "id"
      )
  ) %>%
  select(-synapse_id) %>%
  unnest(canonical)

write_csv(
  compounds_canonical,
  file.path(dir_release, "all_compounds_canonical.csv.gz")
)


# Create molecular fingerprints ------------------------------------------------
###############################################################################T

fingerprinting_args <- tribble(
  ~fp_name, ~fp_type, ~fp_args,
  "morgan_chiral", "morgan", list(useChirality = TRUE),
  "morgan_normal", "morgan", list(useChirality = FALSE),
  "topological_normal", "topological", NULL
)

plan(multicore, workers = 10)
cmpd_fingerprints <- compounds_canonical %>%
  distinct(name = id, compound = inchi) %>%
  # Some HMSL compounds don't report inchi or smiles or anything, should have removed
  # them earlier, removing them here now
  drop_na() %>%
  # slice(1:100) %>%
  chunk_df(30) %>%
  enframe("index", "compounds") %>%
  mutate(compounds = map(compounds, ~set_names(.x[["compound"]], .x[["name"]]))) %>%
  crossing(fingerprinting_args) %>%
  mutate(
    fp_res = future_pmap(
      select(., compounds, fp_type, fp_args),
      ~safely(calculate_fingerprints)(..1, fingerprint_type = ..2, fingerprint_args = ..3),
      .progress = TRUE
    )
  )

write_rds(cmpd_fingerprints %>% select(-compounds), file.path(dir_release, "all_compounds_fingerprints_raw.rds"))
# cmpd_fingerprints <- read_rds(file.path(dir_release, "all_compounds_fingerprints_raw.rds"))

# Check if any errors occured
map(cmpd_fingerprints[["fp_res"]], "error") %>% map_lgl(is.null) %>% all()
# TRUE
# None occured

cmpd_fingerprints_all <- cmpd_fingerprints %>%
  transmute(
    fp_name,
    fp_type,
    fp_res = map(fp_res, "result")
  ) %>%
  unnest(fp_res) %>%
  rename(id = names, fingerprint = fingerprints) %>%
  group_nest(fp_name, fp_type)


skipped_cmpds <- compounds_canonical %>%
  anti_join(
    cmpd_fingerprints_all %>%
      unnest(data),
    by = "id"
  )

write_csv(skipped_cmpds, file.path(dir_release, "all_compounds_fingerprints_skipped.csv"))

write_rds(
  cmpd_fingerprints_all,
  file.path(dir_release, "all_compounds_fingerprints.rds")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

fingerprint_activity <- Activity(
  name = "Calculate molecular fingerprints",
  used = c(
    "syn20692439",
    "syn20692442"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/03_fingerprints.R"
)

fp_folder <- Folder("fingerprints", syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "all_compounds_fingerprints.rds"),
  file.path(dir_release, "all_compounds_canonical.csv.gz"),
  file.path(dir_release, "all_compounds_fingerprints_skipped.csv")
) %>%
  synStoreMany(parent = fp_folder, activity = fingerprint_activity)
