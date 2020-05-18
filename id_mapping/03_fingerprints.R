library(tidyverse)
library(furrr)
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
  ~source, ~canonical, ~raw,
  "chembl", "syn20692439", "syn20692440",
  "hmsl", "syn20692442", "syn20692443"
)

compounds <- compound_sources %>%
  mutate_at(
    vars(-source),
    map,
    function(s) {
      f <- syn(s)
      {
        if (str_ends(f, "gz"))
          read_csv(f, col_types = "ccc")
        else
          read_rds(f)
      } %>%
        rename_all(
          str_replace,
          fixed(if ("hms_id" %in% names(.)) "hms_id" else "chembl_id"),
          "id"
        ) %>%
        rename_all(
          str_replace,
          fixed("standard_inchi"),
          "inchi"
        )
    }
  )

compounds_all <- compounds %>%
  transmute(
    source,
    data = map2(
      canonical, raw,
      ~.y %>%
        select(id, raw_inchi = inchi) %>%
        left_join(
          select(.x, id, canonical_inchi = inchi),
          by = "id"
        )
    )
  ) %>%
  unnest(data) %>%
  # Compute chemical formula
  drop_na(raw_inchi) %>%
  mutate(
    formula = raw_inchi %>%
      str_split_fixed(fixed("/"), 3) %>%
      {.[, 2]},
    formula_vector = formula %>%
      str_match_all("([A-Z][a-z]?)([0-9]*)") %>%
      map(~set_names(replace_na(as.integer(.x[, 3]), 1L), .x[, 2]))
  )

compounds_selected <- compounds_all %>%
  mutate(
    inchi = case_when(
      is.na(canonical_inchi) ~ raw_inchi,
      # Only use canonicalized InChI for compounds that are "drug-like" with
      # at least 3 carbons
      map_int(formula_vector, ~if ("C" %in% names(.x)) .x[["C"]] else 0L) >= 3 ~ canonical_inchi,
      TRUE ~ raw_inchi
    )
  )

write_rds(
  compounds_selected,
  file.path(dir_release, "all_compounds_canonical.rds"),
  compress = "gz"
)
# compounds_selected <- read_rds(
#   file.path(dir_release, "all_compounds_canonical.rds")
# )

# Create molecular fingerprints ------------------------------------------------
###############################################################################T

fingerprinting_args <- tribble(
  ~fp_name, ~fp_type, ~fp_args,
  "morgan_chiral", "morgan", list(useChirality = TRUE),
  "morgan_normal", "morgan", list(useChirality = FALSE),
  "topological_normal", "topological", NULL
)

plan(multicore, workers = 10)
cmpd_fingerprints <- compounds_selected %>%
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

write_rds(
  cmpd_fingerprints %>% select(-compounds),
  file.path(dir_release, "all_compounds_fingerprints_raw.rds"),
  compress = "gz"
)
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


skipped_cmpds <- compounds_selected %>%
  anti_join(
    cmpd_fingerprints_all %>%
      unnest(data),
    by = "id"
  )

write_rds(
  skipped_cmpds,
  file.path(dir_release, "all_compounds_fingerprints_skipped.rds"),
  compress = "gz"
)

write_rds(
  cmpd_fingerprints_all,
  file.path(dir_release, "all_compounds_fingerprints.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

fingerprint_activity <- Activity(
  name = "Calculate molecular fingerprints",
  used = c(
    "syn20692439",
    "syn20692442",
    "syn20692440",
    "syn20692443"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/03_fingerprints.R"
)

fp_folder <- Folder("fingerprints", syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "all_compounds_fingerprints.rds"),
  file.path(dir_release, "all_compounds_canonical.rds"),
  file.path(dir_release, "all_compounds_fingerprints_skipped.rds")
) %>%
  synStoreMany(parent = fp_folder, activity = fingerprint_activity)
