library(tidyverse)
library(furrr)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Loading files ----------------------------------------------------------------
###############################################################################T

all_compounds_fingerprints <- syn("syn20692501") %>%
  read_rds() %>%
  mutate(fn = map_chr(1:n(), ~tempfile(fileext = ".fps")))

cmpds_canonical <- syn("syn20692514") %>%
  read_csv()

# Calculate similarity between all compounds -----------------------------------
###############################################################################T

all_compounds_fingerprints %>%
  # Adding a newline char at the end
  {pwalk(list(.$fingerprint_db, .$fn), write_lines)}

plan(multisession(workers = 4))
similarity_res <- all_compounds_fingerprints %>%
  filter(fp_name == "topological_normal") %>%
  mutate(
    fp_matches = future_map(fn, scan_fingerprint_matches, threshold = 0.9999)
  )

write_rds(similarity_res, file.path(dir_release, "all_compounds_similarity_raw.rds"))
# similarity_res <- read_rds(file.path(dir_release, "all_compounds_similarity_raw.rds"))

similarity_df <- similarity_res %>%
  mutate_at(vars(fp_matches), map, as_tibble) %>%
  select(-skipped, -fingerprint_db, -fn) %>%
  unnest(fp_matches) %>%
  mutate(score = as.double(score))

write_csv(similarity_df, file.path(dir_release, "all_compounds_similarity.csv.gz"))

# Calculate compound mass -----------------------------------
###############################################################################T

# To establish identity between compounds, require that the topological FP
# matches and one of the Morgan FPs and that their molecular mass is identical

plan(multisession(workers = 4))
cmpd_mass_raw <- cmpds_canonical %>%
  distinct(inchi) %>%
  drop_na(inchi) %>%
  chunk_df(12) %>%
  future_map(
    ~molecular_mass(
      set_names(.x$inchi)
    ),
    .progress = TRUE
  ) %>%
  bind_rows() %>%
  distinct()

cmpd_mass <- cmpd_mass_raw %>%
  left_join(
    cmpds_canonical %>%
      select(source_name, id, inchi),
    by = c("compound" = "inchi")
  )

cmpd_mass_map <- cmpd_mass %>%
  distinct(source_name, id, mass)

write_csv(
  cmpd_mass_map,
  file.path(dir_release, "compound_masses.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Calculate chemical similarity between all compounds",
  used = c(
    "syn20692501"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/04_compound_similarity.R"
)

syn_id_mapping <- Folder("id_mapping", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "all_compounds_similarity.csv.gz"),
  file.path(dir_release, "compound_masses.csv.gz")
) %>%
  synStoreMany(parentId = syn_id_mapping, activity = activity)

