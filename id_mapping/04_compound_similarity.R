library(tidyverse)
library(furrr)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Calculate similarity between all compounds -----------------------------------
###############################################################################T

all_compounds_fingerprints <- syn("syn20692501") %>%
  read_rds() %>%
  mutate(fn = map_chr(1:n(), ~tempfile(fileext = ".fps")))

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
  file.path(dir_release, "all_compounds_similarity.csv.gz")
) %>%
  synStoreMany(parentId = syn_id_mapping, activity = activity)

