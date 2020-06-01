library(tidyverse)
library(furrr)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)
library(future.apply)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Loading files ----------------------------------------------------------------
###############################################################################T

all_compounds_fingerprints <- syn("syn20692501") %>%
  read_rds()

compounds_canonical <- syn("syn22080194") %>%
  read_rds()

# Calculate similarity between all compounds -----------------------------------
###############################################################################T

plan(multicore(workers = 3))
options(future.globals.maxSize = 2*2**30)
similarity_res <- all_compounds_fingerprints %>%
  mutate(
    data = data %>%
      map(
        ~set_names(.x[["fingerprint"]], .x[["id"]])
      ) %>%
      future_lapply(
        function(fps) {
          chemical_similarity_threshold(
            fps,
            threshold = 0.9999,
            precalculated = TRUE,
            n_threads = 6
          )
        },
        future.scheduling = 1,
        future.packages = "lspcheminf"
      )
  )

# Reorder query and match such that query is always smaller than match
reorder_matches <- function(df) {
  df %>%
    mutate(
      ordered = map2(
        query, target,
        c
      ) %>%
        map(sort)
    ) %>%
    transmute(
      query = map_chr(ordered, 1),
      target = map_chr(ordered, 2),
      score
    )
}

similarity_res_ordered <- similarity_res %>%
  mutate(
    data = map(
      data,
      reorder_matches
    )
  )

write_rds(
  similarity_res_ordered,
  file.path(dir_release, "all_compounds_similarity.rds"),
  compress = "gz"
)
# similarity_res_ordered <- read_rds(file.path(dir_release, "all_compounds_similarity.rds"))


# Calculate compound mass -----------------------------------
###############################################################################T

# To establish identity between compounds, require that the topological FP
# matches and one of the Morgan FPs and that their molecular mass is identical

plan(multicore(workers = 18))
cmpd_mass_raw <- compounds_canonical %>%
  distinct(inchi) %>%
  drop_na(inchi) %>%
  chunk_df(18*3) %>%
  future_map(
    ~molecular_mass(
      set_names(.x[["inchi"]])
    ),
    .progress = TRUE
  ) %>%
  bind_rows() %>%
  distinct()

cmpd_mass <- cmpd_mass_raw %>%
  left_join(
    compounds_canonical %>%
      select(source, id, inchi),
    by = c("compound" = "inchi")
  )

cmpd_mass_map <- cmpd_mass %>%
  distinct(source, id, mass)

write_csv(
  cmpd_mass_map,
  file.path(dir_release, "compound_masses.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Calculate chemical similarity between all compounds",
  used = c(
    "syn20692501",
    "syn22080194"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/04_compound_similarity.R"
)

syn_id_mapping <- Folder("id_mapping", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "all_compounds_similarity.rds"),
  file.path(dir_release, "compound_masses.csv.gz")
) %>%
  synStoreMany(parentId = syn_id_mapping, activity = activity)

