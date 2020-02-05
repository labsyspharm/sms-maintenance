library(tidyverse)
library(furrr)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(igraph)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Loading files ----------------------------------------------------------------
###############################################################################T

similarity_df <- syn("syn20692550") %>%
  read_csv()

hmsl_cmpds <- syn("syn20692443") %>%
  read_rds()

chembl_cmpds <- syn("syn20692440") %>%
  read_rds()

cmpd_mass <- syn("syn21572844") %>%
  read_csv()

# Map HMSL by name -------------------------------------------------------------
###############################################################################T

# Only done for compounds where for some reason no comparison based on molecular
# fingerprint was possible. Eg. no Inchi was provided in Reagent Tracker

hmsl_name_matches <- hmsl_cmpds %>%
  filter(
    !(hms_id %in% (
      similarity_df %>%
        filter(str_starts(query, "HMSL"), str_starts(match, "CHEMBL")) %>%
        pull(query) %>%
        unique()
    )
    )
  ) %>%
  mutate(
    name = str_to_lower(name)
  ) %>%
  select(-chembl_id) %>%
  left_join(
    chembl_cmpds %>%
      transmute(chembl_id, pref_name = str_to_lower(pref_name)),
    by = c("name" = "pref_name")
  ) %>%
  drop_na(chembl_id)

identity_df <- bind_rows(
  similarity_df %>%
    select(-score),
  hmsl_name_matches %>%
    select(query = hms_id, match = chembl_id) %>%
    tidyr::crossing(
      distinct(similarity_df, fp_name, fp_type)
    )
)

# Defining equivalence classes -------------------------------------------------
###############################################################################T

# Making graph of compounds where edges represent identical compounds
# Extracting the "components" of the graph, isolated subgraphs, they
# represent equivalence clusters of identical compounds

calc_identity_classes <- function(df) {
  cmpds_identical_graph <- df %>%
    mutate(
      id1 = ifelse(match > query, query, match),
      id2 = ifelse(match > query, match, query)
    ) %>%
    distinct(id1, id2) %>%
    as.matrix() %>%
    graph_from_edgelist(directed = FALSE)

  components(cmpds_identical_graph) %>%
    pluck("membership") %>%
    enframe(value = "eq_class", name = "id") %>%
    mutate(eq_class = as.integer(eq_class))
}

assign_func <- function(x, y) mutate(x, !!sym(y) := TRUE)

identity_df_combined <- identity_df %>%
  group_nest(fp_name, fp_type) %>%
  # Making a column for every combination that just contains true
  # for joining later
  mutate(
    data = map2(data, fp_name, assign_func)
  ) %>%
  pull("data") %>%
  reduce(full_join, by = c("match", "query")) %>%
  mutate_if(is.logical, replace_na, FALSE) %>%
  left_join(
    cmpd_mass %>%
      distinct(id, mass),
    by = c("query" = "id")
  ) %>%
  left_join(
    cmpd_mass %>%
      distinct(id, mass),
    by = c("match" = "id")
  ) %>%
  # If masses couldn't be calculated, assume they match...
  mutate(mass_identical = replace_na(mass.x == mass.y, TRUE)) %>%
  # Here checking if either of the morgan fingerprints is identical AND topological AND mass
  mutate_at(
    vars(starts_with("morgan"),  topological_normal),
    function(...) pmap_lgl(list(...), all),
    .$topological_normal, .$mass_identical
  ) %>%
  # Correcting the 6 crazy cases where morgan_normal is FALSE and morgan_chiral is TRUE
  # This should be impossible so setting the morgan_normal to TRUE whenever
  # morgan_chiral is TRUE
  mutate(
    morgan_normal = morgan_chiral | morgan_normal
  )

identity_df_combined_nested <- identity_df_combined %>%
  select(-starts_with("mass")) %>%
  gather("fp_name", "is_identical", everything(), -match, -query) %>%
  mutate(fp_type = str_split_fixed(fp_name, fixed("_"), 2)[, 1]) %>%
  filter(is_identical) %>%
  group_nest(fp_type, fp_name)

cmpd_eq_classes <- identity_df_combined_nested %>%
  mutate(
    data = map(data, calc_identity_classes)
  )

# Sanity check equivalence classes ---------------------------------------------
###############################################################################T

# Checking if the parent_molregno annotation in Chembl is comparable to what we
# find using our canonicalization followed by fingerprint matching approach.

# Augmenting identity map with data from the Chembl parent annotation
chembl_cmpds_with_parent <- chembl_cmpds %>%
  filter(molregno != parent_molregno) %>%
  left_join(
    chembl_cmpds %>%
      select(molregno, parent_chembl_id = chembl_id, parent_standard_inchi = standard_inchi),
    by = c("parent_molregno" = "molregno")
  ) %>%
  select(molregno, chembl_id, parent_chembl_id, standard_inchi, parent_standard_inchi)
#
# identity_df_augmented <- identity_df %>%
#   bind_rows(
#     chembl_cmpds_with_parent %>%
#       select(query = chembl_id, match = parent_chembl_id) %>%
#       distinct() %>%
#       tidyr::crossing(select(similarity_df, fp_name, fp_type))
#   )

chembl_cmpds_with_parent_eq_class <- chembl_cmpds_with_parent %>%
  left_join(
    cmpd_eq_classes %>%
      unnest(data),
    by = c("chembl_id" = "id")
  ) %>%
  left_join(
    cmpd_eq_classes %>%
      unnest(data) %>%
      rename(parent_eq_class = eq_class),
    by = c("parent_chembl_id" = "id", "fp_name", "fp_type")
  )

# Compounds where the eq_class of the Chembl parent compound is different
# from the eq_class of the compound itself
chembl_cmpds_with_parent_disagree <- chembl_cmpds_with_parent_eq_class %>%
  filter(eq_class != parent_eq_class)
# Should be zero, some (~200) disagree but that shouldn't be a huge problem


# Adding equivalence class for all compounds -----------------------------------
###############################################################################T

# Add eq_class for compounds for which no inchi is known or whose inchi is not parseable
add_missing_eq_class <- function(eq_class_df, compound_df) {
  compound_df %>%
    left_join(eq_class_df, by = "id") %>%
    mutate(
      eq_class = if_else(
        is.na(eq_class),
        cumsum(is.na(eq_class)) + max(eq_class, na.rm = TRUE),
        eq_class
      )
    )
}

all_eq_class <- cmpd_eq_classes %>%
  mutate(
    data = map(
      data,
      add_missing_eq_class,
      compound_df = bind_rows(
        chembl_cmpds %>%
          distinct(id = chembl_id) %>%
          mutate(source = "chembl"),
        hmsl_cmpds %>%
          distinct(id = hms_id) %>%
          mutate(source = "hmsl")
      )
    )
  )

write_rds(
  all_eq_class,
  file.path(dir_release, "all_compounds_equivalence_classes.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Calculate compound equivalence classes",
  used = c(
    "syn20692550",
    "syn20692443",
    "syn20692440",
    "syn21572844"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/04_compound_equivalence.R"
)

syn_id_mapping <- Folder("id_mapping", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "all_compounds_equivalence_classes.rds")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = activity)

