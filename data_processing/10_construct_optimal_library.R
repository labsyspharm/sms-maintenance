library(tidyverse)
library(here)
library(furrr)
library(morgancpp)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

source(here("utils", "load_save.R"))

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

inputs <- list(
  compound_dictionary = c("compounds_processed", "compound_dictionary.csv.gz"),
  fingerprints = c("compounds_processed", "lspci_id_fingerprints.csv.gz"),
  selectivity = c("selectivity", "selectivity.csv.gz"),
  affinity_q1 = c("aggregate_data", "dose_response_q1.csv.gz"),
  approval = c("clinical_info", "lspci_id_approval.csv.gz"),
  kinases = c("raw_data", "kinome", "kinase_data.csv")
) %>%
  pluck_inputs(syn_parent = syn_release)

input_data <- inputs %>%
  load_input_data(syn = syn)

# Prepare chemical similarity  files -------------------------------------------
###############################################################################T

fps_selective <- input_data[["fingerprints"]][
  lspci_id %in% input_data[["selectivity"]][
    selectivity_class %in% c("most_selective", "semi_selective")
  ][["lspci_id"]]
]

calc_chemical_sim <- function(df) {
  # browser()
  fps <- MorganFPS$new(
    with(df, set_names(fingerprints, lspci_id))
  )
  set_names(df[["lspci_id"]]) %>%
    map(~fps$tanimoto_all(.x)) %>%
    bind_rows(.id = "query") %>%
    select(lspci_id_1 = query, lspci_id_2 = id, tanimoto_similarity = structural_similarity) %>%
    mutate_at(vars(starts_with("lspci_id")), as.integer) %>%
    filter(lspci_id_1 < lspci_id_2, tanimoto_similarity > 0.2) %>%
    arrange(lspci_id_1, lspci_id_2) %>%
    as_tibble()
}

chemical_sim <- fps_selective %>%
  group_nest(fp_name, fp_type) %>%
  mutate(
    data = map(data, calc_chemical_sim)
  ) %>%
  unnest(data)

fwrite(
  chemical_sim,
  file.path(dir_release, "chemical_sim_selective.csv.gz")
)

# chemical_sim_selective <- read_rds(file.path(dir_release, "chemical_sim_selective.rds"))

# Construct optimal library ----------------------------------------------------
###############################################################################T

SELECTIVITY_CLASSES <- c("most_selective", "semi_selective", "poly_selective", "unknown_selective")

commercially_available_compounds <- input_data[["compound_dictionary"]][
  commercially_available == TRUE
]

find_optimal_compounds <- function(
  selectivity, chemical_similarity
) {
  # If we only have a single compound for this target skip optimal pair
  # computation
  if (nrow(selectivity) < 2) {
    return(
      selectivity %>%
        rename_all(paste0, "_1")
    )
  }
  ranking <- selectivity %>%
    arrange(
      desc(as.integer(commercially_available)),
      as.integer(selectivity_class),
      desc(tool_score),
      ontarget_IC50_Q1
    )
  best_lspci_id <- ranking$lspci_id[[1]]
  # Check if there's another compound with <threshold similarity to best one
  combinations <- tibble(
    lspci_id_1 = if_else(best_lspci_id < ranking$lspci_id, best_lspci_id, ranking$lspci_id),
    lspci_id_2 = if_else(best_lspci_id > ranking$lspci_id, best_lspci_id, ranking$lspci_id)
  ) %>%
    filter(lspci_id_1 != lspci_id_2) %>%
    # Get rid of all compound combinations that have similarity of >threshold,
    # Those are listed in the chemical_similarity data frame
    anti_join(chemical_similarity, by = c("lspci_id_1", "lspci_id_2"))
  bind_cols(
    ranking[1, ] %>%
      rename_all(paste0, "_1"),
    if (nrow(combinations) > 0)
      ranking %>%
        filter(
          lspci_id == combinations[[1, "lspci_id_1"]] |
            lspci_id == combinations[[1, "lspci_id_2"]],
          lspci_id != best_lspci_id
        ) %>%
        rename_all(paste0, "_2")
    else
      NULL
  )
}

find_optimal_compounds_all_targets <- function(
  selectivity, chemical_sim, commercial_info, maximum_tanimoto_similarity = 0.2,
  selectivity_classes = SELECTIVITY_CLASSES
) {
  sel_by_gene <- selectivity %>%
    filter(selectivity_class %in% selectivity_classes) %>%
    mutate(
      selectivity_class = factor(selectivity_class, levels = selectivity_classes),
      commercially_available = lspci_id %in% commercial_info$lspci_id
    ) %>%
    select(entrez_gene_id, lspci_id, tool_score, selectivity_class, ontarget_IC50_Q1, commercially_available) %>%
    group_nest(entrez_gene_id)
  chemical_sim_filtered <- chemical_sim %>%
    filter(tanimoto_similarity >= maximum_tanimoto_similarity)
  cmpds <- sel_by_gene %>%
    # Only keep targets with at least 2 compounds
    # filter(map_lgl(data, ~nrow(.x) >= 2)) %>%
    transmute(
      entrez_gene_id,
      best_pair = map(
        data,
        find_optimal_compounds, chemical_sim_filtered
      )
    ) %>%
    unnest(best_pair) %>%
    mutate(reason_included = "selectivity")
  cmpds
}

optimal_libraries_all <- chemical_sim %>%
  group_nest(fp_name, fp_type) %>%
  mutate(
    data = pmap(
      list(
        list(input_data[["selectivity"]]),
        data,
        list(commercially_available_compounds)
      ),
      find_optimal_compounds_all_targets,
      maximum_tanimoto_similarity = 0.2,
      selectivity_classes = c(
        "most_selective",
        "semi_selective",
        "poly_selective"
      )
    )
  )

qsave(
  optimal_libraries_all,
  file.path(dir_release, "optimal_library.qs"),
  preset = "fast"
)

optimal_libraries <- optimal_libraries_all %>%
  filter(fp_name == "morgan_normal") %>%
  chuck("data", 1)

# Only save library based on morgan fingerprints as csv
fwrite(
  optimal_libraries,
  file.path(dir_release, "optimal_library.csv.gz")
)

# optimal_libraries <- read_rds(file.path(dir_release, "optimal_library.rds"))

# Criteria clinical phase >=1 and affinity IC50_Q1 <= 1 uM
clinical_compounds <- input_data[["approval"]][
  ,
  .(lspci_id, max_phase)
] %>%
  unique() %>%
  inner_join(
    bind_rows(
      input_data[["selectivity"]],
      input_data[["affinity_q1"]] %>%
        rename(
          ontarget_IC50_N = n_measurement,
          ontarget_IC50_Q1 = Q1
        ) %>%
        anti_join(
          input_data[["selectivity"]],
          by = c("lspci_id", "entrez_gene_id")
        )
    ),
    by = "lspci_id"
  ) %>%
  filter(
    ontarget_IC50_Q1 < 1000
  ) %>%
  mutate(
    reason_included = "clinical",
    commercially_available = lspci_id %in% commercially_available_compounds[["lspci_id"]]
  )


fwrite(
  clinical_compounds,
  file.path(dir_release, "clinical_library.csv.gz")
)

# Compare with old version
#
# old_clinical_compounds <- syn("syn22089541") %>%
#   read_rds() %>%
#   chuck("data", 2)
#
# old_optimal_libraries <- syn("syn21092727") %>%
#   read_rds() %>%
#   chuck("data", 2)
#
# old_affinities <- syn("syn20830834") %>%
#   read_rds() %>%
#   chuck("data", 2)
#
# old_clinical_info <- syn("syn21064122") %>%
#   read_rds() %>%
#   chuck("data", 2)
#
# clinical_old_not_new <- old_clinical_compounds %>%
#   anti_join(
#     clinical_compounds,
#     by = c("lspci_id", "gene_id" = "entrez_gene_id")
#   )
#
# clinical_overlap <- full_join(
#   old_clinical_compounds %>%
#     distinct(lspci_id, entrez_gene_id = gene_id, old = TRUE),
#   clinical_compounds %>%
#     distinct(lspci_id, entrez_gene_id, new = TRUE)
# ) %>%
#   mutate(across(c(old, new), replace_na, replace = FALSE))
#
# clinical_overlap %>%
#   count(old, new)
#
# clinical_overlap_source_data <- clinical_overlap %>%
#   left_join(
#     input_data[["affinity_q1"]] %>%
#       group_by(lspci_id, entrez_gene_id) %>%
#       summarize(
#         affinity_data = TRUE,
#         affinity_low = any(Q1 < 1000)
#       )
#   ) %>%
#   left_join(
#     input_data[["approval"]] %>%
#       distinct(lspci_id, approval_data = TRUE)
#   ) %>%
#   mutate(across(where(is.logical), replace_na, replace = FALSE))
#
# clinical_overlap_source_data %>%
#   count(old, new, affinity_data, affinity_low, approval_data)
#
#
# clinical_overlap_source_data_old <- clinical_overlap %>%
#   left_join(
#     old_affinities %>%
#       group_by(lspci_id, entrez_gene_id) %>%
#       summarize(
#         affinity_data = TRUE,
#         affinity_low = any(Q1 < 1000)
#       )
#   ) %>%
#   left_join(
#     old_clinical_info %>%
#       distinct(lspci_id, approval_data = TRUE)
#   ) %>%
#   mutate(across(where(is.logical), replace_na, replace = FALSE))
#
# clinical_overlap_source_data_old %>%
#   count(old, new, affinity_data, affinity_low, approval_data)
#
#
# affinity_comparison <- full_join(
#   input_data[["affinity_q1"]] %>%
#     group_by(lspci_id, entrez_gene_id) %>%
#     summarize(
#       Q1_new = Q1,
#       affinity_new = TRUE,
#       affinity_new_low = any(Q1 < 1000),
#       .groups = "drop"
#     ),
#   old_affinities %>%
#     group_by(lspci_id, entrez_gene_id) %>%
#     summarize(
#       Q1_old = Q1,
#       affinity_old = TRUE,
#       affinity_old_low = any(Q1 < 1000),
#       .groups = "drop"
#     )
# ) %>%
#   mutate(across(where(is.logical), replace_na, replace = FALSE))
#
# affinity_comparison %>%
#   count(affinity_new, affinity_new_low, affinity_old, affinity_old_low)
#
# input_data[["approval"]] %>%
#   filter(lspci_id %in% clinical_old_not_new[["lspci_id"]])

# Combine libraries ------------------------------------------------------------
###############################################################################T

optimal_libraries_combined <- bind_rows(
    optimal_libraries %>%
      select(entrez_gene_id, reason_included, ends_with("_1")) %>%
      rename_all(str_replace, fixed("_1"), ""),
    optimal_libraries %>%
      select(entrez_gene_id, reason_included, ends_with("_2")) %>%
      rename_all(str_replace, fixed("_2"), "") %>%
      drop_na(lspci_id)
  ) %>%
  bind_rows(
    anti_join(
      clinical_compounds,
      .,
      by = c("lspci_id", "entrez_gene_id")
    ) %>%
      select(entrez_gene_id, lspci_id, reason_included, selectivity_class, tool_score, ontarget_IC50_Q1, commercially_available)
  ) %>%
  mutate(
    selectivity_class = factor(selectivity_class, levels = SELECTIVITY_CLASSES)
  ) %>%
  arrange(
    entrez_gene_id,
    desc(as.integer(commercially_available)),
    as.integer(selectivity_class),
    desc(tool_score),
    ontarget_IC50_Q1
  ) %>%
  group_by(entrez_gene_id) %>%
  mutate(rank = 1:n()) %>%
  ungroup()

fwrite(
  optimal_libraries_combined,
  file.path(dir_release, "optimal_library_combined.csv.gz")
)

# Calculate liganded genome ----------------------------------------------------
###############################################################################T

# Liganded genome genes with at least 3 cmpds < 10 uM
liganded_genome <- input_data[["affinity_q1"]] %>%
  filter(Q1 < 10000) %>%
  group_by(entrez_gene_id) %>%
  summarize(Q1 = min(Q1), .groups = "drop")

fwrite(
  liganded_genome,
  file.path(dir_release, "liganded_genome.csv.gz")
)

# Subset optimal library to contain only kinase targets ------------------------
###############################################################################T

optimal_kinase_libraries <- optimal_libraries_combined %>%
  filter(
    entrez_gene_id %in% input_data[["kinases"]][["gene_id"]]
  )

write_csv(
  optimal_kinase_libraries,
  file.path(dir_release, "optimal_kinase_library.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

library_activity <- Activity(
  name = "Construct optimal compound library",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/10_construct_optimal_library.R"
)

syn_library_folder <- synMkdir(syn_release, "compound_library")

c(
  file.path(dir_release, "optimal_library.csv.gz"),
  file.path(dir_release, "clinical_library.csv.gz"),
  file.path(dir_release, "optimal_library_combined.csv.gz"),
  file.path(dir_release, "optimal_kinase_library.csv.gz"),
  file.path(dir_release, "liganded_genome.csv.gz")
) %>%
  synStoreMany(parentId = syn_library_folder, activity = library_activity)

