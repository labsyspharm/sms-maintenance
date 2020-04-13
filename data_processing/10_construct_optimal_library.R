library(tidyverse)
library(here)
library(furrr)
library(lspcheminf)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Load files -------------------------------------------------------------------
###############################################################################T

selectivity <- syn("syn20836653") %>%
  read_rds()

canonical_fps <- syn("syn21042105") %>%
  read_rds()

commercial_info <- syn("syn21049601") %>%
  read_rds()

kinases <- syn("syn12617467") %>%
  read_csv(col_types = "iccllllic")

affinity <- syn("syn20830834") %>%
  read_rds()

clinical_info <- syn("syn21064122") %>%
  read_rds()

# Prepare chemical similarity  files -------------------------------------------
###############################################################################T

fps_selective <- canonical_fps %>%
  left_join(rename(selectivity, selectivity_df = data)) %>%
  mutate(
    data = pmap(
      list(data, fp_name, selectivity_df),
      ~filter(
        ..1,
        fp_name == ..2,
        lspci_id %in% filter(..3, selectivity_class %in% c("most_selective", "semi_selective"))$lspci_id
      ) %>%
        select(fingerprint, lspci_id)
    )
  )


plan(multicore(workers = 8))
chemical_sim_selective <- fps_selective %>%
  mutate(
    data = future_map(
      data,
      ~chemical_similarity_threshold(
        set_names(.x[["fingerprint"]], .x[["lspci_id"]]),
        precalculated = TRUE,
        threshold = 0.2
      ) %>%
        mutate_at(
          vars(query, target),
          as.integer
        ) %>%
        distinct(
          lspci_id_1 = if_else(query > target, target, query),
          lspci_id_2 = if_else(query > target, query, target),
          tanimoto_similarity = score
        ) %>%
        arrange(lspci_id_1, lspci_id_2),
      progress = TRUE
    )
  )

write_rds(
  chemical_sim_selective,
  file.path(dir_release, "chemical_sim_selective.rds"),
  compress = "gz"
)

# chemical_sim_selective <- read_rds(file.path(dir_release, "chemical_sim_selective.rds"))

# Construct optimal library ----------------------------------------------------
###############################################################################T

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
  # Check all possible combinations of compounds
  combn(
    sort(unique(selectivity$lspci_id)), 2, simplify = TRUE
  ) %>%
    t() %>%
    `colnames<-`(c("lspci_id_1", "lspci_id_2")) %>%
    as_tibble() %>%
    # Get rid of all compound combinations that have similarity of >threshold,
    # Those are listed in the chemical_similarity data frame
    anti_join(chemical_similarity, by = c("lspci_id_1", "lspci_id_2")) %>%
    left_join(
      rename_all(selectivity, paste0, "_1"),
      by = "lspci_id_1"
    ) %>%
    left_join(
      rename_all(selectivity, paste0, "_2"),
      by = "lspci_id_2"
    ) %>%
    mutate(tool_score_sum = tool_score_1 + tool_score_2) %>%
    arrange(
      desc(as.integer(commercially_available_1) + as.integer(commercially_available_2)),
      as.integer(selectivity_class_1) + as.integer(selectivity_class_2),
      desc(tool_score_sum),
      ontarget_IC50_Q1_1 + ontarget_IC50_Q1_2
    )
}

find_optimal_compounds_all_targets <- function(
  selectivity, chemical_sim, commercial_info, maximum_tanimoto_similarity = 0.2,
  selectivity_classes = c("most_selective", "semi_selective")
) {
  sel_by_gene <- selectivity %>%
    filter(selectivity_class %in% selectivity_classes) %>%
    mutate(
      selectivity_class = factor(selectivity_class, levels = selectivity_classes),
      commercially_available = lspci_id %in% commercial_info$lspci_id
    ) %>%
    select(gene_id, lspci_id, tool_score, selectivity_class, ontarget_IC50_Q1, commercially_available) %>%
    group_nest(gene_id)
  chemical_sim_filtered <- chemical_sim %>%
    filter(tanimoto_similarity >= maximum_tanimoto_similarity)
  cmpds <- sel_by_gene %>%
    # Only keep targets with at least 2 compounds
    # filter(map_lgl(data, ~nrow(.x) >= 2)) %>%
    transmute(
      gene_id,
      all_pairs = map(
        data,
        find_optimal_compounds, chemical_sim_filtered
      ),
      best_pair = map(
        all_pairs,
        slice, 1
      )
    ) %>%
    unnest(best_pair)
  cmpds
}

optimal_libraries <- chemical_sim_selective %>%
  left_join(rename(commercial_info, commercial_info_df = data), by = c("fp_type", "fp_name")) %>%
  rename(chemical_similarity_df = data) %>%
  mutate(
    data = pmap(
      list(selectivity_df, chemical_similarity_df, commercial_info_df),
      find_optimal_compounds_all_targets
    )
  )

write_rds(
  optimal_libraries %>%
    select(fp_name, fp_type, data),
  file.path(dir_release, "optimal_library.rds"),
  compress = "gz"
)

# optimal_libraries <- read_rds(file.path(dir_release, "optimal_library.rds"))

# Criteria clinical phase >=1 and affinity IC50_Q1 <= 1 uM
clinical_compounds <- clinical_info %>%
  left_join(rename(affinity, affinity_df = data), by = c("fp_type", "fp_name")) %>%
  mutate(
    data = map2(
      data, affinity_df,
      ~inner_join(
        distinct(.x, lspci_id, max_phase),
        filter(.y, Q1 <= 1000),
        by = "lspci_id"
      )
    )
  )

# Calculate liganded genome ----------------------------------------------------
###############################################################################T

# Liganded genome 3 cmpds < 10 uM
# liganded_genome <- affinity %>%
#   mutate(
#     data = map(
#       data,
#
#     )
#   )

# Subset optimal library to contain only kinase targets ------------------------
###############################################################################T

optimal_kinase_libraries <- optimal_libraries %>%
  mutate(
    data = map(
      data,
      filter,
      # Not taking into account kinases annotated in uniprot by Nienke's advice
      gene_id %in% filter(kinases, in_manning | in_kinmap | in_IDG_darkkinases)$gene_id
    )
  )

write_rds(
  optimal_kinase_libraries %>%
    select(fp_type, fp_name, data),
  file.path(dir_release, "optimal_kinase_library.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

library_activity <- Activity(
  name = "Construct optimal compound library",
  used = c(
    "syn12617467",
    "syn20836653",
    "syn21042105",
    "syn21049601",
    "syn20830834",
    "syn21064122"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/10_construct_optimal_library.R"
)

syn_library_folder <- Folder("compound_library", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "optimal_library.rds"),
  file.path(dir_release, "optimal_kinase_library.rds")
) %>%
  synStoreMany(parentId = syn_library_folder, activity = library_activity)

