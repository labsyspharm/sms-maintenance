library(tidyverse)
library(here)
library(furrr)
library(morgancpp)
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
        lspci_id %in% {
          filter(..3, selectivity_class %in% c("most_selective", "semi_selective")) %>%
            pull(lspci_id)
        }
      ) %>%
        select(fingerprint, lspci_id)
    )
  )

chemical_sim_selective <- fps_selective %>%
  mutate(
    data = map(
      data,
      function(df) {
        fps <- MorganFPS$new(
          with(df, set_names(fingerprint, lspci_id))
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
  selectivity_classes = c("most_selective", "semi_selective", "poly_selective", "unknown_selective")
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
      best_pair = map(
        data,
        find_optimal_compounds, chemical_sim_filtered
      )
    ) %>%
    unnest(best_pair) %>%
    mutate(reason_included = "selectivity")
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
  left_join(rename(selectivity, selectivity_df = data), by = c("fp_type", "fp_name")) %>%
  left_join(rename(commercial_info, commercial_info_df = data), by = c("fp_type", "fp_name")) %>%
  mutate(
    data = pmap(
      list(data, affinity_df, selectivity_df, commercial_info_df),
      ~inner_join(
        distinct(..1, lspci_id, max_phase) %>%
          drop_na() %>%
          filter(max_phase >= 1),
        bind_rows(
          ..3,
          anti_join(
              ..2 %>%
                select(lspci_id, gene_id = entrez_gene_id, ontarget_IC50_Q1 = Q1, ontarget_IC50_N = n_measurement),
              ..3,
              by = c("lspci_id", "gene_id")
            )
        ) %>%
          filter(ontarget_IC50_Q1 <= 1000),
        by = "lspci_id"
      ) %>%
        mutate(reason_included = "clinical", commercially_available = lspci_id %in% ..4[["lspci_id"]])
    )
  )

write_rds(
  clinical_compounds %>%
    select(fp_name, fp_type, data),
  file.path(dir_release, "clinical_library.rds"),
  compress = "gz"
)

# Combine libraries ------------------------------------------------------------
###############################################################################T

optimal_libraries_combined <- optimal_libraries %>%
  inner_join(
    select(clinical_compounds, fp_type, fp_name, clinical = data)
  ) %>%
  mutate(
    data = pmap(
      list(data, clinical),
      ~bind_rows(
        ..1 %>%
          select(gene_id, reason_included, ends_with("_1")) %>%
          rename_all(str_replace, fixed("_1"), ""),
        ..1 %>%
          select(gene_id, reason_included, ends_with("_2")) %>%
          rename_all(str_replace, fixed("_2"), "") %>%
          drop_na(lspci_id)
      ) %>%
        mutate(selectivity_class = factor(selectivity_class, levels = levels(..2[["selectivity_class"]]))) %>%
        bind_rows(
          anti_join(..2, ., by = c("lspci_id", "gene_id")) %>%
            select(gene_id, lspci_id, reason_included, selectivity_class, tool_score, ontarget_IC50_Q1, commercially_available)
        ) %>%
        arrange(
          gene_id,
          desc(as.integer(commercially_available)),
          as.integer(selectivity_class),
          desc(tool_score),
          ontarget_IC50_Q1
        ) %>%
        group_by(gene_id) %>%
        mutate(rank = 1:n()) %>%
        ungroup()
    )
  )

write_rds(
  optimal_libraries_combined %>%
    select(fp_name, fp_type, data),
  file.path(dir_release, "optimal_library_combined.rds"),
  compress = "gz"
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

optimal_kinase_libraries <- optimal_libraries_combined %>%
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
  file.path(dir_release, "optimal_kinase_library_combined.rds"),
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
  file.path(dir_release, "optimal_kinase_library.rds"),
  file.path(dir_release, "clinical_library.rds"),
  file.path(dir_release, "optimal_library_combined.rds"),
  file.path(dir_release, "optimal_kinase_library_combined.rds")
) %>%
  synStoreMany(parentId = syn_library_folder, activity = library_activity)

