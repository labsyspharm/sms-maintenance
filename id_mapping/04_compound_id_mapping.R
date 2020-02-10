library(tidyverse)
library(furrr)
library(data.table)
library(bit64)
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

hmsl_cmpds <- syn("syn20692443") %>%
  read_rds()

chembl_cmpds <- syn("syn20692440") %>%
  read_rds()

eq_classes <- syn("syn20830516") %>%
  read_rds()

pubchem_names <- syn("syn21572728") %>%
  read_csv()

chembl_canonical <- syn("syn20692439") %>%
  read_csv()

hmsl_canonical <- syn("syn20692442") %>%
  read_csv()

# Finding canonical member of equivalence class --------------------------------
###############################################################################T

# Find canonical member of equivalence class
# 1. Choose member with highest clinical phase
# 2. Choose member with annotated pref_name
# 3. If present, choose member annotated as parent by Chembl
# 4. Choose member with highest n_assays
# 5. Choose member with lowest id (chembl or HMS)
find_canonical_member <- function(eq_class_df, compound_df) {
  eq_class_df %>%
    left_join(
      compound_df,
      by = "id"
    ) %>%
    mutate(
      has_parent = if_else(parent_molregno != molregno, parent_molregno, NA_integer_)
    ) %>%
    as.data.table() %>%
    .[
      ,
      `:=`(
        annotated_as_parent = as.integer(molregno) %in% has_parent,
        annotated_pref_name = !is.na(pref_name),
        annotated_assays = replace_na(n_assays, 0),
        id_number = as.integer(str_extract(id, "\\d+"))
      ),
      keyby = eq_class
    ] %>%
    .[
      order(
        eq_class,
        -max_phase,
        -annotated_pref_name,
        -annotated_as_parent,
        -annotated_assays,
        -id_number
      ),
      head(.SD, 1),
      keyby = .(eq_class, source)
      ]
}

canonical_members <- eq_classes %>%
  mutate(
    data = map(
      data,
      find_canonical_member,
      compound_df = bind_rows(
        chembl_cmpds %>%
          select(id = chembl_id, pref_name, max_phase, molregno, parent_molregno, n_assays) %>%
          mutate_if(is.integer64, as.integer),
        hmsl_cmpds %>%
          select(id = hms_id, pref_name = name, n_assays = n_batches),
      )
    )
  )

write_rds(
  canonical_members,
  file.path(dir_release, "canonical_members.rds")
)
# canonical_members <- read_rds(
#   file.path(dir_release, "canonical_members.rds")
# )

# Finding canonical names ------------------------------------------------------
###############################################################################T

all_names <- list(
  "chembl" = chembl_cmpds %>%
    transmute(
      id = chembl_id,
      name = map2(
        pref_name, synonyms,
        ~c(if (!is.na(.x)) .x, if (!is.null(.y)) .y)
      )
    ) %>%
    unchop(name),
  "pubchem" = pubchem_names %>%
    select(id = chembl_id, name) %>%
    filter(str_starts(name, "S?CHEMBL", negate = TRUE)),
  "hmsl" = hmsl_cmpds %>%
    transmute(
      id = hms_id,
      name = map2(
        name, alternate_names,
        ~c(if (!is.na(.x)) .x, if (!is.null(.y)) .y)
      )
    ) %>%
    unchop(name)
) %>%
  bind_rows(.id = "source") %>%
  # Arrange sources by most to least trust-worthy
  mutate(
    source = factor(source, levels = c("chembl", "hmsl", "pubchem"))
  ) %>%
  arrange(id, source) %>%
  distinct()

find_canonical_names <- function(eq_class_df, name_df) {
  as.data.table(eq_class_df) %>%
    .[, list(eq_class, id)] %>%
    .[
      name_df,
      on = "id",
      nomatch = NULL
    ] %>%
    .[
      order(eq_class, source)
    ] %>%
    .[
      ,
      .(
        pref_name = name[[1]],
        alt_names = list(name[-1])
      ),
      by = eq_class
    ]
}

canonical_names <- eq_classes %>%
  mutate(
    data = map(
      data,
      find_canonical_names,
      name_df = as.data.table(all_names)
    ) %>%
      map(as_tibble)
  )

# Finding canonical Inchis -----------------------------------------------------
###############################################################################T

all_inchis <- list(
  "chembl" = chembl_canonical %>%
    distinct(id = chembl_id, inchi),
  "hmsl" = hmsl_canonical %>%
    distinct(id = hms_id, inchi)
) %>%
  bind_rows(.id = "source") %>%
  # Arrange sources by most to least trust-worthy
  mutate(
    source = factor(source, levels = c("chembl", "hmsl"))
  ) %>%
  arrange(id, source) %>%
  distinct() %>%
  drop_na()

find_canonical_inchi <- function(eq_class_df, raw_inchis) {
  as.data.table(eq_class_df) %>%
    .[, list(eq_class, id)] %>%
    .[
      raw_inchis,
      on = "id",
      nomatch = NULL
    ] %>%
    .[
      order(eq_class, source)
    ] %>%
    .[
      ,
      .(
        inchi = inchi[[1]]
      ),
      by = "eq_class"
    ]
}

canonical_inchis <- eq_classes %>%
  mutate(
    data = map(
      data,
      find_canonical_inchi,
      raw_inchis = all_inchis %>%
        as.data.table()
    ) %>%
      map(as_tibble)
  )

# Find canonical SMILES representation -----------------------------------------
###############################################################################T

plan(multisession(workers = 6))
canonical_smiles <- canonical_inchis %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        {set_names(.[["inchi"]], .[["eq_class"]])} %>%
        split(sample(1:18, length(.), replace = TRUE)) %>%
        future_map(
          convert_compound_identifier,
          identifier = "inchi",
          target_identifier = "smiles",
          .progress = TRUE
        ) %>%
        bind_rows() %>%
        select(eq_class = names, smiles = compounds) %>%
        mutate_at(vars(eq_class), as.integer)
    )
  )

# Create table of canonical compounds ------------------------------------------
###############################################################################T

create_canonical_table <- function(members, names, inchis, smiles) {
  members %>%
    select(eq_class, source, id) %>%
    mutate(source = paste0(source, "_id")) %>%
    spread(source, id) %>%
    left_join(
      names,
      by = "eq_class"
    ) %>%
    left_join(
      inchis,
      by = "eq_class"
    ) %>%
    left_join(
      smiles,
      by = "eq_class"
    ) %>%
    as_tibble() %>%
    rename(
      lspci_id = eq_class,
      hms_id = hmsl_id
    ) %>%
    # Remove empty lists and replace with NULL
    mutate(
      alt_names = map(
        alt_names,
        ~if(length(.x) == 0) NULL else .x
      )
    )
}

canonical_table <- canonical_members %>%
  rename(members = data) %>%
  left_join(
    canonical_names %>%
      rename(names = data),
    by = c("fp_name", "fp_type")
  ) %>%
  left_join(
    canonical_inchis %>%
      rename(inchis = data),
    by = c("fp_name", "fp_type")
  ) %>%
  left_join(
    canonical_smiles %>%
      rename(smiles = data),
    by = c("fp_name", "fp_type")
  ) %>%
  mutate(
    data = pmap(
      list(members, names, inchis, smiles),
      create_canonical_table
    )
  )

write_rds(
  canonical_table %>%
    select(fp_name, fp_type, data),
  file.path(dir_release, "canonical_table.rds"),
  compress = "gz"
)

# canonical_table <- read_rds(file.path(dir_release, "canonical_table.rds"))


# Checking proportion of mapped HMSL ids ---------------------------------------
###############################################################################T

canonical_table %>%
  select(fp_name, fp_type, data) %>%
  unnest(data) %>%
  transmute(fp_name, fp_type, lspci_id, chembl_id = !is.na(chembl_id), hms_id = !is.na(hms_id)) %>%
  group_by(fp_name, fp_type, chembl_id, hms_id) %>%
  summarize(n = n())
# # A tibble: 9 x 5
# # Groups:   fp_name, fp_type, chembl_id [6]
# fp_name            fp_type     chembl_id hms_id       n
# <chr>              <chr>       <lgl>     <lgl>    <int>
# 1 morgan_chiral      morgan      FALSE     TRUE       183
# 2 morgan_chiral      morgan      TRUE      FALSE  1757877
# 3 morgan_chiral      morgan      TRUE      TRUE      1799
# 4 morgan_normal      morgan      FALSE     TRUE        94
# 5 morgan_normal      morgan      TRUE      FALSE  1696169
# 6 morgan_normal      morgan      TRUE      TRUE      1871
# 7 topological_normal topological FALSE     TRUE        94
# 8 topological_normal topological TRUE      FALSE  1694947
# 9 topological_normal topological TRUE      TRUE      1871

rt_map <- read_csv(
  here("rt_chembl_matches_20190617.txt")
) %>%
  rename_all(paste0, "_legacy") %>%
  group_by(hmsl_id_legacy) %>%
  summarize_at(vars(chembl_id_legacy), list)

canonical_table %>%
  select(fp_name, fp_type, data) %>%
  unnest(data) %>%
  left_join(rt_map, by = c("hms_id" = "hmsl_id_legacy")) %>%
  mutate_at(vars(hms_id, chembl_id), negate(is.na)) %>%
  mutate(chembl_id_legacy = map_lgl(chembl_id_legacy, negate(is.null))) %>%
  group_by(fp_name, fp_type, hms_id, chembl_id, chembl_id_legacy) %>%
  summarize(n = n())
# # A tibble: 15 x 6
# # Groups:   fp_name, fp_type, hms_id, chembl_id [9]
# fp_name            fp_type     hms_id chembl_id chembl_id_legacy       n
# <chr>              <chr>       <lgl>  <lgl>     <lgl>              <int>
#   1 morgan_chiral      morgan      FALSE  TRUE      FALSE            1757877
# 2 morgan_chiral      morgan      TRUE   FALSE     FALSE                166
# 3 morgan_chiral      morgan      TRUE   FALSE     TRUE                  17
# 4 morgan_chiral      morgan      TRUE   TRUE      FALSE               1510
# 5 morgan_chiral      morgan      TRUE   TRUE      TRUE                 289
# 6 morgan_normal      morgan      FALSE  TRUE      FALSE            1696169
# 7 morgan_normal      morgan      TRUE   FALSE     FALSE                 93
# 8 morgan_normal      morgan      TRUE   FALSE     TRUE                   1
# 9 morgan_normal      morgan      TRUE   TRUE      FALSE               1567
# 10 morgan_normal      morgan      TRUE   TRUE      TRUE                 304
# 11 topological_normal topological FALSE  TRUE      FALSE            1694947
# 12 topological_normal topological TRUE   FALSE     FALSE                 93
# 13 topological_normal topological TRUE   FALSE     TRUE                   1
# 14 topological_normal topological TRUE   TRUE      FALSE               1567
# 15 topological_normal topological TRUE   TRUE      TRUE                 304

# Only 1 compounds frorm HMSL not mapped to Chembl that where found in previous
# version. But 1567 additional compounds in new version (in non-chiral mode)

# Making non-canonical compound table ------------------------------------------
###############################################################################T

create_non_canonical_table <- function(all_eq_class, cmpd_inchis) {
  all_eq_class %>%
    inner_join(
      cmpd_inchis %>%
        distinct(id, inchi),
      by = "id"
    ) %>%
    rename(lspci_id = eq_class)
}

alt_table <- eq_classes %>%
  mutate(
    data = map(
      data,
      create_non_canonical_table,
      cmpd_inchis = all_inchis
    )
  )

write_rds(
  alt_table,
  file.path(dir_release, "lspci_id_compound_inchi_mapping.rds"),
  compress =  "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

cmpd_wrangling_activity <- Activity(
  name = "Map and wrangle compound IDs",
  used = c(
    "syn20692501",
    "syn20692443",
    "syn20692440",
    "syn20692514",
    "syn20830516",
    "syn21572728",
    "syn20692439",
    "syn20692442"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/04_compound_id_mapping.R"
)

syn_id_mapping <- Folder("id_mapping", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "canonical_table.rds"),
  file.path(dir_release, "lspci_id_compound_inchi_mapping.rds")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = cmpd_wrangling_activity)

