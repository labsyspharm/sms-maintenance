library(tidyverse)
library(furrr)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(igraph)
library(parallel)
library(qs)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Loading files ----------------------------------------------------------------
###############################################################################T

inputs <- c(
  masses =  synPluck(syn_release, "id_mapping", "compound_masses.csv.gz"),
  inchi_id_vendor_map = synPluck(syn_release, "canonicalization", "inchi_id_vendor_map.csv.gz"),
  similarities = synPluck(syn_release, "id_mapping", "all_compounds_similarity.csv.gz"),
  chembl_raw = synPluck(syn_release, "raw_data", "chembl_compounds_raw.rds"),
  hmsl_raw = synPluck(syn_release, "raw_data", "hmsl_compounds_raw.rds"),
  emolecules_suppliers = synPluck(syn_release, "id_mapping", "emolecules", "suppliers.tsv.gz"),
  emolecules_vendor_info = synPluck(syn_release, "id_mapping", "emolecules", "emolecules_vendor_info.csv.gz"),
  old_sms = "syn21094266"
)

input_data <- inputs %>%
  map(syn) %>%
  map(
    function(x)
      list(
        `.csv` = partial(fread, colClasses = c(inchi_id = "integer")),
        `.tsv` = fread,
        `.rds` = read_rds
      ) %>%
      magrittr::extract2(which(str_detect(x, fixed(names(.))))) %>%
      {.(x)}
  )

# https://stackoverflow.com/questions/15280472/in-r-how-do-i-create-consecutive-id-numbers-for-each-repetition-in-a-separate-v

# Generate compound name map ---------------------------------------------------
###############################################################################T

trusted_suppliers <- c(
  729, # Enamine
  1296,
  417,
  681,
  40,
  1362, #Target Mol
  1054, # MCE
  400, # Selleck
  1806, # Tocris
  1805
)

cmpd_names <- tribble(
  ~source, ~data,
  "chembl", input_data[["chembl_raw"]] %>%
    drop_na(pref_name) %>%
    transmute(id = chembl_id, name = pref_name, name_preference = "primary") %>%
    bind_rows(
      input_data[["chembl_raw"]] %>%
        transmute(id = chembl_id, name = synonyms, name_preference = "secondary") %>%
        filter(map_lgl(name, negate(is.null))) %>%
        unchop(name, ptype = tibble(name = character()))
    ),
  "hmsl", input_data[["hmsl_raw"]] %>%
    drop_na(name) %>%
    transmute(id = hms_id, name, name_preference = "primary") %>%
    bind_rows(
      input_data[["hmsl_raw"]] %>%
        transmute(id = hms_id, name = alternate_names, name_preference = "secondary") %>%
        filter(map_lgl(name, negate(is.null))) %>%
        unchop(name, ptype = tibble(name = character()))
    ),
  "emolecules", input_data[["emolecules_vendor_info"]] %>%
    filter(supplier_id %in% trusted_suppliers) %>%
    transmute(id = as.character(parent_id), name = chemical_name, name_preference = "primary") %>%
    drop_na()
) %>%
  unnest(data) %>%
  inner_join(
    input_data[["inchi_id_vendor_map"]],
    by = c("source", "id" = "vendor_id")
  ) %>%
  drop_na() %>%
  distinct()

fwrite(
  cmpd_names,
  file.path(dir_release, "inchi_id_compound_name_map.csv.gz")
)

# Map HMSL by name -------------------------------------------------------------
###############################################################################T

# Only done for compounds where for some reason no comparison based on molecular
# fingerprint was possible. Eg. no Inchi was provided in Reagent Tracker

norm_drug <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^a-zA-Z0-9]", "")
}

hmsl_name_matches <- cmpd_names %>%
  transmute(name = norm_drug(name), inchi_id) %>%
  inner_join(
    cmpd_names %>%
      filter(source == "hmsl") %>%
      transmute(inchi_id, name = norm_drug(name)),
    by = "name"
  ) %>%
  distinct(inchi_id.x, inchi_id.y) %>%
  # Getting rid of one particular false positive match between
  filter(
    # Nobiletin and Hexamethoxyflavone
    !(inchi_id.x == 14475558L & inchi_id.y == 14622253L),
    # G2-02 and G-202
    !(inchi_id.x == 9415742L & inchi_id.y == 1309010L),
    # MI-2
    !(inchi_id.x == 14112358L & inchi_id.y == 16973627L),
    # Dabrafenib and GSK2118436
    !(inchi_id.x == 260420L & inchi_id.y == 9465114L),
    # INO-1001 and 3-AMINOBENZAMIDE
    !(inchi_id.x %in% c(6381976L, 8073517L) & inchi_id.y == 16043830L)
  )


hmsl_name_matches %>%
  filter(inchi_id.x != inchi_id.y)


# Defining equivalence classes -------------------------------------------------
###############################################################################T

identity_combined <- merge(
  input_data[["similarities"]][fp_name == "morgan_normal"][
    , .(identity_group_morgan_normal = identity_group, inchi_id)
  ],
  input_data[["similarities"]][fp_name == "topological_normal"][
    , .(identity_group_topological_normal = identity_group, inchi_id)
  ],
  by = "inchi_id",
  all = TRUE
) %>%
  merge(
    input_data[["masses"]],
    by = "inchi_id",
    all = TRUE
  ) %>%
  merge(
    input_data[["inchi_id_vendor_map"]][
      source == "old_sms", .(inchi_id, lspci_id = as.integer(vendor_id))
    ],
    by = "inchi_id",
    all = TRUE
  )

# In the graph connect all inchi_ids that share the same fingerprint and mass
identity_graph <- identity_combined[
  ,
  if (.N > 1) .(
    inchi_id_1 = head(inchi_id, n = -1L),
    inchi_id_2 = tail(inchi_id, n = -1L),
    identity_type = "fingerprint"
  ) else .(
    inchi_id_1 = inchi_id,
    inchi_id_2 = inchi_id,
    identity_type = "fingerprint"
  ),
  keyby = .(
    # Replacing NA with unique pseudo ID's so that NAs don't match
    identity_group_morgan_normal = identity_group_morgan_normal %>% {
      na_mask <- is.na(.)
      .[na_mask] <- -1L * seq_len(sum(na_mask))
      .
    },
    identity_group_topological_normal = identity_group_topological_normal %>% {
      na_mask <- is.na(.)
      .[na_mask] <- -1L * seq_len(sum(na_mask))
      .
    },
    # However, we do want to match if mass couldn't be calculated so leaving intact
    mass
  )
] %>%
  # Adding HMSL name match equivalences
  bind_rows(
    hmsl_name_matches %>%
      filter(inchi_id.x != inchi_id.y) %>%
      transmute(
        inchi_id_1 = inchi_id.x,
        inchi_id_2 = inchi_id.y,
        identity_type = "hmsl_name_match"
      )
  )

fwrite(
  identity_graph,
  file.path(dir_release, "identity_graph.csv.gz")
)

# identity_graph <- fread(file.path(dir_release, "identity_graph.csv.gz"))

identity_groups <- identity_graph %>%
  distinct(inchi_id_1, inchi_id_2) %>%
  mutate(across(.fns = as.character)) %>%
  as.matrix() %>%
  graph_from_edgelist(directed = FALSE) %>%
  components() %>%
  pluck("membership") %>%
  enframe(value = "identity_group", name = "inchi_id") %>%
  mutate(across(.fns = as.integer))

identity_mapped <- identity_combined[
  setDT(identity_groups),
  on = "inchi_id",
  nomatch = NULL
][
  ,
  .(
    lspci_id = {
      lspci_ids <- lspci_id %>%
        na.omit() %>%
        unique()
      if (length(lspci_ids) != 1L) NA_integer_ else lspci_ids
    },
    inchi_id = {
      inchi_ids <- inchi_id %>%
        na.omit() %>%
        unique()
      if (length(inchi_ids) == 0L) NA_integer_ else inchi_ids
    }
  ),
  keyby = identity_group
]

fwrite(
  identity_mapped,
  file.path(dir_release, "identity_mapped_raw.csv.gz")
)

# identity_mapped <- fread(file.path(dir_release, "identity_mapped_raw.csv.gz"))

#
# inchis <- fread("chembl_v27/canonical_inchi_ids.csv.gz")
#
# identity_groups_wo_hmsl <- identity_graph %>%
#   filter(identity_type != "hmsl_name_match") %>%
#   distinct(inchi_id_1, inchi_id_2) %>%
#   mutate(across(.fns = as.character)) %>%
#   as.matrix() %>%
#   graph_from_edgelist(directed = FALSE) %>%
#   components() %>%
#   pluck("membership") %>%
#   enframe(value = "identity_group", name = "inchi_id") %>%
#   mutate(across(.fns = as.integer))
#
# identity_groups_conflicts <- identity_groups %>%
#   inner_join(
#     identity_groups_wo_hmsl %>%
#       rename(identity_group_wo_hmsl = identity_group),
#     by = "inchi_id"
#   ) %>%
#   group_by(identity_group) %>%
#   filter(length(unique(identity_group_wo_hmsl)) > 1) %>%
#   ungroup() %>%
#   left_join(
#     input_data[["inchi_id_vendor_map"]],
#     by = "inchi_id"
#   )
#
# identity_groups_conflicts_inchis <- identity_groups_conflicts %>%
#   distinct(identity_group, identity_group_wo_hmsl, inchi_id) %>%
#   left_join(
#     inchis,
#     by = "inchi_id"
#   ) %>%
#   arrange(identity_group)
#
# library(lspcheminf)
#
# group_walk(
#   group_by(identity_groups_conflicts_inchis, identity_group),
#   function(df, key) {
#     draw_compound_grid(
#       compounds(
#         with(
#           df,
#           set_names(
#             canonical_inchi,
#             inchi_id
#           )
#         ),
#         descriptor = "inchi"
#       ),
#       here(paste0("hmsl_name_match_compounds_", key[["identity_group"]], ".svg")),
#       draw_args = list(molsPerRow = 2, subImgSize = c(400, 400))
#     )
#   }
# )



# => Matches seem reasonable, adding them. Except for a bunch where
# there are multiple compounds with the same name and the structure
# doesn't match at all. Excluding them manually above


# Checking strange cases with multiple old lspci_ids ---------------------------
###############################################################################T
#
# identity_mapped_multiple_lspci_id <- copy(identity_combined)[
#   ,
#   if (length(unique(na.omit(lspci_id))) > 1) cbind(.SD, identity_group = .GRP),
#   keyby = .(identity_group_morgan_normal, identity_group_topological_normal, mass)
# ]
#
# fwrite(
#   identity_mapped_multiple_lspci_id,
#   file.path(dir_release, "identity_mapped_multiple_lspci_id.csv.gz")
# )
#
# # identity_mapped_multiple_lspci_id <- qread(file.path(dir_release, "identity_mapped_multiple_lspci_id.qs"))
#
# # Checking which compounds map to multiple lspci_ids
# library(lspcheminf)
#
# inchis <- fread("chembl_v27/canonical_inchi_ids.csv.gz")
# #
# old_inchis <- syn("syn21094266") %>%
#   fread()
#
# pwalk(
#   identity_mapped_multiple_lspci_id[
#     , .(data = list(.SD)),
#     keyby = "identity_group"
#   ],
#   function(identity_group, data, ...) {
#     cmpds <- c(
#       set_names(
#         data[["inchi_id"]]
#       )
#     ) %>%
#       na.omit()
#     cmpd_inchis <- cmpds %>%
#       enframe("inchi_name", "inchi_id") %>%
#       inner_join(
#         inchis,
#         by = "inchi_id"
#       ) %>%
#       bind_rows(
#         old_inchis %>%
#           filter(lspci_id %in% data[["lspci_id"]]) %>%
#           transmute(canonical_inchi = inchi, inchi_name = paste0("lspci_id_", lspci_id))
#       )
#     # browser()
#     draw_compound_grid(
#       compounds(
#         with(
#           cmpd_inchis,
#           set_names(canonical_inchi, inchi_name)
#         ),
#         descriptor = "inchi"
#       ),
#       here(paste0(identity_group, "_compounds.svg"))
#     )
#   }
# )


# ==> There are only a couple such cases (58), it is unclear why some of these
# had distinct lspci_ids in the old release to begin with. Fingerprint and mass
# is identical

# Sanity check equivalence classes ---------------------------------------------
###############################################################################T

# Checking if the parent_molregno annotation in Chembl is comparable to what we
# find using our canonicalization followed by fingerprint matching approach.

# Augmenting identity map with data from the Chembl parent annotation
chembl_cmpds_with_parent <- input_data[["chembl_raw"]] %>%
  filter(molregno != parent_molregno) %>%
  left_join(
    input_data[["chembl_raw"]] %>%
      select(molregno, parent_chembl_id = chembl_id),
    by = c("parent_molregno" = "molregno")
  ) %>%
  select(molregno, chembl_id, parent_chembl_id) %>%
  left_join(
    input_data[["inchi_id_vendor_map"]][
      , .(vendor_id, inchi_id)
    ],
    by = c("chembl_id" = "vendor_id")
  ) %>%
  left_join(
    input_data[["inchi_id_vendor_map"]][
      , .(vendor_id, inchi_id_parent = inchi_id)
    ],
    by = c("chembl_id" = "vendor_id")
  )

chembl_cmpds_with_parent %>% filter(inchi_id != inchi_id_parent)

# => 0 rows. That means we actually capture all compound-parent relationships
# annotated in ChEMBL with our pipeline

# Adding equivalence class for all compounds -----------------------------------
###############################################################################T

# Add eq_class for compounds for which no inchi is known or whose inchi is not parseable

identity_mapped_additional <- identity_mapped %>%
  bind_rows(
    {
      additional_inchis <- setdiff(
        input_data[["inchi_id_vendor_map"]][["inchi_id"]],
        .[["inchi_id"]]
      )
      tibble(
        inchi_id = additional_inchis,
        identity_group = seq_along(additional_inchis) + max(.[["identity_group"]])
      )
    }
  ) %>%
  # Assign lspci_ids where non has been carried over from the previous version
  mutate(
    lspci_id = {
      na_mask <- is.na(lspci_id)
      na_identity_groups <- identity_group[na_mask]
      new_lspci_ids <- c(
        0L,
        cumsum(
          head(na_identity_groups, -1L) != tail(na_identity_groups, -1L)
        )
      ) + as.integer(max(input_data[["old_sms"]][["lspci_id"]])) + 1L
      lspci_id[na_mask] <- new_lspci_ids
      lspci_id
    }
  )

fwrite(
  identity_mapped_additional,
  file.path(dir_release, "inchi_id_lspci_id_map.csv.gz")
)

cmpd_names_lspci_id <- cmpd_names %>%
  left_join(
    identity_mapped_additional %>%
      select(inchi_id, lspci_id),
    by = "inchi_id"
  )

fwrite(
  cmpd_names_lspci_id,
  file.path(dir_release, "lspci_id_compound_name_map.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Calculate compound equivalence classes",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/04_compound_equivalence.R"
)

syn_id_mapping <- Folder("id_mapping", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "lspci_id_compound_name_map.csv.gz"),
  file.path(dir_release, "inchi_id_lspci_id_map.csv.gz")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = activity)

