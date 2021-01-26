library(tidyverse)
library(data.table)
library(bit64)
library(here)
library(synapser)
library(synExtra)
library(lspcheminf)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Loading files ----------------------------------------------------------------
###############################################################################T

inputs <- list(
  commercial_info = c("id_mapping", "emolecules", "emolecules_vendor_info.csv.gz"),
  inchi_id_lspci_id_map = c("id_mapping", "inchi_id_lspci_id_map.csv.gz"),
  lspci_id_compound_name_map = c("id_mapping", "lspci_id_compound_name_map.csv.gz"),
  canonical_inchis = c("canonicalization", "canonical_inchi_ids.csv.gz"),
  inchi_id_vendor_map = c("canonicalization", "inchi_id_vendor_map.csv.gz"),
  chembl_raw = c("raw_data", "chembl_compounds_raw.rds"),
  inchis = c("canonicalization", "canonical_inchi_ids.csv.gz"),
  fda_approval = c("raw_data", "lsp_FDA_first_approval_table.csv")
) %>%
  map(~exec(synPluck, !!!c(syn_release, .x)))

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

# Finding canonical member of equivalence class --------------------------------
###############################################################################T

# Find canonical member of equivalence class
# 1. Choose member with highest clinical phase
# 2. Choose member that is commercially available
# 3. Choose member with annotated pref_name
# 4. If present, choose member annotated as parent by Chembl
# 5. Choose member with highest n_assays
# 6. Choose member with lowest id (chembl or HMS)

VENDOR_RANKING <- c(
  "chembl", "hmsl", "emolecules", "old_sms"
)

canonical_members_ranked <- input_data[["inchi_id_lspci_id_map"]] %>%
  left_join(
    input_data[["inchi_id_vendor_map"]],
    by = "inchi_id"
  ) %>%
  left_join(
    input_data[["chembl_raw"]] %>%
      select(chembl_id, molregno, parent_molregno, n_assays, parental_flag, max_phase),
    by = c("vendor_id" = "chembl_id")
  ) %>%
  setDT() %>%
  copy() %>% {
    .[
      ,
      `:=`(
        source = factor(source, levels = VENDOR_RANKING),
        has_parent = if_else(parent_molregno != molregno, parent_molregno, NA_integer64_),
        annotated_pref_name = vendor_id %in% input_data[["lspci_id_compound_name_map"]][
          name_preference == "primary"
        ][["id"]],
        id_number = as.integer(str_extract(vendor_id, "\\d+")),
        commercially_available = inchi_id %in% input_data[["inchi_id_vendor_map"]][
          source == "emolecules"
        ][["inchi_id"]]
      )
    ][
      ,
      `:=`(
        annotated_as_parent = molregno %in% has_parent
      )
    ][
      order(
        lspci_id,
        source,
        -max_phase,
        -commercially_available,
        -annotated_pref_name,
        -annotated_as_parent,
        -n_assays,
        id_number
      )
    ][
      ,
      rank := seq_len(.N),
      keyby = "lspci_id"
    ]
  }


fwrite(
  canonical_members_ranked,
  file.path(dir_release, "lspci_id_canonical_members_ranked.csv.gz")
)
# canonical_members_ranked <- fread(file.path(dir_release, "lspci_id_canonical_members_ranked.csv.gz"))

canonical_members <- canonical_members_ranked[
  rank == 1L
]

fwrite(
  canonical_members,
  file.path(dir_release, "lspci_id_canonical_members.csv.gz")
)

# Finding canonical names ------------------------------------------------------
###############################################################################T

all_names_ranked <- copy(input_data[["lspci_id_compound_name_map"]])[
  ,
  `:=`(
    source = factor(source, levels = VENDOR_RANKING),
    canonical_member = lspci_id %in% canonical_members[["lspci_id"]],
    inchi_id = NULL
  )
] %>%
  unique() %>% {
    .[
      order(
        -canonical_member,
        source,
        match(name_preference, c("primary", "secondary")),
        str_count(name, fixed(" ")),
        str_length(name)
      )
    ][
      ,
      name_rank := seq_len(.N),
      keyby = lspci_id
    ]
  }

fwrite(
  all_names_ranked,
  file.path(dir_release, "lspci_id_compound_compound_names_ranked.csv.gz")
)

# all_names_ranked <- fread(file.path(dir_release, "lspci_id_compound_compound_names_ranked.csv.gz"))

canonical_names <- all_names_ranked[
  name_rank == 1L
][, name_rank := NULL]

fwrite(
  canonical_names,
  file.path(dir_release, "lspci_id_canonical_compound_names.csv.gz")
)

# Finding canonical Inchis -----------------------------------------------------
###############################################################################T

canonical_inchis_ranked <- copy(canonical_members_ranked)[
  ,
  source := factor(source, levels = VENDOR_RANKING)
][
  ,
  {
    inchi_ids <- unique(inchi_id)
    ranks <- seq_along(inchi_ids)
    list(
      inchi_id = inchi_ids,
      rank = ranks
    )
  },
  keyby = .(lspci_id, source)
][
  input_data[["inchis"]],
  on = "inchi_id"
] %>%
  setkey(lspci_id, source, rank)

fwrite(
  canonical_inchis_ranked,
  file.path(dir_release, "lspci_id_canonical_inchis_ranked.csv.gz")
)

# Create table of approval data ------------------------------------------------
###############################################################################T

approval_table <- copy(canonical_members_ranked)[
  ,
  `:=`(
    max_phase = if_else(
      lspci_id %in% input_data[["fda_approval"]][["lspci_id"]],
      4L,
      nafill(max_phase, fill = 0L)
    )
  )
][
  ,
  .(
    max_phase = max(max_phase)
  ),
  keyby = "lspci_id"
]


# Create table of canonical compounds ------------------------------------------
###############################################################################T

pb <- txtProgressBar(
  max = uniqueN(canonical_members_ranked[["lspci_id"]]),
  style = 3
)
compound_dictionary <- copy(canonical_members_ranked)[
  rank == 1L
][
  ,
  c(
    {
      setTxtProgressBar(pb, .GRP)
      list(
        inchi_id = head(inchi_id, n = 1)
      )
    },
    vendor_id[
      match(VENDOR_RANKING, source)
    ] %>%
      set_names(paste0(VENDOR_RANKING, "_id"))
  ),
  keyby = "lspci_id"
][
  ,
  commercially_available := !is.na(emolecules_id)
] %>%
  merge(
    input_data[["inchis"]],
    by = "inchi_id",
    all.x = TRUE
  ) %>%
  merge(
    all_names_ranked[
      name_rank == 1L,
      .(
        lspci_id,
        pref_name = name
      )
    ],
    by = "lspci_id",
    all.x = TRUE
  ) %>%
  merge(
    approval_table,
    by = "lspci_id",
    all.x = TRUE
  )
close(pb)

fwrite(
  compound_dictionary,
  file.path(dir_release, "compound_dictionary.csv.gz")
)

# compound_dictionary <- fread(file.path(dir_release, "compound_dictionary.csv.gz"))

# Create map of vendor IDs to lspci_ids ----------------------------------------
###############################################################################T

lspci_id_vendor_id_map <- merge(
  input_data[["inchi_id_lspci_id_map"]][
    , .(lspci_id, inchi_id)
  ] %>%
    na.omit(),
  input_data[["inchi_id_vendor_map"]] %>%
    na.omit(),
  by = "inchi_id",
  all = FALSE
)[
  ,
  inchi_id := NULL
] %>%
  setkey(lspci_id, source)

fwrite(
  lspci_id_vendor_id_map,
  file.path(dir_release, "lspci_id_vendor_id_map.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

cmpd_wrangling_activity <- Activity(
  name = "Map and wrangle canonical compound data.",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/05_process_mapped_compounds.R"
)

syn_id_mapping <- synMkdir(syn_release, "compounds_processed")

c(
  file.path(dir_release, "lspci_id_canonical_members_ranked.csv.gz"),
  file.path(dir_release, "lspci_id_canonical_members.csv.gz"),
  file.path(dir_release, "lspci_id_compound_compound_names_ranked.csv.gz"),
  file.path(dir_release, "lspci_id_canonical_compound_names.csv.gz"),
  file.path(dir_release, "lspci_id_canonical_inchis_ranked.csv.gz"),
  file.path(dir_release, "compound_dictionary.csv.gz"),
  file.path(dir_release, "lspci_id_vendor_id_map.csv.gz")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = cmpd_wrangling_activity)

