library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(vroom)
library(httr)
library(lspcheminf)
library(furrr)
library(qs)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v29"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

dir_emolecules <- file.path(dir_release, "emolecules")

dir.create(dir_emolecules, showWarnings = FALSE)

# Download Emolecules tables -------------------------------------------------
###############################################################################T

emolecules_base_url <- "https://downloads.emolecules.com/harvard/plus_data/with_prohibited/2021-07-01/"
emolecules_file_names <- c(
  "sample.tsv.gz",
  "suppliers.tsv.gz",
  "parent.smi.gz",
  "catalog_categories.tsv.gz",
  "catalog_tiers.tsv.gz",
  "category.tsv.gz"
)

walk(
  emolecules_file_names,
  ~GET(
    paste0(emolecules_base_url, .x),
    authenticate(Sys.getenv("EMOLECULES_USERNAME"), Sys.getenv("EMOLECULES_PASSWORD")),
    write_disk(file.path(dir_emolecules, .x))
  )
)

emolecules_files <- emolecules_file_names %>%
  keep(str_ends, fixed("tsv.gz")) %>%
  set_names(str_split_fixed(., fixed("."), 2)[, 1]) %>%
  map_chr(~file.path(dir_emolecules, .x)) %>%
  map(fread)


# Process structures -----------------------------------------------------------
###############################################################################T

parent_molecules <- fread(
  file.path(dir_emolecules, "parent.smi.gz"), sep = " "
) %>%
  setnames("isosmiles", "smiles")

plan(multicore(workers = 6))
parent_molecules_inchi <- parent_molecules %>%
  chunk_df(100, seed = 1) %>%
  future_map(
    ~convert_compound_descriptor(
      compounds(
        set_names(.x[["smiles"]], .x[["parent_id"]]),
        descriptor = "smiles"
      ),
      target_descriptor = "inchi"
    ) %>%
      transmute(parent_id = as.double(names), inchi = compounds),
    .progress = TRUE
  ) %>%
  bind_rows()

fwrite(
  parent_molecules_inchi,
  file.path(dir_emolecules, "parent_molecules_inchi.csv.gz")
)

parent_molecules_both <- parent_molecules %>%
  left_join(
    parent_molecules_inchi,
    by = "parent_id"
  )

# Filter compounds -------------------------------------------------------------
###############################################################################T

# Only keep building blocks and commercial compounds, avoid screening compounds
# because they are often only for purchase in large plate format

emolecules_files[["sample"]][
  , .(n = .N),
  keyby = .(supplier_id, catalog_id)
][
  emolecules_files[["catalog_tiers"]][
    , .(supplier_id, catalog_id, tier_name)
  ]
][
  , .(n = sum(n)),
  keyby = .(tier_name)
]

# tier_name        n
# 1:    Tier 1 10225528
# 2:    Tier 2  7360614
# 3:    Tier 3 15971147
# 4:    Tier 4  8804148
# 5:    Tier 5   400161

emolecules_files[["sample"]][
  , .(n = .N),
  keyby = .(catalog_id)
][
  emolecules_files[["catalog_categories"]][
    , .(catalog_id, category_id)
  ]
][
  , .(n = sum(n)),
  keyby = .(category_id)
][
  emolecules_files[["category"]][
    , .(category_id, name, description)
  ], nomatch = NULL
]

# category_id        n                 name                                 description
# 1:           2 12677651      Building blocks                  Catalog of building blocks
# 2:           3 29877623  Screening compounds                 Catalog screening compounds
# 3:           4 42499962 Commercial compounds Catalog of commercially-available compounds

filtered_compounds <- emolecules_files[["sample"]] %>%
  filter(
    catalog_id %in% intersect(
      emolecules_files[["catalog_tiers"]] %>%
        filter(tier_num %in% c(1, 2, 3)) %>%
        pull(catalog_id),
      emolecules_files[["category"]] %>%
        filter(name %in% c("Building blocks", "Commercial compounds")) %>%
        inner_join(
          emolecules_files[["catalog_categories"]],
          by = "category_id"
        ) %>%
        pull(catalog_id)
    )
  )

filtered_compounds$parent_id %>% unique() %>% length()
# [1] 25172847

fwrite(
  filtered_compounds,
  file.path(dir_emolecules, "emolecules_vendor_info.csv.gz")
)

fwrite(
  parent_molecules_both %>%
    filter(parent_id %in% filtered_compounds[["parent_id"]]),
  file.path(dir_emolecules, "emolecules_compounds.csv.gz")
)


# Store to synapse -------------------------------------------------------------
###############################################################################T

emolecules_activity <- Activity(
  name = "Fetch emolecules compound data",
  used = emolecules_base_url,
  executed = "https://github.com/labsyspharm/small-molecule-suite-maintenance/blob/master/id_mapping/01_emolecules.R"
)

syn_emolecules <- synMkdir(syn_release, "id_mapping", "emolecules", .recursive = TRUE)

c(
  file.path(dir_emolecules, "emolecules_vendor_info.csv.gz"),
  file.path(dir_emolecules, "emolecules_compounds.csv.gz"),
  file.path(dir_emolecules, emolecules_file_names)
) %>%
  synStoreMany(syn_emolecules, activity = emolecules_activity, forceVersion = FALSE)
