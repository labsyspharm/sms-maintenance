library(tidyverse)
library(dtplyr)
library(here)
library(synapser)
library(synExtra)
library(vroom)
library(httr)
library(lspcheminf)


synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

dir_emolecules <- file.path(dir_release, "emolecules")

dir.create(dir_emolecules, showWarnings = FALSE)

# Download Emolecules tables -------------------------------------------------
###############################################################################T

emolecules_base_url <- "https://downloads.emolecules.com/harvard/plus_data/with_prohibited/2020-11-01/"
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
  rename(smiles = isosmiles)

parent_molecules_inchi <- parent_molecules %>%
  chunk_df(10, seed = 1) %>%
  map(
    ~convert_compound_descriptor(
      compounds(
        set_names(.x[["smiles"]], .x[["parent_id"]]),
        descriptor = "smiles",
        keep_invalid = TRUE
      ),
      target_descriptor = "inchi"
    ) %>%
      transmute(parent_id = as.double(names), inchi = compounds)
  ) %>%
  bind_rows()

# write_rds(
#   parent_molecules_inchi,
#   file.path(dir_emolecules, "parent_molecules_inchi.csv.gz")
# )

parent_molecules_both <- parent_molecules %>%
  left_join(
    parent_molecules_inchi,
    by = "parent_id"
  )

# Filter compounds -------------------------------------------------------------
###############################################################################T

# Only keep building blocks and commercial compounds, avoid screening compounds
# because they are often only for purchase in large plate format

emolecules_files[["sample"]] %>%
  lazy_dt(key_by = c(supplier_id, catalog_id)) %>%
  group_by(supplier_id, catalog_id) %>%
  summarize(n = n()) %>%
  inner_join(
    emolecules_files[["catalog_tiers"]] %>%
      lazy_dt(key_by = c(supplier_id, catalogue_id)) %>%
      select(supplier_id, catalogue_id, tier_name),
    by = c("supplier_id", "catalog_id" = "catalogue_id")
  ) %>%
  group_by(tier_name) %>%
  summarize(n = sum(n)) %>%
  as_tibble()

# tier_name        n .groups
# <chr>        <int> <chr>
# 1 Tier 1     9358467 drop
# 2 Tier 2     7343877 drop
# 3 Tier 3    10875927 drop
# 4 Tier 4     8950202 drop
# 5 Tier 5      665524 drop

emolecules_files[["sample"]] %>%
  lazy_dt(key_by = c(catalog_id)) %>%
  group_by(catalog_id) %>%
  summarize(n = n()) %>%
  inner_join(
    emolecules_files[["catalog_categories"]] %>%
      lazy_dt(key_by = c(catalog_id)) %>%
      select(catalog_id, category_id),
    by = c("catalog_id")
  ) %>%
  group_by(category_id) %>%
  summarize(n = sum(n)) %>%
  inner_join(
    emolecules_files[["category"]] %>%
      lazy_dt(key_by = c(category_id)) %>%
      select(category_id, name, description),
    by = c("category_id")
  ) %>%
  as_tibble()

# # A tibble: 3 x 4
# category_id        n name                 description
# <dbl>    <int> <chr>                <chr>
# 1           2 12400298 Building blocks      Catalog of building blocks
# 2           3 24620788 Screening compounds  Catalog screening compounds
# 3           4 36965134 Commercial compounds Catalog of commercially-available compounds

filtered_compounds <- emolecules_files[["sample"]] %>%
  filter(
    catalog_id %in% intersect(
      emolecules_files[["catalog_tiers"]] %>%
        filter(tier_num %in% c(1, 2, 3)) %>%
        pull(catalogue_id),
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
#[1] 19692658

write_csv(
  filtered_compounds,
  file.path(dir_emolecules, "emolecules_vendor_info.csv.gz")
)

write_csv(
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

syn_emolecules <- synMkdir(syn_release, "id_mapping", "emolecules")

c(
  file.path(dir_emolecules, "emolecules_vendor_info.csv.gz"),
  file.path(dir_emolecules, "emolecules_compounds.csv.gz"),
  file.path(dir_emolecules, emolecules_file_names)
) %>%
  synStoreMany(syn_emolecules, activity = emolecules_activity)
