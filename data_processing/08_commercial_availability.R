library(tidyverse)
library(synapser)
library(synExtra)
library(data.table)
library(here)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

inputs <- list(
  compound_dictionary = c("compounds_processed", "compound_dictionary.csv.gz"),
  lspci_id_vendor_id_map = c("compounds_processed", "lspci_id_vendor_id_map.csv.gz"),
  emolecules_vendor_info = c("id_mapping", "emolecules", "emolecules_vendor_info.csv.gz"),
  emolecules_supplier_info = c("id_mapping", "emolecules", "suppliers.tsv.gz")
) %>%
  map(~exec(synPluck, !!!c(syn_release, .x)))

input_data <- inputs %>%
  map(syn) %>%
  map(
    function(x)
      list(
        `.csv` = partial(
          fread,
          colClasses = c(
            lspci_id = "integer"
          )
        ),
        `.tsv` = fread,
        `.rds` = read_rds
      ) %>%
      magrittr::extract2(which(str_detect(x, fixed(names(.))))) %>%
      {.(x)}
  )


# Map Emolecules compounds to lspci_id -----------------------------------------
###############################################################################T

compound_commercial_info <- input_data[["lspci_id_vendor_id_map"]] %>%
  filter(source == "emolecules") %>%
  select(-source) %>%
  mutate(emolecules_id = as.integer(vendor_id)) %>%
  select(-vendor_id) %>%
  inner_join(
    input_data[["emolecules_vendor_info"]] %>%
      distinct(
        parent_id,
        supplier_id,
        catalog_number = compound_id
      ),
    by = c("emolecules_id" = "parent_id")
  ) %>%
  inner_join(
    input_data[["emolecules_supplier_info"]] %>%
      distinct(
        supplier_id,
        supplier_name,
        tier
      ),
    by = "supplier_id"
  )

fwrite(
  compound_commercial_info,
  file.path(dir_release, "lspci_id_compound_commercial_info.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

wrangle_activity <- Activity(
  name = "Wrangle commercial availability of compounds from Emolecules",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/08_commercial_availability.R"
)

syn_vendors <- Folder("vendors", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "lspci_id_compound_commercial_info.csv.gz")
) %>%
  synStoreMany(parentId = syn_vendors, activity = wrangle_activity)
