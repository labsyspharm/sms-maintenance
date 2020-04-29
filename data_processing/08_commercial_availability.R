library(tidyverse)
library(synapser)
library(synExtra)
library(here)

synLogin()
syn <- synDownloader(here("tempdl"))


# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

compound_mapping <- syn("syn20934414") %>%
  read_rds()

zinc_vendor_names <- syn("syn21901782") %>%
  read_csv()

zinc_availability <- syn("syn21901780") %>%
  read_csv()

# Combine vendor names with commercial data ------------------------------------
###############################################################################T

zinc_commercial <- zinc_availability %>%
  left_join(
    zinc_vendor_names %>%
      distinct(vendor, vendor_id, vendor_name),
    by = c("vendor", "vendor_id")
  )

# Map Zinc compounds to lspci_id -----------------------------------------------
###############################################################################T

compound_commercial_info <- compound_mapping %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        select(id, lspci_id) %>%
        inner_join(
          zinc_commercial,
          by = c("id" = "chembl_id")
        ) %>%
        distinct(lspci_id, vendor, vendor_id, vendor_name)
    )
  )

write_rds(
  compound_commercial_info,
  file.path(dir_release, "compound_commercial_info_zinc.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

wrangle_activity <- Activity(
  name = "Wrangle commercial availability of compounds from ZINC",
  used = c(
    "syn20934414",
    "syn21901782",
    "syn21901780"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/08_commercial_availability.R"
)

syn_vendors <- Folder("vendors", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "compound_commercial_info_zinc.rds")
) %>%
  synStoreMany(parentId = syn_vendors, activity = wrangle_activity)
