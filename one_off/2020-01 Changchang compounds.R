library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

wd <- here("one_off")
dir.create(wd, showWarnings = FALSE)

# Set directories, import files ------------------------------------------------
###############################################################################T

inputs <- list(
  compound_dictionary = c("compounds_processed", "compound_dictionary.csv.gz"),
  lspci_id_vendor_id_map = c("vendors", "lspci_id_compound_commercial_info.csv.gz")
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

# Create table of commercially available compounds -----------------------------
###############################################################################T

available_compounds <- input_data[["compound_dictionary"]][
  commercially_available == TRUE,
  .(
    lspci_id, pref_name, chembl_id, hmsl_id, emolecules_id, inchi_id, canonical_inchi
  )
][
  input_data[["lspci_id_vendor_id_map"]][
    ,
    .(
      lspci_id, supplier_name, catalog_number, tier
    )
  ],
  on = "lspci_id",
  nomatch = NULL
]

fwrite(
  available_compounds,
  file.path(wd, "commercial_compounds_changchang.csv.gz")
)


# Store to synapse -------------------------------------------------------------
###############################################################################T

aggregation_activity <- Activity(
  name = "Wrangle compound commercial info for Changchang",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/one_off/2020-01 Changchang compounds.R"
)

syn_aggregate <- synMkdir(syn_release, "one_off")

c(
  file.path(wd, "commercial_compounds_changchang.csv.gz")
) %>%
  synStoreMany(parent = syn_aggregate, activity = aggregation_activity)

