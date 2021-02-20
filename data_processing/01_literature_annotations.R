library(tidyverse)
library(synapser)
library(synExtra)
library(here)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

source(here("utils", "load_save.R"))

# Set directories, import files ------------------------------------------------
###############################################################################T

inputs <- list(
  lspci_id_vendor_id_map = c("compounds_processed", "lspci_id_vendor_id_map.csv.gz"),
  target_map = c("id_mapping", "target_mapping.csv.gz")
) %>%
  pluck_inputs(syn_parent = syn_release) %>%
  c(
    literature_annotations_raw = "syn20694521"
  )

input_data <- inputs %>%
  load_input_data(syn = syn)

# Wrangle literature annotations -----------------------------------------------
###############################################################################T

literature_annotations <- input_data[["literature_annotations_raw"]] %>%
  mutate(
    hmsl_id = paste0("HMSL", hms_id)
  ) %>%
  inner_join(
    input_data[["lspci_id_vendor_id_map"]][
      source == "hmsl",
      .(lspci_id, hmsl_id = vendor_id)
    ],
    by = c("hmsl_id")
  ) %>%
  select(lspci_id, entrez_gene_id = gene_id, symbol = gene_symbol) %>%
  left_join(
    input_data[["target_map"]] %>%
      distinct(lspci_target_id, entrez_gene_id),
    by = "entrez_gene_id"
  ) %>%
  mutate(
    tas = 2L,
    evidence = "literature_annotation",
    reference_type = "synapse_id",
    reference_value = "syn20694521"
  ) %>%
  # There are some duplicates in the literature annotation file
  # plus, in one case, two HMSL ids map to the same eq_class
  distinct() %>%
  mutate(
    measurement_id = seq_len(n())
  )

fwrite(
  literature_annotations,
  file.path(dir_release, "literature_annotations.csv.gz")
)

# Storing references -----------------------------------------------------------
###############################################################################T

REFERENCE_PRIORITY <- c(
  "pubmed_id" = "https://pubmed.ncbi.nlm.nih.gov/",
  "doi" = "http://doi.org/",
  "patent_id" = "https://patents.google.com/?q=",
  "synapse_id" = "https://www.synapse.org/#!Synapse:",
  "chembl_id" = "https://www.ebi.ac.uk/chembl/document_report_card/",
  "hmsl_id" = "https://lincs.hms.harvard.edu/db/datasets/"
)

syn_reference_table <- synPluck(syn_release, "reference_table")

current_references <- synTableQuery(sprintf("SELECT * FROM %s", syn_reference_table)) %>%
  as.data.frame()

new_references <- literature_annotations %>%
  anti_join(
    current_references,
    by = c("reference_type", "reference_value")
  ) %>%
  transmute(
    reference_type,
    reference_value,
    url = paste0(REFERENCE_PRIORITY[reference_type], reference_value)
  ) %>%
  distinct()

synStore(
  Table(
    syn_reference_table,
    new_references
  )
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Wrangle literature annotations",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/01_literature_annotations.R"
)

syn_folder <- synMkdir(syn_release, "raw_data", "literature_annotations")

c(
  file.path(dir_release, "literature_annotations.csv.gz")
) %>%
  synStoreMany(parentId = syn_folder, activity = activity, forceVersion = FALSE)

