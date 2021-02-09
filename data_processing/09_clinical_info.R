library(tidyverse)
library(synapser)
library(synExtra)
library(here)

synLogin()
syn <- synDownloader(here("tempdl"))

source(here("utils", "load_save.R"))

# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

inputs <- list(
  compound_dictionary = c("compounds_processed", "compound_dictionary.csv.gz"),
  approval_chembl = c("raw_data", "chembl", "chembl_approval_info_raw.csv.gz"),
  approval_fda = c("raw_data", "lsp_FDA_first_approval_table.csv"),
  lspci_id_vendor_id_map = c("compounds_processed", "lspci_id_vendor_id_map.csv.gz")
) %>%
  pluck_inputs(syn_parent = syn_release)

input_data <- inputs %>%
  load_input_data(syn = syn)

# Map IDs ----------------------------------------------------------------------
###############################################################################T

approval_table <- copy(input_data[["lspci_id_vendor_id_map"]])[
  source == "chembl",
  .(lspci_id, chembl_id = vendor_id)
][
  input_data[["approval_chembl"]][
    ,
    .(
      chembl_id = chembl_id_compound,
      max_phase,
      first_approval,
      oral, parenteral, topical,
      black_box_warning, first_in_class, prodrug, indication_class,
      withdrawn_flag, withdrawn_year, withdrawn_country, withdrawn_reason,
      max_phase_for_indication = max_phase_for_ind, mesh_heading, efo_id, efo_term, ref_type,
      ref_id, ref_url
    )
  ] %>%
    unique(),
  on = "chembl_id",
  nomatch = NULL
][
#   input_data[["approval_fda"]][
#     ,
#     .(
#       lspci_id, fda_application_number = application_number,
#       active_ingredient, current_marketing_status, application_type,
#       submission_status_year, submission_type, submission_class,
#       submission_class_description
#     )
#   ],
#   on = "lspci_id",
#   nomatch = NULL
# ][
  ,
  `:=`(
    max_phase = if_else(
      lspci_id %in% input_data[["fda_approval"]][["lspci_id"]],
      4L,
      max(nafill(max_phase, fill = 0L))
    ),
    max_phase_for_indication = max(nafill(max_phase_for_indication, fill = 0L))
  ),
  keyby = lspci_id
]

fwrite(
  approval_table,
  file.path(dir_release, "lspci_id_approval.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

clinical_info_activity <- Activity(
  name = "Map IDs of clinical information from Chembl",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/09_clinical_info.R"
)

syn_clinical_info <- synMkdir(syn_release, "clinical_info")

c(
  file.path(dir_release, "lspci_id_approval.csv.gz")
) %>%
  synStoreMany(parentId = syn_clinical_info, activity = clinical_info_activity)

