library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(qs)

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
  chembl_biochemical = c("raw_data", "chembl", "chembl_biochemical_raw.csv.gz"),
  chembl_phenotypic = c("raw_data", "chembl", "chembl_phenotypic_raw.csv.gz"),
  chembl_biochemical = c("raw_data", "chembl", "chembl_biochemical_raw.csv.gz"),
  inhouse_dose_response = c("raw_data", "hmsl", "hmsl_doseresponse.csv.gz"),
  inhouse_single_dose = c("raw_data", "hmsl", "hmsl_singledose.csv.gz"),
  chembl_references_best = c("raw_data", "chembl", "chembl_ref_info_best_source.csv.gz"),
  references = c("references", "references.csv.gz")
) %>%
  pluck_inputs(syn_parent = syn_release)

input_data <- inputs %>%
  load_input_data(syn = syn)

# Clean up raw ChEMBL data -----------------------------------------------------
###############################################################################T

biochem_neat <- copy(input_data[["chembl_biochemical"]]) %>%
  mutate(standard_value = as.numeric(standard_value)) %>%
  rename(pref_name_target = pref_name) %>%
  inner_join(
   copy(input_data[["lspci_id_vendor_id_map"]])[
     source == "chembl"
   ][
     ,
     source := NULL
   ],
   by = c("chembl_id_compound" = "vendor_id")
  ) %>%
  drop_na(lspci_id, lspci_target_id)

dose_response_inhouse_neat <- copy(input_data[["inhouse_dose_response"]]) %>%
  inner_join(
    copy(input_data[["lspci_id_vendor_id_map"]])[
      source == "hmsl"
    ][
      ,
      source := NULL
    ],
    by = c("hmsl_id" = "vendor_id")
  ) %>%
  drop_na(lspci_id, lspci_target_id)

# Aggregate dose-response data -------------------------------------------------
###############################################################################T

biochem_rowbind <- biochem_neat %>%
  left_join(
    input_data[["chembl_references_best"]],
    by = "chembl_id_doc"
  ) %>%
  transmute(
    lspci_id, lspci_target_id,
    entrez_gene_id, symbol,
    value = standard_value,
    value_unit = standard_units,
    value_type = standard_type,
    value_relation = standard_relation,
    description_assay = description,
    reference_type, reference_id = reference_value,
    measurement_source = "chembl_activity",
    external_measurement_id = as.character(activity_id)
  )

doseresponse_inhouse_rowbind <- dose_response_inhouse_neat %>%
  transmute(
    lspci_id, lspci_target_id,
    entrez_gene_id, symbol,
    value, value_unit, value_type,
    value_relation, description_assay = description,
    reference_type, reference_id = reference_value,
    measurement_source = "inhouse_doseresponse",
    external_measurement_id = synapse_id
  )

biochem_complete <- bind_rows(
  biochem_rowbind,
  doseresponse_inhouse_rowbind
) %>%
  # Call to distinct important, since some assays can be recorded multiple times
  # for the same eq_class now, when multiple forms of the same drug where mapped
  # to the same eq_class and an assay was stored in the db for all forms
  distinct() %>%
  # Needed to link aggregated data back to specific measurements in DB
  mutate(
    measurement_id = seq_len(n())
  )

fwrite(
  biochem_complete,
  file.path(dir_release, "dose_response_measurements.csv.gz")
)

# biochem_complete <- fread(file.path(dir_release, "dose_response_measurements.csv.gz"))

# Using data.table here for speed

biochem_complete_q1 <- biochem_complete[
  ,
  .(
    Q1 = round(quantile(value, 0.25, names = FALSE), 2),
    n_measurement = .N
  ),
  by = .(lspci_id, lspci_target_id)
]

biochem_complete_q1_refs <- biochem_complete[
  ,
  .(
    lspci_id,
    lspci_target_id,
    reference_type,
    reference_value = reference_id
  )
] %>%
  left_join(
    input_data[["references"]][
      , .(reference_id, reference_type, reference_value)
    ],
    by = c("reference_type", "reference_value")
  ) %>%
  setkey(lspci_id, lspci_target_id) %>%
  unique()

biochem_complete_q1_measurements <- biochem_complete[
  ,
  .(
    lspci_id,
    lspci_target_id,
    measurement_source,
    external_measurement_id,
    measurement_id
  )
] %>%
  setkey() %>%
  unique()

fwrite(
  biochem_complete_q1,
  file.path(dir_release, "dose_response_q1.csv.gz")
)

fwrite(
  biochem_complete_q1_refs,
  file.path(dir_release, "dose_response_q1_references.csv.gz")
)

fwrite(
  biochem_complete_q1_measurements,
  file.path(dir_release, "dose_response_q1_measurements.csv.gz")
)

# Aggregate single-dose data ---------------------------------------------------
###############################################################################T

single_dose_data <- input_data[["inhouse_single_dose"]] %>%
  distinct(
    lspci_id, lspci_target_id, entrez_gene_id, symbol, cmpd_conc_nM,
    percent_control, reference_type, reference_id = reference_value,
  ) %>%
  mutate(
    measurement_source = "inhouse_single_dose",
    measurement_id = seq_len(n())
  )

fwrite(
  single_dose_data,
  file.path(dir_release, "single_dose_measurements.csv.gz")
)

# Also calculate Q1 values for kinomescan data from HMS LINCS for which no
# complete dose response curve is available
hmsl_kinomescan_q1 <- single_dose_data[
  ,
  .(
    percent_control_Q1 = round(quantile(percent_control, 0.25, names = FALSE), 2),
    n_measurement = .N
  ),
  by = .(lspci_id, lspci_target_id, cmpd_conc_nM)
]

hmsl_kinomescan_q1_refs <- single_dose_data[
  ,
  .(
    lspci_id,
    entrez_gene_id,
    symbol,
    cmpd_conc_nM,
    reference_type,
    reference_value = reference_id
  )
] %>%
  left_join(
    input_data[["references"]][
      , .(reference_id, reference_type, reference_value)
    ],
    by = c("reference_type", "reference_value")
  ) %>%
  setkey() %>%
  unique()

hmsl_kinomescan_q1_measurements <- single_dose_data[
  ,
  .(
    lspci_id,
    lspci_target_id,
    cmpd_conc_nM,
    measurement_source,
    measurement_id
  )
] %>%
  setkey() %>%
  unique()

fwrite(
  hmsl_kinomescan_q1,
  file.path(dir_release, "single_dose_q1.csv.gz")
)

fwrite(
  hmsl_kinomescan_q1_refs,
  file.path(dir_release, "single_dose_q1_references.csv.gz")
)

fwrite(
  hmsl_kinomescan_q1_measurements,
  file.path(dir_release, "single_dose_q1_measurements.csv.gz")
)

# Aggregate phenotypic data ----------------------------------------------------
###############################################################################T

pheno_data_neat <- input_data[["chembl_phenotypic"]] %>%
  mutate(standard_value = as.numeric(standard_value)) %>%
  inner_join(
    copy(input_data[["lspci_id_vendor_id_map"]])[
      source == "chembl"
    ][
      ,
      source := NULL
    ],
    by = c("chembl_id_compound" = "vendor_id")
  )  %>%
  inner_join(
    input_data[["chembl_references_best"]],
    by = "chembl_id_doc"
  ) %>%
  transmute(
    lspci_id,
    assay_id,
    value = standard_value,
    value_unit = standard_units,
    value_type = standard_type,
    value_relation = standard_relation,
    description_assay = description,
    reference_type, reference_id = reference_value,
    measurement_source = "chembl_activity",
    external_measurement_id = as.character(activity_id)
  ) %>%
  distinct() %>%
  mutate(
    measurement_id = seq_len(n())
  )

fwrite(
  pheno_data_neat,
  file.path(dir_release, "phenotypic_measurements.csv.gz")
)

pheno_data_q1 <- pheno_data_neat %>%
  setDT() %>%
  {
    n_groups <- uniqueN(., by = c("lspci_id", "assay_id"))
    pb <- txtProgressBar(min = 1, max = n_groups, style = 3)
    .[
      , {
        setTxtProgressBar(pb, .GRP)
        .(
          value_Q1 = round(quantile(value, 0.25, names = FALSE), 2),
          n_measurement = .N
        )
      },
      by = .(lspci_id, assay_id)
    ]
  }

pheno_data_q1_refs <- pheno_data_neat[
  ,
  .(
    lspci_id,
    assay_id,
    reference_type,
    reference_value = reference_id
  )
] %>%
  left_join(
    input_data[["references"]][
      , .(reference_id, reference_type, reference_value)
    ],
    by = c("reference_type", "reference_value")
  ) %>%
  setkey() %>%
  unique()

pheno_data_q1_measurements <- pheno_data_neat[
  ,
  .(
    lspci_id,
    assay_id,
    measurement_source,
    external_measurement_id,
    measurement_id
  )
] %>%
  setkey() %>%
  unique()


fwrite(
  pheno_data_q1,
  file.path(dir_release, "phenotypic_q1.csv.gz")
)

fwrite(
  pheno_data_q1_refs,
  file.path(dir_release, "phenotypic_q1_references.csv.gz")
)

fwrite(
  pheno_data_q1_measurements,
  file.path(dir_release, "phenotypic_q1_measurements.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

aggregation_activity <- Activity(
  name = "Aggregate affinity data",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/02_aggregate_drug_affinities.R"
)

syn_aggregate <- synMkdir(syn_release, "aggregate_data")

c(
  "dose_response",
  "single_dose",
  "phenotypic"
) %>%
  interaction(
    c(
      "_measurements.csv.gz",
      "_q1.csv.gz",
      "_q1_references.csv.gz",
      "_q1_measurements.csv.gz"
    ),
    sep = ""
  ) %>%
  levels() %>%
  {file.path(dir_release, .)} %>%
  synStoreMany(parent = syn_aggregate, activity = aggregation_activity)
