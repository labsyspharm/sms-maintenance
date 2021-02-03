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


# Set directories, import files ------------------------------------------------
###############################################################################T

inputs <- list(
  lspci_id_vendor_id_map = c("compounds_processed", "lspci_id_vendor_id_map.csv.gz"),
  chembl_biochemical = c("raw_data", "chembl", "chembl_biochemical_raw.csv.gz"),
  chembl_phenotypic = c("raw_data", "chembl", "chembl_phenotypic_raw.csv.gz"),
  chembl_biochemical = c("raw_data", "chembl", "chembl_biochemical_raw.csv.gz"),
  target_dictionary = c("id_mapping", "target_dictionary_wide.csv.gz"),
  inhouse_dose_response = c("raw_data", "hmsl", "hmsl_doseresponse.csv.gz"),
  inhouse_single_dose = c("raw_data", "hmsl", "hmsl_singledose.csv.gz")
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

# Clean up raw ChEMBL data -----------------------------------------------------
###############################################################################T

biochem_neat <- copy(input_data[["chembl_biochemical"]]) %>%
  mutate(standard_value = as.numeric(standard_value)) %>%
  rename(pref_name_target = pref_name) %>%
  filter(!is.na(entrez_gene_id) | !is.na(symbol)) %>%
  inner_join(
   copy(input_data[["lspci_id_vendor_id_map"]])[
     source == "chembl"
   ][
     ,
     source := NULL
   ],
   by = c("chembl_id_compound" = "vendor_id")
  ) %>%
  drop_na(lspci_id, entrez_gene_id)

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
  drop_na(lspci_id, entrez_gene_id)

pheno_data_neat <- copy(input_data[["chembl_phenotypic"]]) %>%
  mutate(standard_value = as.numeric(standard_value)) %>%
  inner_join(
    copy(input_data[["lspci_id_vendor_id_map"]])[
      source == "chembl"
    ][
      ,
      source := NULL
    ],
    by = c("chembl_id_compound" = "vendor_id")
  ) %>%
  drop_na(lspci_id, assay_id)

# Aggregate dose-response data -------------------------------------------------
###############################################################################T

biochem_rowbind <- biochem_neat %>%
  mutate(
    value = standard_value,
    value_unit = standard_units,
    references = map2(
      references, chembl_id_doc,
      ~{
        if (!is.null(.x))
          .x
        else
          data.table(
            reference_id = .y,
            reference_type = factor("chembl_id", levels = REFERENCE_PRIORITY)
          )
      }
    ),
    source = "chembl",
    source_id = chembl_id_assay
  ) %>%
  select(
    lspci_id,
    value, value_unit = standard_units, value_type = standard_type,
    value_relation = standard_relation, description_assay = description,
    uniprot_id, entrez_gene_id, symbol,
    references, source, source_id
  )

doseresponse_inhouse_rowbind <- dose_response_inhouse_neat %>%
  mutate(
    references = map2(
      reference_type, reference_id,
      ~data.table(
        reference_id = .y,
        reference_type = factor(.x, levels = REFERENCE_PRIORITY)
      )
    )
  ) %>%
  select(
    lspci_id,
    value, value_unit, value_type,
    value_relation, description_assay = description,
    uniprot_id, entrez_gene_id, symbol,
    references, file_url
  )

biochem_complete <- bind_rows(
  biochem_rowbind,
  doseresponse_inhouse_rowbind
) %>%
  # Remap obsolete entrez_ids
  # 645840 -> 114112
  # 348738 -> 6241
  mutate(
    entrez_gene_id = recode(
      entrez_gene_id, `645840` = 114112L, `348738` = 6241L
    )
  ) %>%
  # Call to distinct important, since some assays can be recorded multiple times
  # for the same eq_class now, when multiple forms of the same drug where mapped
  # to the same eq_class and an assay was stored in the db for all forms
  distinct()

qsave(
  biochem_complete,
  file.path(dir_release, "biochemical_data_complete.qs"),
  preset = "fast"
)

# Using data.table here for speed
calculate_q1 <- function(data) {
  n_groups <- uniqueN(data, by = c("lspci_id", "entrez_gene_id"))
  pb <- txtProgressBar(min = 1, max = n_groups, style = 3)
  on.exit(close(pb))
  data %>%
    as.data.table() %>% {
      .[
        ,
        .(
          Q1 = {
            setTxtProgressBar(pb, .GRP)
            round(quantile(value, 0.25, names = FALSE), 2)
          },
          n_measurement = .N,
          references = references %>%
            rbindlist(use.names = TRUE) %>%
            unique() %>% {
              paste(.[["reference_type"]], .[["reference_id"]], sep = ":")
            } %>%
            paste(collapse = "|")
        ),
        by = .(lspci_id, aggregate_column, aggregate_value)
      ]
    }
}

biochem_complete_q1 <- biochem_complete %>%
  mutate(
    aggregate_column = if_else(
      !is.na(entrez_gene_id) | is.na(symbol) | symbol == "-",
      "entrez_gene_id",
      "symbol"
    ),
    aggregate_value = if_else(
      aggregate_column == "entrez_gene_id",
      as.character(entrez_gene_id),
      symbol
    )
  ) %>%
  setkey(lspci_id, aggregate_column, aggregate_value) %>%
  calculate_q1()

fwrite(
  biochem_complete_q1,
  file.path(dir_release, "biochemical_data_complete_q1.csv.gz")
)

# Aggregate single-dose data ---------------------------------------------------
###############################################################################T


# Also calculate Q1 values for kinomescan data from HMS LINCS for which no
# complete dose response curve is available

hmsl_kinomescan_q1_func <- function(data) {
    n_groups <- uniqueN(data, by = c("lspci_id", "entrez_gene_id", "cmpd_conc_nM"))
    pb <- txtProgressBar(min = 1, max = n_groups, style = 3)
    on.exit(close(pb))
    data %>%
      as.data.table() %>% {
        .[
          ,
          .(
            percent_control_Q1 = {
              setTxtProgressBar(pb, .GRP)
              round(quantile(percent_control, 0.25, names = FALSE), 2)
            },
            n_measurement = .N,
            references = paste(
              reference_type,
              reference_id,
              sep = ":"
            ) %>%
              unique() %>%
              paste(collapse = "|")
          ),
          by = .(lspci_id, entrez_gene_id, cmpd_conc_nM)
        ]
      }
  }

hmsl_kinomescan_q1 <- input_data[["inhouse_single_dose"]] %>%
  distinct(
    lspci_id, entrez_gene_id, cmpd_conc_nM,
    percent_control, reference_type, reference_id
  ) %>%
  hmsl_kinomescan_q1_func()

fwrite(
  hmsl_kinomescan_q1,
  file.path(dir_release, "inhouse_single_dose_data_complete_q1.csv.gz")
)

# Aggregate phenotypic data ----------------------------------------------------
###############################################################################T

pheno_data_formatted <- pheno_data_neat %>%
  inner_join(
    chembl_ref_info_best,
    by = "chembl_id_doc"
  ) %>%
  distinct(
    lspci_id, assay_id,
    value = standard_value, value_unit = standard_units, value_type = standard_type,
    value_relation = standard_relation, description_assay = description,
    references
  )

qsave(
  pheno_data_formatted,
  file.path(dir_release, "pheno_data.qs"),
  preset = "fast"
)

pheno_data_q1 <- pheno_data_formatted %>%
  setDT() %>%
  {
    .[
      ,
      .(
        value_Q1 = quantile(value, 0.25, names = FALSE),
        n_measurement = .N,
        references = references %>%
          rbindlist(use.names = TRUE) %>%
          unique() %>% {
            paste(.[["reference_type"]], .[["reference_id"]], sep = ":")
          } %>%
          paste(collapse = "|")
      ),
      by = .(lspci_id, assay_id)
    ]
  }

fwrite(
  pheno_data_q1,
  file.path(dir_release, "pheno_data_q1.csv.gz")
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
  file.path(dir_release, "biochemical_data_complete_q1.csv.gz"),
  file.path(dir_release, "inhouse_single_dose_data_complete_q1.csv.gz"),
  file.path(dir_release, "pheno_data_q1.csv.gz")
) %>%
  synStoreMany(parent = syn_aggregate, activity = aggregation_activity)
