# Script to calculate TAS vectors from all combined data sources

library(tidyverse)
library(data.table)
library(here)
library(vroom)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"), ifcollision = "overwrite.local")

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

source(here("utils", "load_save.R"))

# Set directories, import files ------------------------------------------------
###############################################################################T

inputs <- c(
  "dose_response",
  "single_dose"
) %>%
  interaction(
    c(
      "_q1.csv.gz",
      "_q1_measurements.csv.gz"
    ),
    sep = ""
  ) %>%
  levels() %>%
  set_names(
    str_replace(., fixed(".csv.gz"), "")
  ) %>%
  map(~c("aggregate_data", .x)) %>%
  c(
    list(
      literature_annotations = c("raw_data", "literature_annotations", "literature_annotations.csv.gz"),
      target_dictionary = c("id_mapping", "target_dictionary_wide.csv.gz")
    )
  ) %>%
    pluck_inputs(syn_parent = syn_release)

input_data <- inputs %>%
  load_input_data(syn = syn)

# calculate TAS ----------------------------------------------------------------
###############################################################################T

# https://stackoverflow.com/a/3443955/4603385
sigfig <- function(vec, n = 3){
  as.numeric(gsub("\\.$", "", formatC(signif(vec, digits = n), digits = n,format = "fg", flag = "#")))
}

dose_response_tas <- input_data[["dose_response_q1"]] %>%
  mutate(
    Q1 = as.numeric(Q1),
    tas = case_when(
      Q1 < 100 ~ 1L,
      Q1 < 1000 ~ 2L,
      Q1 < 10000 ~ 3L,
      Q1 >= 10000 ~ 10L,
      TRUE ~ NA_integer_
    ),
    unit = "nM",
    measurement = as.numeric(sigfig(Q1))
  )

single_dose_tas <- input_data[["single_dose_q1"]] %>%
  mutate(
    tas = case_when(
      cmpd_conc_nM == 10000 & percent_control_Q1 >= 50 ~ 10L,
      cmpd_conc_nM == 10000 & percent_control_Q1 < 0.1 ~ 2L,
      cmpd_conc_nM == 1000 & percent_control_Q1 >= 90 ~ 10L,
      cmpd_conc_nM == 1000 & percent_control_Q1 < 1 ~ 2L,
      cmpd_conc_nM == 100 & percent_control_Q1 >= 75 ~ 10L,
      cmpd_conc_nM == 100 & percent_control_Q1 < 25 ~ 2L,
      TRUE ~ NA_integer_
    )
  )


# Check for contradicting TAS at different compound concentrations
single_dose_tas[
    ,
    .(n = length(unique(na.omit(tas)))),
    by = .(lspci_id, lspci_target_id)
  ][n > 1]

# 43 cases
# Thre are a few cases where there is contradicting info.
# In this case, take the minimum TAS. According to Nienke false positives
# are less likely than false negatives so this seems like the prudent approach.
# Only retain a single measurement as source for TAS. In case of ties
# measurement at lowest concentration.

single_dose_tas_agg <- single_dose_tas[
  order(tas, cmpd_conc_nM),
  .SD[
    1L,
    .(
      tas,
      unit = paste0("% at ", cmpd_conc_nM, " nM"),
      measurement = as.numeric(sigfig(percent_control_Q1)),
      cmpd_conc_nM
    )
  ],
  by = .(lspci_id, lspci_target_id)
]

combined_tas <- rbindlist(
  list(
    dose_response = dose_response_tas,
    single_dose = single_dose_tas_agg,
    literature_annotation = input_data[["literature_annotations"]]
  ),
  idcol = "source",
  use.names = TRUE,
  fill = TRUE
)[
  ,
  .(lspci_id, lspci_target_id, tas, unit, measurement, n_measurement, source, cmpd_conc_nM)
] %>%
  unique()

TAS_SOURCE_ORDER <- c(
  "dose_response",
  "single_dose",
  "literature_annotation"
)

# Prefer results from full affinity measurements over the percent inhibition
# When incorporating Verena's manual annotatations, prefer affinity data.
# Only if affinity data not present, take minimum of of Verena + percent control
combined_tas_agg <- copy(combined_tas)[
  ,
  `:=`(
    source_collapsed = factor(
      source,
      TAS_SOURCE_ORDER
    ) %>%
      fct_collapse(single_dose_literature = c("single_dose", "literature_annotation")) %>%
      fct_relevel("dose_response")
  )
][
  order(source_collapsed, tas),
  .SD[1],
  keyby = .(lspci_id, lspci_target_id)
][
  ,
  `:=`(
    source_collapsed = NULL,
    tas_id = seq_len(.N)
  )
][
  ,
  .(
    tas_id, lspci_id, lspci_target_id, tas, unit, measurement, n_measurement, source, cmpd_conc_nM
  )
] %>%
  left_join(
    input_data[["target_dictionary"]][
      , .(lspci_target_id, entrez_gene_id, symbol)
    ],
    by = "lspci_target_id"
  )

# Check there's only a single value for every compound / target combo
combined_tas_agg %>%
  count(lspci_id, lspci_target_id) %>%
  count(n)

fwrite(
  combined_tas_agg %>%
    select(-cmpd_conc_nM),
  file.path(dir_release, "tas_vector.csv.gz")
)

# TAS measurement and reference mapping ----------------------------------------
###############################################################################T

all_q1_measurements <-   c(
  dose_response = "dose_response_q1_measurements",
  single_dose = "single_dose_q1_measurements",
  literature_annotation = "literature_annotations"
) %>%
  map(~input_data[[.x]]) %>%
  rbindlist(fill = TRUE, idcol = "source")

tas_measurement_map <- all_q1_measurements[
  combined_tas_agg[
    , .(lspci_id, lspci_target_id, source, tas_id, cmpd_conc_nM)
  ],
  nomatch = NULL, on = c("lspci_id", "lspci_target_id", "source", "cmpd_conc_nM")
][, .(source, tas_id, measurement_id)]

# Check that all tas values have a measurement associated with it
combined_tas_agg[
  !tas_id %in% tas_measurement_map$tas_id
]

# Check that there's no NAs
drop_na(tas_measurement_map)

fwrite(
  tas_measurement_map,
  file.path(dir_release, "tas_measurement_map.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

tas_activity <- Activity(
  name = "Calculate target affinity spectrum (TAS) vectors",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/03_calculate_tas_vectors.R"
)

syn_tas_folder <- synMkdir(syn_release, "tas")

c(
  file.path(dir_release, "tas_vector.csv.gz"),
  file.path(dir_release, "tas_measurement_map.csv.gz")
) %>%
  synStoreMany(parentId = syn_tas_folder, activity = tas_activity, forceVersion = FALSE)
