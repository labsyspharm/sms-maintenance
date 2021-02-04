# Script to calculate TAS vectors from all combined data sources

library(tidyverse)
library(data.table)
library(here)
library(vroom)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Set directories, import files ------------------------------------------------
###############################################################################T

inputs <- c(
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
  set_names(
    str_replace(., fixed(".csv.gz"), "")
  ) %>%
  map(~c("aggregate_data", .x)) %>%
  c(
    list(
      lspci_id_vendor_id_map = c("compounds_processed", "lspci_id_vendor_id_map.csv.gz")
    )
  ) %>%
    map(~exec(synPluck, !!!c(syn_release, .x))) %>%
    c(
      literature_annotations = "syn20694521"
    )

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

# Wrangle literature annotations -----------------------------------------------
###############################################################################T

literature_annotations <- input_data[["literature_annotations"]] %>%
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
  mutate(
    tas = 2L,
    evidence = "manual_curation",
    references = "synapse:syn20694521"
  ) %>%
  # There are some duplicates in the literature annotation file
  # plus, in one case, two HMSL ids map to the same eq_class
  distinct()

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
    measurement = as.numeric(sigfig(Q1)),
    tas_id = seq_len(n())
  )

single_dose_tas <- input_data[["single_dose_q1"]] %>%
  mutate(
    tas = case_when(
      cmpd_conc_nM == 10000 & percent_control_Q1 >= 75 ~ 10L,
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
    by = .(lspci_id, entrez_gene_id)
  ][n > 1]

# 43 cases
# Thre are a few cases where there is contradicting info.
# In this case, take the minimum TAS. According to Nienke false positives
# are less likely than false negatives so this seems like the prudent approach.

single_dose_tas_agg <- single_dose_tas[
  order(tas, cmpd_conc_nM),
  .SD[
    tas == tas[1],
    .(
      tas,
      unit = paste0("% at ", cmpd_conc_nM, " nM"),
      measurement = as.numeric(sigfig(percent_control_Q1)),
      cmpd_conc_nM
    )
  ],
  by = .(lspci_id, entrez_gene_id)
][
  ,
  tas_id := .GRP,
  by = .(lspci_id, entrez_gene_id)
]

combined_tas <- rbindlist(
  list(
    dose_response = dose_response_tas,
    single_dose = single_dose_tas_agg,
    literature_annotation = copy(literature_annotations)[
      ,
      tas_id := seq_len(.N)
    ]
  ),
  idcol = "source",
  use.names = TRUE,
  fill = TRUE
)[
  ,
  .(lspci_id, entrez_gene_id, symbol, tas, unit, measurement, n_measurement, source, tas_id)
] %>%
  unique()

TAS_SOURCE_ORDER <- c(
  "dose_response",
  "single_dose",
  "literature_annotation"
)

# Prefer results from full affinity measuremennts over the percent inhibition
# When incorporating Verena's manual annotatations, prefer affinity data.
# If affinity data not present, take minimum of of Verena + percent control
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
  keyby = .(lspci_id, entrez_gene_id, symbol)
][
  ,
  source_collapsed := NULL
]

# TAS measurement and reference mapping ----------------------------------------
###############################################################################T


# Store to synapse -------------------------------------------------------------
###############################################################################T

tas_activity <- Activity(
  name = "Calculate target affinity spectrum (TAS) vectors",
  used = c(
    "syn20830834",
    "syn20830836",
    "syn20830516",
    "syn20694521",
    "syn20693721",
    "syn20835543"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/03_calculate_tas_vectors.R"
)

syn_tas_folder <- Folder("tas", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "tas_vector_sources.rds"),
  file.path(dir_release, "tas_vector.rds"),
  file.path(dir_release, "tas_vector_annotated.rds"),
  file.path(dir_release, "tas_vector_annotated_long.rds"),
  file.path(dir_release, "tas_vector_annotated_long.csv.gz")
) %>%
  synStoreMany(parentId = syn_tas_folder, activity = tas_activity)
