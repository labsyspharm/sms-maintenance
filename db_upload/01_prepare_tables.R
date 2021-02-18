library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(safejoin)

synLogin()
syn <- synDownloader(here("tempdl"))

source(here("utils", "load_save.R"))

# set directories, import files ------------------------------------------------
###############################################################################T

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# read tables ------------------------------------------------------------------
###############################################################################T

# Function to replace empty strings with NA values
replace_empty_string_na <- function(df) {
  mutate(
    df,
    across(
      where(is.character),
      ~magrittr::inset(.x, .x == "", NA_character_)
    )
  )
}

inputs <- list(
  target_dictionary = c("id_mapping", "target_dictionary_wide.csv.gz"),
  compound_dictionary = c("compounds_processed", "compound_dictionary.csv.gz"),
  fingerprints = c("compounds_processed", "lspci_id_fingerprints.csv.gz"),
  selectivity = c("selectivity", "selectivity.csv.gz"),
  approval = c("clinical_info", "lspci_id_approval.csv.gz"),
  compound_names = c("compounds_processed", "lspci_id_compound_compound_names_ranked.csv.gz"),
  vendor_ids = c("compounds_processed", "lspci_id_vendor_id_map.csv.gz"),
  optimal_library = c("compound_library", "optimal_library_combined.csv.gz"),
  inchis = c("compounds_processed", "lspci_id_canonical_inchis_ranked.csv.gz"),
  tas = c("tas", "tas_vector.csv.gz"),
  tas_measurement_map = c("tas", "tas_measurement_map.csv.gz"),
  commercial_info = c("vendors", "lspci_id_compound_commercial_info.csv.gz"),
  manual_curation = c("raw_data", "literature_annotations", "literature_annotations.csv.gz"),
  references = c("reference_table")
) %>%
  c(
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
      set_names(str_replace(., fixed(".csv.gz"), "")) %>%
      map(~exec(c, "aggregate_data", !!!.x))
  ) %>%
  pluck_inputs(syn_parent = syn_release)

input_data <- inputs %>%
  magrittr::extract(names(.) != "references") %>%
  load_input_data(syn = syn) %>%
  map(replace_empty_string_na)

# Make tables ------------------------------------------------------------------
###############################################################################T

left_join_by_target_id <- function(x, y, by, check = "cbvmn") {
  x %>%
    mutate(jc = if_else(is.na(entrez_gene_id), symbol, as.character(entrez_gene_id))) %>%
    safe_left_join(
      y %>%
        mutate(jc = if_else(is.na(entrez_gene_id), symbol, as.character(entrez_gene_id))) %>%
        select(-symbol, -entrez_gene_id),
      by = union(by, "jc"),
      check = check
    )
}

lsp_compound_dictionary <- input_data[["compound_dictionary"]] %>%
  select(
    lspci_id, hmsl_id, chembl_id, emolecules_id, pref_name,
    inchi = canonical_inchi, commercially_available, highest_approval = max_phase
  ) %>%
  arrange(lspci_id)

lsp_structures <- input_data[["inchis"]] %>%
  filter(rank != 1) %>%
  select(
    lspci_id, source, rank, inchi = canonical_inchi
  ) %>%
  arrange(lspci_id, rank)

lsp_compound_names <- input_data[["compound_names"]] %>%
  select(
    lspci_id, source, priority = name_preference, name
  ) %>%
  distinct() %>%
  arrange(lspci_id)

lsp_compound_mapping <- input_data[["vendor_ids"]] %>%
  select(
    lspci_id, source, external_id = vendor_id
  ) %>%
  distinct() %>%
  arrange(lspci_id, external_id)

lsp_target_dictionary <- input_data[["target_dictionary"]] %>%
  select(
    gene_id = entrez_gene_id, symbol, pref_name, tax_id, organism
  ) %>%
  distinct() %>%
  arrange(gene_id)

lsp_target_mapping <- input_data[["target_dictionary"]] %>%
  select(
    gene_id = entrez_gene_id, chembl_id, uniprot_id, target_type
  ) %>%
  distinct() %>%
  arrange(gene_id) %>%
  drop_na(gene_id)

lsp_references <- synPluck(syn_release, "reference_table") %>%
  {synTableQuery(sprintf("SELECT * FROM %s", .))} %>%
  as.data.frame() %>%
  select(
    reference_id = ROW_ID, reference_type, reference_value, url
  ) %>%
  distinct() %>%
  setDT()

dose_response_q1_measurements <- input_data[["dose_response_q1_measurements"]][
  ,
  .(biochem_agg_id = .GRP),
  keyby = .(lspci_id, entrez_gene_id, symbol)
] %>%
  unique()

lsp_biochem <- input_data[["dose_response_measurements"]] %>%
  left_join_by_target_id(
    dose_response_q1_measurements,
    c("lspci_id")
  ) %>%
  rename(reference_value = reference_id) %>%
  safe_left_join(
    lsp_references %>%
      select(reference_type, reference_value, reference_id) %>%
      distinct(),
    by = c("reference_type", "reference_value"),
    check = "cbm"
  ) %>%
  transmute(
    biochem_id = measurement_id,
    biochem_agg_id,
    lspci_id,
    gene_id = entrez_gene_id,
    symbol,
    source = recode(measurement_source, chembl_activity = "chembl", inhouse_doseresponse = "lsp"),
    description_assay,
    value,
    value_type,
    value_unit,
    value_relation,
    reference_id
  ) %>%
  distinct() %>%
  arrange(lspci_id, gene_id, symbol)

lsp_biochem_agg <- input_data[["dose_response_q1"]] %>%
  left_join_by_target_id(
    dose_response_q1_measurements,
    by = c("lspci_id")
  ) %>%
  left_join_by_target_id(
    input_data[["tas"]] %>%
      filter(source == "dose_response") %>%
      distinct(tas_id, lspci_id, entrez_gene_id, symbol),
    by = c("lspci_id")
  ) %>%
  transmute(
    biochem_agg_id,
    lspci_id, gene_id = entrez_gene_id, symbol,
    value = Q1,
    value_unit = "nM",
    tas_id
  ) %>%
  distinct() %>%
  arrange(lspci_id, gene_id, symbol)

phenotypic_q1_measurements <- input_data[["phenotypic_q1_measurements"]][
  ,
  .(phenotypic_agg_id = .GRP),
  keyby = .(lspci_id, assay_id)
] %>%
  unique()

lsp_phenotypic <- input_data[["phenotypic_measurements"]] %>%
  safe_left_join(
    phenotypic_q1_measurements,
    by = c("lspci_id", "assay_id"),
    check = "bcvmn"
  ) %>%
  select(
    phenotypic_id = measurement_id,
    lspci_id,
    assay_id,
    value,
    value_type,
    value_unit,
    description_assay,
    reference_id,
    phenotypic_agg_id
  ) %>%
  distinct() %>%
  arrange(lspci_id, assay_id)

lsp_phenotypic_agg <- input_data[["phenotypic_q1"]] %>%
  safe_left_join(
    phenotypic_q1_measurements,
    by = c("lspci_id", "assay_id"),
    check = "bcuvmn"
  ) %>%
  transmute(
    phenotypic_agg_id,
    lspci_id,
    assay_id,
    value = value_Q1,
    value_unit = "nM"
  ) %>%
  distinct() %>%
  arrange(lspci_id, assay_id)

lsp_manual_curation <- input_data[["manual_curation"]] %>%
  # Unmatched values OK because some manual curations were
  # superceded by actual measurements
  left_join_by_target_id(
    input_data[["tas"]] %>%
      filter(source == "literature_annotation") %>%
      distinct(tas_id, lspci_id, entrez_gene_id, symbol),
    by = c("lspci_id")
  ) %>%
  safe_left_join(
    lsp_references %>%
      select(reference_type, reference_value, reference_id) %>%
      distinct(),
    by = c("reference_type", "reference_value"),
    check = "bcm"
  ) %>%
  transmute(
    lspci_id,
    gene_id = entrez_gene_id,
    symbol,
    reference_id,
    tas_id
  ) %>%
  distinct() %>%
  arrange(lspci_id, gene_id, symbol)

single_dose_q1_measurements <- input_data[["single_dose_q1_measurements"]][
  ,
  .(one_dose_scan_agg_id = .GRP),
  keyby = .(lspci_id, entrez_gene_id, symbol, cmpd_conc_nM)
] %>%
  unique()

lsp_one_dose_scans <- input_data[["single_dose_measurements"]] %>%
  left_join_by_target_id(
    single_dose_q1_measurements,
    by = c("lspci_id", "cmpd_conc_nM")
  ) %>%
  rename(reference_value = reference_id) %>%
  safe_left_join(
    lsp_references %>%
      select(reference_type, reference_value, reference_id) %>%
      distinct(),
    by = c("reference_type", "reference_value"),
    check = "bcm"
  ) %>%
  transmute(
    one_dose_scan_id = measurement_id,
    one_dose_scan_agg_id,
    lspci_id,
    gene_id = entrez_gene_id,
    symbol,
    source = recode(measurement_source, chembl_activity = "chembl", inhouse_single_dose = "lsp"),
    percent_control,
    concentration = cmpd_conc_nM,
    reference_id
  ) %>%
  distinct() %>%
  arrange(lspci_id, gene_id, symbol, concentration)

lsp_one_dose_scan_agg <- input_data[["single_dose_q1"]] %>%
  left_join_by_target_id(
    single_dose_q1_measurements,
    by = c("lspci_id", "cmpd_conc_nM")
  ) %>%
  # Unmatched values OK because some one-dose measurements
  # are superceded by full dose-response assays
  left_join_by_target_id(
    input_data[["tas"]] %>%
      filter(source == "single_dose") %>%
      distinct(tas_id, lspci_id, entrez_gene_id, symbol),
    by = c("lspci_id")
  ) %>%
  transmute(
    one_dose_scan_agg_id,
    lspci_id, gene_id = entrez_gene_id, symbol,
    percent_control = percent_control_Q1,
    concentration = cmpd_conc_nM,
    value_unit = "nM",
    tas_id
  ) %>%
  distinct() %>%
  arrange(lspci_id, gene_id, symbol, concentration)

lsp_clinical_info <- input_data[["approval"]] %>%
  select(
    lspci_id,
    max_phase,
    first_approval,
    oral,
    parenteral,
    topical,
    black_box_warning,
    first_in_class,
    prodrug,
    indication_class,
    withdrawn_flag,
    withdrawn_year,
    withdrawn_country,
    withdrawn_reason
  ) %>%
  distinct() %>%
  arrange(lspci_id)

lsp_commercial_availability <- input_data[["commercial_info"]] %>%
  select(
    lspci_id,
    emolecules_id,
    vendor = supplier_name,
    catalog_number,
    tier
  ) %>%
  distinct() %>%
  arrange(lspci_id)

lsp_fingerprints <- input_data[["fingerprints"]] %>%
  select(
    lspci_id,
    fingerprint_type = fp_name,
    fingerprint = fingerprints
  ) %>%
  distinct() %>%
  arrange(lspci_id, fingerprint_type)

lsp_compound_library <- input_data[["optimal_library"]] %>%
  select(
    lspci_id,
    gene_id = entrez_gene_id,
    rank,
    reason_included
  ) %>%
  distinct() %>%
  arrange(gene_id, rank)

lsp_tas <- input_data[["tas"]] %>%
  select(
    tas_id,
    lspci_id,
    gene_id = entrez_gene_id,
    symbol,
    tas,
    derived_from = source
  ) %>%
  distinct() %>%
  arrange(lspci_id, gene_id, symbol)

lsp_tas_references <- input_data[["tas"]] %>%
  select(tas_id, lspci_id, entrez_gene_id, symbol, source) %>%
  left_join_by_target_id(
    bind_rows(
      single_dose = input_data[["single_dose_q1_references"]],
      dose_response = input_data[["dose_response_q1_references"]],
      literature_annotation = copy(input_data[["manual_curation"]])[
        ,  `:=`(
          reference_id = reference_value
        )
      ],
      .id = "source"
    ) %>%
      select(lspci_id, entrez_gene_id, symbol, source, reference_type, reference_id),
    by = c("lspci_id", "source"),
    check = "bcum"
  ) %>%
  rename(reference_value = reference_id) %>%
  safe_left_join(
    lsp_references %>%
      distinct(reference_id, reference_type, reference_value),
    by = c("reference_type", "reference_value"),
    check = "bcvm"
  ) %>%
  select(
    tas_id,
    reference_id
  ) %>%
  distinct() %>%
  arrange(tas_id)

lsp_selectivity <- input_data[["selectivity"]] %>%
  select(
    lspci_id,
    gene_id = entrez_gene_id,
    selectivity_class,
    investigation_bias,
    strength,
    wilcox_pval,
    selectivity,
    tool_score,
    ic50_difference = IC50_diff,
    ontarget_ic50_q1 = ontarget_IC50_Q1,
    offtarget_ic50_q1 = offtarget_IC50_Q1,
    ontarget_n = ontarget_IC50_N,
    offtarget_n = offtarget_IC50_N
  ) %>%
  distinct() %>%
  arrange(lspci_id, gene_id)

# write csv files --------------------------------------------------------------
###############################################################################T

tables <- tribble(
  ~name, ~table_inputs,
  "lsp_compound_dictionary", c("compound_dictionary"),
  "lsp_structures", c("inchis"),
  "lsp_compound_names", c("compound_names"),
  "lsp_compound_mapping", c("vendor_ids"),
  "lsp_target_dictionary", c("target_dictionary"),
  "lsp_target_mapping", c("target_dictionary"),
  "lsp_references", c("references"),
  "lsp_biochem", c("dose_response_measurements", "dose_response_q1_measurements", "references"),
  "lsp_biochem_agg", c("dose_response_q1", "dose_response_q1_measurements", "tas"),
  "lsp_phenotypic", c("phenotypic_measurements", "phenotypic_q1_measurements"),
  "lsp_phenotypic_agg", c("phenotypic_q1", "phenotypic_q1_measurements"),
  "lsp_tas", c("tas"),
  "lsp_tas_references", c("tas", "single_dose_q1_references", "dose_response_q1_references", "manual_curation", "references"),
  "lsp_manual_curation", c("manual_curation", "tas", "references"),
  "lsp_selectivity", c("selectivity"),
  "lsp_one_dose_scans", c("single_dose_measurements", "single_dose_q1_measurements", "references"),
  "lsp_one_dose_scan_agg", c("single_dose_q1", "single_dose_q1_measurements", "tas"),
  "lsp_clinical_info", c("approval"),
  "lsp_commercial_availability", c("commercial_info"),
  "lsp_fingerprints", c("fingerprints"),
  "lsp_compound_library", c("optimal_library")
) %>%
  mutate(
    syn_used = map(
      table_inputs,
      ~inputs[.x]
    ),
    table = map(
      name,
      get
    ),
    path = file.path(
      dir_release,
      paste0(name, ".csv.gz")
    )
  )

pwalk(
  tables,
  function(name, table, path, ...) {
    message("Writing ", name)
    fwrite(
      table,
      path,
      na = "NULL"
    )
  }
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

syn_db_tables_parent <- synMkdir(syn_release, "db_tables")

pwalk(
  tables,
  function(name, path, syn_used, ...) {
    message("Uploading ", name)
    activity <- Activity(
      name = "Wrangle tables for postgresql import",
      used = unname(syn_used),
      executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/db_upload/01_prepare_tables.R"
    )
    synStoreMany(
      path, parentId = syn_db_tables_parent,
      activity = activity
    )
  }
)

# Import to PostgreSQL ---------------------------------------------------------
###############################################################################T

dir.create(file.path(dir_release, "db_upload"), showWarnings = FALSE)

synChildren("syn20981961") %>%
  magrittr::extract(str_ends(names(.), fixed(".csv.gz"))) %>%
  walk(synGet, downloadLocation = file.path(dir_release, "db_upload"), ifcollision = "overwrite.local")

list.files(file.path(dir_release, "db_upload"), full.names = TRUE) %>%
  walk(
    function(path, ...) {
      name <- str_sub(basename(path), end = -8L)
      colnames <- read.csv(path, nrows = 1)
      paste0(
        "gunzip -cd ", path, " | psql --command=\"COPY ", name, " (", paste("\"", names(colnames), "\"", sep = "", collapse = ", "), ")",
        " FROM STDIN CSV HEADER NULL 'NA';\" sms_db"
      ) %>%
        message()
    }
  )

