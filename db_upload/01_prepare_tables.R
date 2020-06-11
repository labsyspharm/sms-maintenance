library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# read tables ------------------------------------------------------------------
###############################################################################T

# Dealing with target dictionary seperately because it does not depend
# on the fingerprinting algo

tables <- tribble(
  ~name, ~synapse_id, ~used, ~fun,
  "lsp_compound_dictionary", "syn20835543", NULL, function(x)
    select(x, lspci_id, hmsl_id = hms_id, chembl_id, pref_name, inchi, smiles, commercially_available),
  "lsp_biochem", "syn20830825", NULL, function(x)
    select(x, lspci_id, gene_id = entrez_gene_id, description_assay, value, value_unit, value_type, value_relation, reference_id, reference_type, url = file_url),
  "lsp_compound_names", "syn22035396", NULL, function(x)
    transmute(x, lspci_id, source = str_split_fixed(source, fixed("_"), 2)[, 1], priority = source_collapsed, name),
  "lsp_compound_mapping", "syn20830516", NULL, function(x)
    transmute(x, lspci_id = eq_class, source = if_else(str_starts(id, fixed("CHEMBL")), "chembl", "hmsl"), id),
  "lsp_phenotypic_chembl", "syn20976900", NULL, function(x)
    select(x, lspci_id, assay_id, description_assay, value, value_unit, value_type, value_relation, reference_id, reference_type, url = file_url),
  "lsp_tas", "syn20830939", NULL, function(x) select(x, lspci_id, gene_id = entrez_gene_id, tas),
  "lsp_specificity", "syn20836653", NULL, function(x)
    distinct(
      x, lspci_id, gene_id, selectivity_class, investigation_bias, strength,
      wilcox_pval, selectivity, tool_score, IC50_difference = IC50_diff,
      ontarget_IC50_Q1, offtarget_IC50_Q1, ontarget_N = ontarget_IC50_N, offtarget_N = offtarget_IC50_N
    ),
  "lsp_one_dose_scans", "syn20830835", NULL, function(x)
    transmute(
      x, lspci_id, gene_id = entrez_gene_id, percent_control, description, cmpd_conc_nM,
      reference_id, reference_type = recode(reference_type, synapse = "synapse_id", hms_lincs = "hmsl_id"),
      url = file_url
    ) %>%
    distinct(),
  # "lsp_tas_similarity", "syn21052803", NULL,
  "lsp_clinical_info", "syn21064122", NULL, function(x)
    mutate_at(
      x,
      vars(oral, parenteral, topical, black_box_warning, first_in_class, prodrug, withdrawn_flag),
      magrittr::equals, 1
    ) %>%
    distinct(
      lspci_id, max_phase, first_approval, oral, parenteral, topical, black_box_warning,
      first_in_class, prodrug, indication_class, withdrawn_flag, withdrawn_year, withdrawn_country,
      withdrawn_reason, max_phase_for_indication, mesh_id, mesh_heading, efo_id, efo_term,
      reference_type = ref_type, reference_id = ref_id, url = ref_url
    ),
  "lsp_commercial_availability", "syn21049601", NULL, function(x)
    distinct(x, lspci_id, vendor, id = vendor_id, name = vendor_name),
  "lsp_fingerprints", "syn21042105", NULL, function(x)
    select(x, lspci_id, fingerprint_type = fp_name, fingerprint),
  "lsp_compound_library", "syn22153734", NULL, function(x)
    select(x, lspci_id, gene_id, rank, reason_included)
) %>%
  mutate(
    data = map(
      synapse_id,
      . %>%
        syn() %>%
        read_rds()
    )
  )

tables_formatted <- tables %>%
  mutate(
    data = map2(
      data, fun,
      ~mutate_at(.x, vars(data), map, .y)
    )
  )

# write csv files --------------------------------------------------------------
###############################################################################T

dirs_tables <- tables$data[[1]]$fp_name %>%
  unique() %>%
  {file.path(dir_release, .)}
dirs_tables %>%
  walk(dir.create, showWarnings = FALSE)

tables_long <- tables_formatted %>%
  unnest(data) %>%
  mutate(
    fn = file.path(
      dir_release, fp_name, paste0(name, ".csv.gz")
    ),
    # Concatenate list columns for csv export
    data = map(
      data,
      mutate_if, is.list,
      ~if_else(map_lgl(.x, is.null), NA_character_, map_chr(.x, paste, collapse = "|"))
    )
  )

walk2(
  tables_long$data, tables_long$fn,
  write_csv
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

syn_db_tables_parent <- Folder("db_tables", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

syn_db_table_map <- tables_long %>%
  distinct(fp_name) %>%
  mutate(
    parent_dir = map_chr(
      fp_name,
      ~Folder(.x, parent = syn_db_tables_parent) %>%
        synStore() %>%
        chuck("properties", "id")
    )
  )

syn_tables <- tables_long %>%
  rename(source_syn = synapse_id) %>%
  mutate(
    used = map2(used, source_syn, c),
    activity = map(
      used,
      ~Activity(
        name = "Wrangle tables for postgresql import",
        used = .x,
        executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/db_upload/01_prepare_tables.R"
      )
    )
  ) %>%
  left_join(syn_db_table_map, by = "fp_name")


pwalk(
  syn_tables,
  function(fn, parent_dir, activity, ...) {
    synStoreMany(fn, parent = parent_dir, activity = activity)
  }
)

# Link target mapping for each fp_type
syn_db_table_map %>%
  pull(parent_dir) %>%
  walk(
    function(x) {
      link <- Link("syn20693721", parent = x) %>%
        synStore()
      link$properties$name = "lsp_target_dictionary.csv.gz"
      synStore(link)
    }
  )

# Link PFP correlation for each fp_type
x <- tribble(
  ~fp_name, ~syn_source,
  "morgan_normal", "syn21116196",
  "morgan_chiral", "syn21116195",
  "topological_normal", "syn21116206"
) %>%
  left_join(syn_db_table_map, by = "fp_name") %>%
  pwalk(
    .,
    function(fp_name, syn_source, parent_dir, ...) {
      fn <- syn(syn_source)
      link <- Link(syn_source, parent = parent_dir) %>%
        synStore()
      link$properties$name = "lsp_pfp_correlation.csv.gz"
      synStore(link)
    }
  )
