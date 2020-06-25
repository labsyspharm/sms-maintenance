library(tidyverse)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"), ifcollision = "overwrite.local")

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
  ~name, ~synapse_id, ~used, ~sort_by, ~fun,
  "lsp_compound_dictionary", "syn20835543", NULL, c("lspci_id"), function(x)
    select(x, lspci_id, hmsl_id = hms_id, chembl_id, pref_name, inchi, smiles, commercially_available),
  "lsp_biochem", "syn20830825", NULL, c("lspci_id", "gene_id"), function(x)
    select(x, lspci_id, gene_id = entrez_gene_id, description_assay, value, value_unit, value_type, value_relation, reference_id, reference_type, url = file_url),
  "lsp_compound_names", "syn22035396", NULL, NULL, function(x)
    transmute(x, lspci_id, source = str_split_fixed(source, fixed("_"), 2)[, 1], priority = source_collapsed, name),
  "lsp_compound_mapping", "syn20830516", NULL, c("lspci_id"), function(x)
    transmute(x, lspci_id = eq_class, source = if_else(str_starts(id, fixed("CHEMBL")), "chembl", "hmsl"), id),
  "lsp_phenotypic_chembl", "syn20976900", NULL, c("lspci_id", "assay_id"), function(x)
    select(x, lspci_id, assay_id, description_assay, value, value_unit, value_type, value_relation, reference_id, reference_type, url = file_url),
  "lsp_tas", "syn20830939", NULL, c("lspci_id", "gene_id"), function(x) select(x, lspci_id, gene_id = entrez_gene_id, tas),
  "lsp_specificity", "syn20836653", NULL, c("lspci_id", "gene_id"), function(x)
    distinct(
      x, lspci_id, gene_id, selectivity_class, investigation_bias, strength,
      wilcox_pval, selectivity, tool_score, IC50_difference = IC50_diff,
      ontarget_IC50_Q1, offtarget_IC50_Q1, ontarget_N = ontarget_IC50_N, offtarget_N = offtarget_IC50_N
    ) %>%
    rename_all(str_to_lower),
  "lsp_one_dose_scans", "syn20830835", NULL, c("lspci_id", "gene_id"), function(x)
    transmute(
      x, lspci_id, gene_id = entrez_gene_id, percent_control, description, concentration = cmpd_conc_nM,
      reference_id, reference_type = recode(reference_type, synapse = "synapse_id", hms_lincs = "hmsl_id"),
      url = file_url
    ) %>%
    distinct(),
  # "lsp_tas_similarity", "syn21052803", NULL,
  "lsp_clinical_info", "syn21064122", NULL, c("lspci_id"), function(x)
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
  "lsp_commercial_availability", "syn21049601", NULL, c("lspci_id"), function(x)
    distinct(x, lspci_id, vendor, id = vendor_id, name = vendor_name),
  "lsp_fingerprints", "syn21042105", NULL, c("lspci_id", "fingerprint_type"), function(x)
    select(x, lspci_id, fingerprint_type = fp_name, fingerprint),
  "lsp_compound_library", "syn22153734", NULL, c("lspci_id", "gene_id"), function(x)
    select(x, lspci_id, gene_id, rank, reason_included)
)

tables_dl <- tables %>%
  mutate(
    data = map(
      synapse_id,
      . %>%
        syn() %>%
        read_rds()
    )
  )

process_table <- function(data, fun, sort_by, ...) {
  data %>%
    mutate_at(vars(data), map, fun) %>%
    mutate_at(vars(data), map, arrange, !!!rlang::syms(sort_by))
}

tables_formatted <- tables_dl %>%
  mutate(
    data = pmap(., process_table)
  )

# write csv files --------------------------------------------------------------
###############################################################################T

dirs_tables <- tables_formatted$data[[1]]$fp_name %>%
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

# Make target mapping table ----------------------------------------------------
###############################################################################T

target_dict <- syn("syn20693721") %>%
  read_csv()

target_table <- target_dict %>%
  drop_na(entrez_gene_id) %>%
  # Replace nonsense symbol "-" with uniprot ID
  mutate(symbol = recode(symbol, `-` = paste0("uniprot_", uniprot_id))) %>%
  distinct(gene_id = entrez_gene_id, symbol, tax_id, organism, description = entrez_description) %>%
  arrange(gene_id)

target_mapping <- target_dict %>%
  drop_na(entrez_gene_id) %>%
  distinct(gene_id = entrez_gene_id, chembl_id, target_type, uniprot_id, pref_name) %>%
  filter(!pmap_lgl(list(chembl_id, target_type, uniprot_id, pref_name), ~map_lgl(.x, is.na) %>% all())) %>%
  arrange(gene_id)

write_csv(target_table, file.path(dir_release, "lsp_target_dictionary.csv.gz"))
write_csv(target_mapping, file.path(dir_release, "lsp_target_mapping.csv.gz"))

target_activity <- Activity(
  name = "Wrangle target table for postgresql import",
  used = "syn20693721",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/db_upload/01_prepare_tables.R"
)

walk(
  unname(synChildren(syn_db_tables_parent)),
  function(syn_id) {
    file.path(dir_release, c("lsp_target_dictionary.csv.gz", "lsp_target_mapping.csv.gz")) %>%
      synStoreMany(
        parentId = syn_id,
        activity = target_activity
      )
  }
)


# # Link PFP correlation for each fp_type
# x <- tribble(
#   ~fp_name, ~syn_source,
#   "morgan_normal", "syn21116196",
#   "morgan_chiral", "syn21116195",
#   "topological_normal", "syn21116206"
# ) %>%
#   left_join(syn_db_table_map, by = "fp_name") %>%
#   pwalk(
#     .,
#     function(fp_name, syn_source, parent_dir, ...) {
#       fn <- syn(syn_source)
#       link <- Link(syn_source, parent = parent_dir) %>%
#         synStore()
#       link$properties$name = "lsp_pfp_correlation.csv.gz"
#       synStore(link)
#     }
#   )

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

