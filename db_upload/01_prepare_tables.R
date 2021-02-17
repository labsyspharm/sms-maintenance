library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)

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

# Dealing with target dictionary seperately because it does not depend
# on the fingerprinting algo

lsp_compound_dictionary <- synPluck(syn_release, "compounds_processed", "compound_dictionary.csv.gz") %>%
  syn() %>%
  fread() %>%
  replace_empty_string_na() %>%
  select(
    lspci_id, hmsl_id, chembl_id, emolecules_id, pref_name,
    inchi = canonical_inchi, commercially_available, highest_approval = max_phase
  ) %>%
  arrange(lspci_id)

tables <- tribble(
  ~name, ~synapse_path, ~used, ~sort_by, ~fun,
  "lsp_compound_dictionary",
    c("compounds_processed", "compound_dictionary.csv.gz"),
    NULL,
    c("lspci_id"),
    function(x) select(
      x, lspci_id, hmsl_id, chembl_id, emolecules_id, pref_name,
      inchi = canonical_inchi, commercially_available, highest_approval = max_phase
    ),
  "lsp_compound_names",
    c("compounds_processed", "")
)

tables_dl <- tables %>%
  mutate(
    synapse_id = pluck_inputs(synapse_path, syn_parent = syn_release),
    data = load_input_data(synapse_id, syn = syn)
  )

process_table <- function(data, fun, sort_by, ...) {
  browser()
  data %>%
    fun() %>%
    arrange(!!!rlang::syms(sort_by))
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

