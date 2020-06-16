library(tidyverse)
library(here)
library(fst)
library(data.table)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

syn_tables <- "syn20981852"

selectivity_classes <- syn("syn20836653") %>%
  read_rds()

dose_response <- syn("syn20830834") %>%
  read_rds()

chemical_probes <- syn("syn21627808") %>%
  read_rds()

# Selectivity table ------------------------------------------------------------
###############################################################################T

combine_tables <- function(selectivity, biochemical, ...) {
  df <- selectivity %>%
    distinct(
      lspci_id,
      gene_id,
      selectivity,
      selectivity_class = recode(
        selectivity_class,
        most_selective = "Most selective",
        semi_selective = "Semi-selective",
        poly_selective = "Poly-selective",
        unknown_selective = "Unknown",
        other_selective = "Other"
      ),
      toolscore = tool_score,
      affinity_Q1 = ontarget_IC50_Q1,
      affinity_N = ontarget_IC50_N,
      offtarget_affinity_Q1 = offtarget_IC50_Q1,
      offtarget_affinity_N = offtarget_IC50_N,
      affinity_Q1_diff = IC50_diff,
      investigation_bias,
      wilcox_pval,
      strength
    ) %>%
    left_join(
      biochemical %>%
        distinct(lspci_id, gene_id = entrez_gene_id, references),
      by = c("lspci_id", "gene_id")
    ) %>%
    as.data.table()
  setkey(df, lspci_id, gene_id)
  df
}

selectivity_table <- selectivity_classes %>%
  rename(selectivity = data) %>%
  inner_join(
    dose_response %>%
      rename(biochemical = data),
    by = c("fp_type", "fp_name")
  ) %>%
  transmute(
    fp_type,
    fp_name,
    data = pmap(
      .,
      combine_tables
    )
  )

activity <- Activity(
  name = "Wrangle selectivity table",
  used = c(
    "syn20836653",
    "syn20830834"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/03_shiny_selectivity.R"
)

pwalk(
  selectivity_table,
  function(fp_name, data, ...) {
    write_fst(
      data,
      file.path(dir_release, paste0("shiny_selectivity_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_selectivity_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_selectivity.fst"
    ) %>%
      synStore(activity = activity)
  }
)


# Binding data -----------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Wrangle biochemical table",
  used = c(
    "syn20830834"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/03_shiny_selectivity.R"
)

pwalk(
  dose_response,
  function(fp_name, data, ...) {
    data <- data %>%
      distinct(lspci_id, gene_id = entrez_gene_id, affinity_Q1 = Q1, affinity_N = n_measurement, references) %>%
      as.data.table()
    setkey(data, lspci_id, gene_id)
    write_fst(
      data,
      file.path(dir_release, paste0("shiny_biochemical_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_biochemical_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_biochemical.fst"
    ) %>%
      synStore(activity = activity)
  }
)


# Chemical probes --------------------------------------------------------------
###############################################################################T


pwalk(
  chemical_probes,
  function(fp_name, data, ...) {
    df <- data %>%
      distinct(lspci_id, name, gene_id = entrez_id, n_reviews, avg_rating) %>%
      as.data.table()
    setkey(df, lspci_id, gene_id)
    write_fst(
      df,
      file.path(dir_release, paste0("shiny_chemical_probes_", fp_name, ".fst"))
    )
    syn_fp <- synPluck(syn_tables, fp_name)
    syn_parent <- Folder("website", parent = syn_fp) %>%
      synStore() %>%
      chuck("properties", "id")
    File(
      file.path(dir_release, paste0("shiny_chemical_probes_", fp_name, ".fst")),
      parent = syn_parent,
      name = "shiny_chemical_probes.fst"
    ) %>%
      synStore(activity = activity)
  }
)
