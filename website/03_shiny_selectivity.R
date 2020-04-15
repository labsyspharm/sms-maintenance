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

biochemical <- syn("syn20830834") %>%
  read_rds()

clinical <- syn("syn21064122") %>%
  read_rds()

# chemical_probes <- syn("syn21627808") %>%
#   read_rds()

# Selectivity table ------------------------------------------------------------
###############################################################################T

combine_tables <- function(selectivity, biochemical, clinical, ...) {
  df <- selectivity %>%
    distinct(
      lspci_id,
      gene_id,
      symbol = gene_symbol,
      chembl_id,
      hms_id,
      toolscore = tool_score,
      ontarget_IC50_Q1,
      ontarget_IC50_N,
      offtarget_IC50_Q1,
      offtarget_IC50_N,
      IC50_diff,
      selectivity_class = recode(
        selectivity_class,
        most_selective = "Most selective",
        semi_selective = "Semi-selective",
        poly_selective = "Poly-selective",
        unknown_selective = "Unknown",
        other_selective = "Other"
      )
    ) %>%
    left_join(
      biochemical %>%
        distinct(lspci_id, gene_id = entrez_gene_id, Kd_Q1 = Q1, n_measurement_kd = n_measurement),
      by = c("lspci_id", "gene_id")
    ) %>%
    left_join(
      clinical %>%
        group_by(lspci_id) %>%
        summarize(
          max_phase = max(max_phase, na.rm = TRUE) %>% {if_else(. > 4, NA_real_, .)},
          first_approval = min(first_approval, na.rm = TRUE) %>% {if_else(. < 0, NA_real_, .)}
        ) %>%
        ungroup(),
      by = "lspci_id"
    ) %>%
    as.data.table()
  setkey(df, lspci_id, gene_id)
  df
}

selectivity_table <- selectivity_classes %>%
  rename(selectivity = data) %>%
  inner_join(
    biochemical %>%
      rename(biochemical = data),
    by = c("fp_type", "fp_name")
  ) %>%
  inner_join(
    clinical %>%
      rename(clinical = data),
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
    "syn20830834",
    "syn21064122"
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


activity <- Activity(
  name = "Wrangle biochemical table",
  used = c(
    "syn20830834"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/03_shiny_selectivity.R"
)

pwalk(
  biochemical,
  function(fp_name, data, ...) {
    data <- data %>%
      rename(gene_id = entrez_gene_id, Kd_Q1 = Q1, n_measurement_kd = n_measurement) %>%
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
