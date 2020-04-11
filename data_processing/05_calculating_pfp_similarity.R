
library(tidyverse)
library(data.table)
library(here)
library(Matrix)
library(furrr)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))



`%nin%` <- Negate(`%in%`)

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

pheno_data <- syn("syn20841032") %>%
  read_rds()

# convert assay results to r-scores --------------------------------------------
###############################################################################T

calculate_r_score_per_assay <- function(df) {
  if(nrow(df) <= 2)
    return(NULL)
  c.med <- median(df$log10_value_Q1)
  c.mad <- mad(df$log10_value_Q1)
  c.result <- tibble(
    lspci_id = df$lspci_id,
    rscore = (df$log10_value_Q1 - c.med)/c.mad
  ) %>%
    mutate(
      rscore_tr = ifelse(
        rscore < 0,
        5 * (1 / (1+ exp(-(rscore + 2.5)*2)) - 1) + 0.0335,
        5 * (1 / (1+ exp(-(rscore - 2.5)*2))) - 0.0335
      )
    )
  c.result
}

calculate_r_score <- function(df) {
  plan(multicore(workers = 8))
  df %>%
    group_nest(assay_id) %>%
    mutate(
      data = future_map(
        data, calculate_r_score_per_assay,
        .progress = TRUE
      )
    )
}

rscores_raw <- pheno_data %>%
  mutate(
    data = map(
      data, calculate_r_score
    )
  )

rscores <- rscores_raw %>%
  mutate(
    data = map(
      data,
      ~set_names(.x$data, .x$assay_id) %>%
        rbindlist(idcol = "assay_id") %>%
        as_tibble()
    )
  )

write_rds(
  rscores,
  file.path(dir_release, "pheno_data_rscores.rds"),
  compress = "gz"
)

# Upload to synapse ------------------------------------------------------------
###############################################################################T

pfp_activity <- Activity(
  name = "Calculate phenotypic correlation",
  used = c(
    "syn20841032"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/05_calculating_pfp_similarity.R"
)

syn_pfp_sim <- Folder("pfp_similarity", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "pheno_data_rscores.rds")
) %>%
  synStoreMany(parent = syn_pfp_sim, activity = pfp_activity)

