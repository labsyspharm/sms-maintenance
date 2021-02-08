
library(tidyverse)
library(data.table)
library(here)
library(furrr)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))



`%nin%` <- Negate(`%in%`)

# set directories, import files ------------------------------------------------
###############################################################################T

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

inputs <- list(
  pheno_q1 = c("aggregate_data", "phenotypic_q1.csv.gz")
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

plan(multicore(workers = 10))

rscores_raw <- copy(input_data[["pheno_q1"]])[
  ,
  `:=`(
    value_Q1 = as.numeric(value_Q1),
    log10_value_Q1 = log10(as.numeric(value_Q1))
  )
][
  ,
  .(data = list(.SD)),
  keyby = .(assay_id)
][
  ,
  data := future_map(
    data,
    calculate_r_score_per_assay,
    .progress = TRUE,
    .options = furrr_options(
      scheduling = 10L
    )
  )
]

rscores <- rscores_raw %>%
  unnest(data) %>%
  setDT() %>%
  setkey(assay_id, lspci_id)

fwrite(
  rscores,
  file.path(dir_release, "phenotypic_rscores.csv.gz")
)

# Upload to synapse ------------------------------------------------------------
###############################################################################T

pfp_activity <- Activity(
  name = "Calculate phenotypic R scores",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/05_calculating_pfp_similarity.R"
)

syn_pfp_sim <- synMkdir(syn_release, "pfp_similarity")

c(
  file.path(dir_release, "phenotypic_rscores.csv.gz")
) %>%
  synStoreMany(parent = syn_pfp_sim, activity = pfp_activity)

