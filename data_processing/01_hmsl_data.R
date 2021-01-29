library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Loading files ----------------------------------------------------------------
###############################################################################T

inputs <- list(
  hmsl_doseresponse_20 = "syn20692433",
  hmsl_doseresponse_21 = "syn24210560"
)

input_data <- inputs %>%
  map(syn) %>%
  map(
    function(x)
      list(
        `.csv` = partial(fread, colClasses = c(inchi_id = "integer")),
        `.tsv` = fread,
        `.rds` = read_rds,
        `.xlsx` = openxlsx::read.xlsx
      ) %>%
      magrittr::extract2(which(str_detect(x, fixed(names(.))))) %>%
      {.(x)}
  )

# Compile data -----------------------------------------------------------------
###############################################################################T

doseresponse <- input_data[["hmsl_doseresponse_20"]] %>%
  rename(hmsl_id = hms_id) %>%
  mutate(
    hmsl_id = if_else(
      str_starts(hmsl_id, fixed("HMSL")),
      as.character(hmsl_id),
      paste0("HMSL", hmsl_id)
    )
  ) %>%
  bind_rows(
    input_data[["hmsl_doseresponse_21"]] %>%
      # Remove compounds where no assay was available
      drop_na(Km) %>%
      transmute(
        hmsl_id,
        value = nafill(
          Ki, fill = 10000
        ),
        value_type = "Ki",
        value_unit = Ki_unit,
        value_relation = if_else(
          value != 10000,
          "=", ">"
        ),
        cmpd_name = compound,
        symbol = target,
        synapse_id = "syn24210560",
        file_url = "https://www.synapse.org/#!Synapse:syn24210560"
      ) %>%
      genebabel::join_hgnc(
        "symbol", c("symbol", "alias_symbol"), c("uniprot_ids", "name", "entrez_id")
      ) %>%
      rename(
        uniprot_id = uniprot_ids,
        description = name,
        gene_id = entrez_id
      ) %>%
      mutate(
        gene_id = as.integer(gene_id)
      ) %>%
      unchop(uniprot_id)
  )

fwrite(
  doseresponse,
  here(release, "hmsl_doseresponse.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

hmsl_activity <- Activity(
  name = "Wrangle old and new internal HMSL data.",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/01_hmsl_data.R"
)

syn_id_mapping <- synMkdir(syn_release, "raw_data", "hmsl")

c(
  here(release, "hmsl_doseresponse.csv.gz")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = hmsl_activity)

