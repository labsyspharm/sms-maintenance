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

inputs <- inputs <- list(
  target_dictionary = synPluck(syn_release, "id_mapping", "target_dictionary_wide.csv.gz"),
  lspci_id_vendor_id_map = synPluck(syn_release, "compounds_processed", "lspci_id_vendor_id_map.csv.gz")
) %>%
  c(
    hmsl_doseresponse_20 = "syn20692433",
    hmsl_doseresponse_21 = "syn24210560",
    inhouse_single_dose = "syn20692432"
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
  rename(entrez_gene_id = gene_id) %>%
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
      inner_join(
        input_data[["target_dictionary"]][
          ,
          .(
            symbol,
            entrez_gene_id,
            uniprot_id,
            description = pref_name
          )
        ],
        by = "symbol"
      )
  ) %>%
  mutate(
    reference_value = synapse_id,
    reference_type = "synapse_id"
  ) %>%
  distinct()

fwrite(
  doseresponse,
  here(release, "hmsl_doseresponse.csv.gz")
)

single_dose <- input_data[["inhouse_single_dose"]] %>%
  mutate(
    hmsl_id = if_else(
      str_starts(hms_id, fixed("HMSL")),
      hms_id,
      paste0("HMSL", hms_id)
    )
  ) %>%
  select(-hms_id) %>%
  inner_join(
    input_data[["lspci_id_vendor_id_map"]][
      source == "hmsl",
      .(lspci_id, hmsl_id = vendor_id)
    ],
    by = "hmsl_id"
  ) %>%
  rename(symbol = gene_symbol) %>%
  inner_join(
    input_data[["target_dictionary"]][
      ,
      .(
        symbol,
        entrez_gene_id,
        uniprot_id
      )
    ],
    by = "symbol"
  ) %>%
  mutate(
    reference_type = if_else(
      str_starts(url, fixed("http://lincs.hms.harvard.edu")),
      "hmsl_id",
      "synapse_id"
    ),
    reference_value = source_assay_id
  ) %>%
  distinct()

fwrite(
  single_dose,
  here(release, "hmsl_singledose.csv.gz")
)

# Storing references -----------------------------------------------------------
###############################################################################T

REFERENCE_PRIORITY <- c(
  "pubmed_id" = "https://pubmed.ncbi.nlm.nih.gov/",
  "doi" = "http://doi.org/",
  "patent_id" = "https://patents.google.com/?q=",
  "synapse_id" = "https://www.synapse.org/#!Synapse:",
  "chembl_id" = "https://www.ebi.ac.uk/chembl/document_report_card/",
  "hmsl_id" = "https://lincs.hms.harvard.edu/db/datasets/"
)

syn_reference_table <- synPluck(syn_release, "reference_table")

current_references <- synTableQuery(sprintf("SELECT * FROM %s", syn_reference_table)) %>%
  as.data.frame()

new_references <- bind_rows(
  doseresponse,
  single_dose
) %>%
  anti_join(
    current_references,
    by = c("reference_type", "reference_value")
  ) %>%
  transmute(
    reference_type,
    reference_value,
    url = paste0(REFERENCE_PRIORITY[reference_type], reference_value)
  ) %>%
  distinct()

synStore(
  Table(
    syn_reference_table,
    new_references
  )
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
  here(release, "hmsl_doseresponse.csv.gz"),
  here(release, "hmsl_singledose.csv.gz")
) %>%
  synStoreMany(parent = syn_id_mapping, activity = hmsl_activity)

