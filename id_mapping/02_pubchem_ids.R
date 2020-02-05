library(tidyverse)
library(synapser)
library(synExtra)
library(vroom)
library(here)

synLogin()
syn <- synDownloader(here("tempdl"))

# Pubchem offers a very comprehensive list of compound names which is lacking
# in Chembl, fetching them from Pubchem in addition
# Mapping from Pubchem to Chembl IDs is provided by Unichem

# Set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

chembl_pubchem_mapping <- syn("syn21572661") %>%
  read_csv()


# Fetch compound synonyms from Pubchem -----------------------------------------
###############################################################################T

download.file(
  "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/2020-02-01/Extras/CID-Synonym-filtered.gz",
  file.path(dir_release, "pubchem_synonyms.csv.gz")
)

pubchem_synonyms_raw <- vroom(
  file.path(dir_release, "pubchem_synonyms.csv.gz"),
  col_names = c("pubchem_id", "name"),
  col_types = "ic"
)

pubchem_synonyms_chembl <- pubchem_synonyms_raw %>%
  inner_join(
    chembl_pubchem_mapping %>%
      distinct(pubchem_id, chembl_id),
    by = "pubchem_id"
  )

write_csv(
  pubchem_synonyms_chembl,
  file.path(dir_release, "pubchem_synonyms_chembl.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

wrangle_activity <- Activity(
  name = "Wrangle compound names and synonyms from Pubchem",
  used = c(
    "syn21572661",
    "ftp://ftp.ncbi.nlm.nih.gov/pubchem/Compound/Monthly/2020-02-01/Extras/CID-Synonym-filtered.gz"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/02_pubchem_ids.R"
)

syn_pubchem <- Folder("pubchem", parent = "syn20830877") %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "pubchem_synonyms_chembl.csv.gz")
) %>%
  synStoreMany(parentId = syn_pubchem, activity = wrangle_activity)
