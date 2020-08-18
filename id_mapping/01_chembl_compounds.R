library(tidyverse)
library(httr)
library(furrr)
library(RPostgres)
library(bit64)
library(jsonlite)
library(data.table)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synMkdir("syn18457321", release)

## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("Postgres")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "chembl_27",
                 host = "localhost", port = 5432,
                 user = "chug")

# Fetch Chembl compound data ---------------------------------------------------
###############################################################################T

all_cmpds <- dbGetQuery(
  con,
  "SELECT DISTINCT dict.molregno, dict.pref_name, dict.chembl_id, dict.max_phase,
    hi.parent_molregno, hi.active_molregno, struct.standard_inchi, struct.standard_inchi_key,
    struct.canonical_smiles, dict.inorganic_flag, ass.n_assays
  FROM molecule_dictionary AS dict
  LEFT JOIN molecule_hierarchy AS hi ON dict.molregno = hi.molregno
  LEFT JOIN compound_structures AS struct ON dict.molregno=struct.molregno
  LEFT JOIN (
    SELECT COUNT(activity_id) AS n_assays, molregno FROM activities GROUP BY molregno
  ) AS ass ON dict.molregno = ass.molregno"
) %>%
  as_tibble() %>%
  mutate(parental_flag = ifelse(molregno != parent_molregno, 1L, 0L))

all_synonyms <- dbGetQuery(
  con,
  "SELECT DISTINCT dict.molregno, syn.synonyms
  FROM molecule_dictionary AS dict
  LEFT JOIN molecule_synonyms AS syn ON DICT.molregno = syn.molregno"
) %>%
  drop_na() %>%
  as.data.table() %>%
  .[
    ,
    .(synonyms = list(synonyms)),
    keyby = molregno
  ] %>%
  as_tibble()

all_cmpds <- all_cmpds %>%
  left_join(
    all_synonyms,
    by = "molregno"
  )

write_rds(
  all_cmpds,
  file.path(dir_release, "chembl_compounds_raw.rds"),
  compress = "gz"
)

# all_cmpds <- read_rds(file.path(dir_release, "chembl_compounds_raw.rds"))

# Store to synapse -------------------------------------------------------------
###############################################################################T

fetch_chembl_activity <- Activity(
  name = "Fetch ChEMBL compound data",
  used = "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_27/chembl_27_postgresql.tar.gz",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/01_chembl_compounds.R"
)

syn_raw <- synMkdir(syn_release, "raw_data")

list(
  file.path(dir_release, "chembl_compounds_raw.rds")
) %>%
  map(
    . %>%
      File(parent = syn_raw) %>%
      synStore(activity = fetch_chembl_activity)
  )
