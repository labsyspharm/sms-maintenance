library(tidyverse)
library(vroom)
library(data.table)
# library(biomaRt)
library(bit64)
library(here)
library(synapser)
library(synExtra)
library(RPostgres)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# connect to chembl ------------------------------------------------------------
###############################################################################T


## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("Postgres")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "chembl_27",
                 host = "localhost", port = 5432,
                 user = "chug")

# create target conversion table -----------------------------------------------
###############################################################################T

download.file(
  "ftp://ftp.ebi.ac.uk/pub/databases/chembl/ChEMBLdb/releases/chembl_27/chembl_uniprot_mapping.txt",
  file.path(dir_release, "chembl_uniprot_mapping.txt")
)

uniprot_id_mapping_urls <- c(
  "Homo sapiens" = "HUMAN_9606_idmapping.dat.gz",
  "Mus musculus" = "MOUSE_10090_idmapping.dat.gz",
  "Rattus norvegicus" = "RAT_10116_idmapping.dat.gz"
) %>% map_chr(
    ~paste0(
      "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/by_organism/",
      .x
    )
  )

uniprot_id_mapping_files <- uniprot_id_mapping_urls %>%
  imap_chr(
    ~file.path(dir_release, paste0("uniprot_id_mapping_", str_replace_all(.y, fixed(" "), "_"), ".tsv.gz"))
  )

download.file(
  uniprot_id_mapping_urls,
  uniprot_id_mapping_files,
  method = "libcurl"
)

download.file(
  "ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/idmapping/idmapping_selected.tab.gz",
  file.path(dir_release, "uniprot_id_mapping.tsv.gz")
)

# step 1 --> import chembl generated file & convert to gene_ID

organism_biomart_mapping <- c(
  "Homo sapiens" = "hsapiens_gene_ensembl",
  "Rattus norvegicus" = "rnorvegicus_gene_ensembl",
  "Mus musculus" = "mmusculus_gene_ensembl"
)

allowed_target_types = c(
  "SELECTIVITY GROUP",
  "PROTEIN NUCLEIC-ACID COMPLEX",
  "PROTEIN FAMILY",
  "CHIMERIC PROTEIN",
  "PROTEIN COMPLEX",
  "SINGLE PROTEIN",
  "PROTEIN COMPLEX GROUP",
  "PROTEIN-PROTEIN INTERACTION",
  "UNKNOWN"
)

chembl_info_all_targets <- dbGetQuery(
  con,
  paste0(
    "SELECT dict.tid, dict.pref_name, dict.tax_id, dict.organism, dict.chembl_id, dict.target_type
    FROM target_dictionary AS dict"
  )
) %>%
  as_tibble() %>%
  filter(
    organism %in% names(organism_biomart_mapping),
    target_type %in% allowed_target_types
  )

chembl_info_all_targets %>%
  count(chembl_id, target_type) %>%
  count(n)
# # A tibble: 1 x 2
# n    nn
# <int> <int>
#   1     1  5842
# Only a single target type per chembl_id


map_uniprot_chembl <- read_tsv(
  file.path(dir_release, "chembl_uniprot_mapping.txt"), skip = 1,
  col_names = c("uniprot_id", "chembl_id", "pref_name", "target_type"),
  col_types = "cccc"
)

chembl_map_with_uniprot <- chembl_info_all_targets %>%
  left_join(
    map_uniprot_chembl %>%
      dplyr::select(chembl_id, uniprot_id),
    by = "chembl_id"
  )

uniprot_mapping_official <- uniprot_id_mapping_files %>%
  map(fread, sep = "\t", col.names = c("uniprot_id", "external_db", "external_id")) %>%
  rbindlist(use.names = TRUE, fill = TRUE, idcol = "organism")


download.file(
  "ftp://ftp.ncbi.nih.gov/gene/DATA/gene_info.gz",
  file.path(dir_release, "gene_info_20200113.gz")
)

# Using vroom here instead of loading the entire csv because it is downright massive
# and vroom is much faster
gene_info <- vroom(
  file.path(dir_release, "gene_info_20200113.gz"),
  delim = "\t",
  col_names = c(
    "tax_id", "entrez_gene_id", "entrez_symbol", "locus_tag", "entrez_synonyms", "db_xrefs", "chromosome",
    "map_location", "entrez_description", "entrez_type_of_gene", "symbol", "entrez_name", "nomenclature_status",
    "other_designations", "modification_date", "feature_type"
  ),
  col_types = "iic_c___cccc____",
  skip = 1
)

gene_info_relevant <- gene_info %>%
  distinct(
    tax_id, entrez_gene_id, symbol, entrez_symbol, entrez_synonyms, entrez_type_of_gene, entrez_name
  ) %>%
  mutate(
    across(where(is.character), ~if_else(.x == "-", NA_character_, .x))
  )

# Using gene symbols now for inital mapping, recovering more than with entrez IDs directly

tax_id_organism_map <- chembl_map_with_uniprot %>%
  distinct(tax_id, organism)

uniprot_entrez_map <- bind_rows(
  uniprot_mapping_official %>%
    filter(external_db == "GeneID") %>%
    inner_join(
      tax_id_organism_map,
      by = "organism"
    ) %>%
    transmute(
      uniprot_id, entrez_gene_id = as.integer(external_id), tax_id
    ),
  uniprot_mapping_official %>%
    filter(external_db == "Gene_Name") %>%
    inner_join(
      tax_id_organism_map,
      by = "organism"
    ) %>%
    select(
      uniprot_id, symbol = external_id, tax_id
    ) %>%
    inner_join(
      gene_info_relevant %>%
        select(entrez_gene_id, symbol = entrez_symbol, tax_id),
      by = c("tax_id", "symbol")
    ) %>%
    select(uniprot_id, entrez_gene_id, tax_id)
) %>%
  drop_na(entrez_gene_id) %>%
  left_join(
    gene_info_relevant %>%
      select(entrez_gene_id, symbol = entrez_symbol, entrez_synonyms),
    by = "entrez_gene_id"
  ) %>%
  bind_rows(
      filter(
        uniprot_mapping_official,
        !uniprot_id %in% .[["uniprot_id"]],
        external_db == "Gene_Name"
      ) %>%
        inner_join(tax_id_organism_map, by = "organism") %>%
        select(
          uniprot_id, symbol = external_id, tax_id
        ) %>%
        drop_na(symbol)
  ) %>%
  distinct(tax_id, entrez_gene_id, symbol, uniprot_id, entrez_synonyms)

# add back matches where only symbol is available

chembl_map_with_symbol <- chembl_map_with_uniprot %>%
  left_join(
    uniprot_entrez_map,
    by = c("tax_id", "uniprot_id")
  ) %>%
  distinct()


chembl_map_with_symbol %>%
  filter(is.na(entrez_gene_id), is.na(symbol)) %>%
  dplyr::count(organism)
# A tibble: 3 x 2
# organism              n
# <chr>             <int>
# 1 Homo sapiens         24
# 2 Mus musculus          5
# 3 Rattus norvegicus    16


# Only 24 human uniprot IDs left without a matching Entrez ID / Symbol, good enough...

chembl_map_with_symbol %>%
  drop_na(entrez_gene_id) %>%
  dplyr::count(entrez_gene_id) %>%
  dplyr::count(n)
# A tibble: 17 x 2
# n    nn
# <int> <int>
# 1     1  4737
# 2     2   792
# 3     3   237
# 4     4   110
# 5     5    67
# 6     6    27
# 7     7    14
# 8     8     9
# 9     9     8
# 10    10     7
# 11    11     2
# 12    12     4
# 13    13     2
# 14    14     2
# 15    15     1
# 16    19     1
# 17    23     1

# step 4 --> match with gene_info

# 645840 -> 114112
# 348738 -> 6241

map_chemblID_geneID <- chembl_map_with_symbol %>%
  # Append all human gene's not already in to cover genes that are not covered
  # in chembl but only in the LSP single doses etc
  bind_rows(
    filter(
      gene_info_relevant,
      !entrez_gene_id %in% .$entrez_gene_id,
      entrez_type_of_gene == "protein-coding"
    ) %>%
      inner_join(tax_id_organism_map, by = "tax_id") %>%
      transmute(
        tax_id,
        entrez_gene_id,
        symbol = entrez_symbol,
        entrez_synonyms,
        organism,
        pref_name = entrez_name,
        target_type = "SINGLE PROTEIN"
      )
      # left_join(
      #   uniprot_mapping_official %>%
      #     filter(external_db == "GeneID") %>%
      #     transmute(organism, uniprot_id, entrez_gene_id = as.integer(external_id)),
      #   by = c("organism", "entrez_gene_id")
      # )
  ) %>%
  distinct() %>%
  arrange(entrez_gene_id) %>%
  filter(!(is.na(entrez_gene_id) & is.na(symbol)))

# Deduplicated list of targets
map_chemblID_geneID_table <- map_chemblID_geneID %>%
  setDT() %>% {
    .[
      ,
      .(
        lspci_target_id = .GRP,
        pref_name = pref_name %>%
          na.omit() %>%
          magrittr::extract(1)
      ),
      keyby = .(entrez_gene_id, symbol, tax_id, organism)
    ][
      order(entrez_gene_id, symbol)
    ]
  }

# Mapping of all targets to deduplicated table
map_chemblID_geneID_map <- map_chemblID_geneID %>%
  inner_join(
    map_chemblID_geneID_table %>%
      select(lspci_target_id, entrez_gene_id, symbol),
    by = c("entrez_gene_id", "symbol")
  )

fwrite(map_chemblID_geneID_table, file.path(dir_release, "target_dictionary_wide.csv.gz"))

fwrite(map_chemblID_geneID_map, file.path(dir_release, "target_mapping.csv.gz"))

# map_chemblID_geneID <- syn("syn20693721") %>% read_csv()



# Store to synapse -------------------------------------------------------------
###############################################################################T

target_wrangling_activity <- Activity(
  name = "Map and wrangle drug target data",
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/id_mapping/06_target_id_mapping.R"
)

syn_id_mapping <- synMkdir(syn_release, "id_mapping")

list(
  file.path(dir_release, "target_dictionary_wide.csv.gz"),
  file.path(dir_release, "target_mapping.csv.gz")
) %>%
  synStoreMany(syn_id_mapping, activity = target_wrangling_activity, forceVersion = FALSE)

