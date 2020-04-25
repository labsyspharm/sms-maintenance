library(tidyverse)
library(here)
library(rvest)
library(furrr)
library(lspcheminf)
library(data.table)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v25"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# set directories & import files -----------------------------------------------
###############################################################################T

compound_table <- syn("syn20835543") %>%
  read_rds()

# Retrieve probes --------------------------------------------------------------
###############################################################################T

probes_html <- xml2::read_html("http://www.chemicalprobes.org/browse_probes")

probe_links <- probes_html %>%
  html_nodes("td.views-field a") %>%
  {tibble(link = paste0("http://www.chemicalprobes.org", html_attr(., "href")), compound = html_text(.))}

extract_smiles <- function(link) {
  probe_html <- xml2::read_html(link)
  smiles_node <- probe_html %>%
    html_node(xpath = "//div[contains(text(), 'SMILES')]/../node()[not(self::div)]")
  smiles_node %>%
    html_text(trim = TRUE)
}

probes_raw <- probes_html %>%
  html_node(".views-table") %>%
  html_table() %>%
  inner_join(
    probe_links,
    by = c("Probe Name" = "compound")
  ) %>%
  mutate(
    smiles = map_chr(
      link, extract_smiles
    )
  ) %>%
  # Add SMILES for S63845 which is not available on the website
  mutate(
    smiles = if_else(
      `Probe Name` == "S63845",
      "CC1=C(C=CC(=C1Cl)OCCN2CCN(CC2)C)C3=C(SC4=NC=NC(=C34)OC(CC5=CC=CC=C5OCC6=CC=NN6CC(F)(F)F)C(=O)O)C7=CC=C(O7)F",
      smiles
    )
  )

write_csv(
  probes_raw,
  file.path(dir_release, "chemical_probes_raw.csv.gz")
)
# probes_raw <- read_csv(file.path(dir_release, "chemical_probes_raw.csv.gz"))

# Map compound and target IDs --------------------------------------------------
###############################################################################T


probes_canonicalized <- probes_raw %>%
  {set_names(.$smiles, .$`Probe Name`)} %>%
  convert_compound_identifier(identifier = "smiles", target_identifier = "inchi") %>%
  {set_names(.$compounds, .$names)} %>%
  canonicalize_compound()

find_matches <- function(query, targets) {
  targets %>%
    drop_na(lspci_id, inchi) %>%
    chunk_df(12) %>%
    map(~set_names(.x$inchi, .x$lspci_id)) %>%
    future_map(
      ~compound_identity(set_names(query$inchi, query$compound), .x),
      .progress = TRUE
    ) %>%
    bind_rows()
}

plan(multicore(workers = 4))
probe_matches <- compound_table %>%
  mutate(
    data = map(
      data,
      ~find_matches(query = probes_canonicalized, targets = .x)
    )
  )

probes_table <- probes_raw %>%
  transmute(
    name = `Probe Name`,
    target = map(`Protein target`, str_split, fixed(", ")) %>%
      map(1),
    n_reviews = `Number of SAB Reviews`,
    avg_rating = `Avg Rating (in cells)`,
    link,
    smiles
  ) %>%
  unchop(target) %>%
  genebabel::join_hgnc("target", c("symbol", "alias_symbol"), "entrez_id")


chemical_probes <- probe_matches %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        transmute(name = query, lspci_id = as.integer(target)) %>%
        right_join(probes_table, by = "name")
    )
  )

write_rds(
  chemical_probes,
  file.path(dir_release, "chemical_probes.rds")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Wrangle chemical probe data",
  used = c(
    "syn20835543",
    "http://www.chemicalprobes.org/browse_probes"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/11_chemical_probes.R"
)

syn_probes <- Folder("chemical_probes", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "chemical_probes_raw.csv.gz"),
  file.path(dir_release, "chemical_probes.rds")
) %>%
  synStoreMany(parent = syn_probes, activity = activity)
