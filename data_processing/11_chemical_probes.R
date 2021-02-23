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

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

source(here("utils", "load_save.R"))

# set directories & import files -----------------------------------------------
###############################################################################T

inputs <- list(
  compound_dictionary = c("compounds_processed", "compound_dictionary.csv.gz"),
  compound_name_map = c("id_mapping", "lspci_id_compound_name_map.csv.gz")
) %>%
  pluck_inputs(syn_parent = syn_release)

input_data <- inputs %>%
  load_input_data(syn = syn)

old_probes <- syn("syn21627808") %>%
  read_rds() %>%
  chuck("data", 2)

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

fwrite(
  probes_raw,
  file.path(dir_release, "chemical_probes_raw.csv.gz")
)
# probes_raw <- read_csv(file.path(dir_release, "chemical_probes_raw.csv.gz"))

# Map compound and target IDs --------------------------------------------------
###############################################################################T


probes_canonicalized <- probes_raw %>%
  {compounds(set_names(.$smiles, .$`Probe Name`), descriptor = "smiles")} %>%
  canonicalize_compound()

manual_compound_matches <- list(
  "Bafetinib" = 103222L,
  "A-770041" = 103629L
)

sanitize_compund_name <- function(x) {
  x %>%
    str_to_lower() %>%
    str_replace_all("[^\\w]", "")
}

probes_matches <-bind_rows(
  "new" = probes_canonicalized %>%
    left_join(
      input_data[["compound_dictionary"]][
        , .(lspci_id, inchi = canonical_inchi)
      ],
      by = "inchi"
    ),
  "old" = old_probes %>%
    select(compound = name, lspci_id),
  "name" = input_data[["compound_name_map"]] %>%
    transmute(lspci_id, compound_norm = sanitize_compund_name(name)) %>%
    inner_join(
      probes_raw %>%
        transmute(compound = `Probe Name`, compound_norm = sanitize_compund_name(`Probe Name`)),
      by = "compound_norm"
    ),
  .id = "source"
  ) %>%
  distinct(compound, lspci_id, source) %>%
  drop_na()

probes_matches_unique <- probes_matches %>%
  filter(lspci_id %in% input_data[["compound_dictionary"]][["lspci_id"]]) %>%
  mutate(
    source = factor(source, levels = c("new", "name", "old"))
  ) %>%
  arrange(compound, source) %>%
  group_nest(
    compound
  ) %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        filter(source == source[1])
    )
  ) %>%
  unnest(data)

manual_target_matches <- list(
  "HTR2C" = 3358L,
  "BET BD2" = c(6046L, 8019L, 23476L),
  "PDGFRα" = 5156L,
  "PDGFRβ" = 5159L,
  "IDH1" = 3417L,
  "IDH2" = 3418L,
  "KRAS G12C" = 3845L,
  "BCR-ABL" = c(613L, 25L),
  "HSP90" = 3320L,
  "Actin" = 58L,
  "SEC61 complex" = c(29927L, 55176L),
  "Tubulin" = c(7846L),
  "BRAF V600E" = 673L,
  "LXR-alpha" = 10062L,
  "LXR-beta" = 7376L,
  "PCAF" = 8850L,
  "IDH1 R132C/R132H/R132G" = 3417L,
  "PDGFR" = c(5156L, 5159L),
  "Menin" = 4221L,
  "PARPs" = 142L,
  "HIF-2α" = 2034L,
  "MLK1" = 4293L,
  "MLK2" = 4294L,
  "MLK3" = 4296L
)

target_matches <- probes_raw %>%
  transmute(
    target = map(`Target name`, str_split, fixed(", ")) %>%
      map(1)
  ) %>%
  unchop(target) %>%
  distinct() %>%
  genebabel::join_hgnc("target", c("symbol", "alias_symbol", "prev_symbol"), "entrez_id") %>%
  mutate(across(entrez_id, as.integer)) %>%
  bind_rows(
    enframe(
      manual_target_matches, "target", "entrez_id"
    ) %>%
      unchop(entrez_id)
  ) %>%
  rename(entrez_gene_id = entrez_id) %>%
  drop_na()

probes_table <- probes_raw %>%
  transmute(
    name = `Probe Name`,
    target = map(`Target name`, str_split, fixed(", ")) %>%
      map(1),
    n_reviews = `Reviews`,
    rating_cells = `Avg. Rating in Cells`,
    rating_in_vitro = `Avg. Rating in Vivo`,
    max_rating = pmax(`Avg. Rating in Cells`, `Avg. Rating in Vivo`, na.rm = TRUE),
    link
  ) %>%
  unchop(target) %>%
  left_join(
    target_matches,
    by = "target"
  ) %>%
  left_join(
    probes_matches_unique %>%
      distinct(name = compound, lspci_id),
    by = "name"
  )

fwrite(
  probes_table,
  file.path(dir_release, "chemical_probes.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Wrangle chemical probe data",
  used = c(
    unname(inputs),
    "http://www.chemicalprobes.org/browse_probes"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/11_chemical_probes.R"
)

syn_probes <- synMkdir(syn_release, "chemical_probes")

c(
  file.path(dir_release, "chemical_probes_raw.csv.gz"),
  file.path(dir_release, "chemical_probes.csv.gz")
) %>%
  synStoreMany(parent = syn_probes, activity = activity)
