library(tidyverse)
library(data.table)
library(fst)
library(here)
library(synapser)
library(synExtra)
library(morgancpp)

synLogin()
syn <- synDownloader(here("tempdl"), followLink = TRUE)

source(here("utils", "load_save.R"))

# set directories, import files ------------------------------------------------
###############################################################################T
release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")

# Function to replace empty strings with NA values
replace_empty_string_na <- function(df) {
  mutate(
    df,
    across(
      where(is.character),
      ~magrittr::inset(.x, .x == "", NA_character_)
    )
  )
}

inputs <- synPluck(syn_release, "db_tables") %>%
  synGetChildren() %>%
  as.list() %>% {
    set_names(
      map_chr(., "id"),
      map_chr(., "name") %>%
        str_replace(fixed(".csv.gz"), "")
    )
  } %>%
  c(
    pluck_inputs(
      list(
        chemical_probes = c("chemical_probes", "chemical_probes.csv.gz")
      ),
      syn_parent = syn_release
    )
  )

input_data <- inputs %>%
  # magrittr::extract(names(.) != "lsp_fingerprints") %>%
  load_input_data(syn = syn) %>%
  map(replace_empty_string_na)

# Create tables ----------------------------------------------------------------
###############################################################################T

compounds <- input_data[["lsp_compound_dictionary"]][
  ,
  .(lspci_id, chembl_id, emolecules_id, pref_name, commercially_available, highest_approval)
] %>%
  setkey(lspci_id)

inchis <- input_data[["lsp_compound_dictionary"]][
  ,
  .(lspci_id, inchi)
] %>%
  setkey(lspci_id)

compound_names <- copy(input_data[["lsp_compound_names"]])[
  ,
  `:=`(
    source = factor(source, levels = c("chembl", "hsml", "emolecules")),
    priority = NULL
  )
] %>% unique() %>% {
  .[
    order(lspci_id, name, source)
  ][
    ,
    .SD[source == source[1]],
    keyby = .(lspci_id, name)
  ][
    ,
    name_id := paste(lspci_id, seq_len(.N), sep = "_"),
    keyby = .(lspci_id)
  ][
    ,
    # Make compound names unique
    unique_name := if (.N > 1) paste0(name, "{{", seq_len(.N), "}}") else name,
    keyby = .(name)
  ]
} %>%
  select(lspci_id, source, name_id, name = unique_name) %>%
  setkey(name)

targets <- input_data[["lsp_target_dictionary"]][
  ,
  .(
    lspci_target_id,
    gene_id,
    symbol,
    tax_id
  )
] %>%
  setkey(lspci_target_id)

target_map <- targets %>%
  melt(
    measure.vars = c("gene_id", "symbol"),
    variable.name = "type",
    value.name = "name",
    na.rm = TRUE
  ) %>%
  setkey(name)

pfp <- input_data[["lsp_phenotypic_agg"]][
  !is.na(rscore_tr) & is.finite(rscore_tr),
  .(lspci_id, assay_id, rscore_tr)
][
  # Remove any assay with less than 6 compounds
  ,
  if (.N >= 6) .SD,
  keyby = .(assay_id)
] %>%
  setkey(assay_id, lspci_id)

tas <- input_data[["lsp_tas"]][
  ,
  .(tas_id, lspci_id, lspci_target_id, tas)
][
  input_data[["lsp_tas_references"]][
    input_data[["lsp_references"]][
      ,
      .(reference_id, reference_type, reference_value)
    ],
    on = .(reference_id),
    nomatch = NULL
  ][
    ,
    .(
      tas_id,
      references = paste(
        reference_type, reference_value, sep = ":", collapse = "|"
      )
    ),
    keyby = .(tas_id)
  ],
  on = .(tas_id),
  nomatch = NULL
][
  ,
  .(lspci_id, lspci_target_id, tas, references)
] %>%
  setkey(lspci_id, lspci_target_id)

selectivity <-  input_data[["lsp_selectivity"]][
  ,
  .(
    lspci_id, lspci_target_id, selectivity_class, investigation_bias,
    strength, wilcox_pval, selectivity, tool_score, ic50_difference,
    ontarget_ic50_q1, offtarget_ic50_q1, ontarget_n, offtarget_n
  )
][
  input_data[["lsp_biochem_agg"]][
    ,
    .(biochem_agg_id, lspci_id, lspci_target_id)
  ][
    input_data[["lsp_biochem"]][
      ,
      .(biochem_agg_id, reference_id)
    ],
    on = .(biochem_agg_id),
    nomatch = NULL
  ][
    input_data[["lsp_references"]],
    on = .(reference_id),
    nomatch = NULL
  ][
    ,
    .(
      references = paste(
        reference_type, reference_value, sep = ":", collapse = "|"
      )
    ),
    keyby = .(lspci_id, lspci_target_id)
  ],
  on = .(lspci_id, lspci_target_id),
  nomatch = NULL
] %>%
  setkey(lspci_id, lspci_target_id)

library <- input_data[["lsp_compound_library"]][
  ,
  .(lspci_id, lspci_target_id, rank, reason_included)
][
  compounds[
    ,
    .(lspci_id, highest_approval)
  ],
  on = .(lspci_id),
  nomatch = NULL
][
  rbindlist(list(
    selectivity[
      ,
      .(
        lspci_id, lspci_target_id, selectivity_class,
        ontarget_ic50_q1, offtarget_ic50_q1,
        ontarget_n
      )
    ],
    # Add affinities where no selectivity could be calculated
    anti_join(
      input_data[["lsp_biochem_agg"]],
      selectivity,
      by = c("lspci_id", "lspci_target_id")
    ) %>% {
      .[
        input_data[["lsp_biochem"]][
          ,
          .(ontarget_n = .N),
          keyby = .(lspci_id, lspci_target_id)
        ],
        on = .(lspci_id, lspci_target_id),
        nomatch = NULL
      ][
        ,
        .(
          lspci_id, lspci_target_id,
          ontarget_ic50_q1 = value,
          ontarget_n
        )
      ]
    }
  ), use.names = TRUE, fill = TRUE),
  on = .(lspci_id, lspci_target_id),
  nomatch = NULL
] %>%
  setkey(lspci_target_id, rank)

chemical_probes <- input_data[["chemical_probes"]][
  input_data[["lsp_target_dictionary"]][
    !is.na(gene_id),
    .(lspci_target_id, gene_id)
  ],
  on = c("entrez_gene_id" = "gene_id"),
  nomatch = NULL
][
  ,
  .(lspci_id, lspci_target_id, n_reviews, max_rating, link)
] %>%
  drop_na() %>%
  setkey(lspci_target_id, lspci_id)

fps <- syn("syn24874143") %>%
  fread()

fingerprints_morgan_normal <- fps[
  fingerprint_type == "morgan_normal"
] %>% {
  set_names(.[["fingerprint"]], .[["lspci_id"]])
}

qsave(
  fingerprints_morgan_normal,
  file.path(dir_release, "fps.qs"),
  preset = "fast"
)

fingerprints <- MorganFPS$new(
  fingerprints_morgan_normal
)

fingerprints_small <- MorganFPS$new(
  fingerprints_morgan_normal[1:100000]
)

fingerprints_small$save_file(
  file.path(dir_release, "website_tables", paste0("shiny_fingerprints_small.bin")),
  compression_level = 9
)


library(microbenchmark)
microbenchmark(
  fingerprints_small$save_file(
    file.path(dir_release, "website_tables", paste0("shiny_fingerprints_small.bin")),
    compression_level = 10
  ), times = 3L
)


fingerprints$save_file(
  file.path(dir_release, "website_tables", paste0("shiny_fingerprints.bin")),
  compression_level = 9
)

x <- MorganFPS$new(
  file.path(dir_release, "website_tables", paste0("shiny_fingerprints_small.bin")),
  from_file = TRUE
)

library(microbenchmark)
microbenchmark(
  x <- MorganFPS$new(
    file.path(dir_release, "website_tables", paste0("shiny_fingerprints_small.bin")),
    from_file = TRUE
  ), times = 3L
)

x <- MorganFPS$new(
  file.path(dir_release, "website_tables", paste0("shiny_fingerprints.bin")),
  from_file = TRUE
)

# Assemble tables --------------------------------------------------------------
###############################################################################T

all_table_names <- tribble(
  ~name,
  "compounds",
  "inchis",
  "compound_names",
  "targets",
  "pfp",
  "tas",
  "selectivity",
  "library",
  "chemical_probes",
  "target_map"
) %>%
  mutate(
    path = file.path(dir_release, "website_tables", paste0("shiny_", name, ".fst"))
  )

all_tables <- all_table_names %>%
  mutate(
    table = map(
      name,
      get,
      env = .GlobalEnv
    )
  )

all_tables <- all_table_names %>%
  mutate(
    table = map(
      path,
      read_fst,
      as.data.table = TRUE
    )
  )

dir.create(file.path(dir_release, "website_tables"))

pwalk(
  all_tables,
  function(name, path, table, ...) {
    message(name)
    write_fst(
      table,
      path,
      compress = 100
    )
  }
)

# Upload to synapse ------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Wrangle tables for SMS website",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/01_shiny_tables.R"
)

syn_dir <- synMkdir(syn_release, "website_tables")

synStoreMany(
  all_tables[["path"]],
  syn_dir,
  activity = activity,
  forceVersion = FALSE
)
