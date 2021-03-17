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
  as.list() %>%
  keep(~str_detect(.x[["name"]], fixed(".fst"))) %>% {
    set_names(
      map_chr(., "id"),
      map_chr(., "name") %>%
        str_replace(fixed(".fst"), "")
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

sigfig <- function(vec, n = 3){
  as.numeric(gsub("\\.$", "", formatC(signif(vec, digits = n), digits = n,format = "fg", flag = "#")))
}

compounds <- input_data[["lsp_compound_dictionary"]][
  ,
  .(
    lspci_id,
    chembl_id,
    emolecules_id,
    pref_name = fcoalesce(
      pref_name,
      chembl_id,
      hmsl_id,
      paste0("Emolecules: ", emolecules_id)
    ),
    commercially_available,
    max_phase
  )
] %>%
  setkey(lspci_id)

inchis <- input_data[["lsp_compound_dictionary"]][
  ,
  .(lspci_id, inchi)
] %>%
  setkey(lspci_id)

compound_names <- rbindlist(list(
  input_data[["lsp_compound_names"]],
  input_data[["lsp_compound_mapping"]][
    source %chin% c("emolecules", "chembl", "hmsl"),
    .(
      lspci_id,
      source = paste0(source, "_id"),
      priority = "secondary",
      name = external_id
    )
  ]
), fill = TRUE)[
  ,
  `:=`(
    source = factor(source, levels = c("chembl", "hmsl", "emolecules", "chembl_id", "hmsl_id", "emolecules_id"))
  )
] %>% unique() %>% {
  .[
    order(lspci_id, name, priority, source)
  ]
} %>% {
  # Only keep identical names from the highest source and priority
  .[
    # Self join. V1 is logical vector of the desired entries
    .[
      ,
      source == source[1] & priority == priority[1],
      by = .(lspci_id, name)
    ]$V1
  ][
    order(lspci_id, priority, source, name)
  ][
    ,
    name_id := paste(lspci_id, seq_len(.N), sep = "-"),
    keyby = .(lspci_id)
  ][
    ,
    # Make compound names unique
    unique_name := if (.N > 1) paste0(name, "{{", seq_len(.N), "}}") else name,
    keyby = .(name)
  ][
    order(
      # Reorder data.table so that some nice normal compounds are at the front
      !name_id %in% c("16241-1", "78621-1", "90319-1", "96316-1", "76418-1", "78036-1", "83706-1", "81903-1", "72090-1", "97590-1"),
      priority, source, name
    )
  ]
} %>%
  select(lspci_id, priority, source, name_id, name = unique_name)

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
  ) %>% {
    .[
      ,
      type := factor(type, levels = c("symbol", "gene_id"))
    ][
      order(lspci_target_id, type, name)
    ][
      ,
      # Make unique target id
      lspci_target_id_unique := paste(lspci_target_id, seq_len(.N), sep = "-"),
      keyby = .(lspci_target_id)
    ][
      # Reorder data.table so that some nice normal targets are at the front
      order(
        type,
        !lspci_target_id %in% c(2525, 487, 34, 1373, 4620, 160, 4950, 1479, 5123, 2407),
        lspci_target_id,
        name
      )
    ]
  }

pfp <- input_data[["lsp_phenotypic_agg"]][
  !is.na(rscore_tr) & is.finite(rscore_tr),
  .(
    lspci_id,
    assay_id,
    rscore_tr
  )
][
  # Remove any assay with less than 6 compounds
  ,
  if (.N >= 6) .SD,
  keyby = .(assay_id)
] %>%
  setkey(assay_id, lspci_id)

tas_measurements <- rbindlist(list(
  input_data[["lsp_one_dose_scan_agg"]][
    !is.na(tas_id),
    .(
      tas_id,
      measurement = sigfig(percent_control, 2),
      unit = paste0("% control at ", concentration, " nM")
    )
  ],
  input_data[["lsp_biochem_agg"]][
    !is.na(tas_id),
    .(
      tas_id,
      measurement = sigfig(as.numeric(value), 2),
      unit = paste0("nM")
    )
  ],
  input_data[["lsp_manual_curation"]][
    !is.na(tas_id),
    .(
      tas_id,
      measurement = NA_real_,
      unit = NA_character_
    )
  ]
))

tas <- input_data[["lsp_tas"]][
  ,
  .(tas_id, lspci_id, lspci_target_id, tas, derived_from)
][
  tas_measurements,
  on = .(tas_id),
  nomatch = NULL
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
  .(lspci_id, lspci_target_id, tas, measurement, unit, derived_from, references)
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
  mutate(
    across(
      c(investigation_bias, wilcox_pval, selectivity, tool_score,
        ic50_difference, ontarget_ic50_q1, offtarget_ic50_q1),
      ~as.numeric(.x) %>%
        sigfig()
    )
  ) %>%
  setkey(lspci_id, lspci_target_id)

library <- input_data[["lsp_compound_library"]][
  ,
  .(lspci_id, lspci_target_id, rank, reason_included)
][
  compounds[
    ,
    .(lspci_id, max_phase)
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
  mutate(
    across(
      c(ontarget_ic50_q1, offtarget_ic50_q1),
      ~as.numeric(.x) %>%
        sigfig()
    )
  ) %>%
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
  merge(
    input_data[["lsp_selectivity"]][
      ,
      .(lspci_id, lspci_target_id, selectivity_class, ontarget_ic50_q1, ontarget_n)
    ],
    by = c("lspci_id", "lspci_target_id"),
    all.x = TRUE
  ) %>%
  merge(
    input_data[["lsp_compound_dictionary"]][
      ,
      .(lspci_id, max_phase)
    ],
    by = c("lspci_id"),
    all.x = TRUE
  ) %>%
  drop_na() %>%
  setkey(lspci_target_id, lspci_id)

# qsave(
#   fingerprints_morgan_normal,
#   file.path(dir_release, "fps.qs"),
#   preset = "fast"
# )
#
# fingerprints <- MorganFPS$new(
#   fingerprints_morgan_normal
# )
#
# fingerprints_small <- MorganFPS$new(
#   fingerprints_morgan_normal[1:100000]
# )
#
# fingerprints_small$save_file(
#   file.path(dir_release, "website_tables", paste0("shiny_fingerprints_small.bin"))
# )
#
#
# library(microbenchmark)
# microbenchmark(
#   fingerprints_small$save_file(
#     file.path(dir_release, "website_tables", paste0("shiny_fingerprints_small.bin")),
#     compression_level = 10
#   ), times = 3L
# )
#
#
# fingerprints$save_file(
#   file.path(dir_release, "website_tables", paste0("shiny_fingerprints.bin"))
# )
#
# x <- MorganFPS$new(
#   file.path(dir_release, "website_tables", paste0("shiny_fingerprints_small.bin")),
#   from_file = TRUE
# )
#
# library(microbenchmark)
# microbenchmark(
#   x <- MorganFPS$new(
#     file.path(dir_release, "website_tables", paste0("shiny_fingerprints_small.bin")),
#     from_file = TRUE
#   ), times = 3L
# )
#
# x <- MorganFPS$new(
#   file.path(dir_release, "website_tables", paste0("shiny_fingerprints.bin")),
#   from_file = TRUE
# )


# Chemfp FBP fingerprint file --------------------------------------------------
###############################################################################T
#
# unloadNamespace("synExtra")
# unloadNamespace("synapser")
# unloadNamespace("PythonEmbedInR")
#
# library(reticulate)
# use_python("~/miniconda/envs/lspcheminf/bin/python3.7", required = TRUE)
#
# main <- import_main()
# builtins <- import_builtins()
# chemfp <- import("chemfp")
#
# arena <- chemfp$load_fingerprints(py$, format = "fpb.gz", allow_mmap = TRUE)

# Restrict compound space to compounds with annotation -------------------------
###############################################################################T

eligible_lspci_ids <- c(
  compounds[
    !is.na(chembl_id)
  ][["lspci_id"]],
  pfp[["lspci_id"]],
  selectivity[["lspci_id"]],
  library[["lspci_id"]],
  tas[["lspci_id"]],
  chemical_probes[["lspci_id"]]
) %>%
  unique()

qsave(
  eligible_lspci_ids,
  file.path(dir_release, "eligible_lspci_ids.qs")
)

# eligible_lspci_ids <- qread(file.path(dir_release, "eligible_lspci_ids.qs"))

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
    ) %>%
      # Only take coompounds with annotation
      map(
        ~{
          if ("lspci_id" %in% colnames(.x))
            .x[
              lspci_id %in% eligible_lspci_ids
            ]
          else
            .x
        }
      )
  )

# all_tables <- all_table_names %>%
#   mutate(
#     table = map(
#       path,
#       read_fst,
#       as.data.table = TRUE
#     )
#   )

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

# Fingerprints -----------------------------------------------------------------
###############################################################################T

fps <- input_data[["lsp_fingerprints"]]

fingerprints_morgan_normal <- fps[
  fingerprint_type == "morgan_normal" & lspci_id %in% eligible_lspci_ids
] %>% {
  set_names(.[["fingerprint"]], .[["lspci_id"]])
}

fingerprints <- MorganFPS$new(
  fingerprints_morgan_normal
)

fingerprints$save_file(
  file.path(dir_release, "website_tables", paste0("shiny_fingerprints.bin")),
  compression_level = 22
)

# f <- MorganFPS$new(
#   file.path(dir_release, "website_tables", paste0("shiny_fingerprints.bin")),
#   from_file = TRUE
# )

# Upload to synapse ------------------------------------------------------------
###############################################################################T

activity <- Activity(
  name = "Wrangle tables for SMS website",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/website/01_shiny_tables.R"
)

syn_dir <- synMkdir(syn_release, "website_tables")

synStoreMany(
  c(
    all_tables[["path"]],
    file.path(dir_release, "website_tables", paste0("shiny_fingerprints.bin"))
  ),
  syn_dir,
  activity = activity,
  forceVersion = FALSE
)
