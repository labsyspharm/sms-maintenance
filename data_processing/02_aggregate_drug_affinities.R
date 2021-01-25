library(tidyverse)
library(data.table)
library(here)
library(synapser)
library(synExtra)
library(qs)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")


# Set directories, import files ------------------------------------------------
###############################################################################T

inputs <- list(
  lspci_id_vendor_id_map = c("compounds_processed", "lspci_id_vendor_id_map.csv.gz"),
  chembl_ref_info = c("raw_data", "chembl", "chembl_ref_info_raw.csv.gz"),
  chembl_biochemical = c("raw_data", "chembl", "chembl_biochemical_raw.csv.gz"),
  chembl_phenotypic = c("raw_data", "chembl", "chembl_phenotypic_raw.csv.gz"),
  chembl_biochemical = c("raw_data", "chembl", "chembl_biochemical_raw.csv.gz")
) %>%
  map(~exec(synPluck, !!!c(syn_release, .x))) %>%
  c(
    inhouse_dose_response = "syn20692433",
    inhouse_single_dose = "syn20692432"
  )

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

# Clean up raw ChEMBL data -----------------------------------------------------
###############################################################################T

REFERENCE_PRIORITY <- c(
  "pubmed_id",
  "doi",
  "patent_id",
  "synapse_id",
  "chembl_id"
)

chembl_ref_info_best <- copy(input_data[["chembl_ref_info"]])[
  ,
  reference_type := factor(reference_type, levels = REFERENCE_PRIORITY)
][
  order(reference_type),
  .(
    references = .SD[
      # Remove DOI if Pubmed ID is present
      if (any(reference_type == "pubmed_id"))
        !reference_type == "doi"
      else
        TRUE,
      .(reference_type, reference_id)
    ] %>%
      unique() %>%
      list()
  ),
  keyby = .(chembl_id_doc)
]

biochem_neat <- copy(input_data[["chembl_biochemical"]]) %>%
  mutate(standard_value = as.numeric(standard_value)) %>%
  rename(pref_name_target = pref_name) %>%
  filter(tax_id == 9606, !is.na(entrez_gene_id)) %>%
  inner_join(
   copy(input_data[["lspci_id_vendor_id_map"]])[
     source == "chembl"
   ][
     ,
     source := NULL
   ],
   by = c("chembl_id_compound" = "vendor_id")
  )

dose_response_inhouse_neat <- copy(input_data[["inhouse_dose_response"]]) %>%
  mutate(hms_id = paste0("HMSL", hms_id)) %>%
  rename(entrez_gene_id = gene_id, hmsl_id = hms_id) %>%
  inner_join(
    copy(input_data[["lspci_id_vendor_id_map"]])[
      source == "hmsl"
    ][
      ,
      source := NULL
    ],
    by = c("hmsl_id" = "vendor_id")
  )

pheno_data_neat <- copy(input_data[["chembl_phenotypic"]]) %>%
  mutate(standard_value = as.numeric(standard_value)) %>%
  inner_join(
    copy(input_data[["lspci_id_vendor_id_map"]])[
      source == "chembl"
    ][
      ,
      source := NULL
    ],
    by = c("chembl_id_compound" = "vendor_id")
  )

# Aggregate dose-response data -------------------------------------------------
###############################################################################T

biochem_rowbind <- biochem_neat %>%
  left_join(
    chembl_ref_info_best,
    by = "chembl_id_doc"
  ) %>%
  mutate(
    file_url = paste0("https://www.ebi.ac.uk/chembl/document_report_card/", chembl_id_doc),
    value = standard_value,
    value_unit = standard_units,
    references = map2(
      references, chembl_id_doc,
      ~{
        if (!is.null(.x))
          .x
        else
          data.table(
            reference_id = .y,
            reference_type = factor("chembl_id", levels = REFERENCE_PRIORITY)
          )
      }
    )
  ) %>%
  select(
    lspci_id,
    value, value_unit = standard_units, value_type = standard_type,
    value_relation = standard_relation, description_assay = description,
    uniprot_id, entrez_gene_id,
    references, file_url
  )

doseresponse_inhouse_rowbind <- input_data[["inhouse_dose_response"]] %>%
  mutate(
    hms_id = paste0("HMSL", hms_id)
  ) %>%
  inner_join(
    input_data[["lspci_id_vendor_id_map"]] %>%
      filter(source == "hmsl") %>%
      select(lspci_id, vendor_id),
    by = c("hms_id" = "vendor_id")
  ) %>%
  mutate(
    references = map(
      synapse_id,
      ~data.table(
        reference_id = .x,
        reference_type = factor("synapse_id", levels = REFERENCE_PRIORITY)
      )
    ),
    reference_id = synapse_id,
    reference_type = "synapse_id"
  ) %>%
  select(
    lspci_id,
    value, value_unit, value_type,
    value_relation, description_assay = description,
    uniprot_id, entrez_gene_id = gene_id,
    references, file_url
  )

biochem_complete <- bind_rows(
  biochem_rowbind,
  doseresponse_inhouse_rowbind
) %>%
  # Call to distinct important, since some assays can be recorded multiple times
  # for the same eq_class now, when multiple forms of the same drug where mapped
  # to the same eq_class and an assay was stored in the db for all forms
  distinct() %>%
  # Remap obsolete entrez_ids
  # 645840 -> 114112
  # 348738 -> 6241
  mutate(
    entrez_gene_id = recode(
      entrez_gene_id, `645840` = 114112L, `348738` = 6241L
    )
  )

qsave(
  biochem_complete,
  file.path(dir_release, "biochemical_data_complete.qs"),
  preset = "fast"
)

# Using data.table here for speed
calculate_q1 <- function(data) {
  n_groups <- uniqueN(data, by = c("lspci_id", "entrez_gene_id"))
  pb <- txtProgressBar(min = 1, max = n_groups, style = 3)
  on.exit(close(pb))
  data %>%
    as.data.table() %>% {
      .[
        ,
        .(
          Q1 = {
            setTxtProgressBar(pb, .GRP)
            round(quantile(value, 0.25, names = FALSE), 2)
          },
          n_measurement = .N,
          references = references %>%
            rbindlist(use.names = TRUE) %>%
            unique() %>% {
              paste(.[["reference_type"]], .[["reference_id"]], sep = ":")
            } %>%
            paste(collapse = "|")
        ),
        by = .(lspci_id, entrez_gene_id)
      ]
    }
}

biochem_complete_q1 <- biochem_complete %>%
  setkey(lspci_id, entrez_gene_id) %>%
  head(n = 10000) %>%
  calculate_q1()

fwrite(
  biochem_complete_q1,
  file.path(dir_release, "biochemical_data_complete_q1.csv.gz")
)

# Aggregate single-dose data ---------------------------------------------------
###############################################################################T


# Also calculate Q1 values for kinomescan data from HMS LINCS for which no
# complete dose response curve is available
hmsl_kinomescan <- syn("syn20692432") %>%
  read_csv(col_types = "cccdcdcc") %>%
  # Some of the HMSL IDs don't start with "HMSL", fix that
  mutate(hms_id = if_else(str_starts(hms_id, "HMSL"), hms_id, paste0("HMSL", hms_id)))

# check how common the situation is where a single gene was tested in different
# variants (mutants, post-translational modification, etc.)
hmsl_kinomescan %>%
  count(hms_id, gene_symbol, cmpd_conc_nM) %>%
  count(n)
# # A tibble: 21 x 2
# n    nn
# <int> <int>
#   1     1 60448
# 2     2  4467
# 3     3   283
# 4     4   250
# 5     5    32
# 6     6    57
# 7     7    56
# 8     8    88
# 9     9   122
# 10    10   161
# # … with 11 more rows
# It's very common, so we have to decide how to aggregate info for every
# target/compound combo.
# I think it makes sense to try to filter out mutant data, but keep post-translational
# modified targets (especially phosphorylation), because those are likely to be
# broadly relevant. if somebody is interested in specific mutations they will
# have to query the Kinomescan directly and shoulnd't rely on TAS.

hmsl_kinomescan_cleaned <- hmsl_kinomescan %>%
  filter(
    # Filter mutant annotation e.g. S154A
    !grepl("[^A-Za-z]+[A-Z][0-9]+([A-Z]|del|Del)[^A-Za-z]", description),
    !grepl("domain", description),
    !grepl("Dom", description),
    !grepl("inhibited", description)
  )

hmsl_kinomescan_mapped <- all_cmpds_eq_classes %>%
  mutate(
    data = map(
      data,
      ~hmsl_kinomescan_cleaned %>%
        genebabel::join_hgnc(
          "gene_symbol",
          c("symbol", "alias_symbol", "prev_symbol"),
          c("entrez_id", "name", "uniprot_ids")
        ) %>%
        # I checked, no gene_symbol maps to multiple uniprot, so this is safe
        mutate(
          uniprot_id = as.character(uniprot_ids),
          # Either synapse ID or HMSL ID
          reference_type = if_else(str_starts(source_assay_id, fixed("syn")), "synapse", "hms_lincs")
        ) %>%
        select(-uniprot_ids) %>%
        left_join(
          .x %>%
            select(id, eq_class),
          by = c("hms_id" = "id")
        ) %>%
        rename(
          entrez_gene_id = entrez_id,
          pref_name_cmpd = pref_name,
          pref_name_target = name,
          lspci_id = eq_class,
          reference_id = source_assay_id,
          file_url = url
        ) %>%
        drop_na(entrez_gene_id) %>%
        # Remap obsolete entrez_ids
        # 645840 -> 114112
        # 348738 -> 6241
        mutate(entrez_gene_id = recode(as.integer(entrez_gene_id), `645840` = 114112L, `348738` = 6241L))
    )
  )

write_rds(
  hmsl_kinomescan_mapped,
  file.path(dir_release, "biochemicaldata_single_dose_inhouse.rds"),
  compress = "gz"
)

hmsl_kinomescan_q1 <- hmsl_kinomescan_mapped %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        as.data.table() %>%
        {
          .[
            ,
            reference_type := recode(
              reference_type,
              pubmed_id = "pubmed", patent_id = "patent", chembl_id = "chembl", synapse_id = "synapse"
            )
          ][
            ,
            .(
              percent_control_Q1 = quantile(percent_control, 0.25, names = FALSE),
              n_measurement = .N,
              references = .SD[, .(reference_type, reference_id)] %>%
                unique() %>%
                with(
                  paste(
                    reference_type,
                    reference_id,
                    sep = ":",
                    collapse = "|"
                  )
                )
            ),
            by = .(lspci_id, entrez_gene_id, cmpd_conc_nM)
          ]
        } %>%
        as_tibble()
    )
  )


write_rds(
  hmsl_kinomescan_q1,
  file.path(dir_release, "biochemicaldata_single_dose_inhouse_Q1.rds"),
  compress = "gz"
)

# Aggregate phenotypic data ----------------------------------------------------
###############################################################################T

pheno_data_formatted <- pheno_data_neat %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        left_join(
          select(chembl_ref_info_best, chembl_id_doc, reference_type, reference_id),
          by = "chembl_id_doc"
        ) %>%
        mutate(
          reference_id = if_else(is.na(reference_id), chembl_id_doc, reference_id),
          reference_type = if_else(is.na(reference_type), "chembl_id", as.character(reference_type)),
          file_url = paste0("https://www.ebi.ac.uk/chembl/document_report_card/", chembl_id_doc)
        ) %>%
        select(
          lspci_id, assay_id,
          value = standard_value, value_unit = standard_units, value_type = standard_type,
          value_relation = standard_relation, description_assay = description,
          reference_id, reference_type, file_url
        )
    )
  )

write_rds(
  pheno_data_formatted,
  file.path(dir_release, "pheno_data.rds"),
  compress = "gz"
)

pheno_data_q1 <- pheno_data_formatted %>%
  mutate(
    data = map(
      data,
      ~.x %>%
        as.data.table() %>%
        {
          .[
            ,
            reference_type := recode(
              reference_type,
              pubmed_id = "pubmed", patent_id = "patent", chembl_id = "chembl", synapse_id = "synapse"
            )
            ][
            ,
            .(
              standard_value_Q1 = quantile(value, 0.25, names = FALSE),
              log10_value_Q1 = quantile(log10(value), 0.25, names = FALSE),
              n_measurement = .N,
              references = .SD[, .(reference_type, reference_id)] %>%
                unique() %>%
                with(
                  paste(
                    reference_type,
                    reference_id,
                    sep = ":",
                    collapse = "|"
                  )
                )
            ),
            by = .(lspci_id, assay_id)
          ]
        } %>%
        as_tibble()
    )
  )

write_rds(
  pheno_data_q1,
  file.path(dir_release, "pheno_data_Q1.rds"),
  compress = "gz"
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

aggregation_activity <- Activity(
  name = "Aggregate affinity data",
  used = c(
    "syn20692432",
    "syn20692433",
    "syn20693825",
    "syn20693827",
    "syn20830516",
    "syn21652213"
  ),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/02_aggregate_drug_affinities.R"
)

syn_aggregate <- Folder("aggregate_data", parent = syn_release) %>%
  synStore() %>%
  chuck("properties", "id")

c(
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl.rds"),
  file.path(dir_release, "biochemicaldata_complete_inhouse_chembl_Q1.rds"),
  file.path(dir_release, "biochemicaldata_single_dose_inhouse.rds"),
  file.path(dir_release, "biochemicaldata_single_dose_inhouse_Q1.rds"),
  file.path(dir_release, "pheno_data.rds"),
  file.path(dir_release, "pheno_data_Q1.rds")
) %>%
  synStoreMany(parent = syn_aggregate, activity = aggregation_activity)
