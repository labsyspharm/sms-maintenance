library(tidyverse)
library(vroom)
library(data.table)
library(qs)
library(biomaRt)
library(bit64)
library(RPostgres)
library(here)
library(synapser)
library(synExtra)

synLogin()
syn <- synDownloader(here("tempdl"))

release <- "chembl_v27"
dir_release <- here(release)
syn_release <- synFindEntityId(release, "syn18457321")


# Set directories, import files ------------------------------------------------
###############################################################################T

inputs <- list(
  target_map = synPluck(syn_release, "id_mapping", "target_mapping.csv.gz")
)

input_data <- inputs %>%
  map(syn) %>%
  map(
    function(x)
      list(
        `.csv` = partial(fread, colClasses = c(inchi_id = "integer")),
        `.tsv` = fread,
        `.rds` = read_rds
      ) %>%
      magrittr::extract2(which(str_detect(x, fixed(names(.))))) %>%
      {.(x)}
  )


# connect to chembl v. 24_1  ---------------------------------------------------
###############################################################################T


## in terminal: ssh -L 5433:pgsql96.orchestra:5432 nm192@transfer.rc.hms.harvard.edu
# first portnumber can change
# loads the PostgreSQL driver
drv <- dbDriver("Postgres")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = str_replace(release, fixed("v"), ""),
                 host = "localhost", port = 5432,
                 user = "chug")


# get biochemical data all compounds -------------------------------------------
###############################################################################T

BAO_format <- dbGetQuery(con, paste0("select *
                        from bioassay_ontology
                                   "))
#where bao_id in ('",paste(biochem_test$bao_format%>%unique,collapse="','"),"')
View(BAO_format)

# Normalize all units to nM, saving the factors here
standard_unit_map <- c(
  'M' = 1,
  'mol/L' = 1,
  'nM' = 10^-9,
  'nmol/L' = 10^-9,
  'nmol.L-1' = 10^-9,
  'pM' = 10^-12,
  'pmol/L' = 10^-12,
  'pmol/ml' = 10^-9,
  'um' = 10^-6,
  'uM' = 10^-6,
  'umol/L' = 10^-6,
  'umol/ml' = 10^-3,
  'umol/uL' = 1
) %>%
  magrittr::multiply_by(10^9)


activities_biochem_1 <- dbGetQuery(
  con,
  paste0(
    "select A.doc_id, ACT.activity_id, A.assay_id, ACT.molregno, MOL_DICT.chembl_id as chembl_id_compound, ACT.standard_relation, ACT.standard_type,
     ACT.standard_value, ACT.standard_units,
     A.tid,
     A.description,A.chembl_id as chembl_id_assay, BAO.label, DOCS.chembl_id as chembl_id_doc, DOCS.pubmed_id as pubmed_id, TARGET_DICT.organism AS organism
     from activities as ACT
     left join assays as A
     on ACT.assay_id = A.assay_id
     left join DOCS
     on DOCS.doc_id=A.doc_id
     left join bioassay_ontology as BAO
     on A.bao_format=BAO.bao_id
     LEFT JOIN molecule_dictionary AS MOL_DICT
     on ACT.molregno = MOL_DICT.molregno
     LEFT JOIN target_dictionary AS TARGET_DICT
     on A.tid = TARGET_DICT.tid
     WHERE ACT.standard_value is not null
     AND TARGET_DICT.organism in ('Homo sapiens', 'Rattus norvegicus', 'Mus musculus')
     and A.assay_type = 'B'
     and A.relationship_type in ('D', 'H', 'M', 'U')
     and A.bao_format not in ('BAO_0000221', 'BAO_0000219','BAO_0000218')
     and ACT.standard_units in (", paste(paste0("'", names(standard_unit_map), "'"), collapse = ","), ")
     and ACT.standard_type in ('IC50','Ki','EC50','Kd','IC90','CC50','ID50','AC50','Inhibition','MIC','Potency','Activity','ED50')
     and A.assay_cell_type is NULL"
  )
)



# Data from "Navigating the Kinome" paper is annotated using F (functional) assay type
# whereas most other data B (binding)
activities_biochem_2<-dbGetQuery(con, paste0("select A.doc_id, ACT.activity_id, A.assay_id, ACT.molregno, MOL_DICT.chembl_id as chembl_id_compound, ACT.standard_relation, ACT.standard_type,
                                             ACT.standard_value,ACT.standard_units,
                                             A.tid,
                                             A.description,A.chembl_id as chembl_id_assay, BAO.label, DOCS.chembl_id as chembl_id_doc, DOCS.pubmed_id as pubmed_id, TARGET_DICT.organism AS organism
                                             from activities as ACT
                                             left join assays as A
                                             on ACT.assay_id = A.assay_id
                                             left join DOCS
                                             on DOCS.doc_id=A.doc_id
                                             left join bioassay_ontology as BAO
                                             on A.bao_format=BAO.bao_id
                                             LEFT JOIN molecule_dictionary AS MOL_DICT
                                             on ACT.molregno = MOL_DICT.molregno
                                             LEFT JOIN target_dictionary AS TARGET_DICT
                                             ON A.tid = TARGET_DICT.tid
                                             WHERE ACT.standard_value is not null
                                             AND TARGET_DICT.organism in ('Homo sapiens', 'Rattus norvegicus', 'Mus musculus')
                                             and A.assay_type = 'F'
                                             and A.description like '%Navigating the Kinome%'
                                             "))

activities_biochem <- data.table::rbindlist(
  list(
    activities_biochem_1 %>%
      filter(standard_value > 0) %>%
      mutate(
        standard_value = standard_value / standard_unit_map[standard_units],
        # log10_value = log10(standard_value * standard_unit_map[standard_units]),
        standard_units = "nM"
      ),
    activities_biochem_2
  )
)

target_dict_relevant <- input_data[["target_map"]] %>%
  distinct(
    lspci_target_id,
    tid,
    chembl_id_target = chembl_id,
    target_type,
    pref_name,
    uniprot_id,
    symbol,
    entrez_gene_id,
    tax_id,
    organism
  )

activities_biochem_geneid <- activities_biochem %>%
  # Only two targets are not present in dictionary, both are DNA. Not relevant
  inner_join(
    target_dict_relevant,
    by = c("organism", "tid")
  )


fwrite(
  activities_biochem_geneid,
  file.path(dir_release, "chembl_biochemical_raw.csv.gz")
)

# get phenotypic data all compounds --------------------------------------------
###############################################################################T


# assay_freq<-dbGetQuery(con, paste0(" select assay_id, count(molregno)
#                                    from activities
#                                    where standard_units in ('",paste(standard_units_ok,collapse="','"),"')
#                                    group by assay_id
#                                    "))

activities_1<-dbGetQuery(con, paste0("select A.doc_id, ACT.activity_id, A.assay_id, ACT.molregno, MOL_DICT.chembl_id as chembl_id_compound, ACT.standard_relation, ACT.standard_type,
                                             ACT.standard_value,ACT.standard_units,
                                             A.tid,
                                             A.description,A.chembl_id as chembl_id_assay, BAO.label, DOCS.chembl_id as chembl_id_doc, TARGET_DICT.organism AS organism
                                             from activities as ACT
                                             left join assays as A
                                             on ACT.assay_id = A.assay_id
                                             left join DOCS
                                             on DOCS.doc_id=A.doc_id
                                             left join bioassay_ontology as BAO
                                             on A.bao_format=BAO.bao_id
                                             LEFT JOIN molecule_dictionary AS MOL_DICT
                                             on ACT.molregno = MOL_DICT.molregno
                                             LEFT JOIN target_dictionary AS TARGET_DICT
                                             ON A.tid = TARGET_DICT.tid
                                     where ACT.standard_value is not null
                                     and A.relationship_type like 'N'
                                     and ACT.standard_units in (", paste(paste0("'", names(standard_unit_map), "'"), collapse = ","), ")
                                     "))


pheno_activities <- activities_1 %>%
  filter(standard_value > 0) %>%
  mutate(
    standard_value = standard_value / standard_unit_map[standard_units],
    # log10_value = log10(standard_value * standard_unit_map[standard_units]),
    standard_units = "nM"
  )

fwrite(
  pheno_activities,
  file.path(dir_release, "chembl_phenotypic_raw.csv.gz")
)

# get clinical &  info ---------------------------------------------------------
###############################################################################T

approval_info<-dbGetQuery(con, paste0("select pref_name, chembl_id as chembl_id_compound, MOLDICT.molregno, max_phase, first_approval, oral, parenteral,
                              topical, black_box_warning,first_in_class, prodrug, indication_class, withdrawn_flag, withdrawn_year, withdrawn_country,
                              withdrawn_reason, max_phase_for_ind, mesh_id, mesh_heading, efo_id, efo_term, ref_type, ref_id, ref_url
                              from molecule_dictionary as MOLDICT
                              left join drug_indication as DRUGIND on MOLDICT.molregno = DRUGIND. molregno
                              left join indication_refs as INDREF on DRUGIND.drugind_id=INDREF.drugind_id
                              where max_phase >0"))
View(approval_info)

fwrite(
  approval_info,
  file.path(dir_release, "chembl_approval_info_raw.csv.gz")
)


# get document data ------------------------------------------------------------
###############################################################################T

REFERENCE_PRIORITY <- c(
  "pubmed_id" = "https://pubmed.ncbi.nlm.nih.gov/",
  "doi" = "http://doi.org/",
  "patent_id" = "https://patents.google.com/?q=",
  "synapse_id" = "https://www.synapse.org/#!Synapse:",
  "chembl_id" = "https://www.ebi.ac.uk/chembl/document_report_card/",
  "hmsl_id" = "https://lincs.hms.harvard.edu/db/datasets/"
)


doc_info <- dbGetQuery(
  con,
  "SELECT chembl_id AS chembl_id_doc, title, pubmed_id, doi, patent_id
  FROM docs"
)

doc_info_long <- doc_info %>%
  as_tibble() %>%
  dplyr::select(chembl_id_doc, pubmed_id, doi, patent_id) %>%
  mutate_at(vars(pubmed_id, doi, patent_id), as.character) %>%
  gather("reference_type", "reference_value", pubmed_id, doi, patent_id, na.rm = TRUE) %>%
  # Add ChEMBL doc ID if nothing else is present
  bind_rows(
    doc_info %>%
      filter(!chembl_id_doc %in% .[["chembl_doc_id"]]) %>%
      transmute(
        chembl_id_doc,
        reference_type = "chembl_id",
        reference_value = chembl_id_doc
      )
  ) %>%
  mutate(
    url = paste0(REFERENCE_PRIORITY[reference_type], reference_value)
  )

fwrite(
  doc_info_long,
  file.path(dir_release, "chembl_ref_info_raw.csv.gz")
)

syn_reference_table <- synPluck(syn_release, "reference_table")

current_references <- synTableQuery(sprintf("SELECT * FROM %s", syn_reference_table)) %>%
  as.data.frame()

synStore(
  Table(
    syn_reference_table,
    anti_join(
      doc_info_long,
      current_references,
      by = c("reference_type", "reference_value")
    )
  )
)

# wrangle best reference source ------------------------------------------------
###############################################################################T

chembl_ref_info_best <- as.data.table(doc_info_long)[
  ,
  reference_type := factor(reference_type, levels = names(REFERENCE_PRIORITY))
][
  order(reference_type),
  .(
    reference_type = reference_type[1],
    reference_value = reference_value[1]
  ),
  keyby = .(chembl_id_doc)
]

fwrite(
  chembl_ref_info_best,
  file.path(dir_release, "chembl_ref_info_best_source.csv.gz")
)

# Store to synapse -------------------------------------------------------------
###############################################################################T

fetch_chembl_activity <- Activity(
  name = "Fetch ChEMBL actitivy data",
  used = unname(inputs),
  executed = "https://github.com/clemenshug/small-molecule-suite-maintenance/blob/master/data_processing/01_chembl_data.R"
)

chembl_raw_syn <- synMkdir(syn_release, "raw_data", "chembl")

c(
  file.path(dir_release, "chembl_biochemical_raw.csv.gz"),
  file.path(dir_release, "chembl_phenotypic_raw.csv.gz"),
  file.path(dir_release, "chembl_ref_info_raw.csv.gz"),
  file.path(dir_release, "chembl_approval_info_raw.csv.gz"),
  file.path(dir_release, "chembl_ref_info_best_source.csv.gz")
) %>%
  synStoreMany(parent = chembl_raw_syn, activity = fetch_chembl_activity)
