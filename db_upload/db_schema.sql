-- SQL dump generated using DBML (dbml-lang.org)
-- Database: PostgreSQL
-- Generated at: 2021-03-09T14:29:54.189Z

CREATE TYPE "approval_tiers" AS ENUM (
  '0',
  '1',
  '2',
  '3',
  '4'
);

CREATE TYPE "name_priorities" AS ENUM (
  'primary',
  'secondary'
);

CREATE TYPE "compound_sources" AS ENUM (
  'chembl',
  'hmsl',
  'emolecules',
  'old_sms'
);

CREATE TYPE "target_types" AS ENUM (
  'SINGLE PROTEIN',
  'CHIMERIC PROTEIN',
  'PROTEIN FAMILY',
  'SELECTIVITY GROUP',
  'PROTEIN COMPLEX GROUP',
  'PROTEIN COMPLEX',
  'PROTEIN NUCLEIC-ACID COMPLEX',
  'PROTEIN-PROTEIN INTERACTION'
);

CREATE TYPE "biochem_value_types" AS ENUM (
  'IC50',
  'Ki',
  'Kd',
  'MIC',
  'EC50',
  'Inhibition',
  'AC50',
  'ED50',
  'IC90',
  'Activity',
  'Potency',
  'ID50',
  'CC50'
);

CREATE TYPE "biochem_value_relations" AS ENUM (
  '=',
  '>',
  '<',
  '>>',
  '<<',
  '~',
  '<=',
  '>='
);

CREATE TYPE "biochem_value_units" AS ENUM (
  'nM'
);

CREATE TYPE "reference_types" AS ENUM (
  'pubmed_id',
  'doi',
  'chembl_id',
  'patent_id',
  'synapse_id',
  'hmsl_id'
);

CREATE TYPE "binding_data_types" AS ENUM (
  'dose_response',
  'single_dose',
  'literature_annotation'
);

CREATE TYPE "selectivity_classes" AS ENUM (
  'most_selective',
  'semi_selective',
  'poly_selective',
  'other_selective',
  'unknown_selective'
);

CREATE TYPE "commercial_tiers" AS ENUM (
  'Tier 1',
  'Tier 2',
  'Tier 3',
  'Tier NA'
);

CREATE TYPE "fingerprint_types" AS ENUM (
  'morgan_normal',
  'morgan_chiral',
  'topological_normal'
);

CREATE TYPE "include_reasons" AS ENUM (
  'selectivity',
  'clinical'
);

CREATE TABLE "lsp_compound_dictionary" (
  "lspci_id" int PRIMARY KEY,
  "hmsl_id" varchar,
  "chembl_id" varchar,
  "emolecules_id" int,
  "pref_name" varchar,
  "inchi" varchar,
  "commercially_available" boolean,
  "max_phase" approval_tiers
);

CREATE TABLE "lsp_structures" (
  "lspci_id" int,
  "source" compound_sources,
  "rank" int,
  "inchi" varchar
);

CREATE TABLE "lsp_compound_names" (
  "lspci_id" int,
  "source" compound_sources,
  "priority" name_priorities,
  "name" varchar
);

CREATE TABLE "lsp_compound_mapping" (
  "lspci_id" int,
  "source" compound_sources,
  "external_id" varchar
);

CREATE TABLE "lsp_target_dictionary" (
  "lspci_target_id" int PRIMARY KEY,
  "gene_id" int,
  "symbol" varchar,
  "pref_name" varchar,
  "tax_id" int,
  "organism" varchar
);

CREATE TABLE "lsp_target_mapping" (
  "lspci_target_id" int,
  "chembl_id" varchar,
  "uniprot_id" varchar,
  "target_type" target_types
);

CREATE TABLE "lsp_references" (
  "reference_id" int PRIMARY KEY,
  "reference_type" reference_types,
  "reference_value" varchar,
  "url" varchar
);

CREATE TABLE "lsp_biochem" (
  "biochem_id" int PRIMARY KEY,
  "biochem_agg_id" int,
  "lspci_id" int,
  "lspci_target_id" int,
  "gene_id" int,
  "symbol" varchar,
  "source" compound_sources,
  "description_assay" varchar,
  "value" float,
  "value_type" biochem_value_types,
  "value_unit" biochem_value_units,
  "value_relation" biochem_value_relations,
  "reference_id" int
);

CREATE TABLE "lsp_biochem_agg" (
  "biochem_agg_id" int PRIMARY KEY,
  "lspci_id" int,
  "lspci_target_id" int,
  "gene_id" int,
  "symbol" varchar,
  "value" float,
  "value_unit" biochem_value_units,
  "tas_id" int
);

CREATE TABLE "lsp_phenotypic" (
  "phenotypic_id" int PRIMARY KEY,
  "lspci_id" int,
  "assay_id" int,
  "value" float,
  "value_type" varchar,
  "value_unit" biochem_value_units,
  "description_assay" varchar,
  "reference_id" int,
  "phenotypic_agg_id" int
);

CREATE TABLE "lsp_phenotypic_agg" (
  "phenotypic_agg_id" int PRIMARY KEY,
  "lspci_id" int,
  "assay_id" int,
  "value" float,
  "value_unit" biochem_value_units,
  "rscore" float,
  "rscore_tr" float
);

CREATE TABLE "lsp_tas" (
  "tas_id" int PRIMARY KEY,
  "lspci_id" int,
  "lspci_target_id" int,
  "gene_id" int,
  "symbol" varchar,
  "tas" int,
  "derived_from" binding_data_types
);

CREATE TABLE "lsp_tas_references" (
  "tas_id" int,
  "reference_id" int
);

CREATE TABLE "lsp_manual_curation" (
  "lspci_id" int,
  "lspci_target_id" int,
  "gene_id" int,
  "symbol" varchar,
  "reference_id" int,
  "tas_id" int
);

CREATE TABLE "lsp_selectivity" (
  "lspci_id" int,
  "lspci_target_id" int,
  "gene_id" int,
  "symbol" varchar,
  "selectivity_class" selectivity_classes,
  "investigation_bias" float,
  "strength" int,
  "wilcox_pval" float,
  "selectivity" float,
  "tool_score" float,
  "ic50_difference" float,
  "ontarget_ic50_q1" float,
  "offtarget_ic50_q1" float,
  "ontarget_n" int,
  "offtarget_n" int
);

CREATE TABLE "lsp_one_dose_scans" (
  "one_dose_scan_id" int PRIMARY KEY,
  "lspci_id" int,
  "lspci_target_id" int,
  "source" compound_sources,
  "gene_id" int,
  "symbol" varchar,
  "percent_control" float,
  "concentration" float,
  "reference_id" varchar,
  "one_dose_scan_agg_id" int
);

CREATE TABLE "lsp_one_dose_scan_agg" (
  "one_dose_scan_agg_id" int PRIMARY KEY,
  "lspci_id" int,
  "lspci_target_id" int,
  "gene_id" int,
  "symbol" varchar,
  "percent_control" float,
  "concentration" float,
  "tas_id" int
);

CREATE TABLE "lsp_clinical_info" (
  "lspci_id" int,
  "max_phase" int,
  "first_approval" int,
  "oral" boolean,
  "parenteral" boolean,
  "topical" boolean,
  "black_box_warning" boolean,
  "first_in_class" boolean,
  "prodrug" boolean,
  "indication_class" varchar,
  "withdrawn_flag" boolean,
  "withdrawn_year" int,
  "withdrawn_country" varchar,
  "withdrawn_reason" varchar
);

CREATE TABLE "lsp_commercial_availability" (
  "lspci_id" int,
  "emolecules_id" int,
  "vendor" varchar,
  "catalog_number" varchar,
  "tier" commercial_tiers
);

CREATE TABLE "lsp_fingerprints" (
  "lspci_id" int,
  "fingerprint_type" fingerprint_types,
  "fingerprint" char(512)
);

CREATE TABLE "lsp_compound_library" (
  "lspci_id" int,
  "lspci_target_id" int,
  "gene_id" int,
  "symbol" varchar,
  "rank" int,
  "reason_included" include_reasons
);

ALTER TABLE "lsp_structures" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_compound_names" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_compound_mapping" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_target_mapping" ADD FOREIGN KEY ("lspci_target_id") REFERENCES "lsp_target_dictionary" ("lspci_target_id");

ALTER TABLE "lsp_biochem" ADD FOREIGN KEY ("biochem_agg_id") REFERENCES "lsp_biochem_agg" ("biochem_agg_id");

ALTER TABLE "lsp_biochem" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_biochem" ADD FOREIGN KEY ("lspci_target_id") REFERENCES "lsp_target_dictionary" ("lspci_target_id");

ALTER TABLE "lsp_biochem" ADD FOREIGN KEY ("reference_id") REFERENCES "lsp_references" ("reference_id");

ALTER TABLE "lsp_biochem_agg" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_biochem_agg" ADD FOREIGN KEY ("lspci_target_id") REFERENCES "lsp_target_dictionary" ("lspci_target_id");

ALTER TABLE "lsp_biochem_agg" ADD FOREIGN KEY ("tas_id") REFERENCES "lsp_tas" ("tas_id");

ALTER TABLE "lsp_phenotypic" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_phenotypic" ADD FOREIGN KEY ("reference_id") REFERENCES "lsp_references" ("reference_id");

ALTER TABLE "lsp_phenotypic" ADD FOREIGN KEY ("phenotypic_agg_id") REFERENCES "lsp_phenotypic_agg" ("phenotypic_agg_id");

ALTER TABLE "lsp_phenotypic_agg" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_tas" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_tas" ADD FOREIGN KEY ("lspci_target_id") REFERENCES "lsp_target_dictionary" ("lspci_target_id");

ALTER TABLE "lsp_tas_references" ADD FOREIGN KEY ("tas_id") REFERENCES "lsp_tas" ("tas_id");

ALTER TABLE "lsp_tas_references" ADD FOREIGN KEY ("reference_id") REFERENCES "lsp_references" ("reference_id");

ALTER TABLE "lsp_manual_curation" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_manual_curation" ADD FOREIGN KEY ("lspci_target_id") REFERENCES "lsp_target_dictionary" ("lspci_target_id");

ALTER TABLE "lsp_manual_curation" ADD FOREIGN KEY ("reference_id") REFERENCES "lsp_references" ("reference_id");

ALTER TABLE "lsp_manual_curation" ADD FOREIGN KEY ("tas_id") REFERENCES "lsp_tas" ("tas_id");

ALTER TABLE "lsp_selectivity" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_selectivity" ADD FOREIGN KEY ("lspci_target_id") REFERENCES "lsp_target_dictionary" ("lspci_target_id");

ALTER TABLE "lsp_one_dose_scans" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_one_dose_scans" ADD FOREIGN KEY ("lspci_target_id") REFERENCES "lsp_target_dictionary" ("lspci_target_id");

ALTER TABLE "lsp_one_dose_scans" ADD FOREIGN KEY ("one_dose_scan_agg_id") REFERENCES "lsp_one_dose_scan_agg" ("one_dose_scan_agg_id");

ALTER TABLE "lsp_one_dose_scan_agg" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_one_dose_scan_agg" ADD FOREIGN KEY ("lspci_target_id") REFERENCES "lsp_target_dictionary" ("lspci_target_id");

ALTER TABLE "lsp_one_dose_scan_agg" ADD FOREIGN KEY ("tas_id") REFERENCES "lsp_tas" ("tas_id");

ALTER TABLE "lsp_clinical_info" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_commercial_availability" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_fingerprints" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_compound_library" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_compound_library" ADD FOREIGN KEY ("lspci_target_id") REFERENCES "lsp_target_dictionary" ("lspci_target_id");

CREATE INDEX ON "lsp_compound_dictionary" ("hmsl_id");

CREATE INDEX ON "lsp_compound_dictionary" ("chembl_id");

CREATE INDEX ON "lsp_compound_dictionary" ("emolecules_id");

CREATE INDEX ON "lsp_compound_dictionary" ("commercially_available");

CREATE INDEX ON "lsp_structures" ("lspci_id");

CREATE INDEX ON "lsp_compound_names" ("lspci_id");

CREATE INDEX ON "lsp_compound_mapping" ("lspci_id");

CREATE INDEX ON "lsp_compound_mapping" ("external_id");

CREATE INDEX ON "lsp_target_dictionary" ("symbol");

CREATE INDEX ON "lsp_target_dictionary" ("tax_id");

CREATE INDEX ON "lsp_target_dictionary" ("organism");

CREATE INDEX ON "lsp_target_mapping" ("lspci_target_id");

CREATE INDEX ON "lsp_references" ("reference_value");

CREATE INDEX ON "lsp_biochem" ("lspci_id");

CREATE INDEX ON "lsp_biochem" ("lspci_target_id");

CREATE INDEX ON "lsp_biochem" ("lspci_id", "lspci_target_id");

CREATE INDEX ON "lsp_biochem" ("biochem_agg_id");

CREATE INDEX ON "lsp_biochem" ("reference_id");

CREATE INDEX ON "lsp_biochem_agg" ("lspci_id");

CREATE INDEX ON "lsp_biochem_agg" ("lspci_target_id");

CREATE INDEX ON "lsp_biochem_agg" ("lspci_id", "lspci_target_id");

CREATE INDEX ON "lsp_biochem_agg" ("biochem_agg_id");

CREATE INDEX ON "lsp_biochem_agg" ("tas_id");

CREATE INDEX ON "lsp_phenotypic" ("lspci_id");

CREATE INDEX ON "lsp_phenotypic" ("assay_id");

CREATE INDEX ON "lsp_phenotypic" ("lspci_id", "assay_id");

CREATE INDEX ON "lsp_phenotypic" ("phenotypic_agg_id");

CREATE INDEX ON "lsp_phenotypic_agg" ("lspci_id");

CREATE INDEX ON "lsp_phenotypic_agg" ("assay_id");

CREATE INDEX ON "lsp_phenotypic_agg" ("lspci_id", "assay_id");

CREATE INDEX ON "lsp_tas" ("lspci_id");

CREATE INDEX ON "lsp_tas" ("lspci_target_id");

CREATE INDEX ON "lsp_tas" ("lspci_id", "lspci_target_id");

CREATE INDEX ON "lsp_tas" ("tas");

CREATE INDEX ON "lsp_tas_references" ("tas_id");

CREATE INDEX ON "lsp_tas_references" ("reference_id");

CREATE INDEX ON "lsp_manual_curation" ("lspci_id");

CREATE INDEX ON "lsp_manual_curation" ("lspci_target_id");

CREATE INDEX ON "lsp_manual_curation" ("lspci_id", "lspci_target_id");

CREATE INDEX ON "lsp_manual_curation" ("tas_id");

CREATE INDEX ON "lsp_selectivity" ("lspci_id");

CREATE INDEX ON "lsp_selectivity" ("lspci_target_id");

CREATE INDEX ON "lsp_selectivity" ("lspci_id", "lspci_target_id");

CREATE INDEX ON "lsp_one_dose_scans" ("lspci_id");

CREATE INDEX ON "lsp_one_dose_scans" ("lspci_target_id");

CREATE INDEX ON "lsp_one_dose_scans" ("lspci_id", "lspci_target_id");

CREATE INDEX ON "lsp_one_dose_scans" ("one_dose_scan_id");

CREATE INDEX ON "lsp_one_dose_scan_agg" ("lspci_id");

CREATE INDEX ON "lsp_one_dose_scan_agg" ("lspci_target_id");

CREATE INDEX ON "lsp_one_dose_scan_agg" ("tas_id");

CREATE INDEX ON "lsp_one_dose_scan_agg" ("lspci_id", "lspci_target_id");

CREATE INDEX ON "lsp_one_dose_scan_agg" ("one_dose_scan_agg_id");

CREATE INDEX ON "lsp_clinical_info" ("lspci_id");

CREATE INDEX ON "lsp_clinical_info" ("max_phase");

CREATE INDEX ON "lsp_commercial_availability" ("lspci_id");

CREATE INDEX ON "lsp_fingerprints" ("lspci_id");

CREATE INDEX ON "lsp_compound_library" ("lspci_id");

CREATE INDEX ON "lsp_compound_library" ("lspci_target_id");

CREATE INDEX ON "lsp_compound_library" ("lspci_id", "lspci_target_id");

COMMENT ON TABLE "lsp_compound_dictionary" IS 'Primary table listing all compounds in the database.     During compound processing distinct salts of the same compound     are aggregated into a single compound entry in this table.     The constituent compound IDs for each compound in this table     are available in the lsp_compound_mapping table.';

COMMENT ON COLUMN "lsp_compound_dictionary"."lspci_id" IS 'Internal compound ID';

COMMENT ON COLUMN "lsp_compound_dictionary"."hmsl_id" IS 'Primary HMS LINCS compound ID, if available';

COMMENT ON COLUMN "lsp_compound_dictionary"."chembl_id" IS 'Primary ChEMBL compound ID, if available';

COMMENT ON COLUMN "lsp_compound_dictionary"."emolecules_id" IS 'Primary EMolecules compound ID, if available';

COMMENT ON COLUMN "lsp_compound_dictionary"."pref_name" IS 'Preferred name of the compound';

COMMENT ON COLUMN "lsp_compound_dictionary"."inchi" IS 'InChI chemical indentifier for the compound, standardized     using https://github.com/chembl/ChEMBL_Structure_Pipeline';

COMMENT ON COLUMN "lsp_compound_dictionary"."commercially_available" IS 'Indicates that compound is commercially available at a vendor in the eMolecules database';

COMMENT ON COLUMN "lsp_compound_dictionary"."max_phase" IS 'Approval and clinical trial status of the compound.';

COMMENT ON TABLE "lsp_structures" IS 'Additional secondary InChIs for compounds.';

COMMENT ON COLUMN "lsp_structures"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_structures"."source" IS 'Source for the compound name';

COMMENT ON COLUMN "lsp_structures"."rank" IS 'Rank of InCHI according to source reliability and annotation quality';

COMMENT ON COLUMN "lsp_structures"."inchi" IS 'InChI chemical indentifier for the compound, standardized     using https://github.com/chembl/ChEMBL_Structure_Pipeline';

COMMENT ON TABLE "lsp_compound_names" IS 'Table of all annotated names for compounds. The sources     for compound names generally distinguish between primary and     alternative (secondary) names.';

COMMENT ON COLUMN "lsp_compound_names"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_compound_names"."source" IS 'Source for the compound name';

COMMENT ON COLUMN "lsp_compound_names"."priority" IS 'Annotated as primary or secondary name at source';

COMMENT ON COLUMN "lsp_compound_names"."name" IS 'Compound name';

COMMENT ON TABLE "lsp_compound_mapping" IS 'Table of mappings between compound IDs from different sources to     the internal lspci_ids.';

COMMENT ON COLUMN "lsp_compound_mapping"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_compound_mapping"."source" IS 'Source for the compound ID.';

COMMENT ON COLUMN "lsp_compound_mapping"."external_id" IS 'ID of the compound at the external source.';

COMMENT ON TABLE "lsp_target_dictionary" IS 'Table of drug targets. The original drug targets are mostly     annotated as ChEMBL or UniProt IDs. For convenience we converted     these IDs to Entrez gene IDs. The original mapping between     ChEMBL and UniProt target IDs are in the table `lsp_target_mapping`';

COMMENT ON COLUMN "lsp_target_dictionary"."lspci_target_id" IS 'Internal target ID';

COMMENT ON COLUMN "lsp_target_dictionary"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_target_dictionary"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_target_dictionary"."pref_name" IS 'Gene name';

COMMENT ON COLUMN "lsp_target_dictionary"."tax_id" IS 'Entrez taxonomy ID';

COMMENT ON COLUMN "lsp_target_dictionary"."organism" IS 'Organism for which gene is annotated';

COMMENT ON TABLE "lsp_target_mapping" IS 'Mapping between the original ChEMBL target IDs,     their corresponding UniProt IDs and Entrez gene IDs.     A single UniProt or ChEMBL ID can refer to protein complexes,     therefore multiple gene IDs often map to the same UniProt or     ChEMBL ID.';

COMMENT ON COLUMN "lsp_target_mapping"."lspci_target_id" IS 'Foreign key to lsp_target_dictionary table';

COMMENT ON COLUMN "lsp_target_mapping"."chembl_id" IS 'ChEMBL target ID';

COMMENT ON COLUMN "lsp_target_mapping"."uniprot_id" IS 'UniProt ID';

COMMENT ON COLUMN "lsp_target_mapping"."target_type" IS 'The type of the original target before translation to Entrez IDs';

COMMENT ON TABLE "lsp_references" IS 'External references for the data in the database.';

COMMENT ON COLUMN "lsp_references"."reference_id" IS 'Primary key';

COMMENT ON COLUMN "lsp_references"."reference_type" IS 'The source of the measurement.';

COMMENT ON COLUMN "lsp_references"."reference_value" IS 'The reference accesion number at the source.';

COMMENT ON COLUMN "lsp_references"."url" IS 'URL to access the reference.';

COMMENT ON TABLE "lsp_biochem" IS 'Table of biochemical affinity measurements.';

COMMENT ON COLUMN "lsp_biochem"."biochem_id" IS 'Primary key';

COMMENT ON COLUMN "lsp_biochem"."biochem_agg_id" IS 'Foreign key to lsp_biochem_agg table. All measurments with the same ID are aggregated.';

COMMENT ON COLUMN "lsp_biochem"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_biochem"."lspci_target_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_biochem"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_biochem"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_biochem"."source" IS 'Source of the measurement';

COMMENT ON COLUMN "lsp_biochem"."description_assay" IS 'Description of the assay that the measurement is derived from';

COMMENT ON COLUMN "lsp_biochem"."value" IS 'Measurement value';

COMMENT ON COLUMN "lsp_biochem"."value_type" IS 'The type of measurement performed.';

COMMENT ON COLUMN "lsp_biochem"."value_unit" IS 'The unit of the measurement.';

COMMENT ON COLUMN "lsp_biochem"."value_relation" IS 'Some assays can""t determine the measured value precisely. This column gives relationship between the actual and the measured value.';

COMMENT ON COLUMN "lsp_biochem"."reference_id" IS 'Foreing key to the ls_reference table for this measurement.';

COMMENT ON TABLE "lsp_biochem_agg" IS 'Table of aggregated biochemical affinity measurements. All
available data for a single compound target pair were aggregated
by taking the first quartile.';

COMMENT ON COLUMN "lsp_biochem_agg"."biochem_agg_id" IS 'All measurements with this ID in lsp_biochem were aggregated.';

COMMENT ON COLUMN "lsp_biochem_agg"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_biochem_agg"."lspci_target_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_biochem_agg"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_biochem_agg"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_biochem_agg"."value" IS 'Aggregated measurement value';

COMMENT ON COLUMN "lsp_biochem_agg"."value_unit" IS 'The unit of the measurement.';

COMMENT ON COLUMN "lsp_biochem_agg"."tas_id" IS 'Foreign key to lsp_tas table. Indicates that this value was used to calculate the referenced TAS.';

COMMENT ON TABLE "lsp_phenotypic" IS 'Table of phenotypic assays performed on the compounds.';

COMMENT ON COLUMN "lsp_phenotypic"."phenotypic_id" IS 'Primary key for phenotypic data.';

COMMENT ON COLUMN "lsp_phenotypic"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_phenotypic"."assay_id" IS 'Unique ID of the assays in the database. Correspond to ChEMBL assay IDs.';

COMMENT ON COLUMN "lsp_phenotypic"."value" IS 'Measurement value';

COMMENT ON COLUMN "lsp_phenotypic"."value_type" IS 'The type of measurement performed.';

COMMENT ON COLUMN "lsp_phenotypic"."value_unit" IS 'The unit of the measurement.';

COMMENT ON COLUMN "lsp_phenotypic"."description_assay" IS 'Description of the assay.';

COMMENT ON COLUMN "lsp_phenotypic"."reference_id" IS 'Foreing key to the ls_reference table for this measurement.';

COMMENT ON COLUMN "lsp_phenotypic"."phenotypic_agg_id" IS 'Foreign key to lsp_phenotypic_agg table. All measurements with the same ID are aggregated.';

COMMENT ON TABLE "lsp_phenotypic_agg" IS 'Table of aggregated phenotypic assays performed on the compounds.
All available data for a single assay and compound target pair were
aggregated by taking the first quartile.';

COMMENT ON COLUMN "lsp_phenotypic_agg"."phenotypic_agg_id" IS 'All measurements with this ID in lsp_phenotypic were aggregated.';

COMMENT ON COLUMN "lsp_phenotypic_agg"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_phenotypic_agg"."assay_id" IS 'Unique ID of the assays in the database. Correspond to ChEMBL assay IDs.';

COMMENT ON COLUMN "lsp_phenotypic_agg"."value" IS 'Aggregated measurement value';

COMMENT ON COLUMN "lsp_phenotypic_agg"."value_unit" IS 'The unit of the measurement.';

COMMENT ON COLUMN "lsp_phenotypic_agg"."rscore" IS 'Standardized measurement value.';

COMMENT ON COLUMN "lsp_phenotypic_agg"."rscore_tr" IS 'Rscore normalized using a double logistic function.';

COMMENT ON TABLE "lsp_tas" IS 'Table of Target Affinity Spectrum (TAS) values for the affinity between     compound and target. TAS enables aggregation of affinity measurements from     heterogeneous sources and assays into a single value. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "lsp_tas"."tas_id" IS 'ID of this TAS value. The lsp_biochem_agg, lsp_one_dose_scan_agg, and lsp_manual_curation tables reference this ID in order to designate which value was used to calculate the TAS.';

COMMENT ON COLUMN "lsp_tas"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_tas"."lspci_target_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_tas"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_tas"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_tas"."tas" IS 'Target Affinity Spectrum (TAS). A number between 1 and 10 for the affinity between compound and target. 1 being most strongly binding.';

COMMENT ON COLUMN "lsp_tas"."derived_from" IS 'Data type that the TAS value was derived from.';

COMMENT ON TABLE "lsp_tas_references" IS 'Table that makes it easier to link TAS values to the references that were used to compute the TAS values';

COMMENT ON COLUMN "lsp_tas_references"."tas_id" IS 'Foreing key to lsp_tas table. References in this table were used to compute the linked TAS values.';

COMMENT ON COLUMN "lsp_tas_references"."reference_id" IS 'Foreign key to lsp_references table. Points to all references associated with a given TAS value.';

COMMENT ON TABLE "lsp_manual_curation" IS 'Table of manual compund target binding assertions.';

COMMENT ON COLUMN "lsp_manual_curation"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_manual_curation"."lspci_target_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_manual_curation"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_manual_curation"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_manual_curation"."reference_id" IS 'Foreing key to the ls_reference table for this measurement.';

COMMENT ON COLUMN "lsp_manual_curation"."tas_id" IS 'Foreign key to lsp_tas table. Indicates that this value was used to calculate the referenced TAS.';

COMMENT ON TABLE "lsp_selectivity" IS 'Table of selectivity assertions of compounds to their targets. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "lsp_selectivity"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_selectivity"."lspci_target_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_selectivity"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_selectivity"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_selectivity"."selectivity_class" IS 'Assertion for the selectivity of the compound to the given target.';

COMMENT ON COLUMN "lsp_selectivity"."ic50_difference" IS 'Difference between ontarget and offtarget IC50 Q1 measurements';

COMMENT ON COLUMN "lsp_selectivity"."ontarget_ic50_q1" IS 'First quartile of ontarget IC50 measurements';

COMMENT ON COLUMN "lsp_selectivity"."offtarget_ic50_q1" IS 'First quartile of offtarget IC50 measurements';

COMMENT ON COLUMN "lsp_selectivity"."ontarget_n" IS 'Number of ontarget IC50 measurements';

COMMENT ON COLUMN "lsp_selectivity"."offtarget_n" IS 'Number of offtarget IC50 measurements';

COMMENT ON TABLE "lsp_one_dose_scans" IS 'Table of single dose compound activity measurements as     opposed to full dose-response affinity measurements.';

COMMENT ON COLUMN "lsp_one_dose_scans"."one_dose_scan_id" IS 'Primary key for single dose measurements';

COMMENT ON COLUMN "lsp_one_dose_scans"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_one_dose_scans"."lspci_target_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_one_dose_scans"."source" IS 'Source of the measurement';

COMMENT ON COLUMN "lsp_one_dose_scans"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_one_dose_scans"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_one_dose_scans"."percent_control" IS 'Remaining activity of target at the given compound concentration.';

COMMENT ON COLUMN "lsp_one_dose_scans"."concentration" IS 'Concentration of the compound in the assay.';

COMMENT ON COLUMN "lsp_one_dose_scans"."reference_id" IS 'The reference for the measurement.';

COMMENT ON COLUMN "lsp_one_dose_scans"."one_dose_scan_agg_id" IS 'Foreign key to lsp_phenotypic_agg table. All measurements with the same ID are aggregated.';

COMMENT ON TABLE "lsp_one_dose_scan_agg" IS 'Table of single dose compound activity measurements as     opposed to full dose-response affinity measurements.
All available data for a single concentration and compound target pair were
aggregated by taking the first quartile.';

COMMENT ON COLUMN "lsp_one_dose_scan_agg"."one_dose_scan_agg_id" IS 'All measurements with this ID in lsp_one_dose_scans were aggregated.';

COMMENT ON COLUMN "lsp_one_dose_scan_agg"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_one_dose_scan_agg"."lspci_target_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_one_dose_scan_agg"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_one_dose_scan_agg"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_one_dose_scan_agg"."percent_control" IS 'Aggregated remaining activity of target at the given compound concentration.';

COMMENT ON COLUMN "lsp_one_dose_scan_agg"."concentration" IS 'Concentration of the compound in the assay.';

COMMENT ON COLUMN "lsp_one_dose_scan_agg"."tas_id" IS 'Foreign key to lsp_tas table. Indicates that this value was used to calculate the referenced TAS.';

COMMENT ON TABLE "lsp_clinical_info" IS 'Table of the clinical approval status of compounds.     Sourced from ChEMBL';

COMMENT ON COLUMN "lsp_clinical_info"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON TABLE "lsp_commercial_availability" IS 'Table of the commercial availability of compounds.     Sourced from eMolecules (https://www.emolecules.com/).';

COMMENT ON COLUMN "lsp_commercial_availability"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_commercial_availability"."emolecules_id" IS 'ID of compound at eMolecules';

COMMENT ON COLUMN "lsp_commercial_availability"."vendor" IS 'The name of the vendor';

COMMENT ON COLUMN "lsp_commercial_availability"."catalog_number" IS 'Catalog number of the compound.';

COMMENT ON COLUMN "lsp_commercial_availability"."tier" IS 'Shipment timeframe.';

COMMENT ON TABLE "lsp_fingerprints" IS 'Table of specificity assertions of compounds to their targets. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "lsp_fingerprints"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_fingerprints"."fingerprint_type" IS 'Three different fingerprint types are available: Morgan fingerprints (either ignoring or respecting stereochemistry) or RDKit topological fingerprints (ignoring stereochemistry).';

COMMENT ON COLUMN "lsp_fingerprints"."fingerprint" IS '256 byte hex-encoded fingerprint';

COMMENT ON TABLE "lsp_compound_library" IS 'Library of optimal compounds for each target. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "lsp_compound_library"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_compound_library"."lspci_target_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_compound_library"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_compound_library"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_compound_library"."rank" IS 'Ranking of compounds for the same target. First compound is best.';

COMMENT ON COLUMN "lsp_compound_library"."reason_included" IS 'Compounds can be included either because of their selectivity or because they are approved/in clinical trials.';
