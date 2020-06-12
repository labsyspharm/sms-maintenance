CREATE TYPE "name_priorities" AS ENUM (
  'primary',
  'secondary'
);

CREATE TYPE "compound_sources" AS ENUM (
  'chembl',
  'hmsl',
  'vendor'
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

CREATE TYPE "selectivity_classes" AS ENUM (
  'most_selective',
  'semi_selective',
  'poly_selective',
  'other_selective',
  'unknown_selective'
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
  "pref_name" varchar,
  "inchi" varchar,
  "smiles" varchar,
  "commercially_available" boolean
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
  "id" varchar
);

CREATE TABLE "lsp_target_dictionary" (
  "gene_id" int PRIMARY KEY,
  "symbol" varchar,
  "pref_name" varchar,
  "description" varchar,
  "tax_id" int,
  "organism" varchar
);

CREATE TABLE "lsp_target_mapping" (
  "gene_id" int,
  "chembl_id" varchar,
  "uniprot_id" varchar
);

CREATE TABLE "lsp_biochem" (
  "lspci_id" int,
  "gene_id" int,
  "description_assay" varchar,
  "value" float,
  "value_type" biochem_value_types,
  "value_unit" biochem_value_units,
  "value_relation" biochem_value_relations,
  "reference_id" varchar,
  "reference_type" reference_types,
  "url" varchar
);

CREATE TABLE "lsp_phenotypic_chembl" (
  "lspci_id" int,
  "assay_id" int,
  "value" float,
  "value_type" biochem_value_types,
  "value_unit" biochem_value_units,
  "value_relation" biochem_value_relations,
  "reference_id" varchar,
  "reference_type" reference_types,
  "url" varchar
);

CREATE TABLE "lsp_tas" (
  "lspci_id" int,
  "gene_id" int,
  "tas" int
);

CREATE TABLE "lsp_specificity" (
  "lspci_id" int,
  "gene_id" int,
  "selectivity_class" selectivity_classes,
  "investigation_bias" float,
  "strength" int,
  "wilcox_pval" float,
  "selectivity" float,
  "tool_score" float,
  "IC50_difference" float,
  "ontarget_IC50_Q1" float,
  "offtarget_IC50_Q1" float,
  "ontarget_N" int,
  "offtarget_N" int
);

CREATE TABLE "lsp_one_dose_scans" (
  "lspci_id" int,
  "gene_id" int,
  "percent_control" float,
  "description" varchar,
  "cmpd_conc_nM" float,
  "reference_id" varchar,
  "reference_type" reference_types,
  "url" varchar
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
  "withdrawn_reason" varchar,
  "max_phase_for_indication" varchar,
  "mesh_id" varchar,
  "mesh_heading" varchar,
  "efo_id" varchar,
  "efo_term" varchar,
  "reference_type" varchar,
  "reference_id" varchar,
  "url" varchar
);

CREATE TABLE "lsp_commercial_availability" (
  "lspci_id" int,
  "vendor" varchar,
  "id" varchar,
  "name" varchar
);

CREATE TABLE "lsp_fingerprints" (
  "lspci_id" int,
  "fingerprint_type" fingerprint_types,
  "fingerprint" char(512)
);

CREATE TABLE "lsp_compound_library" (
  "lspci_id" int,
  "gene_id" int,
  "rank" int,
  "reason_included" include_reasons
);

ALTER TABLE "lsp_compound_names" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_compound_mapping" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_target_mapping" ADD FOREIGN KEY ("gene_id") REFERENCES "lsp_target_dictionary" ("gene_id");

ALTER TABLE "lsp_biochem" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_biochem" ADD FOREIGN KEY ("gene_id") REFERENCES "lsp_target_dictionary" ("gene_id");

ALTER TABLE "lsp_phenotypic_chembl" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_tas" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_tas" ADD FOREIGN KEY ("gene_id") REFERENCES "lsp_target_dictionary" ("gene_id");

ALTER TABLE "lsp_specificity" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_specificity" ADD FOREIGN KEY ("gene_id") REFERENCES "lsp_target_dictionary" ("gene_id");

ALTER TABLE "lsp_one_dose_scans" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_one_dose_scans" ADD FOREIGN KEY ("gene_id") REFERENCES "lsp_target_dictionary" ("gene_id");

ALTER TABLE "lsp_clinical_info" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_commercial_availability" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_fingerprints" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_compound_library" ADD FOREIGN KEY ("lspci_id") REFERENCES "lsp_compound_dictionary" ("lspci_id");

ALTER TABLE "lsp_compound_library" ADD FOREIGN KEY ("gene_id") REFERENCES "lsp_target_dictionary" ("gene_id");

CREATE UNIQUE INDEX ON "lsp_compound_dictionary" ("lspci_id");

CREATE INDEX ON "lsp_compound_names" ("lspci_id");

CREATE INDEX ON "lsp_compound_mapping" ("lspci_id");

CREATE UNIQUE INDEX ON "lsp_target_dictionary" ("gene_id");

CREATE INDEX ON "lsp_target_dictionary" ("symbol");

CREATE INDEX ON "lsp_target_mapping" ("gene_id");

CREATE INDEX ON "lsp_biochem" ("lspci_id");

CREATE INDEX ON "lsp_biochem" ("gene_id");

CREATE INDEX ON "lsp_biochem" ("lspci_id", "gene_id");

CREATE INDEX ON "lsp_phenotypic_chembl" ("lspci_id");

CREATE INDEX ON "lsp_phenotypic_chembl" ("assay_id");

CREATE INDEX ON "lsp_tas" ("lspci_id");

CREATE INDEX ON "lsp_tas" ("gene_id");

CREATE INDEX ON "lsp_tas" ("lspci_id", "gene_id");

CREATE INDEX ON "lsp_specificity" ("lspci_id");

CREATE INDEX ON "lsp_specificity" ("gene_id");

CREATE INDEX ON "lsp_specificity" ("lspci_id", "gene_id");

CREATE INDEX ON "lsp_one_dose_scans" ("lspci_id");

CREATE INDEX ON "lsp_one_dose_scans" ("gene_id");

CREATE INDEX ON "lsp_one_dose_scans" ("lspci_id", "gene_id");

CREATE INDEX ON "lsp_clinical_info" ("lspci_id");

CREATE INDEX ON "lsp_clinical_info" ("max_phase");

CREATE INDEX ON "lsp_commercial_availability" ("lspci_id");

CREATE INDEX ON "lsp_fingerprints" ("lspci_id");

CREATE INDEX ON "lsp_compound_library" ("lspci_id");

CREATE INDEX ON "lsp_compound_library" ("gene_id");

CREATE INDEX ON "lsp_compound_library" ("lspci_id", "gene_id");

COMMENT ON TABLE "lsp_compound_dictionary" IS 'Primary table listing all compounds in the database.     During compound processing distinct salts of the same compound     are aggregated into a single compound entry in this table.     The constituent compound IDs for each compound in this table     are available in the lsp_compound_mapping table.';

COMMENT ON COLUMN "lsp_compound_dictionary"."lspci_id" IS 'Internal compound ID';

COMMENT ON COLUMN "lsp_compound_dictionary"."hmsl_id" IS 'Primary HMS LINCS compound ID, if available';

COMMENT ON COLUMN "lsp_compound_dictionary"."chembl_id" IS 'Primary ChEMBL compound ID, if available';

COMMENT ON COLUMN "lsp_compound_dictionary"."pref_name" IS 'Preferred name of the compound';

COMMENT ON COLUMN "lsp_compound_dictionary"."inchi" IS 'InChI chemical indentifier for the compound, standardized     using https://molvs.readthedocs.io/en/latest/guide/standardize.html';

COMMENT ON COLUMN "lsp_compound_dictionary"."smiles" IS 'SMILES chemical indentifier for the compound, standardized     using https://molvs.readthedocs.io/en/latest/guide/standardize.html';

COMMENT ON COLUMN "lsp_compound_dictionary"."commercially_available" IS 'Commercially available at a vendor in the ZINC database';

COMMENT ON TABLE "lsp_compound_names" IS 'Table of all annotated names for compounds. The sources     for compound names generally distinguish between primary and     alternative (secondary) names.';

COMMENT ON COLUMN "lsp_compound_names"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_compound_names"."source" IS 'Source for the compound name';

COMMENT ON COLUMN "lsp_compound_names"."priority" IS 'Annotated as primary or secondary name at source';

COMMENT ON COLUMN "lsp_compound_names"."name" IS 'Compound name';

COMMENT ON TABLE "lsp_compound_mapping" IS 'Table of mappings between compound IDs from different sources to     the internal lspci_ids.';

COMMENT ON COLUMN "lsp_compound_mapping"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_compound_mapping"."source" IS 'Source for the compound ID';

COMMENT ON COLUMN "lsp_compound_mapping"."id" IS 'Compound ID associated with the given lspci_id';

COMMENT ON TABLE "lsp_target_dictionary" IS 'Table of drug targets. The original drug targets are mostly     annotated as ChEMBL or UniProt IDs. For convenience we converted     these IDs to Entrez gene IDs.';

COMMENT ON COLUMN "lsp_target_dictionary"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "lsp_target_dictionary"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "lsp_target_dictionary"."pref_name" IS 'Preferred gene name';

COMMENT ON COLUMN "lsp_target_dictionary"."description" IS 'Description of gene function';

COMMENT ON COLUMN "lsp_target_dictionary"."tax_id" IS 'Entrez taxonomy ID';

COMMENT ON COLUMN "lsp_target_dictionary"."organism" IS 'Organism for which gene is annotated';

COMMENT ON TABLE "lsp_target_mapping" IS 'Mapping between the original ChEMBL target IDs,     their corresponding UniProt IDs and Entrez gene IDs.     A signle UniProt or ChEMBL ID can refer to protein complexes,     therefore multiple gene IDs often map to the same UniProt or     ChEMBL ID.';

COMMENT ON COLUMN "lsp_target_mapping"."gene_id" IS 'Foreign key to Entrez gene ID';

COMMENT ON COLUMN "lsp_target_mapping"."chembl_id" IS 'ChEMBL target ID';

COMMENT ON COLUMN "lsp_target_mapping"."uniprot_id" IS 'UniProt ID';

COMMENT ON TABLE "lsp_biochem" IS 'Table of biochemical affinity measurements.';

COMMENT ON COLUMN "lsp_biochem"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_biochem"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_biochem"."description_assay" IS 'Description of the assay that the measurement is derived from';

COMMENT ON COLUMN "lsp_biochem"."value" IS 'Measurement value';

COMMENT ON COLUMN "lsp_biochem"."value_type" IS 'The type of measurement performed.';

COMMENT ON COLUMN "lsp_biochem"."value_unit" IS 'The unit of the measurement.';

COMMENT ON COLUMN "lsp_biochem"."value_relation" IS 'Some assays can't determine the measured value precisely. This column gives relationship between the actual and the measured value.';

COMMENT ON COLUMN "lsp_biochem"."reference_id" IS 'The reference for the measurement.';

COMMENT ON COLUMN "lsp_biochem"."reference_type" IS 'The source of the measurement.';

COMMENT ON COLUMN "lsp_biochem"."url" IS 'URL to access the reference.';

COMMENT ON TABLE "lsp_phenotypic_chembl" IS 'Table of phenotypic assays performed on the compounds.';

COMMENT ON COLUMN "lsp_phenotypic_chembl"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_phenotypic_chembl"."assay_id" IS 'Unique ID of the assays in the database. Correspond to ChEMBL assay IDs.';

COMMENT ON COLUMN "lsp_phenotypic_chembl"."value" IS 'Measurement value';

COMMENT ON COLUMN "lsp_phenotypic_chembl"."value_type" IS 'The type of measurement performed.';

COMMENT ON COLUMN "lsp_phenotypic_chembl"."value_unit" IS 'The unit of the measurement.';

COMMENT ON COLUMN "lsp_phenotypic_chembl"."value_relation" IS 'Some assays can't determine the measured value precisely. This column gives relationship between the actual and the measured value.';

COMMENT ON COLUMN "lsp_phenotypic_chembl"."reference_id" IS 'The reference for the measurement.';

COMMENT ON COLUMN "lsp_phenotypic_chembl"."reference_type" IS 'The source of the measurement.';

COMMENT ON COLUMN "lsp_phenotypic_chembl"."url" IS 'URL to access the reference.';

COMMENT ON TABLE "lsp_tas" IS 'Table of Target Affinity Spectrum (TAS) values for the affinity between     compound and target. TAS enables aggregation of affinity measurements from     heterogeneous sources and assays into a single value. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "lsp_tas"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_tas"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_tas"."tas" IS 'Target Affinity Spectrum (TAS). A number between 1 and 10 for the affinity between compound and target. 1 being most strongly binding.';

COMMENT ON TABLE "lsp_specificity" IS 'Table of specificity assertions of compounds to their targets. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "lsp_specificity"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_specificity"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_specificity"."selectivity_class" IS 'Assertion for the selectivity of the compound to the given target.';

COMMENT ON COLUMN "lsp_specificity"."IC50_difference" IS 'Difference between ontarget and offtarget IC50 Q1 measurements';

COMMENT ON COLUMN "lsp_specificity"."ontarget_IC50_Q1" IS 'First quartile of ontarget IC50 measurements';

COMMENT ON COLUMN "lsp_specificity"."offtarget_IC50_Q1" IS 'First quartile of offtarget IC50 measurements';

COMMENT ON COLUMN "lsp_specificity"."ontarget_N" IS 'Number of ontarget IC50 measurements';

COMMENT ON COLUMN "lsp_specificity"."offtarget_N" IS 'Number of offtarget IC50 measurements';

COMMENT ON TABLE "lsp_one_dose_scans" IS 'Table of single dose compound activity measurements as     opposed to full dose-response affinity measurements.';

COMMENT ON COLUMN "lsp_one_dose_scans"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_one_dose_scans"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_one_dose_scans"."percent_control" IS 'Remaining activity of target at the given compound concentration.';

COMMENT ON COLUMN "lsp_one_dose_scans"."description" IS 'Description of target.';

COMMENT ON COLUMN "lsp_one_dose_scans"."cmpd_conc_nM" IS 'Concentration of the compound.';

COMMENT ON COLUMN "lsp_one_dose_scans"."reference_id" IS 'The reference for the measurement.';

COMMENT ON COLUMN "lsp_one_dose_scans"."reference_type" IS 'The source of the measurement.';

COMMENT ON COLUMN "lsp_one_dose_scans"."url" IS 'URL to access the reference.';

COMMENT ON TABLE "lsp_clinical_info" IS 'Table of the clinical approval status of compounds.     Sourced from ChEMBL';

COMMENT ON COLUMN "lsp_clinical_info"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON TABLE "lsp_commercial_availability" IS 'Table of the commercial availability of compounds.     Sourced from ZINC (https://zinc.docking.org/).';

COMMENT ON COLUMN "lsp_commercial_availability"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_commercial_availability"."vendor" IS 'The name of the vendor';

COMMENT ON COLUMN "lsp_commercial_availability"."id" IS 'Catalog number of the compound.';

COMMENT ON COLUMN "lsp_commercial_availability"."name" IS 'Compound name from the vendor, if available.';

COMMENT ON TABLE "lsp_fingerprints" IS 'Table of specificity assertions of compounds to their targets. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "lsp_fingerprints"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_fingerprints"."fingerprint_type" IS 'Three different fingerprint types are available: Morgan fingerprints (either ignoring or respecting stereochemistry) or RDKit topological fingerprints (ignoring stereochemistry).';

COMMENT ON COLUMN "lsp_fingerprints"."fingerprint" IS '256 byte hex-encoded fingerprint';

COMMENT ON TABLE "lsp_compound_library" IS 'Library of optimal compounds for each target. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "lsp_compound_library"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "lsp_compound_library"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "lsp_compound_library"."rank" IS 'Ranking of compounds for the same target. First compound is best.';

COMMENT ON COLUMN "lsp_compound_library"."reason_included" IS 'Compounds can be included either because of their selectivity or because they are approved/in clinical trials.';
