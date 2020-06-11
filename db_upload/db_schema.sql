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

CREATE TABLE "LSP_COMPOUND_DICTIONARY" (
  "lspci_id" int PRIMARY KEY,
  "hmsl_id" varchar,
  "chembl_id" varchar,
  "pref_name" varchar,
  "inchi" varchar,
  "smiles" varchar,
  "commercially_available" boolean
);

CREATE TABLE "LSP_COMPOUND_NAMES" (
  "lspci_id" int,
  "source" compound_sources,
  "priority" name_priorities,
  "name" varchar
);

CREATE TABLE "LSP_COMPOUND_MAPPING" (
  "lspci_id" int,
  "source" compound_sources,
  "id" varchar
);

CREATE TABLE "LSP_TARGET_DICTIONARY" (
  "gene_id" int PRIMARY KEY,
  "symbol" varchar,
  "pref_name" varchar,
  "description" varchar,
  "tax_id" int,
  "organism" varchar
);

CREATE TABLE "LSP_TARGET_MAPPING" (
  "gene_id" int,
  "chembl_id" varchar,
  "uniprot_id" varchar
);

CREATE TABLE "LSP_BIOCHEM" (
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

CREATE TABLE "LSP_PHENOTYPIC_CHEMBL" (
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

CREATE TABLE "LSP_TAS" (
  "lspci_id" int,
  "gene_id" int,
  "tas" int
);

CREATE TABLE "LSP_SPECIFICITY" (
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

CREATE TABLE "LSP_ONE_DOSE_SCANS" (
  "lspci_id" int,
  "gene_id" int,
  "percent_control" float,
  "description" varchar,
  "cmpd_conc_nM" float,
  "reference_id" varchar,
  "reference_type" reference_types,
  "url" varchar
);

CREATE TABLE "LSP_CLINICAL_INFO" (
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

CREATE TABLE "LSP_COMMERCIAL_AVAILABILITY" (
  "lspci_id" int,
  "vendor" varchar,
  "id" varchar,
  "name" varchar
);

CREATE TABLE "LSP_FINGERPRINTS" (
  "lspci_id" int,
  "fingerprint_type" fingerprint_types,
  "fingerprint" char(512)
);

CREATE TABLE "LSP_COMPOUND_LIBRARY" (
  "lspci_id" int,
  "gene_id" int,
  "rank" int,
  "reason_included" include_reasons
);

ALTER TABLE "LSP_COMPOUND_NAMES" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_COMPOUND_MAPPING" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_TARGET_MAPPING" ADD FOREIGN KEY ("gene_id") REFERENCES "LSP_TARGET_DICTIONARY" ("gene_id");

ALTER TABLE "LSP_BIOCHEM" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_BIOCHEM" ADD FOREIGN KEY ("gene_id") REFERENCES "LSP_TARGET_DICTIONARY" ("gene_id");

ALTER TABLE "LSP_PHENOTYPIC_CHEMBL" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_TAS" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_TAS" ADD FOREIGN KEY ("gene_id") REFERENCES "LSP_TARGET_DICTIONARY" ("gene_id");

ALTER TABLE "LSP_SPECIFICITY" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_SPECIFICITY" ADD FOREIGN KEY ("gene_id") REFERENCES "LSP_TARGET_DICTIONARY" ("gene_id");

ALTER TABLE "LSP_ONE_DOSE_SCANS" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_ONE_DOSE_SCANS" ADD FOREIGN KEY ("gene_id") REFERENCES "LSP_TARGET_DICTIONARY" ("gene_id");

ALTER TABLE "LSP_CLINICAL_INFO" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_COMMERCIAL_AVAILABILITY" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_FINGERPRINTS" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_COMPOUND_LIBRARY" ADD FOREIGN KEY ("lspci_id") REFERENCES "LSP_COMPOUND_DICTIONARY" ("lspci_id");

ALTER TABLE "LSP_COMPOUND_LIBRARY" ADD FOREIGN KEY ("gene_id") REFERENCES "LSP_TARGET_DICTIONARY" ("gene_id");

CREATE UNIQUE INDEX ON "LSP_COMPOUND_DICTIONARY" ("lspci_id");

CREATE INDEX ON "LSP_COMPOUND_NAMES" ("lspci_id");

CREATE INDEX ON "LSP_COMPOUND_MAPPING" ("lspci_id");

CREATE UNIQUE INDEX ON "LSP_TARGET_DICTIONARY" ("gene_id");

CREATE INDEX ON "LSP_TARGET_DICTIONARY" ("symbol");

CREATE INDEX ON "LSP_TARGET_MAPPING" ("gene_id");

CREATE INDEX ON "LSP_BIOCHEM" ("lspci_id");

CREATE INDEX ON "LSP_BIOCHEM" ("gene_id");

CREATE INDEX ON "LSP_BIOCHEM" ("lspci_id", "gene_id");

CREATE INDEX ON "LSP_PHENOTYPIC_CHEMBL" ("lspci_id");

CREATE INDEX ON "LSP_PHENOTYPIC_CHEMBL" ("assay_id");

CREATE INDEX ON "LSP_TAS" ("lspci_id");

CREATE INDEX ON "LSP_TAS" ("gene_id");

CREATE INDEX ON "LSP_TAS" ("lspci_id", "gene_id");

CREATE INDEX ON "LSP_SPECIFICITY" ("lspci_id");

CREATE INDEX ON "LSP_SPECIFICITY" ("gene_id");

CREATE INDEX ON "LSP_SPECIFICITY" ("lspci_id", "gene_id");

CREATE INDEX ON "LSP_ONE_DOSE_SCANS" ("lspci_id");

CREATE INDEX ON "LSP_ONE_DOSE_SCANS" ("gene_id");

CREATE INDEX ON "LSP_ONE_DOSE_SCANS" ("lspci_id", "gene_id");

CREATE INDEX ON "LSP_CLINICAL_INFO" ("lspci_id");

CREATE INDEX ON "LSP_CLINICAL_INFO" ("max_phase");

CREATE INDEX ON "LSP_COMMERCIAL_AVAILABILITY" ("lspci_id");

CREATE INDEX ON "LSP_FINGERPRINTS" ("lspci_id");

CREATE INDEX ON "LSP_COMPOUND_LIBRARY" ("lspci_id");

CREATE INDEX ON "LSP_COMPOUND_LIBRARY" ("gene_id");

CREATE INDEX ON "LSP_COMPOUND_LIBRARY" ("lspci_id", "gene_id");

COMMENT ON TABLE "LSP_COMPOUND_DICTIONARY" IS 'Primary table listing all compounds in the database.     During compound processing distinct salts of the same compound     are aggregated into a single compound entry in this table.     The constituent compound IDs for each compound in this table     are available in the LSP_COMPOUND_MAPPING table.';

COMMENT ON COLUMN "LSP_COMPOUND_DICTIONARY"."lspci_id" IS 'Internal compound ID';

COMMENT ON COLUMN "LSP_COMPOUND_DICTIONARY"."hmsl_id" IS 'Primary HMS LINCS compound ID, if available';

COMMENT ON COLUMN "LSP_COMPOUND_DICTIONARY"."chembl_id" IS 'Primary ChEMBL compound ID, if available';

COMMENT ON COLUMN "LSP_COMPOUND_DICTIONARY"."pref_name" IS 'Preferred name of the compound';

COMMENT ON COLUMN "LSP_COMPOUND_DICTIONARY"."inchi" IS 'InChI chemical indentifier for the compound, standardized     using https://molvs.readthedocs.io/en/latest/guide/standardize.html';

COMMENT ON COLUMN "LSP_COMPOUND_DICTIONARY"."smiles" IS 'SMILES chemical indentifier for the compound, standardized     using https://molvs.readthedocs.io/en/latest/guide/standardize.html';

COMMENT ON COLUMN "LSP_COMPOUND_DICTIONARY"."commercially_available" IS 'Commercially available at a vendor in the ZINC database';

COMMENT ON TABLE "LSP_COMPOUND_NAMES" IS 'Table of all annotated names for compounds. The sources     for compound names generally distinguish between primary and     alternative (secondary) names.';

COMMENT ON COLUMN "LSP_COMPOUND_NAMES"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_COMPOUND_NAMES"."source" IS 'Source for the compound name';

COMMENT ON COLUMN "LSP_COMPOUND_NAMES"."priority" IS 'Annotated as primary or secondary name at source';

COMMENT ON COLUMN "LSP_COMPOUND_NAMES"."name" IS 'Compound name';

COMMENT ON TABLE "LSP_COMPOUND_MAPPING" IS 'Table of mappings between compound IDs from different sources to     the internal lspci_ids.';

COMMENT ON COLUMN "LSP_COMPOUND_MAPPING"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_COMPOUND_MAPPING"."source" IS 'Source for the compound ID';

COMMENT ON COLUMN "LSP_COMPOUND_MAPPING"."id" IS 'Compound ID associated with the given lspci_id';

COMMENT ON TABLE "LSP_TARGET_DICTIONARY" IS 'Table of drug targets. The original drug targets are mostly     annotated as ChEMBL or UniProt IDs. For convenience we converted     these IDs to Entrez gene IDs.';

COMMENT ON COLUMN "LSP_TARGET_DICTIONARY"."gene_id" IS 'Entrez gene ID';

COMMENT ON COLUMN "LSP_TARGET_DICTIONARY"."symbol" IS 'Entrez gene symbol';

COMMENT ON COLUMN "LSP_TARGET_DICTIONARY"."pref_name" IS 'Preferred gene name';

COMMENT ON COLUMN "LSP_TARGET_DICTIONARY"."description" IS 'Description of gene function';

COMMENT ON COLUMN "LSP_TARGET_DICTIONARY"."tax_id" IS 'Entrez taxonomy ID';

COMMENT ON COLUMN "LSP_TARGET_DICTIONARY"."organism" IS 'Organism for which gene is annotated';

COMMENT ON TABLE "LSP_TARGET_MAPPING" IS 'Mapping between the original ChEMBL target IDs,     their corresponding UniProt IDs and Entrez gene IDs.     A signle UniProt or ChEMBL ID can refer to protein complexes,     therefore multiple gene IDs often map to the same UniProt or     ChEMBL ID.';

COMMENT ON COLUMN "LSP_TARGET_MAPPING"."gene_id" IS 'Foreign key to Entrez gene ID';

COMMENT ON COLUMN "LSP_TARGET_MAPPING"."chembl_id" IS 'ChEMBL target ID';

COMMENT ON COLUMN "LSP_TARGET_MAPPING"."uniprot_id" IS 'UniProt ID';

COMMENT ON TABLE "LSP_BIOCHEM" IS 'Table of biochemical affinity measurements.';

COMMENT ON COLUMN "LSP_BIOCHEM"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_BIOCHEM"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "LSP_BIOCHEM"."description_assay" IS 'Description of the assay that the measurement is derived from';

COMMENT ON COLUMN "LSP_BIOCHEM"."value" IS 'Measurement value';

COMMENT ON COLUMN "LSP_BIOCHEM"."value_type" IS 'The type of measurement performed.';

COMMENT ON COLUMN "LSP_BIOCHEM"."value_unit" IS 'The unit of the measurement.';

COMMENT ON COLUMN "LSP_BIOCHEM"."value_relation" IS 'Some assays can't determine the measured value precisely. This column gives relationship between the actual and the measured value.';

COMMENT ON COLUMN "LSP_BIOCHEM"."reference_id" IS 'The reference for the measurement.';

COMMENT ON COLUMN "LSP_BIOCHEM"."reference_type" IS 'The source of the measurement.';

COMMENT ON COLUMN "LSP_BIOCHEM"."url" IS 'URL to access the reference.';

COMMENT ON TABLE "LSP_PHENOTYPIC_CHEMBL" IS 'Table of phenotypic assays performed on the compounds.';

COMMENT ON COLUMN "LSP_PHENOTYPIC_CHEMBL"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_PHENOTYPIC_CHEMBL"."assay_id" IS 'Unique ID of the assays in the database. Correspond to ChEMBL assay IDs.';

COMMENT ON COLUMN "LSP_PHENOTYPIC_CHEMBL"."value" IS 'Measurement value';

COMMENT ON COLUMN "LSP_PHENOTYPIC_CHEMBL"."value_type" IS 'The type of measurement performed.';

COMMENT ON COLUMN "LSP_PHENOTYPIC_CHEMBL"."value_unit" IS 'The unit of the measurement.';

COMMENT ON COLUMN "LSP_PHENOTYPIC_CHEMBL"."value_relation" IS 'Some assays can't determine the measured value precisely. This column gives relationship between the actual and the measured value.';

COMMENT ON COLUMN "LSP_PHENOTYPIC_CHEMBL"."reference_id" IS 'The reference for the measurement.';

COMMENT ON COLUMN "LSP_PHENOTYPIC_CHEMBL"."reference_type" IS 'The source of the measurement.';

COMMENT ON COLUMN "LSP_PHENOTYPIC_CHEMBL"."url" IS 'URL to access the reference.';

COMMENT ON TABLE "LSP_TAS" IS 'Table of Target Affinity Spectrum (TAS) values for the affinity between     compound and target. TAS enables aggregation of affinity measurements from     heterogeneous sources and assays into a single value. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "LSP_TAS"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_TAS"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "LSP_TAS"."tas" IS 'Target Affinity Spectrum (TAS). A number between 1 and 10 for the affinity between compound and target. 1 being most strongly binding.';

COMMENT ON TABLE "LSP_SPECIFICITY" IS 'Table of Target Affinity Spectrum (TAS) values for the affinity between     compound and target. TAS enables aggregation of affinity measurements from     heterogeneous sources and assays into a single value. See     10.1016/j.chembiol.2019.02.018 for details.';

COMMENT ON COLUMN "LSP_SPECIFICITY"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_SPECIFICITY"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "LSP_SPECIFICITY"."selectivity_class" IS 'Assertion for the selectivity of the compound to the given target.';

COMMENT ON COLUMN "LSP_SPECIFICITY"."IC50_difference" IS 'Difference between ontarget and offtarget IC50 Q1 measurements';

COMMENT ON COLUMN "LSP_SPECIFICITY"."ontarget_IC50_Q1" IS 'First quartile of ontarget IC50 measurements';

COMMENT ON COLUMN "LSP_SPECIFICITY"."offtarget_IC50_Q1" IS 'First quartile of offtarget IC50 measurements';

COMMENT ON COLUMN "LSP_SPECIFICITY"."ontarget_N" IS 'Number of ontarget IC50 measurements';

COMMENT ON COLUMN "LSP_SPECIFICITY"."offtarget_N" IS 'Number of offtarget IC50 measurements';

COMMENT ON TABLE "LSP_ONE_DOSE_SCANS" IS 'Table of single dose compound activity measurements as     opposed to full dose-response affinity measurements.';

COMMENT ON COLUMN "LSP_ONE_DOSE_SCANS"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_ONE_DOSE_SCANS"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "LSP_ONE_DOSE_SCANS"."percent_control" IS 'Remaining activity of target at the given compound concentration.';

COMMENT ON COLUMN "LSP_ONE_DOSE_SCANS"."description" IS 'Description of target.';

COMMENT ON COLUMN "LSP_ONE_DOSE_SCANS"."cmpd_conc_nM" IS 'Concentration of the compound.';

COMMENT ON COLUMN "LSP_ONE_DOSE_SCANS"."reference_id" IS 'The reference for the measurement.';

COMMENT ON COLUMN "LSP_ONE_DOSE_SCANS"."reference_type" IS 'The source of the measurement.';

COMMENT ON COLUMN "LSP_ONE_DOSE_SCANS"."url" IS 'URL to access the reference.';

COMMENT ON TABLE "LSP_CLINICAL_INFO" IS 'Table of the clinical approval status of compounds.     Sourced from ChEMBL';

COMMENT ON COLUMN "LSP_CLINICAL_INFO"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON TABLE "LSP_COMMERCIAL_AVAILABILITY" IS 'Table of the commercial availability of compounds.     Sourced from ZINC (https://zinc.docking.org/).';

COMMENT ON COLUMN "LSP_COMMERCIAL_AVAILABILITY"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_COMMERCIAL_AVAILABILITY"."vendor" IS 'The name of the vendor';

COMMENT ON COLUMN "LSP_COMMERCIAL_AVAILABILITY"."id" IS 'Catalog number of the compound.';

COMMENT ON COLUMN "LSP_COMMERCIAL_AVAILABILITY"."name" IS 'Compound name from the vendor, if available.';

COMMENT ON COLUMN "LSP_FINGERPRINTS"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_FINGERPRINTS"."fingerprint_type" IS 'Three different fingerprint types are available: Morgan fingerprints (either ignoring or respecting stereochemistry) or RDKit topological fingerprints (ignoring stereochemistry).';

COMMENT ON COLUMN "LSP_FINGERPRINTS"."fingerprint" IS '512';

COMMENT ON COLUMN "LSP_COMPOUND_LIBRARY"."lspci_id" IS 'Foreign key for compound ID';

COMMENT ON COLUMN "LSP_COMPOUND_LIBRARY"."gene_id" IS 'Foreign key for gene ID';

COMMENT ON COLUMN "LSP_COMPOUND_LIBRARY"."rank" IS 'Ranking of compounds for the same target. First compound is best.';

COMMENT ON COLUMN "LSP_COMPOUND_LIBRARY"."reason_included" IS 'Compounds can be included either because of their selectivity or because they are approved/in clinical trials.';
