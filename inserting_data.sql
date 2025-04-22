LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/combined_genes.csv'
INTO TABLE Genes
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(gene_symbol, Ensembl_ID, Entrez_ID, chromosome, start_position, end_position, strand);

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/cell_type.csv'
INTO TABLE Cell_Type
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(cell);

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/conditions.csv'
INTO TABLE Conditions
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(name, disease_category);

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/merged_cres.csv'
INTO TABLE Merged_CRES
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(chromosome, start_position, end_position);

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/tfs.csv'
INTO TABLE Transcription_Factors
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(name);

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/pathways.csv'
INTO TABLE Biological_Pathways
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(name);

-- insert relationship tables through intermediate tables 
CREATE TABLE import_differential_expression (
    entrez VARCHAR(50) NOT NULL,
    `condition` VARCHAR(100) NOT NULL,
    cell_type VARCHAR(50) NOT NULL,
    baseMean FLOAT,
    log2foldchange FLOAT,
    p_value DOUBLE,
    padj DOUBLE
);

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/combined_de.csv'
INTO TABLE import_differential_expression
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(entrez, `condition`, cell_type, baseMean, log2foldchange, p_value, padj)

INSERT INTO Differential_Expression (gid, cdid, cell_id, baseMean, log2foldchange, p_value, padj)
SELECT 
  g.gid,
  c.cdid,
  ct.cell_id,
  i.baseMean,
  i.log2foldchange,
  i.p_value,
  i.padj
FROM import_differential_expression i
JOIN Genes g ON g.Entrez_ID = i.entrez
JOIN Conditions c ON c.name = i.`condition`
JOIN Cell_Type ct ON ct.cell = i.cell_type;

CREATE TABLE import_cres ( 
    `condition` VARCHAR(100) NOT NULL,
    cell_type VARCHAR(50) NOT NULL,
    chromosome VARCHAR(50),
    start_position BIGINT,
    end_position BIGINT,
    cre_log2foldchange FLOAT,
    merged_chromosome VARCHAR(50),
    merged_start_position BIGINT,
    merged_end_position BIGINT);

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/unique_cres.csv'
INTO TABLE import_cres
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(`condition`, cell_type, chromosome, start_position, end_position, cre_log2foldchange, merged_chromosome, merged_start_position, merged_end_position);

INSERT INTO Cis_Regulatory_Elements (cdid, cell_id, chromosome, start_position, end_position, cre_log2foldchange, mcid)
SELECT 
  c.cdid,
  ct.cell_id,
  i.chromosome,
  i.start_position,
  i.end_position,
  i.cre_log2foldchange,
  mc.mcid
FROM import_cres i
JOIN Conditions c ON c.name = i.`condition`
JOIN Cell_Type ct ON ct.cell = i.cell_type
JOIN Merged_CRES mc ON mc.chromosome = i.merged_chromosome
  AND mc.start_position = i.merged_start_position
  AND mc.end_position = i.merged_end_position;

CREATE TABLE import_cres_gene ( 
    `condition` VARCHAR(100) NOT NULL,
    cell_type VARCHAR(50) NOT NULL,
    chromosome VARCHAR(50),
    start_position BIGINT,
    end_position BIGINT,
    entrez VARCHAR(50),
    distance_to_tss INT);

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/unique_cres_genes.csv'
INTO TABLE import_cres_gene
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(`condition`, cell_type, chromosome, start_position, end_position, entrez, distance_to_tss);

-- indexes to improve loading
CREATE INDEX idx_cre_chromosome_start_end ON Cis_Regulatory_Elements (chromosome, start_position, end_position, cdid, cell_id);
CREATE INDEX idx_conditions_name ON Conditions (name);
CREATE INDEX idx_cell_type_cell ON Cell_Type (cell);
CREATE INDEX idx_genes_entrez_id ON Genes (Entrez_ID);

INSERT INTO CRE_Gene_Interactions (cid, gid, distance_to_TSS)
SELECT 
  cre.cid,          -- Reference to Cis_Regulatory_Elements table's cid
  g.gid,            -- Reference to Genes table's gid
  i.distance_to_tss -- Distance to TSS from import_cres_gene
FROM import_cres_gene i
JOIN Genes g ON g.Entrez_ID = i.entrez
JOIN Conditions c ON c.name = i.`condition`
JOIN Cell_Type ct ON ct.cell = i.cell_type
JOIN Cis_Regulatory_Elements cre 
  ON cre.chromosome = i.chromosome
  AND cre.start_position = i.start_position
  AND cre.end_position = i.end_position
  AND cre.cdid = c.cdid              -- Ensure the correct cxsondition match
  AND cre.cell_id = ct.cell_id; 

CREATE TABLE import_gene_pathways ( 
    entrez VARCHAR(50),
    pathway VARCHAR(500));

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/filtered_gene_pathways.csv'
INTO TABLE import_gene_pathways
FIELDS TERMINATED BY ',' 
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(entrez, pathway);

INSERT INTO Gene_Pathway_Associations (gid, pid)
SELECT
  g.gid,
  p.pid
FROM import_gene_pathways i
JOIN Genes g ON g.Entrez_ID = i.entrez
JOIN Biological_Pathways p ON p.name = i.pathway;

CREATE TABLE import_cre_tfs (
    merged_chromosome VARCHAR(50),
    merged_start_position BIGINT,
    merged_end_position BIGINT,
    transcription_factor VARCHAR(50),
    `condition` VARCHAR(100) NOT NULL,
    cell_type VARCHAR(50) NOT NULL);

LOAD DATA LOCAL INFILE '/Users/nathan/Desktop/Bioinformatics/Projects/AD_Database/Tables/merged_cre_tf.csv'
INTO TABLE import_cre_tfs
FIELDS TERMINATED BY ','
LINES TERMINATED BY '\n'
IGNORE 1 ROWS
(merged_chromosome, merged_start_position, merged_end_position, transcription_factor, `condition`, cell_type);

-- indexes to improve loading
CREATE INDEX idx_tf_name ON Transcription_Factors (name);
CREATE INDEX idx_tf_cre_chromosome_start_end ON Merged_CRES (chromosome, start_position, end_position);
CREATE INDEX idx_import_cre_tfs_lookup 
ON import_cre_tfs (transcription_factor, `condition`, cell_type, merged_chromosome, merged_start_position, merged_end_position);

INSERT INTO TF_CRE_Interactions (tfid, mcid, cdid, cell_id)
SELECT 
  tf.tf_id,          -- Reference to Transcription_Factors table's tfid
  cre.mcid,          -- Reference to Merged_CRES table's mcid
  c.cdid,            -- Reference to Conditions table's cdid
  ct.cell_id         -- Reference to Cell_Type table's cell_id
FROM import_cre_tfs i 
JOIN Transcription_Factors tf ON tf.name = i.transcription_factor
JOIN Conditions c ON c.name = i.`condition`
JOIN Cell_Type ct ON ct.cell = i.cell_type
JOIN Merged_CRES cre 
  ON cre.chromosome = i.merged_chromosome
  AND cre.start_position = i.merged_start_position
  AND cre.end_position = i.merged_end_position;

-- Clean up temporary tables
DROP TABLE IF EXISTS import_differential_expression;
DROP TABLE IF EXISTS import_cres;
DROP TABLE IF EXISTS import_cres_gene;
DROP TABLE IF EXISTS import_gene_pathways;
DROP TABLE IF EXISTS import_cre_tfs;

