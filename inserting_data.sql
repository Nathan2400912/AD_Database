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

