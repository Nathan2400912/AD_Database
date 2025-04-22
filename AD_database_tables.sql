SET FOREIGN_KEY_CHECKS = 0;

-- Drop tables in reverse order of dependency
DROP TABLE IF EXISTS TF_CRE_Interactions;
DROP TABLE IF EXISTS CRE_Gene_Interactions;
DROP TABLE IF EXISTS Gene_Pathway_Associations;
DROP TABLE IF EXISTS Differential_Expression;
DROP TABLE IF EXISTS Cis_Regulatory_Elements;
DROP TABLE IF EXISTS Merged_CRES;
DROP TABLE IF EXISTS Transcription_Factors;
DROP TABLE IF EXISTS Biological_Pathways;
DROP TABLE IF EXISTS Genes;
DROP TABLE IF EXISTS Cell_Type;
DROP TABLE IF EXISTS Conditions;

-- Re-enable foreign key checks
SET FOREIGN_KEY_CHECKS = 1;

CREATE TABLE Genes ( 
    gid INT not null auto_increment,
    gene_symbol VARCHAR(50),
    Ensembl_ID VARCHAR(50),
    Entrez_ID VARCHAR(50),
    chromosome VARCHAR(50),
    start_position BIGINT,
    end_position BIGINT,
    strand enum('+', '-'),
    unique (Entrez_ID),
    primary key (gid));

CREATE Table Cell_Type (
    cell_id INT not null auto_increment,
    cell varchar(30) not null unique,
    Primary key(cell_id));

CREATE TABLE Conditions (
    cdid INT not null auto_increment,
    name VARCHAR(100) not null unique,
    disease_category VARCHAR(100),
    Primary key (cdid));

CREATE TABLE Differential_Expression (
    gid INT not null,
    cdid INT not null,
    cell_id INT not null,
    baseMean FLOAT,
    log2foldchange FLOAT,
    p_value DOUBLE,
    padj DOUBLE,
    FOREIGN KEY (gid) REFERENCES Genes(gid), -- merge based on entrez symbol
    FOREIGN KEY (cdid) REFERENCES Conditions(cdid), -- merge based on condition name
    Foreign key (cell_id) references Cell_Type (cell_id), -- merge based on cell type
    Primary key (gid, cdid, cell_id));

CREATE TABLE Cis_Regulatory_Elements ( 
    cid INT not null auto_increment,
    cdid INT not null,
    cell_id INT not null,
    chromosome VARCHAR(50),
    start_position BIGINT,
    end_position BIGINT,
    cre_log2foldchange FLOAT,
    mcid INT not null,
    FOREIGN KEY (mcid) REFERENCES Merged_CRES(mcid), -- merge based on merge cre chr, start and end
    FOREIGN KEY (cdid) REFERENCES Conditions(cdid), -- same as above
    Foreign key (cell_id) references Cell_Type (cell_id),
    Primary key (cid),
    UNIQUE (cdid, cell_id, chromosome, start_position, end_position));

CREATE TABLE Merged_CRES (
    mcid INT not null auto_increment,
    chromosome VARCHAR(50),
    start_position BIGINT,
    end_position BIGINT,
    Primary key (mcid),
    UNIQUE (chromosome, start_position, end_position));

CREATE TABLE Transcription_Factors (
    tfid INT not null auto_increment,
    name VARCHAR(50) not null unique,
    Primary key (tfid));

CREATE TABLE CRE_Gene_Interactions ( -- many to many  
    cid INT not null,
    gid INT not null,
    distance_to_TSS INT,
    FOREIGN KEY (cid) REFERENCES Cis_Regulatory_Elements(cid), -- same as above
    FOREIGN KEY (gid) REFERENCES Genes(gid), -- merge based on entrez
    Primary key (gid, cid));

CREATE TABLE TF_CRE_Interactions ( -- quartnery relationship
    tfid INT not null,
    mcid INT not null,
    cdid INT not null,
    cell_id INT not null,
    FOREIGN KEY (tfid) REFERENCES Transcription_Factors(tfid), -- merge based on tf name 
    FOREIGN KEY (mcid) REFERENCES Merged_CRES(mcid), -- merge based on merge cre chr, start and end
    FOREIGN KEY (cdid) REFERENCES Conditions(cdid), -- same as above
    Foreign key (cell_id) references Cell_Type (cell_id), -- same as above
    Primary key (tfid, mcid, cdid, cell_id));

CREATE TABLE Biological_Pathways (
    pid INT not null auto_increment,
	name VARCHAR(500) not null unique,
    Primary key (pid));

CREATE TABLE Gene_Pathway_Associations ( 
    gid INT not null,
    pid INT not null,   
    FOREIGN KEY (gid) REFERENCES Genes(gid), -- merge based on entrez
    FOREIGN KEY (pid) REFERENCES Biological_Pathways(pid), -- merge based on pathway name 
    Primary key (gid, pid));








