CREATE TABLE Genes (
    gid INT not null auto_increment,
    gene_symbol VARCHAR(50),
    Ensembl_ID VARCHAR(50),
    Entrez_ID VARCHAR(50),
    chromosome INT,
    start_position INT,
    end_position INT,
    strand enum('+', '-'),
    primary key (gid));

CREATE Table Cell_Type (
    cell_id not null auto_increment,
    Cell_type varchar(30),
    Primary key(cell_id));

CREATE TABLE Conditions (
    cdid INT not null auto_increment,
    name VARCHAR(100),
    disease_category VARCHAR(100),
    Primary key (cdid));

CREATE TABLE Differential_Expression (
    gid INT not null,
    cdid INT not null,
    cell_id INT not null,
    lfcSE FLOAT,
    stat FLOAT,
    log2foldchange FLOAT,
    p_value FLOAT,
    padj FLOAT,
    -- expression_status ENUM('upregulated', 'downregulated', 'not significant'),
    FOREIGN KEY (gid) REFERENCES Genes(gid), -- merge based on entrez symbol
    FOREIGN KEY (cdid) REFERENCES Conditions(cdid), -- merge based on condition name
    Foreign key (cell_id) references CELL_TYPE (cell_id), -- merge based on cell type
    Primary key (gid, cdid, cell_id));

CREATE TABLE Cis_Regulatory_Elements (
    cid INT not null auto_increment,
    chromosome INT,
    start_position INT,
    end_position INT,
    mcid INT not null,
    FOREIGN KEY (mcid) REFERENCES Merged_CRES(mcid), -- merge based on merge cre chr, start and end
    Primary key (cid));

CREATE TABLE CRE_condition ( -- ternary relationship
    cid INT not null,
    cdid INT not null,
    cell_id INT not null,
    cre_log2foldchange FLOAT,
    FOREIGN KEY (cid) REFERENCES Cis_Regulatory_Elements(cid), -- merge based on chr, start and end
    FOREIGN KEY (cdid) REFERENCES Conditions(cdid), -- same as above
    Foreign key (cell_id) references CELL_TYPE (cell_id), -- same as above
    Primary key (cid, cdid, cell_id));

CREATE TABLE Merged_CRES (
    mcid INT not null auto_increment,
    chromosome INT,
    start_position INT,
    end_position INT,
    Primary key (mcid));

CREATE TABLE Transcription_Factors (
    tfid INT not null auto_increment,
    name VARCHAR(50),
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
    cid INT not null,
    cdid INT not null,
    cell_id not null,
    FOREIGN KEY (tfid) REFERENCES Transcription_Factors(tfid), -- merge based on tf name 
    FOREIGN KEY (cid) REFERENCES Cis_Regulatory_Elements(cid), -- same as above
    FOREIGN KEY (cdid) REFERENCES Conditions(cdid), -- same as above
    Foreign key (cell_id) references cell_type (cell_id), -- same as above
    Primary key (tfid, cid, cdid, cell_id));

CREATE TABLE Biological_Pathways (
    pid INT not null auto_increment,
	name VARCHAR(100),
    Primary key (pid));

CREATE TABLE Gene_Pathway_Associations (
    gid INT not null,
    pid INT not null,
    FOREIGN KEY (gid) REFERENCES Genes(gid), -- merge based on entrez
    FOREIGN KEY (pid) REFERENCES Biological_Pathways(pid), -- merge based on pathway name 
    Primary key (gid, pid));








