#!/bin/bash

# Download NCBI taxonomy infomration
wget ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz;
mkdir taxonomy && tar -xxvf taxdump.tar.gz -C taxonomy;

# Make blastDB
makeblastdb -in refseq_NN.fasta -dbtype nucl -out nt;

# Extract FASTA and taxonomy mapping files
blastdbcmd -db nt -entry all > nt.fna;
blastdbcmd -db nt -entry all -outfmt "%a %T" > nt.fna.taxidmapping;

# Create and annotate mmseqs DB
mmseqs createdb nt.fna nt.fnaDB;
mmseqs createtaxdb nt.fnaDB tmp --ncbi-tax-dump taxonomy/ --tax-mapping-file nt.fna.taxidmapping;
