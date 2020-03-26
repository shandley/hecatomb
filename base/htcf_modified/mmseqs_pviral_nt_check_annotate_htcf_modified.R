#!/usr/bin/env Rscript

# Load libraries
library("tidyverse")
library("taxonomizr")

# Read *.m8 file with best hit NCBI accession number
m8 <- read_tsv("./results/mmseqs_nt_checked_out/resultDB.firsthit.m8", col_names = FALSE)
colnames(m8) <- c("query","target","pident","alnlen","mismatch","gapopen","qstart","qend","tstart","tend","evalue","bits")

# Create NCBI accession number vector
acc <- m8$target

# Convert accessions to taxonomic IDs
ids <- accessionToTaxa(acc, "/scratch/ref/hecatomb_databases/taxonomizr/accessionTaxa.sql")

# Get taxonomic lineage
ncbi_tax <- as_tibble(getTaxonomy(ids, "/scratch/ref/hecatomb_databases/taxonomizr/accessionTaxa.sql"))

# Bind m8 file to lineage
seqids <- select(m8, "query")
mmseqs_pviral_nt_lineage <- cbind(seqids, ncbi_tax)

# Write results to table
write_tsv(mmseqs_pviral_nt_lineage, path = "./results/mmseqs_nt_checked_out/mmseqs_pviral_nt_checked_lineage.tsv")
