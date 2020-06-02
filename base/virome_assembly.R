#!/usr/bin/env Rscript

# Load libraries
library("tidyverse")
library("data.table")
library("phyloseq")
library("speedyseq")

#### Combine FPKM tables
files <- list.files(path = "./assembly/contig_dictionary", pattern = "*.fpkm.filt", full.names = TRUE)

# Reduce filtered fpkm tables into a single table
fpkm.all <- files %>%
  map(read_tsv) %>%
  reduce(full_join, by = "ID") %>%
  replace(is.na(.), 0)

fpkmToTpm <- function(fpkm) {

exp(log(fpkm) - log(sum(fpkm)) + log(1e6))

}

###################
#### Load data ####
###################

#### USER INPUT ####
# Load mapping file
MAP <- import_qiime_sample_data("./results/map_phage_attack_virome.txt")

# Load contig taxonomic annotation
tax_table <- fread("./assembly/contig_dictionary/contig.taxonomy", header=TRUE, stringsAsFactors=TRUE, fill=TRUE)

# Import coverage table
covtable <- fread("./assembly/contig_dictionary/covtable.all", header=TRUE, sep = "\t")

# Import per contig CAT scores
CAT.scores <- fread("./assembly/contig_dictionary/CAT.scores", header=TRUE, sep = "\t", fill=TRUE)
