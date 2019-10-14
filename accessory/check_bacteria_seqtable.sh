#!/usr/bin/env Rscript

# Load libraries
library("tidyverse")

# Read seqtables
files <- list.files(path = "./results/bacterial_check", pattern = "*_seqtable.txt", full.names = TRUE)

# Reduce seqtables to a single table
seqtable.all <- files %>%
  map(read_tsv) %>%
  reduce(full_join, by = "sequence")

# Write count table
dir.create(path = "./results/bacterial_check/results", showWarnings = FALSE)
write_tsv(seqtable.all, path = "./results/bacterial_check/results/seqtable.all", col_names=TRUE)

# Write	tabular	fasta (tab2fx)
seqs <- tibble(`sequence` = seqtable.all$sequence)
seqs.df <- seqs %>% rownames_to_column(var = "names")
write_tsv(seqs.df, "./results/bacterial_check/results/seqtable.tab2fx", col_names = FALSE)
