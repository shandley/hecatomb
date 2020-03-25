#!/usr/bin/env Rscript

# Load libraries
library("tidyverse")

# Read seqtables
files <- list.files(path = "./QC/step_8/clustered/", pattern = "*_seqtable.txt", full.names = TRUE)

# Reduce seqtables to a single table
seqtable.all <- files %>%
  map(read_tsv) %>%
  reduce(full_join, by = "sequence") %>%
  mutate(id = row_number() - 1) %>%
  select(id, everything()) %>%
  mutate_if(is.numeric, as.integer)

# Write count table
dir.create(path = "results", showWarnings = FALSE)
write_tsv(seqtable.all, path = "./results/seqtable.all", col_names=TRUE)

# Write	tabular	fasta (tab2fx)
seqs <- tibble(`sequence` = seqtable.all$sequence)
seqs.df <- seqs %>%
	mutate(id = row_number() - 1) %>%
	select(id, everything()) %>%
	mutate_if(is.numeric, as.integer)
write_tsv(seqs.df, "./results/seqtable.tab2fx", col_names = FALSE)
