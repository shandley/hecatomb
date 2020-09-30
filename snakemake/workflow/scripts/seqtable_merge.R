#!/usr/bin/env Rscript

# Load libraries
library("tidyverse")

# Read seqtables
files <- snakemake@input[["files"]]

# Reduce seqtables to a single table
seqtable.all <- files %>%
  map(read_tsv) %>%
  reduce(full_join, by = "sequence") %>%
  mutate(id = row_number() - 1) %>%
  select(id, everything()) %>%
  mutate_if(is.numeric, as.integer)

# Write count table
dir.create(path = snakemake@params[["resultsdir"]], showWarnings = FALSE)
write_tsv(seqtable.all, path = snakemake@output[["seqtable"]], col_names=TRUE)

# Write	tabular	fasta (tab2fx)
seqs <- tibble(`sequence` = seqtable.all$sequence)
seqs.df <- seqs %>%
	mutate(id = row_number() - 1) %>%
	select(id, everything()) %>%
	mutate_if(is.numeric, as.integer)
write_tsv(seqs.df, path = snakemake@output[["tab2fx"]], col_names = FALSE)
