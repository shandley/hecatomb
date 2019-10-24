#!/usr/bin/Rscript

# Load libraries
library("tidyverse")

# Read seqtables
files <- list.files(path = "./QC/step_5/", pattern = "*_edited", full.names = TRUE)

# Reduce to a single table
con.table <- files %>%
	map(read_tsv) %>%
	reduce(full_join, by = "contaminant")

# Write joined table
write_tsv(con.table, path = "./QC/step_5/contaminant.table", col_names=TRUE)
