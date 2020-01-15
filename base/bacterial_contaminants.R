#!/usr/bin/Rscript

# Load libraries
library("tidyverse")

# Read seqtables
files <- list.files(path = "./QC/step_9/", pattern = "*_bac_contaminants_collapsed.txt", full.names = TRUE)

# Reduce to a single table
con.table <- files %>%
	map(read_tsv) %>%
	reduce(full_join, by = "bacteria")

# Write joined table
write_tsv(con.table, path = "./QC/step_9/bacterial_contaminant.table", col_names=TRUE)
