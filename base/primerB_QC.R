#!/usr/bin/env Rscript

# Load libraries
library("tidyverse")

# Read seqtables
files.num_reads <- list.files(path = "./QC/step_1/", pattern = "*_primerB_num_reads.out", full.names = TRUE)
files.perc_reads <- list.files(path = "./QC/step_1/", pattern = "*_primerB_perc_reads.out", full.names = TRUE)

# Join
primerB.num_reads <- files.num_reads %>%
  map(read_tsv) %>%
  reduce(full_join, by = "primerB")

primerB.perc_reads <- files.perc_reads %>%
  map(read_tsv) %>%
  reduce(full_join, by = "primerB")

# Write
write_tsv(primerB.num_reads, path = "./QC/stats/primerB.num_reads.txt", col_names=TRUE)
write_tsv(primerB.perc_reads, path = "./QC/stats/primerB.perc_reads.txt", col_names=TRUE)
