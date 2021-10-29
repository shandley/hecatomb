#!/usr/bin/env Rscript

# Load libraries
library("tidyverse")

# Read
files.num_reads <- list.files(path = "../STATS/primer_stats/", pattern = "*.primer_stats_num_reads.out", full.names = TRUE)
files.perc_reads <- list.files(path = "../STATS/primer_stats/", pattern = "*.primer_stats_percent.out", full.names = TRUE)

# Join
primerB.num_reads <- files.num_reads %>%
  map(read_tsv) %>%
  reduce(full_join, by = "primerB")

primerB.perc_reads <- files.perc_reads %>%
  map(read_tsv) %>%
  reduce(full_join, by = "primerB")

# Write
write_tsv(primerB.num_reads, path = "../STATS/primer_stats/primer_stats.num_reads.tsv", col_names=TRUE)
write_tsv(primerB.perc_reads, path = "../STATS/primer_stats/primer_stats.perc_reads.tsv", col_names=TRUE)
