#!/usr/bin/Rscript

#### What this script does ####
# Generates 3 output files ready for import into R for additional analysis
# 1) A raw count PhyloSeq object
# 2) A data frame on counts standardized to input library size
# 3) A data frame of individual reads and their alingment statistics

#### Instructions ####
# There are several places that require user input. Each section is preceeded with a #### USER INPUT #### tag
# 1) Need to supply mapping file location for MAP
# 2) Need to change directroy and filenames to write (two saveRDS commands below)
# 3) Joining the library size information requires a common key (column name) and should be set in the Standardization section 

# Load libraries
library("tidyverse")
library("data.table")
library("phyloseq")
library("speedyseq")

###################
#### Load data ####
###################

#### USER INPUT ####
# Load mapping file
MAP <- import_qiime_sample_data("./results/mapping_pot.txt")

# Load viral taxonomy and remove any non-viral sequences
viral_table <- fread(file = "./results/phage_tax_table.tsv", header = TRUE, stringsAsFactors=TRUE)
setkey(viral_table, Kingdom)
viral_table <- viral_table[!"Bacteria"]

# Import and adjust count table
seqtable <- fread(file = "./results/seqtable.all", header = TRUE, sep = "\t", drop = 2)
seqtable <- setnafill(seqtable, fill=0)

# Import alignment files
# aa_checked alignment
aln.aa <- fread(file = "./results/aa.aln.m8", header = TRUE, sep = "\t")
aln.aa <- aln.aa %>%
  rename(id = query)

# nt_checked alignment
aln.nt <- fread(file = "./results/nt.aln.m8", header = TRUE, sep = "\t")
aln.nt <- aln.nt %>%
  rename(id = query)

# Bind aa and nt alignment
aln.aa <- aln.aa %>%
  mutate(query_type = "aa")
aln.nt <- aln.nt %>%
  mutate(query_type = "nt")

aln.all <- bind_rows(aln.aa, aln.nt)
aln.all <- data.table(aln.all)

# Standardization data
# Standardize to per sample library size
col.sums <- as_tibble(colSums(seqtable[,-1]))
samples <- as_tibble(colnames(seqtable[,-1]))
library.size <- bind_cols(samples, col.sums)
rm(col.sums, samples); gc()

#### USER INPUT ####
library.size <- library.size %>%
  rename(sample_ID = value) %>% # Rename to whatever your sequence file name / sample name is in your mapping file
  rename(library_size = value1)
library.size <- data.table(library.size)

###############################################
#### 1) GENERATE RAW COUNT PHYLOSEQ OBJECT ####
###############################################
# Merge viral taxonomy table with alignments stats
alnmerge <- merge(viral_table, aln.all, by = "id")

# Merge viral_aln with count (seqtable) table
stmerge <- merge(viral_table, seqtable, all.x = TRUE, by = "id")

# Sum per taxon counts
merged_counts <- stmerge %>%
  select(-id) %>%
  unite("lineage", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "_", remove = FALSE) %>%
  group_by(lineage, Species) %>%
  summarise_if(is.numeric, funs(sum(as.numeric(.)))) %>%
  separate("lineage", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "_") %>%
  ungroup()

# Average aln stats
merged_aln <- alnmerge %>%
  select(-id) %>%
  unite("lineage", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "_", remove = FALSE) %>%
  group_by(lineage, Species) %>%
  summarise_if(is.numeric, funs(mean(as.numeric(.)))) %>%
  separate("lineage", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "_") %>%
  ungroup()

# Remove taxon information to create final count table (called otu_table in phyloseq)
count_table <- merged_counts[,-c(1:7)]

# Convert summarized table to a matrix of just the taxon in preperation for phyloseq import
tax.tbl.m <- merged_aln %>%
  as.matrix()

# Create PhyloSeq object
OTU <- otu_table(count_table, taxa_are_rows = TRUE)
TAX <- tax_table(tax.tbl.m)
ps0 = phyloseq(OTU, TAX, MAP)

#### USER INPUT ####
# Save PhyloSeq object as an RDS file
saveRDS(ps0, file = "./results/ps0_euk_viruses.RDS")

# Remove some files to free up memory
rm(aln.aa, aln.nt, alnmerge, count_table, merged_aln, merged_counts, OTU, TAX, tax.tbl.m, stmerge); gc()

##################################################
#### 2) GENERATE STANDARDIZED PHYLOSEQ OBJECT ####
##################################################
# Melt PhyloSeq Object
# This may take a lot of time and require a lot of RAM
ps.melt <- ps0 %>%
  speedyseq::psmelt() %>%
  as_tibble()

# Adjust values
ps.melt.value <- ps.melt %>%
  select(c((ncol(MAP)+11):(ncol(ps.melt)))) %>%
   mutate_if(is.factor, ~ as.numeric(levels(.x))[.x])

# Adjust factors
ps.melt.fact <- ps.melt %>%
  select(c(1:(ncol(MAP)+10))) %>%
  mutate_if(is.character, as.factor)

ps.melt.fixed <- bind_cols(ps.melt.fact, ps.melt.value)

# Convert to data table. This may take a bit of time as well
ps.melt.fixed <- data.table(ps.melt.fixed)

# Add Baltimore virus classifications
baltimore <- fread("/mnt/data1/databases/hecatomb/2020_03_11_family_Baltimore_classification_lookup.txt", stringsAsFactors = TRUE)
ps.melt.fixed <- left_join(ps.melt.fixed, baltimore, by = "Family")
ps.melt.fixed$Family <- as.factor(ps.melt.fixed$Family) # left_join may have converted to character
ps.melt.fixed <- data.table(ps.melt.fixed)
rm(baltimore)

# Merge with library size
ps.melt.fixed <- merge(ps.melt.fixed, library.size, all.x = TRUE, by = "sample_ID")

# Add proportion, min and mean and median standardizations
min.ls <- min(library.size$library_size)
mean.ls <- mean(library.size$library_size)
median.ls <- median(library.size$library_size)

ps.melt.scaled <- ps.melt.fixed %>%
        mutate(proportional_abundance = Abundance/library_size) %>%
        mutate(scaled_abundance_min = (min.ls*proportional_abundance)) %>%
        mutate(scaled_abundance_mean = (mean.ls*proportional_abundance)) %>%
        mutate(scaled_abundance_median = (median.ls*proportional_abundance))

ps.melt.scaled <- ps.melt.fixed %>%
        mutate(proportional_abundance = Abundance/library_size) %>%
        mutate(scaled_abundance_min = (min.ls*proportional_abundance)) %>%
        mutate(scaled_abundance_mean = (mean.ls*proportional_abundance)) %>%
        mutate(scaled_abundance_median = (median.ls*proportional_abundance))

ps.melt.scaled <- ps.melt.scaled %>%
        rename(avg_percent_id = percent_id) %>%
        rename(avg_alignment_length = alignment_length) %>%
        rename(avg_num_mismatches = num_mismatches) %>%
        rename(avg_number_gaps = number_gaps) %>%
        rename(avg_start_query = start_query) %>%
        rename(avg_end_query = end_query) %>%
        rename(avg_start_target = start_target) %>%
        rename(avg_end_target = end_target) %>%
        rename(avg_e_value = e_value) %>%
        rename(avg_bit_score = bit_score)

#### USER INPUT ####
# Save PhyloSeq object as an RDS file
saveRDS(ps0, file = "./results/ps0_euk_viruses_standardized.RDS")

#################################################
#### 3) GENERATE INDIVIDUAL READ STATS TABLE ####
#################################################
# Merge viral taxonomy and aa/nt alignment table
setkey(viral_table, id)
setkey(aln.all, id)
vir.aln <- merge(viral_table, aln.all, all.x = TRUE)

# Merge virus + alignment table with count table
setkey(seqtable, id)
aln.merge <- merge(vir.aln, seqtable, all.x = TRUE)
rm(seqtable); gc()

# Extract count table
aln.count <- aln.merge[,-c(1:20)]
OTU.ALN <- otu_table(aln.count, taxa_are_rows = TRUE)

# Extract tax table
aln.tax.tbl <- aln.merge[,c(1:20)] %>%
  as_tibble()

# Convert to matrix in preperation for PhyloSeq
aln.tax.tbl.m <- aln.merge[,c(1:20)] %>%
  as.matrix()

TAX.ALN <- tax_table(aln.tax.tbl.m)

# Create Phyloseq Object of alignment results
ps.aln <- phyloseq(OTU.ALN, TAX.ALN, MAP)

# Adjust PhyloSeq alignment object
# This may take a lot of time and require a lot of RAM
ps.aln.adj <- ps.aln %>%
  speedyseq::psmelt() %>%
  as_tibble()

# Adjust values
ps.aln.value <- ps.aln.adj %>%
  select(c((ncol(MAP)+13):(ncol(ps.aln.adj)-1))) %>%
   mutate_if(is.factor, ~ as.numeric(levels(.x))[.x])

# Adjust factors
ps.aln.fact <- ps.aln.adj %>%
  select(c(1:(ncol(MAP)+12)), ncol(ps.aln.adj)) # include query type (last column)

ps.aln.fixed <- bind_cols(ps.aln.fact, ps.aln.value)
ps.aln.fixed <- data.table(ps.aln.fixed)

# Add UniProt target descriptions
descriptions <- fread("/mnt/data1/databases/virus_uniprot/uniprot.description.db", stringsAsFactors = TRUE, col.names = c("target", "protein"))
ps.aln.fixed <- merge(ps.aln.fixed, descriptions, all.x = TRUE, by = "target")
ps.aln.fixed <- ps.aln.fixed %>%
        select(-target, everything())
rm(descriptions); gc()

# Add Baltimore virus classifications
baltimore <- fread("/mnt/data1/databases/hecatomb/2020_03_11_family_Baltimore_classification_lookup.txt", stringsAsFactors = TRUE)
ps.aln.fixed <- left_join(ps.aln.fixed, baltimore, by = "Family")
ps.aln.fixed$Family <- as.factor(ps.aln.fixed$Family) # left_join may have converted to character
ps.aln.fixed <- data.table(ps.aln.fixed)
rm(baltimore); gc()

# Merge with library size
ps.aln.fixed <- merge(ps.aln.fixed, library.size, all.x = TRUE, by = "sample_ID")

ps.aln.scaled <- ps.aln.fixed %>%
        mutate(proportional_abundance = Abundance/library_size) %>%
        mutate(scaled_abundance_min = (min.ls*proportional_abundance)) %>%
        mutate(scaled_abundance_mean = (mean.ls*proportional_abundance)) %>%
        mutate(scaled_abundance_median = (median.ls*proportional_abundance))

# Adjust nt alignment length
ps.aln.scaled <- ps.aln.scaled %>%
        mutate(alignment_length_adjusted = ifelse(query_type == "aa", alignment_length*3, alignment_length*1))

#### USER INPUT ####
# Write full alignment table
# This is the full data set (every sequence) and is likely a very large file
saveRDS(ps.aln.scaled, file = "./results/aln_euk_viruses.RDS")

# Write zero abundance filtered file
ps.aln.fixed.filt <- ps.aln.scaled %>%
        filter(Abundance > 0)
saveRDS(ps.aln.fixed.filt, file = "./results/aln_euk_viruses_filt.RDS")

