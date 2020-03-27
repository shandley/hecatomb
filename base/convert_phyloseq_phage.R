#!/usr/bin/Rscript

#### Instructions ###
# There are several places that require user input. Each section is preceeded with a #### USER INPUT #### tag
# 1) Need to supply mapping file location for MAP
# 2) Need to change directroy and filenames to write (two saveRDS commands below)
# 3) Joining the library size information requires a common key (column name) and should be set in the Standardization section 

# Load libraries
library("tidyverse")
library("data.table")
library("phyloseq")
library("speedyseq")

#### USER INPUT ####
# Import mapping file
MAP <- import_qiime_sample_data("./results/mapping_rotabiome_virome.txt")

# Import viral taxonomy table
viral_table <- fread(file = "./results/phage_tax_table.tsv", header = TRUE, stringsAsFactors=TRUE)

# Remove any unfiltered bacerial sequences
setkey(viral_table, Kingdom)
viral_table <- viral_table[!"Bacteria"]

# Import and adjust count table
seqtable <- fread(file = "./results/seqtable.all", header = TRUE, sep = "\t", drop = 2)
seqtable <- setnafill(seqtable, fill=0)

# Merge count table with viral tax_table
setkey(viral_table, id)
setkey(seqtable, id)
stmerge <- merge(viral_table, seqtable, all.x = TRUE)

# Sum per taxon counts
merged_table <- stmerge %>%
  select(-id) %>%
  unite("lineage", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "_", remove = FALSE) %>%
  group_by(lineage, Species) %>%
  summarise_if(is.numeric, funs(sum(as.numeric(.)))) %>%
  separate("lineage", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus"), sep = "_") %>%
  ungroup()

# Remove taxon information to create final count table (called otu_table in phyloseq)
count_table <- merged_table[,-c(1:8)]

# Convert summarized table to a matrix of just the taxon in preperation for phyloseq import
tax.tbl.m <- merged_table %>%
  select("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %>%
  as.matrix()

# Create Phyloseq Object
OTU = otu_table(count_table, taxa_are_rows = TRUE)
TAX = tax_table(tax.tbl.m)
ps0 = phyloseq(OTU, TAX, MAP)

#### USER INPUT ####
# Save PhyloSeq object as an RDS file
saveRDS(ps0, file = "./results/ps0_phage.RDS")

# Clean up environment
rm(count_table, merged_table, OTU, ps0, seqtable.all, stmerge, TAX, tax.tbl.m); gc()

##### Alignment Files #####
# aa_checked alignment
aln.aa <- fread(file = "./results/aa.aln.m8", header = TRUE, sep = "\t")
aln.aa <- aln.aa %>%
  rename(id = query)

# nt_checked alignment
aln.nt <- fread(file = "./results/nt.aln.m8", header = TRUE, sep = "\t")
aln.nt <- aln.nt %>%
  rename(id = query)

# bind aa and nt alignment
aln.aa <- aln.aa %>%
  mutate(query_type = "aa")
aln.nt <- aln.nt %>%
  mutate(query_type = "nt")

aln.all <- bind_rows(aln.aa, aln.nt)
aln.all <- data.table(aln.all)

# Merge viral taxonomy and aa/nt alignment table
setkey(viral_table, id)
setkey(aln.all, id)
vir.aln <- merge(viral_table, aln.all, all.x = TRUE)

# Merge virus + alignment table with count table
setkey(seqtable, id)
aln.merge <- merge(vir.aln, seqtable, all.x = TRUE)

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

# Convert to data table. This may take a bit of time as well
ps.aln.fixed <- data.table(ps.aln.fixed)

# Free up some memory
rm(aln.aa, aln.all, aln.count, aln.merge, aln.nt, aln.tax.tbl, aln.tax.tbl.m, OTU.ALN, ps.aln, ps.aln.adj, ps.aln.fact, ps.aln.value, TAX.ALN, vir.aln, viral_table); gc()

# Add UniProt target descriptions
descriptions <- fread("/mnt/data1/databases/virus_uniprot/uniprot.description.db", stringsAsFactors = TRUE, col.names = c("target", "protein"))
setkey(ps.aln.fixed, target)
setkey(descriptions, target)
ps.aln.fixed <- merge(ps.aln.fixed, descriptions, all.x = TRUE)
rm(descriptions); gc()

# Add Baltimore virus classifications
baltimore <- fread("/mnt/data1/databases/hecatomb/2020_03_11_family_Baltimore_classification_lookup.txt", stringsAsFactors = TRUE)
setkey(ps.aln.fixed, Family)
setkey(baltimore, Family)
ps.aln.fixed <- merge(ps.aln.fixed, baltimore, all.x = TRUE)
rm(baltimore); gc()

# Standardize to per sample library size
col.sums <- as_tibble(colSums(seqtable[,-1]))
samples <- as_tibble(colnames(seqtable[,-1]))
library.size <- bind_cols(samples, col.sums)
rm(col.sums, samples); gc()

#### USER INPUT ####
library.size <- library.size %>%
  rename(seq_id = value) %>% # Rename to whatever your sequence file name / sample name is in your mapping file
  rename(library_size = value1)
library.size <- data.table(library.size)

setkey(ps.aln.fixed, seq_id)
setkey(library.size, seq_id)
ps.aln.fixed <- ps.aln.fixed[library.size, nomatch = 0]

# Add proportion, min and mean and median standardizations
min.ls <- min(library.size$library_size)
mean.ls <- mean(library.size$library_size)
median.ls <- median(library.size$library_size)

ps.aln.scaled <- ps.aln.fixed %>%
	mutate(proportional_abundance = Abundance/library_size) %>%
	mutate(scaled_abundance_min = (min.ls*proportional_abundance)) %>%
	mutate(scaled_abundance_mean = (mean.ls*proportional_abundance)) %>%
	mutate(scaled_abundance_median = (median.ls*proportional_abundance))

# Adjust nt alignment length
ps.aln.scaled <- ps.aln.scaled %>%
	mutate(alignment_length_adjusted = ifelse(query_type == "aa", alignment_length*3, alignment_length*1))

#### USER INPUT ####
# Write alignment table
saveRDS(ps.aln.scaled, file = "./results/aln_phage.RDS")

# Write zero abundance filtered file
ps.aln.fixed.filt <- ps.aln.scaled %>%
	filter(Abundance > 0)
saveRDS(ps.aln.fixed.filt, file = "./results/aln_phage_filt.RDS")
