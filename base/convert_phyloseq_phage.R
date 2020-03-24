#!/usr/bin/Rscript

# Load libraries
library("tidyverse")
library("data.table")
library("phyloseq")
library("speedyseq")

# Import mapping file
MAP <- import_qiime_sample_data("./results/map_phage_attack_virome.txt")

# Import in viral taxonomy table
phage_table <- read_tsv(file = "./results/phage_tax_table.tsv", col_names = TRUE, col_types = "ifffffff")

# Remove unfiltered bacerial sequences
phage_table <- phage_table %>%
  filter(Kingdom != "Bacteria")

# Import and adjust count table
seqtable <- fread(file = "./results/seqtable.all", header = TRUE, sep = "\t", drop = 2)
seqtable <- setnafill(seqtable, fill=0)

# Merge counts table with tax_table
stmerge <- left_join(phage_table, seqtable, by = "id")

# Sumamrise per taxon counts
merged_table <- stmerge %>%
  select(-id) %>%
  unite("lineage", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "_") %>%
  group_by(lineage) %>%
  summarise_if(is.numeric, funs(sum(as.numeric(.)))) %>%
  separate("lineage", c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"), sep = "_") %>%
  ungroup()

# Remove taxon information to create final counts table (called otu_table in phyloseq)
count_table <- merged_table[,-c(1:8)]

# Convert summarized table to a matrix of just the taxon
tax.tbl.m <- merged_table %>%
  select("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species") %>%
  as.matrix()

# Create Phyloseq Object
OTU = otu_table(count_table, taxa_are_rows = TRUE)
TAX = tax_table(tax.tbl.m)
ps0 = phyloseq(OTU, TAX, MAP)

# Write PhyloSeq object as an RDS
saveRDS(ps0, file = "./results/ps0_phage.RDS")

##### Alignment Files #####

# aa_checked alignment
aln.aa <- read_delim("./results/aa.aln.m8", col_names = TRUE, delim = "\t")
aln.aa <- aln.aa %>%
  rename(id = query)

# nt_checked alignment
aln.nt <- read_delim("./results/nt.aln.m8", col_names = TRUE, delim = "\t")
aln.nt <- aln.nt %>%
  rename(id = query)

# bind aa and nt alignment
aln.aa <- aln.aa %>%
  mutate(query_type = "aa")
aln.nt <- aln.nt %>%
  mutate(query_type = "nt")

aln.all <- bind_rows(aln.aa, aln.nt)

# Merge viral taxonomy with aa alignment table
vir.aln <- left_join(phage_table, aln.all, by = "id")

# Merge virus + alignment table with seqtable
aln.merge <- left_join(vir.aln, seqtable, by = "id")

# Extract count table
aln.count <- aln.merge[,-c(1:20)]
OTU.ALN <- otu_table(aln.count, taxa_are_rows = TRUE)

# Create tax table
aln.tax.tbl <- aln.merge[,c(1:20)] %>%
  as_tibble()

aln.tax.tbl.m <- aln.merge[,c(1:20)] %>%
  as.matrix()

TAX.ALN <- tax_table(aln.tax.tbl.m)

# Create Phyloseq Object of alignment results
ps.aln <- phyloseq(OTU.ALN, TAX.ALN, MAP)

# Adjust PhyloSeq alignment object
ps.aln.adj <- ps.aln %>%
  speedyseq::psmelt() %>%
  as_tibble()

ps.aln.value <- ps.aln.adj %>%
  select(c((ncol(MAP)+13):(ncol(ps.aln.adj)-1))) %>%
   mutate_if(is.factor, ~ as.numeric(levels(.x))[.x])

ps.aln.fact <- ps.aln.adj %>%
  select(c(1:(ncol(MAP)+12)), ncol(ps.aln.adj)) # include query type (last column)

ps.aln.fixed <- bind_cols(ps.aln.fact, ps.aln.value)

# Add UniProt definitions
descriptions <- fread("/mnt/data1/databases/virus_uniprot/uniprot.description.db", stringsAsFactors = TRUE, col.names = c("target", "protein"))
ps.aln.fixed <- left_join(ps.aln.fixed, descriptions, by = "target")
ps.aln.fixed$target <- as.factor(ps.aln.fixed$target) # left_join may have converted to character

# Add Baltimore classifications
baltimore <- fread("/mnt/data1/databases/hecatomb/2020_03_11_family_Baltimore_classification_lookup.txt", stringsAsFactors = TRUE)
ps.aln.fixed <- left_join(ps.aln.fixed, baltimore, by = "Family")
ps.aln.fixed$Family <- as.factor(ps.aln.fixed$Family) # left_join may have converted to character

# Standardize
col.sums <- as_tibble(colSums(seqtable[,-1]))
samples <- as_tibble(colnames(seqtable[,-1]))
library.size <- bind_cols(samples, col.sums)

library.size <- library.size %>%
  rename(SampleID = value) %>%
  rename(library_size = value1)

ps.aln.fixed <- inner_join(ps.aln.fixed, library.size, by = "SampleID")
ps.aln.fixed$SampleID <- as.factor(ps.aln.fixed$SampleID) # join may have converted to character

# Add proportion, min and mean and median standardizations
min.ls <- min(library.size$library_size)
mean.ls <- mean(library.size$library_size)
median.ls <- median(library.size$library_size)

ps.aln.scaled <- ps.aln.fixed %>%
	mutate(proportional_abundance = Abundance/library_size) %>%
	mutate(scaled_abundance_min = (min.ls*proportional_abundance)) %>%
	mutate(scaled_abundance_mean = (mean.ls*proportional_abundance)) %>%
	mutate(scaled_abundance_median = (median.ls*proportional_abundance))

# Write alignment tables as PhyloSeq object
saveRDS(ps.aln.fixed, file = "./results/ps_aln_phage.RDS")

