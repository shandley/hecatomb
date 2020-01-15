#!/usr/bin/Rscript

# Load libraries
library("tidyverse")
library("data.table")
library("phyloseq")
library("speedyseq")

# Import mapping file
MAP <- import_qiime_sample_data("./results/mapping_rotabiome_virome.txt")

# Import in viral taxonomy table
viral_table <- read_tsv(file = "./results/viruses_tax_table.tsv", col_names = TRUE, col_types = "ifffffff")

# Remove unfiltered bacerial sequences
viral_table <- viral_table %>%
  filter(Kingdom != "Bacteria")

# Import and adjust count table
seqtable <- fread(file = "./results/seqtable.all", header = TRUE, sep = "\t", drop = 2)
seqtable <- setnafill(seqtable, fill=0)

# Merge counts table with tax_table
stmerge <- left_join(viral_table, seqtable, by = "id")

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
saveRDS(ps0, file = "./results/rotabiome_ps0.RDS")

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
vir.aln <- left_join(viral_table, aln.all, by = "id")

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
  select(c(1:(ncol(MAP)+12)), ncol(ps.aln.adj)) # include query tupe (last column)

ps.aln.fixed <- bind_cols(ps.aln.fact, ps.aln.value)

# Write alignment tables as PhyloSeq object
saveRDS(ps.aln.fixed, file = "./results/ps.aln.RDS")


