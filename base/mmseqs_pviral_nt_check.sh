#!/bin/bash
# Query probable viral hits vs. UniClust30 proteinDB to remove false-positives
        # Uniclust: https://uniclust.mmseqs.com

# References:
# Heavy reliance on:
        # mmseqs2: https://github.com/soedinglab/MMseqs2
        # pullseq: https://github.com/bcthomas/pullseq
        # SeqKit: https://bioinf.shenwei.me/seqkit/

# REQUIRES that targetDB has already been indexed
# If it has not been index then run the following script in the directory of your choice: uniprot_viral_DB_build.sh (found in /accessory)
# Note: mmseqs2 taxonomy is currently most useful if you have UniProt formatted fasta databases
# more details about database building can be found at: https://github.com/soedinglab/mmseqs2/wiki#taxonomy-assignment-using-mmseqs-taxonomy

# Set targetDB
# This is a targetDB consisting of all UniProtKB entires clustered at 30% ID (UniClust30) concatenated to Virus UniProt entries clustered at $
# More information about the UniClust db's as well as download links to UniClust30 and UniClust90 are available at: https://uniclust.mmseqs.c$
DB=/mnt/data1/databases/refseq_viral_genomes/refseq_virus_nn_nt/masked/bac_virus_masked/nt.fnaDB

# Set phage lineage file path
PHAGE=~/virome_analysis/scripts/base/phage_taxonomic_lineages.txt

# Create output directory
OUT=./results/mmseqs_nt_checked_out;

## Adjust taxonomy table and extract viral lineages
# Extract phage lineages
tail -n+2 $OUT/mmseqs_pviral_nt_lineage.tsv | \
	grep -f $PHAGE | \
	sort -n -k1 > $OUT/phage_nt_table.tsv;
cut -f1 $OUT/phage_nt_table.tsv > $OUT/phage_nt_table.list;
pullseq -i ./results/seqtable.fasta -n $OUT/phage_nt_table.list -l 5000 > $OUT/phage_nt_seqs.fasta;

# Extract non-phage viral lineages
tail -n+2 $OUT/mmseqs_pviral_nt_lineage.tsv | \
        grep -v -f $PHAGE | \
        sort -n -k1 > $OUT/pviral_virus_nt_table.tsv;
cut -f1	$OUT/pviral_virus_nt_table.tsv > $OUT/pviral_virus_nt_table.list;
pullseq	-i ./results/seqtable.fasta -n $OUT/pviral_virus_nt_table.list -l 5000 > $OUT/pviral_virus_nt_seqs.fasta;

# Create Query databases
mmseqs createdb $OUT/pviral_virus_nt_seqs.fasta $OUT/seqtable_queryDB --dbtype 2;

# mmseqs search
mmseqs search $OUT/seqtable_queryDB $DB $OUT/resultDB $OUT/tmp_nt_check -a -e 0.000001 --search-type 3 --cov-mode 2 -c 0.95;

# extract top-hit
mmseqs filterdb $OUT/resultDB $OUT/resultDB.firsthit --extract-lines 1;
mmseqs convertalis $OUT/seqtable_queryDB $DB $OUT/resultDB.firsthit $OUT/resultDB.firsthit.m8;

# Annotate
mmseqs_pviral_nt_check_annotate.R
