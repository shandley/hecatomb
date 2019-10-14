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
# This is a targetDB consisting of all UniProtKB entires clustered at 30% ID (UniClust30) concatenated to Virus UniProt entries clustered at 99%
# More information about the UniClust db's as well as download links to UniClust30 and UniClust90 are available at: https://uniclust.mmseqs.com
DB=/mnt/data1/databases/uniclust/uniclust30/uni_plus_virus/targetDB

# Set phage lineage file path
PHAGE=~/virome_analysis/scripts/base/phage_taxonomic_lineages.txt

# Create output directory
mkdir -p ./results/mmseqs_aa_checked_out;
OUT=./results/mmseqs_aa_checked_out;

# Create Query databases
mmseqs createdb ./results/mmseqs_aa_out/viruses_seqs.fasta $OUT/viral_seqs_queryDB --dont-shuffle 0 --dbtype 0;

## mmseqs2 taxonomy search
mmseqs taxonomy $OUT/viral_seqs_queryDB $DB $OUT/taxonomyResult $OUT/tmp_aa_checked \
        -a -s 7 --search-type 2 --tax-output-mode 1;

mmseqs convertalis $OUT/viral_seqs_queryDB $DB $OUT/taxonomyResult $OUT/aln.m8 \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln";

mmseqs lca $DB $OUT/taxonomyResult $OUT/lcaDB --tax-lineage true \
        --lca-ranks "superkingdom:phylum:class:order:family:genus:species";

# extract top-hit alignment db
mmseqs filterdb $OUT/taxonomyResult $OUT/taxonomyResult.firsthit --extract-lines 1;
mmseqs convertalis $OUT/viral_seqs_queryDB $DB $OUT/taxonomyResult.firsthit $OUT/taxonomyResult.firsthit.m8;

# create taxonomy table (tsv)
mmseqs createtsv $OUT/viral_seqs_queryDB $OUT/lcaDB $OUT/taxonomyResult.tsv;

# create kraken-style report
mmseqs taxonomyreport $DB $OUT/lcaDB $OUT/taxonomyResult.report;

# Extract non-phage viral lineages and generate taxonomy table for import into R as PhyloSeq object
grep -v 'Bacteria:' $OUT/taxonomyResult.tsv | \
	grep 'Viruses:' | \
	grep -v -f $PHAGE | cut -f1,5 | \
	sed 's/:/\t/g' | \
        sort -n -k1 > $OUT/viruses_checked_aa_table.tsv;
cut -f1 $OUT/viruses_checked_aa_table.tsv > $OUT/viruses_checked_aa_seqs.list;
pullseq -i ./results/seqtable.fasta -n $OUT/viruses_checked_aa_seqs.list -l 5000 > $OUT/viruses_checked_aa_seqs.fasta;

seqkit fx2tab $OUT/viruses_checked_aa_seqs.fasta > $OUT/viruses_checked_aa_seqs.fx2tab;
join $OUT/viruses_checked_aa_seqs.fx2tab $OUT/viruses_checked_aa_table.tsv | awk -F ' ' '{ print $2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7,"\t",$8,"\t",$9 }' | \
        sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > $OUT/viruses_checked_aa_tax_table.tsv;

# Extract non-viral lineages
grep -v 'Viruses:' $OUT/taxonomyResult.tsv | cut -f1,5 | \
        sed 's/:/\t/g' | \
        sort -n -k1 > $OUT/unclassified_checked_aa_seqs.list;
pullseq -i ./results/seqtable.fasta -n $OUT/unclassified_checked_aa_seqs.list -l 5000 > $OUT/unclassified_checked_aa_seqs.fasta;
