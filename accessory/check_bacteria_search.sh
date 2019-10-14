#!/bin/bash
# Script to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

# References:
# Heavy reliance on:
        # mmseqs2: https://github.com/soedinglab/MMseqs2
        # pullseq: https://github.com/bcthomas/pullseq
	# SeqKit: https://bioinf.shenwei.me/seqkit/

# REQUIRES that targetDB has already been indexed
# If it has not been index then run the following script in the directory of your choice: uniprot_viral_DB_build.sh (found in /accessory)
# Note: mmseqs2 taxonomy is currently most useful if you have UniProt formatted fasta databases
# more details about database building can be found at: https://github.com/soedinglab/mmseqs2/wiki#taxonomy-assignment-using-mmseqs-taxonomy

# Set targetDB (UniProt Viral Proteins Cluster at 99% identity)
DB=/mnt/data1/databases/uniclust/uniclust30/uni_plus_virus/targetDB

# Create output directory
OUT=./results/bacterial_check/results

# Convert seqtable.tab2fx to fasta
seqkit tab2fx $OUT/seqtable.tab2fx -w 5000 -o $OUT/seqtable.fasta;

# Create Query databases
mmseqs createdb $OUT/seqtable.fasta $OUT/seqtable_queryDB --dont-shuffle 0 --dbtype 0;

## mmseqs2 taxonomy search
mmseqs taxonomy $OUT/seqtable_queryDB $DB $OUT/taxonomyResult $OUT/tmp_aa \
	-a -s 7 --search-type 2 --tax-output-mode 1;

mmseqs convertalis $OUT/seqtable_queryDB $DB $OUT/taxonomyResult $OUT/aln.m8 \
	--format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln";

mmseqs lca $DB $OUT/taxonomyResult $OUT/lcaDB --tax-lineage true \
	--lca-ranks "superkingdom:phylum:class:order:family:genus:species";

# create taxonomy table (tsv)
mmseqs createtsv $OUT/seqtable_queryDB $OUT/lcaDB $OUT/taxonomyResult.tsv;

# create kraken-style report
mmseqs taxonomyreport $DB $OUT/lcaDB $OUT/taxonomyResult.report;


