#!/bin/bash
# Script to deduplicate, count and assign taxonomy to sequences assigned as BACTERIAL
# References:
# Heavy reliance on:
        # BBtools: https://jgi.doe.gov/data-and-tools/bbtools/
	# mmseqs2: https://github.com/soedinglab/MMseqs2
        # pullseq: https://github.com/bcthomas/pullseq
        # SeqKit: https://bioinf.shenwei.me/seqkit/

# Prep output directory
mkdir -p ./results/bacterial_check
IN=./QC/step_9
OUT=./results/bacterial_check

# Set file names
for i in $IN/*_bacterial.fastq; do
        F=`basename $i _bacterial.fastq`;

# Summary:
	# Step 1: Remove exact duplicates (consider to be PCR artifacts)
	# Step 2: Deduplicate
	# Step 3: Reformat and create sequence table

# Remove exact duplicates
dedupe.sh in=$IN/"$F"_bacterial.fastq ow=t \
	out=$OUT/"$F"_bacterial.deduped.out.fastq \
	ac=f \
	-Xmx64g;

# Deduplication 
dedupe.sh in=$OUT/"$F"_bacterial.deduped.out.fastq ow=t \
	s=4 \
	rnc=t pbr=t \
	csf=$OUT/"$F"_stats.txt out=$OUT/"$F"_best.fasta \
	-Xmx64g;

# Convert to fasta
reformat.sh in=$OUT/"$F"_best.fasta out=$OUT/"$F"_reformated.fasta \
	deleteinput=t fastawrap=0 \
	ow=t \
	-Xmx64g;

# Parse and combine stats and contig files
# Extract sequences
grep -v '>' $OUT/"$F"_reformated.fasta | sed '1i sequence' > $OUT/"$F"_seqs.txt;

# Extract sequence IDs
# grep '>' $OUT/"$F"_reformated.fasta | sed 's|>Cluster_||' | awk -F "," '{ print$1 }' | sort -n | sed '1i contig_ids' > $OUT/"$F"_contig_ids.txt;

# Extract counts
cut -f 2 $OUT/"$F"_stats.txt | sed "1s/size/$F/" > $OUT/"$F"_counts.txt;

# Create sequence table
paste $OUT/"$F"_seqs.txt $OUT/"$F"_counts.txt > $OUT/"$F"_seqtable.txt;

done

# Remove zero byte (no reads) files
find $OUT -type f -size 0 | sed '1i##These files contained zero reads following derreplication\n' > $OUT/zero_read_files.txt
find $OUT -type f -size 0 -exec rm -f {} \;

# Create sequence table of clutered results
check_bacteria_seqtable.R

# search bacterial sequences against viral database
# REQUIRES that targetDB has already been indexed
# If it has not been index then run the following script in the directory of your choice: uniprot_viral_DB_build.sh (found in /accessory)
# Note: mmseqs2 taxonomy is currently most useful if you have UniProt formatted fasta databases
# more details about database building can be found at: https://github.com/soedinglab/mmseqs2/wiki#taxonomy-assignment
-using-mmseqs-taxonomy

# Set targetDB (UniProt Viral Proteins Cluster at 99% identity)
DB=/mnt/data1/databases/uniclust/uniclust30/uni_plus_virus/targetDB

# Convert seqtable.tab2fx to fasta
seqkit tab2fx $OUT/results_bacterial_check/seqtable.tab2fx -w 5000 -o $OUT/results_bacterial_check/seqtable.fasta;

# Create Query databases
mmseqs createdb $OUT/results_bacterial_check/seqtable.fasta $OUT/results_bacterial_check/seqtable_queryDB --dont-shuffle 0 --dbtype 0;

## mmseqs2 taxonomy search
mmseqs taxonomy $OUT/results_bacterial_check/seqtable_queryDB $DB $OUT/results_bacterial_check/taxonomyResult $OUT/results_bacterial_check/tmp_aa \
	-a -s 7 --search-type 2 --tax-output-mode 1;

mmseqs convertalis $OUT/results_bacterial_check/seqtable_queryDB $DB $OUT/results_bacterial_check/taxonomyResult $OUT/results_bacterial_check/aln.m8 \
	--format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln";

mmseqs lca $DB $OUT/results_bacterial_check/taxonomyResult $OUT/results_bacterial_check/lcaDB --tax-lineage true \
	--lca-ranks "superkingdom:phylum:class:order:family:genus:species";

# create taxonomy table (tsv)
mmseqs createtsv $OUT/results_bacterial_check/seqtable_queryDB $OUT/results_bacterial_check/lcaDB $OUT/results_bacterial_check/taxonomyResult.tsv;

# create kraken-style report
mmseqs taxonomyreport $DB $OUT/results_bacterial_check/lcaDB $OUT/results_bacterial_check/taxonomyResult.report;
