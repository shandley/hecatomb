#!/bin/bash

## Concatenate and fix eukaryotic virus results
# Correct aa annotation so it jives with nt annotation
sed 's/uc_//g' ./results/mmseqs_aa_checked_out/viruses_checked_aa_table.tsv > ./results/mmseqs_aa_checked_out/viruses_checked_aa_table_edited.tsv

# Correct nt annotation so it jives with aa annotation
tail -n +2 ./results/mmseqs_nt_checked_out/mmseqs_pviral_nt_checked_lineage.tsv | \
	sed 's/NA/unknown/g' > ./results/mmseqs_nt_checked_out/mmseqs_pviral_nt_checked_lineage_edited.tsv

# Happily marry the corrected aa and nt files
cat ./results/mmseqs_aa_checked_out/viruses_checked_aa_table_edited.tsv ./results/mmseqs_nt_checked_out/mmseqs_pviral_nt_checked_lineage_edited.tsv | \
	sort -n -k 1 > ./results/viruses_tax_table.tsv

sed -i '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' ./results/viruses_tax_table.tsv

## Fix some strange taxonomy names
sed 's/uc_//g' ./results/mmseqs_aa_out/phage_table.tsv > ./results/phage_tax_table.tsv
sed -i '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' ./results/phage_tax_table.tsv
cut -f1-8 ./results/viruses_tax_table.tsv > tmp & mv tmp ./results/viruses_tax_table.tsv

# Adjust alignment files

sed -i '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' ./results/mmseqs_aa_checked_out/taxonomyResult.firsthit.m8;
sed -i '1iquery\ttarget\tpercent_id\talignment_length\tnum_mismatches\tnumber_gaps\tstart_query\tend_query\tstart_target\tend_target\te_value\tbit_score' ./results/mmseqs_nt_checked_out/resultDB.firsthit.m8;

cp ./results/mmseqs_aa_checked_out/taxonomyResult.firsthit.m8 ./results/aa.aln.m8;
cp ./results/mmseqs_nt_checked_out/resultDB.firsthit.m8 ./results/nt.aln.m8;

# Create fasta entires for every family in ./results/viruses_tax_table.tsv

# Create directory
mkdir -p ./results/family_reads
OUT=./results/family_reads

# Create list of families present in viruses_tax_table.sh
tail -n +2 ./results/viruses_tax_table.tsv | cut -f6 | sort | uniq > ./results/family.list;

# Create a temporary file list of each family to make for-do-looping simple
cd $OUT;
xargs touch <../family.list;
for f in *; do
	mv "$f" $(echo "$f".id);
done;

# Pull the id for each read within each family
for i in *.id; do
        f=`basename $i .id`;
        grep $f ../viruses_tax_table.tsv > "$f".reads;
done;

# Use pullseq to pull the fasta records for the *.reads ids
for i in *.reads; do
        f=`basename $i .reads`;
        pullseq -i ../seqtable.fasta -n "$f".reads > "$f".fasta;
done;

rm *.id;
rm *.reads;

cd ../../;
