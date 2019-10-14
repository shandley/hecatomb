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

## Fix phage results
sed 's/uc_//g' ./results/mmseqs_aa_out/phage_table.tsv > ./results/phage_tax_table.tsv
sed -i '1iid\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' ./results/phage_tax_table.tsv

# Concatentate alignment files

