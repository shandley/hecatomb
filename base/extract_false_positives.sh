#!/bin/bash
# Script to extract common false-positive viral contaminant sequences for subsequent analysis

##### !!!!! Need to update directory paths for final script. Currently in a testing form

# Create output directories
mkdir -p ./common_contaminant_seqs/
OUT=./common_contaminant_seqs

# Extract Mimiviridae (family level)
grep 'Mimiviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Mimiviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Mimiviridae_seqs.list -l 5000 > $OUT/Mimiviridae_seqs.fasta;

# Extract Poxviridae (family level)
grep 'Poxviridae' viruses_checked_aa_table.tsv | \
	cut -f1 | \
       	sort -n -k1 > $OUT/Poxviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Poxviridae_seqs.list -l 5000 > $OUT/Poxviridae_seqs.fasta;

# Extract Marseilleviridae (family level)
grep 'Marseilleviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Marseilleviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Marseilleviridae_seqs.list -l 5000 > $OUT/Marseilleviridae_seqs.fasta;

# Extract Baculoviridae (family level)
grep 'Baculoviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Baculoviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Baculoviridae_seqs.list -l 5000 > $OUT/Baculoviridae_seqs.fasta;

# Extract Iridoviridae (family level)
grep 'Iridoviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Iridoviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Iridoviridae_seqs.list -l 5000 > $OUT/Iridoviridae_seqs.fasta;

# Extract Asfariviridae (family level)
grep 'Asfariviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Asfariviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Asfariviridae_seqs.list -l 5000 > $OUT/Asfariviridae_seqs.fasta;

# Extract Ascoviridae (family level)
grep 'Ascoviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Ascoviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Ascoviridae_seqs.list -l 5000 > $OUT/Ascoviridae_seqs.fasta;

# Extract Ascoviridae (family level)
grep 'Ascoviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Ascoviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Ascoviridae_seqs.list -l 5000 > $OUT/Ascoviridae_seqs.fasta;

# Extract Megaviridae (family level)
grep 'Megaviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Megaviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Megaviridae_seqs.list -l 5000 > $OUT/Megaviridae_seqs.fasta;

# Extract Pandoravirus (genus level)
grep 'Pandoravirus' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Pandoravirus_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Pandoravirus_seqs.list -l 5000 > $OUT/Pandoravirus_seqs.fasta;

# Extract Phycodnaviridae (family level)
grep 'Phycodnaviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Phycodnaviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Phycodnaviridae_seqs.list -l 5000 > $OUT/Phycodnaviridae_seqs.fasta;

# Extract Pithoviridae (family level)
grep 'Pithoviridae' viruses_checked_aa_table.tsv | \
        cut -f1 | \
        sort -n -k1 > $OUT/Pithoviridae_seqs.list;
pullseq -i ../seqtable.fasta -n $OUT/Pithoviridae_seqs.list -l 5000 > $OUT/Pithoviridae_seqs.fasta;

# Crate contaminant fasta file
cat $OUT/*.fasta > $OUT/common_contaminants.fasta;

# Create Summary Table
grep -c ">" $OUT/*.fasta | sed 's/\:/\t/' > $OUT/common_contaminant.stats;

# Clean up
rm $OUT/*.list;
