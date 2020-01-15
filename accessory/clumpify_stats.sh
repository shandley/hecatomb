#!/bin/bash
# Calculate file size stats pre and post clumping


# Extract file size information
ls -1 *.fastq.gz > filenames
ls -l *.fastq.gz | awk -F ' ' '{ print$5 }' > base.fs
ls -l ./clumped/*.fastq.clumped.gz | awk -F ' ' '{ print$5 }' > clumped.fs

# Create combined file stats table with % reduction
paste -d '\t' filenames base.fs clumped.fs > clumping.stats
awk '$4=100*(1-($3/$2))' clumping.stats > tmp && mv tmp clumping.stats
sed -i '1i sample_id raw clumped perc_shrunk' clumping.stats
rm filenames base.fs clumped.fs

# Ouput spark plot to the command line
cut -f4 -d ' ' clumping.stats | tail -n +2 | spark

