#!/bin/bash

ls -1 *.fastq.gz > filenames
ls -l *.fastq.gz | awk -F ' ' '{ print$5 }' > base.fs
ls -l ./clumped/*.fastq.clumped.gz | awk -F ' ' '{ print$5 }' > clumped.fs

paste -d '\t' filenames base.fs clumped.fs > clumping.stats
awk '$4=100*$3/$2' clumping.stats > tmp && mv tmp clumping.stats
sed -i '1i sample_id raw clumped perc_shrunk' clumping.stats
rm filenames base.fs clumped.fs

cut -f4 -d ' ' clumping.stats | tail -n +2 | spark

