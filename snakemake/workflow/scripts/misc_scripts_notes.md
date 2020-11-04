# Miscellaneous scripts

This directory contains optional scripts useful for pre or post-processing data produced by the Snakemake workflow or managing reference databases.


## primer_stats.sh
A bash script that compiles statistics about the primers detected during step_01 of contaminant_removal.smk.

The script automatically calls primer_stats.R to do the final multi-join. Of note, you will need [tidyverse](https://www.tidyverse.org) installed for this part to work.

The script will output two *.tsv files useful for determining what primers were detected in each sample.
