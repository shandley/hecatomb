These scripts are here for historical reasons and are unlikely to be useful to most users.
By default, they are not installed with the conda build for Hecatomb.


## primer_stats.sh
A bash script that compiles statistics about the primers detected during step_01 of contaminant_removal.smk.

The script automatically calls primer_stats.R to do the final multi-join. 
Of note, you will need [tidyverse](https://www.tidyverse.org) installed for this part to work.

The script will output two *.tsv files useful for determining what primers were detected in each sample.
