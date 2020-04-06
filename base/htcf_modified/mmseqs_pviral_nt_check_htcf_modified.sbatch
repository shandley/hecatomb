#!/bin/bash

# mmseqs_pviral_nt_check_htcf_modified.sbatch

#SBATCH --constraint=cpu_E52650
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G
#SBATCH -J ck_pvrl_nt
#SBATCH --mail-user=rachel.rodgers
#SBATCH --mail-type=ALL

set -ue

module load mmseqs2/10-6d92c
module load pullseq
module load R

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
# This is a targetDB consisting of all UniProtKB entires clustered at 30% ID (UniClust30) concatenated to Virus UniProt entries clustered at $
# More information about the UniClust db's as well as download links to UniClust30 and UniClust90 are available at: https://uniclust.mmseqs.c$
DB=/scratch/ref/hecatomb_databases/bac_virus_masked/nt.fnaDB

# Set phage lineage file path
PHAGE=/scratch/ref/hecatomb_databases/phage_taxonomic_lineages.txt

# Create output directory
OUT=./results/mmseqs_nt_checked_out

## Adjust taxonomy table and extract viral lineages
# Extract phage lineages
echo
echo "Extract phage lineages"
echo
tail -n+2 $OUT/mmseqs_pviral_nt_lineage.tsv | \
	grep -f $PHAGE | \
	sort -n -k1 > $OUT/phage_nt_table.tsv

cut -f1 $OUT/phage_nt_table.tsv > $OUT/phage_nt_table.list

pullseq -i ./results/seqtable.fasta -n $OUT/phage_nt_table.list -l 5000 > $OUT/phage_nt_seqs.fasta

# Extract non-phage viral lineages
echo
echo "Extract non-phage viral lineages"
echo
tail -n+2 $OUT/mmseqs_pviral_nt_lineage.tsv | \
        grep -v -f $PHAGE | \
        sort -n -k1 > $OUT/pviral_virus_nt_table.tsv

cut -f1	$OUT/pviral_virus_nt_table.tsv > $OUT/pviral_virus_nt_table.list

pullseq	-i ./results/seqtable.fasta -n $OUT/pviral_virus_nt_table.list -l 5000 > $OUT/pviral_virus_nt_seqs.fasta

# Create Query databases
echo
echo "Create query databases"
echo
mmseqs createdb $OUT/pviral_virus_nt_seqs.fasta $OUT/seqtable_queryDB --dbtype 2

# mmseqs search
echo
echo "MMSeqs Search"
echo
mmseqs search \
    $OUT/seqtable_queryDB \
    $DB \
    $OUT/resultDB \
    $OUT/tmp_nt_check \
    -a -e 0.000001 --search-type 3 --cov-mode 2 -c 0.95 \
    --threads ${SLURM_CPUS_PER_TASK}

# extract top-hit
echo
echo "MMSeqs filterdb"
echo
mmseqs filterdb \
    $OUT/resultDB $OUT/resultDB.firsthit \
    --extract-lines 1 \
    --threads ${SLURM_CPUS_PER_TASK}

echo
echo "MMSeqs convertalis"
echo
mmseqs convertalis \
    $OUT/seqtable_queryDB \
    $DB \
    $OUT/resultDB.firsthit \
    $OUT/resultDB.firsthit.m8 \
    --threads ${SLURM_CPUS_PER_TASK}

# Annotate
Rscript mmseqs_pviral_nt_check_annotate_htcf_modified.R
