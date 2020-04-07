#!/bin/bash

# mmseqs_pviral_aa_htcf_modified.sbatch

#SBATCH --constraint=cpu_E52650
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G
#SBATCH -J pvrl_aa
#SBATCH --mail-user=rachel.rodgers
#SBATCH --mail-type=ALL

set -ue

module load seqkit
module load mmseqs2/10-6d92c
module load pullseq

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
DB=/scratch/ref/hecatomb_databases/virus_uniprot/targetDB

# Set phage lineage file path
PHAGE=/scratch/ref/hecatomb_databases/phage_taxonomic_lineages.txt

# Create output directory
mkdir -p ./results/mmseqs_aa_out/
OUT=./results/mmseqs_aa_out

#  Convert seqtable.tab2fx to fasta
echo
echo "CONVERTING SEQTABLE.TAB2FX TO FASTA"
echo
seqkit tab2fx ./results/seqtable.tab2fx -w 5000 -o ./results/seqtable.fasta

# Create Query databases
echo
echo "CREATE QUERY DB"
echo
mmseqs createdb ./results/seqtable.fasta $OUT/seqtable_queryDB --dont-shuffle 0 --dbtype 0

## mmseqs2 taxonomy search
echo
echo "MMSEQS TAXONOMY SEARCH"
echo
mmseqs taxonomy \
	$OUT/seqtable_queryDB \
	$DB \
	$OUT/taxonomyResult \
	$OUT/tmp_aa \
	--threads ${SLURM_CPUS_PER_TASK} \
	-a --start-sens 1 --sens-steps 3 -s 7 --search-type 2 --tax-output-mode 1

echo
echo "MMSEQS CONVERTALIS"
echo
mmseqs convertalis $OUT/seqtable_queryDB $DB $OUT/taxonomyResult $OUT/aln.m8 \
	--threads ${SLURM_CPUS_PER_TASK} \
	--format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"

echo
echo "MMSEQS LCA"
echo
mmseqs lca \
	$DB \
	$OUT/taxonomyResult \
	$OUT/lcaDB \
	--tax-lineage 1 \
	--threads ${SLURM_CPUS_PER_TASK} \
	--lca-ranks "superkingdom:phylum:class:order:family:genus:species"

# create taxonomy table (tsv)
echo
echo "CREATE TSV"
echo
mmseqs createtsv $OUT/seqtable_queryDB $OUT/lcaDB $OUT/taxonomyResult.tsv --threads ${SLURM_CPUS_PER_TASK}

# create kraken-style report (optional: unhash if interested in this output)
echo
echo "CREATE TAXONOMY REPORT"
echo
mmseqs taxonomyreport $DB $OUT/lcaDB $OUT/taxonomyResult.report --threads ${SLURM_CPUS_PER_TASK}

## Adjust taxonomy table and extract viral lineages
# Extract all (virus + phage) potential viral sequences
echo
echo "ADJUST TAX TABLE AND EXTRACT VIRAL LINEAGES"
echo
echo "EXTRACT ALL POTENTIAL VIRAL SEQUENCES"
echo
grep 'Viruses:' $OUT/taxonomyResult.tsv | cut -f1,5 | sed 's/phi14:2/phi14_2/g' | \
	sed 's/:/\t/g' | \
	sort -n -k1 > $OUT/all_viruses_table.tsv

# Extract phage lineages and generate taxonomy table for import into R as PhyloSeq object
echo
echo "EXTRACT PHAGE LINEAGES AND GENERATE TAX TABLE FOR R"
echo
grep -f $PHAGE $OUT/all_viruses_table.tsv > $OUT/phage_table.tsv
cut -f1 $OUT/phage_table.tsv > $OUT/phage_seqs.list
pullseq -i ./results/seqtable.fasta -n $OUT/phage_seqs.list -l 5000 > $OUT/phage_seqs.fasta

seqkit fx2tab $OUT/phage_seqs.fasta > $OUT/phage_seqs.fx2tab
join $OUT/phage_seqs.fx2tab $OUT/phage_table.tsv | awk -F ' ' '{ print $2,"\t",$3,"\t",$4,"\t",$5,"\t",$6,"\t",$7,"\t",$8,"\t",$9 }' | \
        sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' > $OUT/phage_tax_table.tsv

# Extract non-phage viral lineages and generate taxonomy table for import into R as PhyloSeq object
echo
echo "EXTRACT NON-PHAGE VIRAL LINEAGES AND GENERATE TAX TABLE FOR R"
echo
grep -v -f $PHAGE $OUT/all_viruses_table.tsv > $OUT/viruses_table.tsv
cut -f1 $OUT/viruses_table.tsv > $OUT/viruses_seqs.list
pullseq -i ./results/seqtable.fasta -n $OUT/viruses_seqs.list -l 5000 > $OUT/viruses_seqs.fasta

# Extract unclassified lineages
echo
echo "EXTRACT UNCLASSIFIED LINEAGES"
echo
grep -v 'Viruses:' $OUT/taxonomyResult.tsv | cut -f1,5 | \
	sed 's/:/\t/g' | \
        sort -n -k1 > $OUT/pviral_aa_unclassified_seqs.list
pullseq -i ./results/seqtable.fasta -n $OUT/pviral_aa_unclassified_seqs.list -l 5000 > $OUT/pviral_aa_unclassified_seqs.fasta
