#!/bin/bash
set -e
set -u

# Script to mask virus and low-entropy (repeats) sequences from host reference genomes

# Dependencies:
# bbtools: https://jgi.doe.gov/data-and-tools/bbtools/

# Location of shredded virus sequences
SHRED=/mnt/data1/databases/hecatomb/nt/ncbi/ncbi_assembly_virus/virus_shred.fasta.gz

# Help function
helpFunction()
{
        echo "Script to mask virus and low-entropy (repeats) sequences from host reference genomes."
        echo ""
        echo "Syntax: scriptTemplate [-h | -r | -e]"
        echo "options:"
	echo "-h Display help"
        echo "-r Reference Genome (FASTA format. Can by gz or not.)"
	echo "-e Entropy (value from 0-1. See bbmask.sh help for more detail)"
        echo ""
        exit 1 # Exit script after printing help
}

# Retrieve and set options/flags
while getopts "hr:e:" opt; do
        case $opt in
                h) # display Help
                helpFunction
                ;;

                r) # reference genome
                REF="$OPTARG"
                echo "The name of the directory containing your reference genome is $OPTARG"
		echo ""
                ;;

		e) # entropy
		ENT="$OPTARG"
		echo "Entropy is set to $OPTARG"
		echo ""
		;;

                \?) # incorrect option
                echo "Usage: cmd [-h] [-r] [-t]"
                exit
                ;;
        esac
done

# Mapping
echo ""
echo The shredded virus genomes are in: ${SHRED}.
echo ""
echo The reference is located in: ${REF}
echo ""

# Set input file name
IN=$(ls ../host/${REF}/*.fa.gz)
echo Input file name: ${IN}
echo ""

bbmap.sh ref=${IN} in=${SHRED} outm=../host/${REF}/${REF}.sam.gz \
	path=../host/${REF} \
	minid=0.90 maxindel=2 \ # See: https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/bbmask-guide/
	ow=t

# Masking
bbmask.sh in=${IN} out=../host/${REF}/${REF}_masked.fa.gz entropy=${ENT} sam=../host/${REF}/${REF}.sam.gz ow=t
