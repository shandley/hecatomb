#!/bin/bash

# This script will download NCBI's UniVec (core) fasta information
# Note 1: wget is not installed by default on Mac OS. See https://brew.sh for wget installation instructions
# Note 2: The file will be timestamped with the donwload date

wget https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core -O vector_contaminants.fa

# Timestamp file
now=`date +"%Y-%m-%d"`

echo The last time this script was run was: $now > UniVec_download.timestamp

echo "#############################################################################################"
echo "#"
echo "# Congatulations!"
echo "# UniVec database downloaded on $now."
echo "# Please move the file "vector_contaminants.fa" to the hecatomb/contaminants directory."
echo "#"
echo "# A file named UniVec_download.timestamp has been saved containing the download date."
echo "#"
echo "#############################################################################################"
echo
echo

