#!/bin/bash

# Simplest way to access the UniProt viral protein resource:
# 1) https://www.uniprot.org/uniprot/?query=taxonomy:%22Viruses%20[10239]%22&fil=reviewed%3Ano&sort=score
# 2) Select Download, Download All, format = FASTA (canonical)
# 3) Download to local computer and FTP to server
# *I know there is likely an easier way to wget this data directly, but I am too lazy to sort that out right now

# RAE:
# Turns out this is really simple. Here is the URL:
# https://www.uniprot.org/uniprot/?query=taxonomy:%22Viruses%20[10239]%22&format=fasta&limit=10&sort=score&fil=reviewed:no
# see: https://www.uniprot.org/help/api_queries for more information on how to get the URL!

# Write	out date/time of DB build
date > uniprotDB_build.timestamp;

# Cluster DB @ 99% ID
# Change -c flag to cluster at a different %ID
cd-hit -i uniprot-viral.fasta -o virus_uniprot_c99.fasta -c 0.99 -M 64000 -T 64;

# Create mmseqs2 DB
mmseqs createdb virus_uniprot_c99.fasta targetDB;
mmseqs createtaxdb targetDB tmp;
