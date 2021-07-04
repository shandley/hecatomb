"""
Snakefile for downloading the databases. Only need to run this once.

Set a custom database install location in hecatomb/snakemake/config/config.yaml

Use the launcher to run:
hecatomb install

Rob Edwards, Feb 2020
Updated: Michael Roach, Q2 2021
"""

import os
import sys


### LOAD DEFAULT CONFIG
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')


# Directory paths
if config['Databases'] is None:
    DBDIR = os.path.join(workflow.basedir,'../../databases')
else:
    DBDIR = config['Databases']

TMPDIR = 'tmp'
LOG = 'logs'
# paths for our databases
BACPATH  = os.path.join(DBDIR, "bac_giant_unique_species")
HOSTPATH = os.path.join(DBDIR, "hosts")
CONPATH  = os.path.join(DBDIR, "contaminants")
PROTPATH = os.path.join(DBDIR, "proteins")
TAXPATH  = os.path.join(DBDIR, "taxonomy")
NUCLPATH = os.path.join(DBDIR, "nucleotides")

UNIVIRDB = os.path.join(PROTPATH, "virus_primary_aa")
#URVPATH = os.path.join(PROTPATH, "uniref_plus_virus") # uniref50 + viruses
UNIREF50VIR = os.path.join(PROTPATH, "virus_secondary_aa")



# create directories
for DIR in [DBDIR, TMPDIR, BACPATH, HOSTPATH, CONPATH, PROTPATH, TAXPATH, NUCLPATH, UNIVIRDB, UNIREF50VIR, LOG]:
    if not os.path.exists(DIR):
        os.mkdir(DIR)


# import the rules
include: "rules/d0_download_dbs.smk"


## bowtie vs bbmap databases
# inputs = []
# if config['Options']['use_bowtie']:
# inputs = [
#     expand(os.path.join(BACPATH, "bac_uniquespecies_giant.masked_Ns_removed.{n}.bt2l"), n=[1,2,3,4]),
#     expand(os.path.join(HOSTPATH, "human_virus_masked.{n}.bt2l"), n=[1,2,3,4]),
#     expand(os.path.join(CONPATH, "line_sine.{n}.bt2"), n=[1,2,3,4]),
#     expand(os.path.join(CONPATH, "line_sine.rev.{m}.bt2"), m=[1,2])
# ]
# else:
#     inputs = [
#         os.path.join(BACPATH, "ref"),
#         os.path.join(HOSTPATH, "ref"),
#         os.path.join(CONPATH, "line_sine.fasta")
#     ]


rule all:
    input:
        # the database directories
        # inputs,
        #os.path.join(PROTPATH, "uniprot_virus.faa"),
        os.path.join(TAXPATH, "uniprot_ncbi_mapping.dat"),
        os.path.join(NUCLPATH, "refseq_virus_nt_UniVec_masked/nt.fnaDB.dbtype"),
        os.path.join(NUCLPATH, "bac_virus_masked/nt.fnaDB.dbtype"),
        os.path.join(NUCLPATH, "refseq_virus_nt_UniVec_masked/nt.fnaDB.index"),
        os.path.join(NUCLPATH, "bac_virus_masked/nt.fnaDB.index"),
        #os.path.join(TAXPATH, "taxonomizr_accessionTaxa.sql"),
        #multiext(os.path.join(UNIVIRDB, "sequenceDB"), ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", ".db_merged.dmp", ".db_delnodes.dmp"),
        #multiext(os.path.join(UNIREF50VIR, "sequenceDB"), ".db_mapping", ".db_names.dmp", ".db_nodes.dmp", ".db_merged.dmp", ".db_delnodes.dmp")
################# TODO: FIX UNIVIRDB / UNIREF50VIR RULES
