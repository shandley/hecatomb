"""
The snakefile that runs hecatomb.

USE THE LAUNCHER:
hecatomb run --reads test_data/ --profile slurm

Manual launch example:
snakemake -s Hecatomb.smk --profile slurm --config Reads=test_data/ Host=human

Rob Edwards, October 2020
Overhauled: Michael Roach, Q2 2021
"""


import os
import sys


"""
Summary:
    # Step 0: Preprocessing (Rule: 01_preprocessing.smk)
    # Step 1: Taxonomic Assignment (Rule: 02_taxonomic_assignment.smk)
    # Step 2: Compile Results (Rule: 03_compile_results.smk)
"""


### DEFAULT CONFIG FILE
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')


### LAUNCHER-CONTROLLED CONFIG--"Reads" MUST BE PASSED TO SNAKEMAKE
READDIR = config['Reads']
HOST = config['Host']
doAssembly = config['Assembly']


### RESOURCE CONFIG
MMSeqsMem = config['MMSeqsMem']
MMSeqsCPU = config['MMSeqsCPU']
MMSeqsMemSplit = str(int(0.75 * int(MMSeqsMem))) + 'M'
BBToolsMem = config['BBToolsMem']
BBToolsCPU = config['BBToolsCPU']
MhitMem = config['MhitMem']
MhitCPU = config['MhitCPU']
MiscMem = config['MiscMem']
MiscCPU = config['MiscCPU']

### DIRECTORIES
include: "rules/00_directories.smk"


### HOST ORGANISM
HOSTFA = os.path.join(HOSTPATH, HOST, "masked_ref.fa.gz")
HOSTINDEX = HOSTFA + '.idx'


# PREFLIGHT CHECKS
fatal_errors = False
fatal_messages = []

###################################################################
#                                                                 #
# REFERENCE DATABASES                                             #
#                                                                 #
###################################################################

# Base path for protein sequence reference databases
# if not os.path.exists(PROTPATH):
#     fatal_messages.append("protein databases")
#     fatal_errors = True

# Primary aa search database:
# The virus protein database, clustered at 99% with cd-hit    
if not os.path.exists(UNIVIRDB):
    fatal_messages.append(UNIVIRDB)
    fatal_errors = True

# Secondary aa search database
# UniRef50 + primary aa search database
if not os.path.exists(UNIREF50VIR):
 fatal_messages.append(UNIREF50VIR)
 fatal_errors = True

# Base path for nucleotide sequence reference databases
#NUCLPATH = os.path.join(DBDIR, "nt")

# The virus nucleotide database, clustered at 100% with linclust
if not os.path.exists(NCBIVIRDB):
 fatal_messages.append(NCBIVIRDB)
 fatal_errors = True
 
# Polymicrobial + plant + virus database
if not os.path.exists(POLYMICRODB):
   fatal_messages.append("nucleotide database {POLYMICRODB}")
   fatal_errors = True

# # output directories for our untranslated (nt-to-nt) searches
# PRIMARY_NT_OUT = os.path.join(RESULTS, "MMSEQS_NT_PRIMARY")
# SECONDARY_NT_OUT = os.path.join(RESULTS, "MMSEQS_NT_SECONDARY")

###################################################################
#                                                                 #
# Taxonomy databases and related information                      #
#                                                                 #
###################################################################

#PHAGE_LINEAGES = os.path.join(DBDIR, "phages", "phage_taxonomic_lineages.txt")
#if not os.path.exists(PHAGE_LINEAGES):
#    fatal_messages.append("phages/phage_taxonomic_lineages.txt")
#    fatal_errors = True

#TAXPATH  = os.path.join(DBDIR, "taxonomy")
#TAXTAX = os.path.join(TAXPATH, "taxonomizr_accessionTaxa.sql")
#if not os.path.exists(TAXTAX):
#    fatal_messages.append(f"taxonomizr database {TAXTAX}")
#    fatal_errors = True

# Bacterial virus masked database for section 07 mmseqs pviral check

#BVMDB = os.path.join(NUCLPATH, "bac_virus_masked", "nt.fnaDB")
#if not os.path.exists(BVMDB):
#    fatal_messages.append(BVMDB)
#    fatal_errors = True


###################################################################
#                                                                 #
# Fatal errors should all be resolved by the download databsaes   #
#                                                                 #
###################################################################

#if fatal_errors:
#    sys.stderr.write("""
#**** FATAL ERRORS ****

#We can't proceed because we can't find one or more of the databases.
#You probably need to download the databases before you can continue.

#Please use the snakefile: 
#   download_databases.smk

#To download and install all the databases.

#Here are a list of the databases that are currently missing:
#""")
#    sys.stderr.write("\n".join(fatal_messages))
#    sys.stderr.write("\n\n")
#    sys.exit(5)

###################################################################
#                                                                 #
# Read the sequence files and parse the file names.               #
#                                                                 #
###################################################################

SAMPLES,EXTENSIONS = glob_wildcards(os.path.join(READDIR, '{sample}_R1{extensions}'))

if not EXTENSIONS:
    sys.stderr.write("""
        FATAL: We could not parse the sequence file names.
        We are expecting {sample}_R1{extension}, and so your files
        should contain the characters '_R1' in the fwd reads
        and '_R2' in the rev reads
        """)
    sys.exit()
# we just get the generic extension. This is changed in Step 1/Volumes/Macintosh HD/Users/shandley/Library/Caches/Nova/42111DAF-02/Volumes/Macintosh HD/Users/shandley/Library/Caches/Nova/42111DAF-0218-485F-908B-6E1034888DEE/10.39.174.207/mnt/data3/shandley/dev/hecatomb_v_2/hecatomb/snakemake/workflow/Snakefile18-485F-908B-6E1034888DEE/10.39.174.207/mnt/data3/shandley/dev/hecatomb_v_2/hecatomb/snakemake/workflow/Hecatomb.smk

file_extension = EXTENSIONS[0]
# a convenience so we don't need to use '{sample}_R1' all the time
PATTERN = '{sample}'
PATTERN_R1 = '{sample}_R1'
PATTERN_R2 = '{sample}_R2'

if len(SAMPLES) == 0:
    sys.stderr.write("FATAL: We could not detect any samples at all.\n")
    sys.stderr.write("You should complain to Rob\n")
    sys.exit()


include: "rules/00_functions.smk"
include: "rules/00_targets.smk"
include: "rules/01_preprocessing.smk"
include: "rules/02_assembly.smk"
include: "rules/02_taxonomic_assignment.smk"
include: "rules/03_mapping.smk"


rule all:
    input:
        #### Output files from 01_preprocessing.smk
        PreprocessingFiles,
        ## Assembly out from 02_assembly.smk
        AssemblyFiles,
        #### Output file from 02_taxonomic_assignment.smk
        ## Primary (virus protein database) translated (nt-to-aa) search
        PrimarySearchFilesAA,
        ## Secondary (likely viral to transkingdom database) translated (nt-to-aa) search
        SecondarySearchFilesAA,
        ## Primary untranslated (nt-to-nt) search
        PrimarySearchFilesNT,
        ## Secondary untranslated (nt-to-nt) search
        SecondarySearchFilesNT,
        ## Contig annotation files from 03_contig_annotation.smk
        ContigAnnotFiles,
        ## Mapping files (read-based contig id)
        MappingFiles
