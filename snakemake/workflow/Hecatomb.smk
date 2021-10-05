"""
The snakefile that runs hecatomb.

USE THE LAUNCHER:
hecatomb run --reads test_data/ --profile slurm

Manual launch example:
snakemake -s Hecatomb.smk --profile slurm --config Reads=test_data/ Host=human

Rob Edwards, October 2020
Overhauled: Michael Roach, Q2 2021
"""


### DEFAULT CONFIG FILE
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')


### LAUNCHER-CONTROLLED CONFIG--"Reads" MUST BE PASSED TO SNAKEMAKE
READDIR = config['Reads']
HOST = config['Host']
doAssembly = config['Assembly']
if config['Fast']:
    MMSeqsSensAA = config['perfAAfast']
    MMSeqsSensNT = config['perfNTfast']
else:
    MMSeqsSensAA = config['perfAA']
    MMSeqsSensNT = config['perfNT']


### RESOURCE CONFIG
MMSeqsMem = config['MMSeqsMem']
MMSeqsCPU = config['MMSeqsCPU']
MMSeqsMemSplit = str(int(0.75 * int(MMSeqsMem))) + 'M'
MMSeqsTimeMin = config['MMSeqsTimeMin']
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


# Check for Database files
dbFail = False
for f in config['dbFiles']:
    dbFile = os.path.join(DBDIR, f)
    if not os.path.isfile(dbFile):
        dbFail = True
        sys.stderr.write("\n"
            f"    ERROR: missing database file {dbFile}"
            "\n")
if dbFail:
    sys.stderr.write("\n"
        "    FATAL: One or more database files is missing.\n"
        "    Please run 'hecatomb install' to download the missing database files.\n"
        "\n")
    sys.exit(1)


# Identify samples
SAMPLES,EXTENSIONS = glob_wildcards(os.path.join(READDIR, '{sample}_R1{extensions}'))

if not EXTENSIONS:
    sys.stderr.write("\n"
        "    FATAL: We could not parse the sequence file names.\n"
        "    We are expecting {sample}_R1{extension}, and so your files\n"
        "    should contain the characters '_R1' in the fwd reads\n"
        "    and '_R2' in the rev reads\n"
        "\n")
    sys.exit(1)

file_extension = EXTENSIONS[0]
# a convenience so we don't need to use '{sample}_R1' all the time
PATTERN = '{sample}'
PATTERN_R1 = '{sample}_R1'
PATTERN_R2 = '{sample}_R2'

if len(SAMPLES) == 0:
    sys.stderr.write("\n"
        "    FATAL: We could not detect any samples at all.\n"
        "    See https://hecatomb.readthedocs.io/en/latest/usage/#read-directory for more info\n"
        "\n")
    sys.exit()


include: "rules/00_functions.smk"
include: "rules/00_targets.smk"
include: "rules/01_preprocessing.smk"
include: "rules/02_assembly.smk"
include: "rules/02_taxonomic_assignment.smk"
include: "rules/03_mapping.smk"
include: "rules/03_contig_annotation.smk"
include: "rules/04_summaries.smk"


rule all:
    input:
        ## Preprocessing
        PreprocessingFiles,
        ## Assembly
        AssemblyFiles,
        ## Translated (nt-to-aa) search
        SecondarySearchFilesAA,
        ## Untranslated (nt-to-nt) search
        SecondarySearchFilesNT,
        ## Contig annotation
        ContigAnnotFiles,
        ## Mapping (read-based contig id)
        MappingFiles,
        ## Summary
        SummaryFiles
