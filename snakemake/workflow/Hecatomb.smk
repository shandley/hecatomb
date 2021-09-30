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
        sys.stderr.write(f'ERROR: missing database file {dbFile}\n')
if dbFail:
    sys.stderr.write('ONE OR MORE DATABASE FILES ARE MISSING\nPlease run "hecatomb install" to download the missing database files\n\n')
    exit(1)


SAMPLES,EXTENSIONS = glob_wildcards(os.path.join(READDIR, '{sample}_R1{extensions}'))

if not EXTENSIONS:
    sys.stderr.write("""
        FATAL: We could not parse the sequence file names.
        We are expecting {sample}_R1{extension}, and so your files
        should contain the characters '_R1' in the fwd reads
        and '_R2' in the rev reads
        """)
    sys.exit()

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
include: "rules/03_contig_annotation.smk"
include: "rules/04_summaries.smk"


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
        MappingFiles,
        ## Summary Files
        SummaryFiles
