"""
The snakefile that runs hecatomb.

USE THE LAUNCHER:
hecatomb run --reads test_data/ --profile slurm

Manual launch example:
snakemake -s Hecatomb.smk --profile slurm --config Reads=test_data/ Host=human

Rob Edwards, October 2020
Overhauled: Michael Roach, Q2 2021
"""


### CONFIG FILES
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, '../', 'config', 'dbFiles.yaml')
configfile: os.path.join(workflow.basedir, '../', 'config', 'immutable.yaml')


### LAUNCHER-CONTROLLED CONFIG--"Reads" MUST BE PASSED TO SNAKEMAKE
READS = config['Reads']
HOST = config['Host']
skipAssembly = config['SkipAssembly']
makeReport = config['Report']
if config['Fast']:
    MMSeqsSensAA = config['perfAAfast']
    MMSeqsSensNT = config['perfNTfast']
else:
    MMSeqsSensAA = config['perfAA']
    MMSeqsSensNT = config['perfNT']


### RESOURCE CONFIG
MMSeqsMem = config['BigJobMem']
MMSeqsCPU = config['BigJobCpu']
MMSeqsMemSplit = str(int(0.75 * int(MMSeqsMem))) + 'M'
MMSeqsTimeMin = config['BigJobTimeMin']
MhitMem = config['MediumJobMem']
MhitCPU = config['MediumJobCpu']
BBToolsMem = config['SmallJobMem']
BBToolsCPU = config['SmallJobCpu']
MiscMem = config['MoreRamMem']
MiscCPU = config['MoreRamCpu']
FastpMem = config['MediumJobMem']
FastpCPU = config['MediumJobCpu']


### DIRECTORIES
include: "rules/00_directories.smk"


### HOST ORGANISM
HOSTFA = os.path.join(HOSTPATH, HOST, "masked_ref.fa.gz")
HOSTINDEX = f"{HOSTFA}.idx"


### PREFLIGHT CHECKS
include: "rules/00_preflight.smk"


# Parse the samples and read files
if config['Sampling'] == 'single':
    include: "rules/00_samples_se.smk"
else:
    include: "rules/00_samples.smk"
sampleReads = parseSamples(READS)
SAMPLES = sampleReads.keys()
# wildcard_constraints:
#     sample="[a-zA-Z0-9._-]+"

# Import rules and functions
include: "rules/00_functions.smk"
include: "rules/00_targets.smk"
if config['QC'] == 'longreads':
    include: "rules/01_preprocessing_longreads.smk"
    include: "rules/02_sample_assembly_longreads.smk"
elif config['QC'] == 'single':
    include: "rules/01_preprocessing_single.smk"
    include: "rules/02_sample_assembly_single.smk"
elif config['QC'] == 'roundAB':
    include: "rules/01_preprocessing_round.smk"
    include: "rules/02_sample_assembly.smk"
else:
    include: "rules/01_preprocessing.smk"
    include: "rules/02_sample_assembly.smk"
include: "rules/02_taxonomic_assignment.smk"
include: "rules/03_population_assembly.smk"
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
