"""
The snakefile that runs hecatomb.

USE THE LAUNCHER:
hecatomb run --reads test_data/ --profile slurm

Manual launch example:
snakemake -s Hecatomb.smk --profile slurm --config Reads=test_data/ Host=human

Rob Edwards, October 2020
Overhauled: Michael Roach, Q2 2021
"""

import attrmap as ap


### CONFIG
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')
configfile: os.path.join(workflow.basedir, '../', 'config', 'dbFiles.yaml')
configfile: os.path.join(workflow.basedir, '../', 'config', 'immutable.yaml')
config = ap.AttrMap(config)


### LAUNCHER-CONTROLLED CONFIG--"Reads" MUST BE PASSED TO SNAKEMAKE
if config.args.search == 'fast':
    config.mmseqs.sensAA = config.mmseqs.perfAAfast
    config.mmseqs.sensNT = config.mmseqs.perfNTfast
else:
    config.mmseqs.sensAA = config.mmseqs.perfAA
    config.mmseqs.sensNT = config.mmseqs.perfNT


### DIRECTORIES
include: "rules/00_directories.smk"


### HOST ORGANISM
dir.dbs.host.fasta = os.path.join(dir.dbs.host.base, config.args.host, "masked_ref.fa.gz")
dir.dbs.host.index = dir.dbs.host.fasta + ".idx"


### PREFLIGHT CHECKS
include: "rules/00_preflight.smk"
include: "rules/00_functions.smk"


### TARGETS
include: "rules/00_targets.smk"


### PARSE SAMPLES
if config.args.preprocess == 'paired':
    include: "rules/00_samples.smk"
elif config.args.preprocess == 'single':
    include: "rules/00_samples_se.smk"
elif config.args.preprocess == 'longread':
    include: "rules/00_samples_se.smk"
else: # config.args.preprocess == 'roundAB'
    include: "rules/00_samples.smk"


samples = ap.AttrMap()
samples.reads = parseSamples(config.args.reads)
samples.names = list(samples.reads.keys())
# wildcard_constraints:
#     sample="[a-zA-Z0-9._-]+"


### PREPROCESSING
if config.args.preprocess == 'paired':
    include: "rules/01_preprocessing.smk"
    include: "rules/02_sample_assembly.smk"
elif config.args.preprocess == 'single':
    include: "rules/01_preprocessing_single.smk"
    include: "rules/02_sample_assembly_single.smk"
elif config.args.preprocess == 'longread':
    include: "rules/01_preprocessing_longreads.smk"
    include: "rules/02_sample_assembly_longreads.smk"
else: # config.args.preprocess == 'roundAB'
    include: "rules/01_preprocessing_round.smk"
    include: "rules/02_sample_assembly.smk"


### REMAINING PIPELINE RULES
include: "rules/02_taxonomic_assignment.smk"
include: "rules/03_population_assembly.smk"
include: "rules/03_mapping.smk"
include: "rules/03_contig_annotation.smk"
include: "rules/04_summaries.smk"


# Mark target rules
target_rules = []
def targetRule(fn):
    assert fn.__name__.startswith('__')
    target_rules.append(fn.__name__[2:])
    return fn


localrules: all, preprocessing, assembly, annotations, ctg_annotations, print_stages


@targetRule
rule all:
    input:
        targets.preprocessing,
        targets.assembly,
        targets.readAnnotations,
        targets.contigAnnotations,
        targets.mapping,
        targets.summary


@targetRule
rule preprocessing:
    input:
        targets.preprocessing


@targetRule
rule assembly:
    input:
        targets.assembly


@targetRule
rule annotations:
    input:
        targets.readAnnotations


@targetRule
rule ctg_annotations:
    input:
        targets.contigAnnotations,
        targets.mapping


@targetRule
rule print_stages:
    run:
        print("\nIndividual Hecatomb stages to run: \n", file=sys.stderr)
        print("* " + "\n* ".join(target_rules) + "\n\n", file=sys.stderr)