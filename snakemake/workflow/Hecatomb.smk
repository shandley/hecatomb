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
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "immutable.yaml")
config = ap.AttrMap(config)


### LAUNCHER-CONTROLLED CONFIG--"Reads" MUST BE PASSED TO SNAKEMAKE
if config.args.search == "fast":
    config.mmseqs.sensAA = config.mmseqs.perfAAfast
    config.mmseqs.sensNT = config.mmseqs.perfNTfast
else:
    config.mmseqs.sensAA = config.mmseqs.perfAA
    config.mmseqs.sensNT = config.mmseqs.perfNT


### DIRECTORIES
include: "rules/preflight/directories.smk"


### HOST ORGANISM
dir.dbs.host.fasta = os.path.join(
    dir.dbs.host.base, config.args.host, "masked_ref.fa.gz"
)
dir.dbs.host.index = dir.dbs.host.fasta + ".idx"


### PREFLIGHT CHECKS
include: "rules/preflight/validate.smk"
include: "rules/preflight/functions.smk"


### TARGETS
include: "rules/preflight/targets.smk"


### PARSE SAMPLES
if config.args.preprocess == "paired":
    include: "rules/preflight/samples.smk"
elif config.args.preprocess == "single":
    include: "rules/preflight/samples_se.smk"
elif config.args.preprocess == "longread":
    include: "rules/preflight/samples_se.smk"
else:  # config.args.preprocess == 'roundAB'
    include: "rules/preflight/samples.smk"


samples = ap.AttrMap()
samples.reads = parseSamples(config.args.reads)
samples.names = list(ap.utils.get_keys(samples.reads))
# wildcard_constraints:
#     sample="[a-zA-Z0-9._-]+"


### PREPROCESSING
if config.args.preprocess == "paired":
    include: "rules/preprocessing/paired_end.smk"
    include: "rules/assembly/paired_end.smk"
elif config.args.preprocess == "single":
    include: "rules/preprocessing/single_end.smk"
    include: "rules/assembly/single_end.smk"
elif config.args.preprocess == "longread":
    include: "rules/preprocessing/longreads.smk"
    include: "rules/assembly/longreads.smk"
else:  # config.args.preprocess == 'roundAB'
    include: "rules/preprocessing/roundAB.smk"
    include: "rules/assembly/paired_end.smk"


### REMAINING PIPELINE RULES
include: "rules/annotation/read_annotation.smk"
include: "rules/assembly/combine_sample_assemblies.smk"
include: "rules/annotation/contig_mapping.smk"
include: "rules/annotation/contig_annotation.smk"
include: "rules/reports/summaries.smk"


# Mark target rules
target_rules = []


def targetRule(fn):
    assert fn.__name__.startswith("__")
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
