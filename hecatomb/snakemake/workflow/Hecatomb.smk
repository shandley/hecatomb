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
import attrmap.utils as au


### CONFIG
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "immutable.yaml")
config = ap.AttrMap(config)


### LAUNCHER-CONTROLLED CONFIG SETTINGS
if config.args.search == "fast":
    config.mmseqs.sensAA = config.mmseqs.perfAAfast
    config.mmseqs.sensNT = config.mmseqs.perfNTfast
else:
    config.mmseqs.sensAA = config.mmseqs.perfAA
    config.mmseqs.sensNT = config.mmseqs.perfNT


### DIRECTORIES
include: os.path.join("rules", "preflight", "directories.smk")


### HOST ORGANISM
dir.dbs.host.fasta = os.path.join(
    dir.dbs.host.base, config.args.host, "masked_ref.fa.gz"
)
dir.dbs.host.index = dir.dbs.host.fasta + ".idx"


### PREFLIGHT CHECKS, PARSE SAMPLES
include: os.path.join("rules", "preflight", "validate.smk")
include: os.path.join("rules", "preflight", "functions.smk")
include: config.modules[config.args.library]["preflight"]


samples = ap.AttrMap()
samples.reads = parseSamples(config.args.reads)
samples.names = list(ap.utils.get_keys(samples.reads))
samples = au.convert_state(samples, read_only=True)

# wildcard_constraints:
#     sample="[a-zA-Z0-9._-]+"


### TARGETS (must be included AFTER parsing samples)
include: os.path.join("rules", "preflight", "targets.smk")


### PREPROCESSING
include: config.modules[config.args.library]["preprocessing"]
include: config.modules[config.args.library]["assembly"]


### REMAINING PIPELINE RULES
include: os.path.join("rules","preprocessing","cluster_seqs.smk")
include: os.path.join("rules","annotation","read_annotation.smk")
include: os.path.join("rules","assembly","combine_sample_assemblies.smk")
include: os.path.join("rules","annotation","contig_mapping.smk")
include: os.path.join("rules","annotation","contig_annotation.smk")
include: os.path.join("rules","reports","summaries.smk")
include: os.path.join("rules","reports","summaries_optional.smk")


# Mark target rules
target_rules = []


def targetRule(fn):
    assert fn.__name__.startswith("__")
    target_rules.append(fn.__name__[2:])
    return fn


localrules: all, preprocess, assemble, annotate, ctg_annotate, print_stages, dumpSamplesTsv


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
rule preprocess:
    input:
        targets.preprocessing


@targetRule
rule assemble:
    input:
        targets.assembly


@targetRule
rule annotate:
    input:
        targets.readAnnotations


@targetRule
rule ctg_annotate:
    input:
        targets.contigAnnotations,
        targets.mapping


@targetRule
rule print_stages:
    run:
        print("\nIndividual Hecatomb stages to run: \n", file=sys.stderr)
        print("* " + "\n* ".join(target_rules) + "\n\n", file=sys.stderr)
