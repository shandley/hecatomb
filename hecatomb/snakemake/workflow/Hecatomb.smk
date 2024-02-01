import os
from metasnek import fastq_finder


### CONFIG
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "immutable.yaml")
resources = config["resources"]

trimnami = config["trimnami"]
trimnami["resources"] = resources
config = config["hecatomb"]


### LAUNCHER-CONTROLLED CONFIG SETTINGS
if config["args"]["search"] == "fast":
    config["mmseqs"]["sens"] = config["mmseqs"]["fast"]
else:
    config["mmseqs"]["sens"] = config["mmseqs"]["sensitive"]


### DIRECTORIES
include: os.path.join("rules", "preflight", "directories.smk")


### HOST ORGANISM
if os.path.isfile(config["args"]["host"]):
    dir["dbs"]["hostFasta"] = config["args"]["host"]
else:
    dir["dbs"]["hostFasta"] = os.path.join(
        dir["dbs"]["hostBase"], config["args"]["host"], "masked_ref.fa.gz"
    )


### PREFLIGHT CHECKS, PARSE SAMPLES
include: os.path.join("rules", "preflight", "validate.smk")
include: os.path.join("rules", "preflight", "functions.smk")


samples = dict()
samples["reads"] = dict(sorted(fastq_finder.parse_samples_to_dictionary(config["args"]["reads"]).items()))
samples["names"] = sorted(list(samples["reads"].keys()))


### TARGETS (must be included AFTER parsing samples)
include: os.path.join("rules", "preflight", "targets.smk")
include: os.path.join("rules", "preprocessing", "preprocessing.smk")


### ASSEMBLY
if config["args"]["trim"] == "filtlong":
    include: os.path.join("rules", "assembly", "longreads.smk")
else:
    include: os.path.join("rules", "assembly", "shortreads.smk")


### REMAINING RULES
include: os.path.join("rules","annotation","read_annotation.smk")
include: os.path.join("rules","assembly","coverage.smk")
include: os.path.join("rules","annotation","contig_mapping.smk")
include: os.path.join("rules","annotation","contig_annotation.smk")
include: os.path.join("rules","reports","summaries.smk")
include: os.path.join("rules","reports","summaries_optional.smk")


target_rules = []

def targetRule(fn):
    assert fn.__name__.startswith("__")
    target_rules.append(fn.__name__[2:])
    return fn


@targetRule
rule all:
    input:
        targets["preprocessing"],
        targets["assembly"],
        targets["readAnnotations"],
        targets["contigAnnotations"],
        targets["mapping"],
        # targets["summary"]


@targetRule
rule preprocessing:
    input:
        targets["preprocessing"]


@targetRule
rule assembly:
    input:
        targets["assembly"]


@targetRule
rule read_annotations:
    input:
        targets["readAnnotations"]


@targetRule
rule contig_annotations:
    input:
        targets["contigAnnotations"]


@targetRule
rule combined_annotations:
    input:
        targets["mapping"]

@targetRule
rule build_envs:
    input:
        targets["envs"]


@targetRule
rule print_stages:
    localrule:
        True
    run:
        print("\nIndividual Hecatomb stages to run: \n", file=sys.stderr)
        print("* " + "\n* ".join(target_rules) + "\n\n", file=sys.stderr)
