import os
import attrmap as ap

from metasnek import fastq_finder


### CONFIG
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "immutable.yaml")
resources = ap.AttrMap(config["resources"])
trimnami = ap.AttrMap(config["trimnami"])
trimnami.resources = resources
config = ap.AttrMap(config["hecatomb"])


### LAUNCHER-CONTROLLED CONFIG SETTINGS
if config.args.search == "fast":
    config.mmseqs.sens = config.mmseqs.fast
else:
    config.mmseqs.sens = config.mmseqs.sensitive


### DIRECTORIES
include: os.path.join("rules", "preflight", "directories.smk")


### HOST ORGANISM
if os.path.isfile(config.args.host):
    dir.dbs.host.fasta = config.args.host
else:
    dir.dbs.host.fasta = os.path.join(
        dir.dbs.host.base, config.args.host, "masked_ref.fa.gz"
    )


### PREFLIGHT CHECKS, PARSE SAMPLES
include: os.path.join("rules", "preflight", "validate.smk")
include: os.path.join("rules", "preflight", "functions.smk")


samples = ap.AttrMap()
samples.reads = fastq_finder.parse_samples_to_dictionary(config.args.reads)
samples.names = sorted(list(ap.utils.get_keys(samples.reads)))


### TARGETS (must be included AFTER parsing samples)
include: os.path.join("rules", "preflight", "targets.smk")
include: os.path.join("rules", "preprocessing", "preprocessing.smk")


### ASSEMBLY
if config.args.trim == "nanopore":
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
        targets.preprocessing,
        targets.assembly,
        targets.readAnnotations,
        targets.contigAnnotations,
        targets.mapping,
        # targets.summary


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
    localrule:
        True
    run:
        print("\nIndividual Hecatomb stages to run: \n", file=sys.stderr)
        print("* " + "\n* ".join(target_rules) + "\n\n", file=sys.stderr)
