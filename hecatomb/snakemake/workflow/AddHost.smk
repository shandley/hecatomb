"""
Snakefile to add a new host genome for use with Hecatomb.

Hecatomb will filter contaminant reads from the host organism for viral metagenomes using nucleotide alignments.
The genomes need to have their viral-like sequences masked before they can be used for filtering in Hecatomb.

Michael Roach, Q2 2021
"""

import attrmap as ap


# load default config
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "immutable.yaml")
config = ap.AttrMap(config)


# directories
include: "rules/preflight/directories.smk"


# host files
dir.dbs.host.fasta = os.path.join(
    dir.dbs.host.base, config.args.hostName, "masked_ref.fa.gz"
)
dir.dbs.host.virShred = os.path.join(dir.dbs.host.base, "virus_shred.fasta.gz")


rule all:
    input:
        dir.dbs.host.fasta


rule map_shred_seq:
    """Add host step 01: Map the virus shred seqs to the input host"""
    input:
        ref = config.args.hostFa,
        shred = dir.dbs.host.virShred
    output:
        temp(os.path.join(dir.out.temp, f'{config.args.hostName}.sam.gz'))
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    resources:
        mem_mb = config.resources.med.mem,
    threads:
        config.resources.med.cpu
    log:
        os.path.join(dir.out.stderr, 'map_shred_seq.log')
    shell:
        """
        bbmap.sh ref={input.ref} in={input.shred} \
            outm={output} path=tmp/ \
            minid=0.90 maxindel=2 ow=t \
            threads={threads} -Xmx{resources.mem_mb}m &> {log}
        """


rule mask_host:
    """Add host step 02: Mask the host"""
    input:
        ref = config.args.hostFa,
        sam = os.path.join(dir.out.temp, f'{config.args.hostName}.sam.gz')
    output:
        fa = temp(os.path.join(dir.out.temp, f'{config.args.hostName}.processed.fasta')),
        gz = dir.dbs.host.fasta
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    resources:
        mem_mb = config.resources.med.mem,
    threads:
        config.resources.med.cpu
    params:
        entropy = config.qc.entropy
    log:
        os.path.join(dir.out.stderr, 'mask_host.log')
    shell:
        """
        bbmask.sh in={input.ref} out={output.fa} \
            entropy={params.entropy} sam={input.sam} ow=t \
            threads={threads} -Xmx{resources.mem_mb}m &> {log}
        gzip -c {output.fa} > {output.gz}
        """
