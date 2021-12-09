"""
Snakefile to add a new host genome for use with Hecatomb.

Hecatomb will filter contaminant reads from the host organism for viral metagenomes using nucleotide alignments.
The genomes need to have their viral-like sequences masked before they can be used for filtering in Hecatomb.

Michael Roach, Q2 2021
"""


# load default config
# configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')


# directories
include: "rules/00_directories.smk"


# parse config
virShred = os.path.join(HOSTPATH, 'virus_shred.fasta.gz')
hostFasta = config['HostFa']
hostName = config['HostName']
entropy = config['ENTROPY']
BBToolsMem = config['BBToolsMem']
BBToolsCPU = config['BBToolsCPU']

# output files
hostOutFasta = os.path.join(HOSTPATH, hostName, 'masked_ref.fa.gz')


rule all:
    input:
        hostOutFasta


rule map_shred_seq:
    """Add host step 01: Map the virus shred seqs to the input host"""
    input:
        ref = hostFasta,
        shred = virShred
    output:
        temp(f'{hostFasta}.sam.gz')
    conda:
        os.path.join("envs", "bbmap.yaml")
    resources:
        mem_mb = BBToolsMem,
        cpus = BBToolsCPU
    log:
        os.path.join(STDERR, 'h01.bbmap.log')
    shell:
        """
        bbmap.sh ref={input.ref} in={input.shred} \
            outm={output} path=tmp/ \
            minid=0.90 maxindel=2 ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m
        """


rule mask_host:
    """Add host step 02: Mask the host"""
    input:
        ref = hostFasta,
        sam = f'{hostFasta}.sam.gz'
    output:
        fa = temp(os.path.join(WORKDIR, 'tmpHost.fa')),
        gz = hostOutFasta
    conda:
        os.path.join("envs", "bbmap.yaml")
    resources:
        mem_mb = BBToolsMem,
        cpus = BBToolsCPU
    log:
        os.path.join(STDERR, 'bbmask.log')
    shell:
        """
        bbmask.sh in={input.ref} out={output.fa} \
            entropy={entropy} sam={input.sam} ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m
        gzip -c {output.fa} > {output.gz}
        """
