"""
Rules for adding a new host genome for use with contaminant removal. Imported by AddHost.smk
"""

rule map_shred_seq:
    """Add host step 01: Map the virus shred seqs to the input host"""
    input:
        ref = hostFasta,
        shred = virShred
    output:
        temp(f'{hostFasta}.sam.gz')
    conda:
        "../envs/bbmap.yaml"
    resources:
        mem_mb=64000,
        cpus=8
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
        "../envs/bbmap.yaml"
    resources:
        mem_mb=64000,
        cpus=8
    log:
        os.path.join(STDERR, 'bbmask.log')
    shell:
        """
        bbmask.sh in={input.ref} out={output.fa} \
            entropy={entropy} sam={input.sam} ow=t \
            threads={resources.cpus} -Xmx{resources.mem_mb}m
        gzip -c {output.fa} > {output.gz}
        """
