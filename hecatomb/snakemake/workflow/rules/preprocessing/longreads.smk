# rules
rule host_removal_mapping:
    """Preprocessing step 07a: Host removal: mapping to host. 

    Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    If your reference is not available you need to add it using 'Hecatomb addHost'
    """
    input:
        r1 = lambda wildcards: samples.reads[wildcards.sample]['R1'],
        host = dir.dbs.host.index,
        # summ = optionalSummary[0]
    output:
        r1=temp(os.path.join(dir.out.temp, "p04", "{sample}_R1.all.fastq"))
    benchmark:
        os.path.join(dir.out.bench,"host_removal_mapping.{sample}.txt")
    log:
        mm=os.path.join(dir.out.stderr,"host_removal_mapping.{sample}.minimap.log"),
        sv=os.path.join(dir.out.stderr,"host_removal_mapping.{sample}.samtoolsView.log"),
        fq=os.path.join(dir.out.stderr,"host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb=config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    group:
        "preprocessing"
    shell:
        """
        minimap2 -ax map-ont -t {threads} --secondary=no {input.host} {input.r1} 2> {log.mm} \
            | samtools view -f 4 -h 2> {log.sv} \
            | samtools fasta > {output.r1} 2> {log.fq}
        rm {log.mm} {log.sv} {log.fq}
        """


rule archive_for_assembly:
    """Copy the files that will be required in the assembly steps; fastq.gz files will be generated from these"""
    input:
        os.path.join(dir.out.temp, "p04", "{sample}_R1.all.fastq")
    output:
        os.path.join(dir.out.assembly,"{sample}_R1.all.fasta.gz")
    params:
        compression = "-" + str(config.qc.compression)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env,"pigz.yaml")
    group:
        "preprocessing"
    shell:
        """
        pigz -p {threads} -c {params.compression} {input} > {output}
        """