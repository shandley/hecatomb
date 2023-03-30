
# rules
rule remove_5prime_primer:
    """Preprocessing step 01: Remove 5' primer.

    Default RdA/B Primer sequences are provided in the file primerB.fa. If your lab uses other primers you will need to
    place them in dir.dbs.contaminants (defined in the Hecatomb.smk) and change the file name from primerB.fa to your file name below.
    """
    input:
        r1=lambda wildcards: samples.reads[wildcards.sample]['R1'],
        r2=lambda wildcards: samples.reads[wildcards.sample]['R2'],
        primers=os.path.join(dir.dbs.contaminants,"primerB.fa")
    output:
        r1=temp(os.path.join(dir.out.temp,"p01","{sample}_R1.s1.out.fastq")),
        r2=temp(os.path.join(dir.out.temp,"p01","{sample}_R2.s1.out.fastq")),
        stats=os.path.join(dir.out.stats,"p01","{sample}.s1.stats.tsv")
    benchmark:
        os.path.join(dir.out.bench,"remove_5prime_primer.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"remove_5prime_primer.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem),
        time = int(0.5 * config.resources.med.time)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    group:
        "preprocessing"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
            removeifeitherbad=f trimpolya=10 ordered=t rcomp=f ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        rm {log}
        """


rule remove_3prime_contaminant:
    """Preprocessing step 02: Remove 3' read through contaminant. 

    This is sequence that occurs if the library fragment is shorter than 250 bases and the sequencer reads through the 
    the 3' end. We use the full length of primerB plus 6 bases of the adapter to detect this event and remove everything
    to the right of that molecule when detected.
    """
    input:
        r1=os.path.join(dir.out.temp,"p01","{sample}_R1.s1.out.fastq"),
        r2=os.path.join(dir.out.temp,"p01","{sample}_R2.s1.out.fastq"),
        primers=os.path.join(dir.dbs.contaminants,"rc_primerB_ad6.fa")
    output:
        r1=temp(os.path.join(dir.out.temp,"p02","{sample}_R1.s2.out.fastq")),
        r2=temp(os.path.join(dir.out.temp,"p02","{sample}_R2.s2.out.fastq")),
        stats=os.path.join(dir.out.stats,"p02","{sample}.s2.stats.tsv")
    benchmark:
        os.path.join(dir.out.bench,"remove_3prime_contaminant.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"remove_3prime_contaminant.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem),
        time = int(0.5 * config.resources.med.time)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    group:
        "preprocessing"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=r removeifeitherbad=f ordered=t rcomp=f ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        rm {log}
        """


rule remove_primer_free_adapter:
    """Preprocessing step 03: Remove primer free adapter (both orientations). 

    Rarely the adapter will be seen in the molecule indpendent of the primer. This removes those instances as well as 
    everything to the right of the detected primer-free adapter.
    """
    input:
        r1=os.path.join(dir.out.temp,"p02","{sample}_R1.s2.out.fastq"),
        r2=os.path.join(dir.out.temp,"p02","{sample}_R2.s2.out.fastq"),
        primers=os.path.join(dir.dbs.contaminants,"nebnext_adapters.fa")
    output:
        r1=temp(os.path.join(dir.out.temp,"p03","{sample}_R1.s3.out.fastq")),
        r2=temp(os.path.join(dir.out.temp,"p03","{sample}_R2.s3.out.fastq")),
        stats=os.path.join(dir.out.stats,"p03","{sample}.s3.stats.tsv")
    benchmark:
        os.path.join(dir.out.bench,"remove_primer_free_adapter.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"remove_primer_free_adapter.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem),
        time = int(0.5 * config.resources.med.time)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    group:
        "preprocessing"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=10 ktrim=r removeifeitherbad=f ordered=t rcomp=t ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        rm {log}
        """


rule remove_adapter_free_primer:
    """Preprocessing step 04: Remove adapter free primer (both orientations). 

    Rarely the primer is detected without the primer. This removes those instances as well as everything to the right 
    of the detected adapter-free primer. 
    """
    input:
        r1=os.path.join(dir.out.temp,"p03","{sample}_R1.s3.out.fastq"),
        r2=os.path.join(dir.out.temp,"p03","{sample}_R2.s3.out.fastq"),
        primers=os.path.join(dir.dbs.contaminants,"rc_primerB_ad6.fa")
    output:
        r1=temp(os.path.join(dir.out.temp,"p04","{sample}_R1.s4.out.fastq")),
        r2=temp(os.path.join(dir.out.temp,"p04","{sample}_R2.s4.out.fastq")),
        stats=os.path.join(dir.out.stats,"p04","{sample}.s4.stats.tsv")
    benchmark:
        os.path.join(dir.out.bench,"remove_adapter_free_primer.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"remove_adapter_free_primer.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem),
        time = int(0.5 * config.resources.med.time)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    group:
        "preprocessing"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=0 removeifeitherbad=f ordered=t rcomp=t ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        rm {log}
        """


rule remove_vector_contamination:
    """Preprocessing step 05: Vector contamination removal (PhiX + NCBI UniVecDB)"""
    input:
        r1=os.path.join(dir.out.temp,"p04","{sample}_R1.s4.out.fastq"),
        r2=os.path.join(dir.out.temp,"p04","{sample}_R2.s4.out.fastq"),
        primers=os.path.join(dir.dbs.contaminants,"vector_contaminants.fa")
    output:
        r1=temp(os.path.join(dir.out.temp,"p05","{sample}_R1.s5.out.fastq")),
        r2=temp(os.path.join(dir.out.temp,"p05","{sample}_R2.s5.out.fastq")),
        stats=os.path.join(dir.out.stats,"p05","{sample}.s5.stats.tsv")
    benchmark:
        os.path.join(dir.out.bench,"remove_vector_contamination.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"remove_vector_contamination.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem),
        time = int(0.5 * config.resources.med.time)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    group:
        "preprocessing"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        rm {log}
        """


rule remove_low_quality:
    """Preprocessing step 06: Remove remaining low-quality bases and short reads. 

    Quality score can be modified in config.yaml (QSCORE).
    """
    input:
        r1=os.path.join(dir.out.temp,"p05","{sample}_R1.s5.out.fastq"),
        r2=os.path.join(dir.out.temp,"p05","{sample}_R2.s5.out.fastq")
    output:
        r1=temp(os.path.join(dir.out.temp,"p06","{sample}_R1.s6.out.fastq")),
        r2=temp(os.path.join(dir.out.temp,"p06","{sample}_R2.s6.out.fastq")),
        stats=os.path.join(dir.out.stats,"p06","{sample}.s6.stats.tsv")
    benchmark:
        os.path.join(dir.out.bench,"remove_low_quality.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"remove_low_quality.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem),
        time = int(0.5 * config.resources.med.time)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    group:
        "preprocessing"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            ordered=t \
            qtrim=r \
            maxns=2 \
            entropy=0.5 \
            entropywindow=25 \
            trimq=15 \
            minlength=90 \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        rm {log}
        """


rule host_removal_mapping:
    """Preprocessing step 02a: Host removal: mapping to host.
    
    Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    If your reference is not available you need to add it using 'Hecatomb addHost'
    """
    input:
        r1 = os.path.join(dir.out.temp,"p06","{sample}_R1.s6.out.fastq"),
        r2 = os.path.join(dir.out.temp,"p06","{sample}_R2.s6.out.fastq"),
        host = dir.dbs.host.index
    output:
        r1 = os.path.join(dir.out.temp, "p02", "{sample}_R1.unmapped.fastq"),
        r2 = os.path.join(dir.out.temp, "p02", "{sample}_R2.unmapped.fastq"),
        s = temp(os.path.join(dir.out.temp, "p02", "{sample}_R1.unmapped.singletons.fastq")),
        o = temp(os.path.join(dir.out.temp, "p02", "{sample}_R1.other.singletons.fastq"))
    benchmark:
        os.path.join(dir.out.bench, "host_removal_mapping.{sample}.txt")
    log:
        mm = os.path.join(dir.out.stderr, "host_removal_mapping.{sample}.minimap.log"),
        sv = os.path.join(dir.out.stderr, "host_removal_mapping.{sample}.samtoolsView.log"),
        fq = os.path.join(dir.out.stderr, "host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb = config.resources.med.mem,
        time = int(0.5 * config.resources.med.time)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    group:
        "preprocessing"
    shell:
        """
        minimap2 -ax sr -t {threads} --secondary=no {input.host} {input.r1} {input.r2} 2> {log.mm} \
            | samtools view -f 4 -h 2> {log.sv} \
            | samtools fastq -NO -1 {output.r1} -2 {output.r2} -0 {output.o} -s {output.s} 2> {log.fq}
        rm {log.mm} {log.sv} {log.fq}
        """


rule nonhost_read_repair:
    """Preprocessing step 03: Parse R1/R2 singletons (if singletons at all)"""
    input:
        s = os.path.join(dir.out.temp, "p02", "{sample}_R1.unmapped.singletons.fastq"),
        o = os.path.join(dir.out.temp, "p02", "{sample}_R1.other.singletons.fastq")
    output:
        sr1 = temp(os.path.join(dir.out.temp, "p03", "{sample}_R1.u.singletons.fastq")),
        sr2 = temp(os.path.join(dir.out.temp, "p03", "{sample}_R2.u.singletons.fastq")),
        or1 = temp(os.path.join(dir.out.temp, "p03", "{sample}_R1.o.singletons.fastq")),
        or2 = temp(os.path.join(dir.out.temp, "p03", "{sample}_R2.o.singletons.fastq"))
    benchmark:
        os.path.join(dir.out.bench, "nonhost_read_repair.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "nonhost_read_repair.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem),
        time = int(0.5 * config.resources.med.time)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    group:
        "preprocessing"
    shell:
        """
        {{ reformat.sh in={input.s} out={output.sr1} out2={output.sr2} \
            -Xmx{resources.javaAlloc}m;
        reformat.sh in={input.o} out={output.or1} out2={output.or2} \
            -Xmx{resources.javaAlloc}m; }} 2>> {log}
        rm {log}
        """


rule nonhost_read_combine:
    """Preprocessing step 04: Combine paired and singleton reads. """
    input:
        r1 = os.path.join(dir.out.temp, "p02", "{PATTERN}_R1.unmapped.fastq"),
        r2 = os.path.join(dir.out.temp, "p02", "{PATTERN}_R2.unmapped.fastq"),
        sr1 = os.path.join(dir.out.temp, "p03", "{PATTERN}_R1.u.singletons.fastq"),
        sr2 = os.path.join(dir.out.temp, "p03", "{PATTERN}_R2.u.singletons.fastq"),
        or1 = os.path.join(dir.out.temp, "p03", "{PATTERN}_R1.o.singletons.fastq"),
        or2 = os.path.join(dir.out.temp, "p03", "{PATTERN}_R2.o.singletons.fastq")
    output:
        t1 = os.path.join(dir.out.temp, "p04", "{PATTERN}_R1.singletons.fastq"),
        t2 = os.path.join(dir.out.temp, "p04", "{PATTERN}_R2.singletons.fastq"),
        r1 = temp(os.path.join(dir.out.temp, "p04", "{PATTERN}_R1.all.fastq")),
        r2 = temp(os.path.join(dir.out.temp, "p04", "{PATTERN}_R2.all.fastq"))
    benchmark:
        os.path.join(dir.out.bench, "nonhost_read_combine.{PATTERN}.txt")
    log:
        os.path.join(dir.out.stderr, "nonhost_read_combine.{PATTERN}.log")
    resources:
        time = config.resources.sml.time
    group:
        "preprocessing"
    shell:
        """
        {{ cat {input.sr1} {input.or1} > {output.t1};
        cat {input.sr2} {input.or2} > {output.t2};
        cat {input.r1} {output.t1} > {output.r1};
        cat {input.r2} {output.t2} > {output.r2}; }} 2> {log}
        rm {log}
        """


rule archive_for_assembly:
    """Copy the files that will be required in the assembly steps; fastq.gz files will be generated from these"""
    input:
        os.path.join(dir.out.temp,"p02","{sample}_R1.unmapped.fastq"),
        os.path.join(dir.out.temp,"p02","{sample}_R2.unmapped.fastq"),
        os.path.join(dir.out.temp,"p04","{sample}_R1.singletons.fastq"),
        os.path.join(dir.out.temp,"p04","{sample}_R2.singletons.fastq"),
        os.path.join(dir.out.temp,"p04","{sample}_R1.all.fastq"),
        os.path.join(dir.out.temp,"p04","{sample}_R2.all.fastq"),
    output:
        os.path.join(dir.out.assembly,"{sample}_R1.unmapped.fastq.gz"),
        os.path.join(dir.out.assembly,"{sample}_R2.unmapped.fastq.gz"),
        os.path.join(dir.out.assembly,"{sample}_R1.singletons.fastq.gz"),
        os.path.join(dir.out.assembly,"{sample}_R2.singletons.fastq.gz"),
        os.path.join(dir.out.assembly,"{sample}_R1.all.fastq.gz"),
        os.path.join(dir.out.assembly,"{sample}_R2.all.fastq.gz"),
    params:
        dir = dir.out.assembly,
        compression = "-" + str(config.qc.compression),
        prefix = os.path.join(dir.out.assembly, "{sample}*.fastq")
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env,"pigz.yaml")
    resources:
        time = config.resources.sml.time
    group:
        "preprocessing"
    shell:
        """
        cp {input} {params.dir}
        pigz -p {threads} {params.compression} {params.prefix}
        """