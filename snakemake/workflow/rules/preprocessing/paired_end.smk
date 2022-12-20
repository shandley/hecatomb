
# Add preprocessing-specific targets
targets.preprocessing += [
        expand(os.path.join(dir.out.assembly,"{sample}_R1.unmapped.fastq.gz"), sample=samples.names),
        expand(os.path.join(dir.out.assembly,"{sample}_R1.singletons.fastq.gz"), sample=samples.names),
        expand(os.path.join(dir.out.assembly,"{sample}_R1.all.fastq.gz"), sample=samples.names),
        expand(os.path.join(dir.out.assembly,"{sample}_R2.singletons.fastq.gz"),sample=samples.names),
        expand(os.path.join(dir.out.assembly,"{sample}_R2.unmapped.fastq.gz"),sample=samples.names),
        expand(os.path.join(dir.out.assembly,"{sample}_R2.all.fastq.gz"),sample=samples.names),
    ]


# rules
rule fastp_preprocessing:
    """Preprocessing step 01: fastp_preprocessing.
    
    Use fastP to remove adaptors, vector contaminants, low quality sequences, poly-A tails and reads shorts than minimum length, plus deduplicate.
    """
    input:
        r1 = lambda wildcards: samples.reads[wildcards.sample]['R1'],
        r2 = lambda wildcards: samples.reads[wildcards.sample]['R2'],
        contaminants = os.path.join(dir.dbs.contaminants, "vector_contaminants.fa"),
        # summ = optionalSummary[0]
    output:
        r1 = temp(os.path.join(dir.out.temp, "p01", "{sample}_R1.s1.out.fastq")),
        r2 = temp(os.path.join(dir.out.temp, "p01", "{sample}_R2.s1.out.fastq")),
        stats = os.path.join(dir.out.stats, "p01", "{sample}.s1.stats.json"),
        html = os.path.join(dir.out.stats, "p01", "{sample}.s1.stats.html")
    benchmark:
        os.path.join(dir.out.bench, "fastp_preprocessing.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "fastp_preprocessing.{sample}.log")
    resources:
        mem_mb = config.resources.sml.mem
    threads:
        config.resources.sml.cpu
    conda:
        os.path.join(dir.env, "fastp.yaml")
    params:
        compression = config.qc.compression,
        qscore = config.qc.qscore,
        readlen = config.qc.readMinLen,
        cuttail = config.qc.cutTailWindow,
        dedupacc = config.qc.dedupAccuracy
    shell:
        """
        fastp -i {input.r1} -I {input.r2} -o {output.r1} -O {output.r2} \
            -z {params.compression} \
            -j {output.stats} -h {output.html} \
            --qualified_quality_phred {params.qscore} \
            --length_required {params.readlen} \
            --detect_adapter_for_pe \
            --cut_tail --cut_tail_window_size {params.cuttail} --cut_tail_mean_quality {params.qscore} \
            --dedup --dup_calc_accuracy {params.dedupacc} \
            --trim_poly_x \
            --thread {threads} 2> {log}
        rm {log}
        """

rule create_host_index:
    """Step 02. Create the minimap2 index for mapping to the host; this will save time."""
    input:
        dir.dbs.host.fasta,
    output:
        dir.dbs.host.index
    benchmark:
        os.path.join(dir.out.bench, "create_host_index.txt")
    log:
        os.path.join(dir.out.stderr, 'create_host_index.log')
    resources:
        mem_mb = config.resources.med.mem
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    shell:
        """
        minimap2 -t {threads} -d {output} <(cat {input}) 2> {log}
        rm {log}
        """

rule host_removal_mapping:
    """Preprocessing step 02a: Host removal: mapping to host.
    
    Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    If your reference is not available you need to add it using 'Hecatomb addHost'
    """
    input:
        r1 = os.path.join(dir.out.temp, "p01", "{sample}_R1.s1.out.fastq"),
        r2 = os.path.join(dir.out.temp, "p01", "{sample}_R2.s1.out.fastq"),
        host = dir.dbs.host.index
    output:
        r1 = temp(os.path.join(dir.out.temp, "p02", "{sample}_R1.unmapped.fastq")),
        r2 = temp(os.path.join(dir.out.temp, "p02", "{sample}_R2.unmapped.fastq")),
        s = temp(os.path.join(dir.out.temp, "p02", "{sample}_R1.unmapped.singletons.fastq")),
        o = temp(os.path.join(dir.out.temp, "p02", "{sample}_R1.other.singletons.fastq"))
    benchmark:
        os.path.join(dir.out.bench, "host_removal_mapping.{sample}.txt")
    log:
        mm = os.path.join(dir.out.stderr, "host_removal_mapping.{sample}.minimap.log"),
        sv = os.path.join(dir.out.stderr, "host_removal_mapping.{sample}.samtoolsView.log"),
        fq = os.path.join(dir.out.stderr, "host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb = config.resources.med.mem
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "minimap2.yaml")
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
        javaAlloc = int(0.9 * config.resources.med.mem)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
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
        t1 = temp(os.path.join(dir.out.temp, "p04", "{PATTERN}_R1.singletons.fastq")),
        t2 = temp(os.path.join(dir.out.temp, "p04", "{PATTERN}_R2.singletons.fastq")),
        r1 = temp(os.path.join(dir.out.temp, "p04", "{PATTERN}_R1.all.fastq")),
        r2 = temp(os.path.join(dir.out.temp, "p04", "{PATTERN}_R2.all.fastq"))
    benchmark:
        os.path.join(dir.out.bench, "nonhost_read_combine.{PATTERN}.txt")
    log:
        os.path.join(dir.out.stderr, "nonhost_read_combine.{PATTERN}.log")
    shell:
        """
        {{ cat {input.sr1} {input.or1} > {output.t1};
        cat {input.sr2} {input.or2} > {output.t2};
        cat {input.r1} {output.t1} > {output.r1};
        cat {input.r2} {output.t2} > {output.r2}; }} 2> {log}
        rm {log}
        """
          
rule cluster_similar_sequences: ### TODO: CHECK IF WE STILL HAVE ANY READS LEFT AT THIS POINT
    """Preprocessing step 05: Cluster similar sequences.
     
     Sequences clustered at CLUSTERID in config.yaml.
    """
    input:
        fq = os.path.join(dir.out.temp, "p04", "{sample}_R1.all.fastq"),
        # summ = optionalSummary[1]
    output:
        temp(os.path.join(dir.out.temp, "p05", "{sample}_R1_rep_seq.fasta")),
        temp(os.path.join(dir.out.temp, "p05", "{sample}_R1_cluster.tsv")),
        temp(os.path.join(dir.out.temp, "p05", "{sample}_R1_all_seqs.fasta"))
    params:
        respath = lambda w, output: os.path.split(output[0])[0],
        tmppath = lambda wildcards, output: os.path.join(os.path.split(output[0])[0], f"{wildcards.sample}_TMP"),
        prefix = '{sample}_R1',
        config = config.mmseqs.linclustParams
    benchmark:
        os.path.join(dir.out.bench, "cluster_similar_sequences.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "cluster_similar_sequences.{sample}.log")
    resources:
        mem_mb = config.resources.big.mem
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    shell:
        """ 
        mmseqs easy-linclust {input.fq} {params.respath}/{params.prefix} {params.tmppath} \
            {params.config} \
            --threads {threads} &> {log}
        rm {log}
        """
        
rule create_individual_seqtables:
    """Preprocessing step 06: Create individual seqtables. 
    
    A seqtable is a count table with each sequence as a row, each column as a sample and each cell the counts of each 
    sequence per sample.
    """
    input:
        seqs = os.path.join(dir.out.temp, "p05", "{sample}_R1_rep_seq.fasta"),
        counts = os.path.join(dir.out.temp, "p05", "{sample}_R1_cluster.tsv"),
        # summ = optionalSummary[2]
    output:
        seqs = temp(os.path.join(dir.out.temp, "p06", "{sample}_R1.seqs")),
        counts = temp(os.path.join(dir.out.temp, "p06", "{sample}_R1.counts")),
        seqtable = temp(os.path.join(dir.out.temp, "p06", "{sample}_R1.seqtable"))
    benchmark:
        os.path.join(dir.out.bench, "individual_seqtables.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "individual_seqtables.{sample}.txt")
    resources:
        mem_mb = config.resources.big.mem
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env, "seqkit.yaml")
    shell:
        """
        {{ seqkit sort {input.seqs} --quiet -j {threads} -w 5000 -t dna \
            | seqkit fx2tab -w 5000 -t dna \
            | sed 's/\\t\\+$//' \
            | cut -f2,3 \
            | sed '1i sequence' > {output.seqs};
        cut -f1 {input.counts} \
            | sort \
            | uniq -c \
            | awk -F ' ' '{{print$2"\\t"$1}}' \
            | cut -f2 \
            | sed "1i {wildcards.sample}" > {output.counts};
        paste {output.seqs} {output.counts} > {output.seqtable}; }} 2> {log}
        rm {log}
        """


rule merge_seq_table:
    """Preprocessing step 07: Merge seq tables
    
    Reads the sequences and counts from each samples' seqtable text file and converts to fasta format for the rest of 
    the pipline.
    """
    input:
        seqtables = expand(os.path.join(dir.out.temp, "p06", "{sample}_R1.seqtable"), sample=samples.names),
    output:
        fa = os.path.join(dir.out.results, "seqtable.fasta"),
        tsv = os.path.join(dir.out.results, "sampleSeqCounts.tsv")
    params:
        samples = samples.names,
        tmpdir = lambda w, input: os.path.split(input[0])[0]
    conda:
        os.path.join(dir.env, 'pysam.yaml')
    benchmark:
        os.path.join(dir.out.bench, "merge_seq_table.txt")
    log:
        os.path.join(dir.out.stderr, 'merge_seq_table.log')
    script:
        os.path.join(dir.scripts,   'mergeSeqTable.py')

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
        temp(os.path.join(dir.out.assembly,"{sample}_R1.unmapped.fastq")),
        temp(os.path.join(dir.out.assembly,"{sample}_R2.unmapped.fastq")),
        temp(os.path.join(dir.out.assembly,"{sample}_R1.singletons.fastq")),
        temp(os.path.join(dir.out.assembly,"{sample}_R2.singletons.fastq")),
        temp(os.path.join(dir.out.assembly,"{sample}_R1.all.fastq")),
        temp(os.path.join(dir.out.assembly,"{sample}_R2.all.fastq")),
    params:
        dir.out.assembly
    shell:
        """cp {input} {params}"""
