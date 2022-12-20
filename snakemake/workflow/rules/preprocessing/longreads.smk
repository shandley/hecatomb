
# Add longread-specific targets
PreprocessingFiles += [
        expand(os.path.join(dir.out.assembly,"{sample}_R1.all.fasta.gz"), sample=samples.names)
    ]


# rules
rule create_host_index:
    """Create the minimap2 index for mapping to the host; this will save time."""
    input:
        dir.dbs.host.fasta,
    output:
        dir.dbs.host.index
    benchmark:
        os.path.join(dir.out.bench,"create_host_index.txt")
    log:
        os.path.join(dir.out.stderr,'create_host_index.log')
    resources:
        mem_mb=config.resources.med.mem
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
    """Preprocessing step 07a: Host removal: mapping to host. 

    Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    If your reference is not available you need to add it using 'Hecatomb addHost'
    """
    input:
        r1 = lambda wildcards: samples.reads[wildcards.sample]['R1'],
        host = dir.dbs.host.index,
        # summ = optionalSummary[0]
    output:
        r1=temp(os.path.join(dir.out.temp,"p01","{sample}_R1.all.fasta")),
    benchmark:
        os.path.join(dir.out.bench,"host_removal_mapping.{sample}.txt")
    log:
        mm=os.path.join(dir.out.stderr,"host_removal_mapping.{sample}.minimap.log"),
        sv=os.path.join(dir.out.stderr,"host_removal_mapping.{sample}.samtoolsView.log"),
        fq=os.path.join(dir.out.stderr,"host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb=config.resources.med.mem
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    shell:
        """
        minimap2 -ax map-ont -t {threads} --secondary=no {input.host} {input.r1} 2> {log.mm} \
            | samtools view -f 4 -h 2> {log.sv} \
            | samtools fasta > {output.r1} 2> {log.fq}
        rm {log.mm} {log.sv} {log.fq}
        """


rule cluster_similar_sequences:  ### TODO: CHECK IF WE STILL HAVE ANY READS LEFT AT THIS POINT
    """Preprocessing step 09: Cluster similar sequences.

     Sequences clustered at CLUSTERID in config.yaml.
    """
    input:
        fa=os.path.join(dir.out.temp,"p01","{sample}_R1.all.fasta"),
        # summ=optionalSummary[1]
    output:
        temp(os.path.join(dir.out.temp,"p05","{sample}_R1_rep_seq.fasta")),
        temp(os.path.join(dir.out.temp,"p05","{sample}_R1_cluster.tsv")),
        temp(os.path.join(dir.out.temp,"p05","{sample}_R1_all_seqs.fasta"))
    params:
        respath=os.path.join(dir.out.temp,"p05"),
        tmppath=os.path.join(dir.out.temp,"p05","{sample}_TMP"),
        prefix='{sample}_R1',
        config=config.mmseqs.linclustParams
    benchmark:
        os.path.join(dir.out.bench,"cluster_similar_sequences.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"cluster_similar_sequences.{sample}.log")
    resources:
        mem_mb=config.resources.big.mem
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    shell:
        """ 
        mmseqs easy-linclust {input.fa} {params.respath}/{params.prefix} {params.tmppath} \
            {params.config} \
            --threads {threads} &> {log}
        rm {log}
        """


rule create_individual_seqtables:
    """Preprocessing step 10: Create individual seqtables. 

    A seqtable is a count table with each sequence as a row, each column as a sample and each cell the counts of each 
    sequence per sample.
    """
    input:
        seqs=os.path.join(dir.out.temp,"p05","{sample}_R1_rep_seq.fasta"),
        counts=os.path.join(dir.out.temp,"p05","{sample}_R1_cluster.tsv"),
        # summ=optionalSummary[2]
    output:
        seqs=temp(os.path.join(dir.out.temp,"p06","{sample}_R1.seqs")),
        counts=temp(os.path.join(dir.out.temp,"p06","{sample}_R1.counts")),
        seqtable=temp(os.path.join(dir.out.temp,"p06","{sample}_R1.seqtable"))
    benchmark:
        os.path.join(dir.out.bench,"individual_seqtables.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"individual_seqtables.{sample}.txt")
    resources:
        mem_mb=config.resources.big.mem
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
    """Preprocessing step 11: Merge seq tables

    Reads the sequences and counts from each samples' seqtable text file and converts to fasta format for the rest of 
    the pipline.
    """
    input:
        seqtables=expand(os.path.join(dir.out.temp,"p06","{sample}_R1.seqtable"),sample=samples.names)
    output:
        fa=os.path.join(dir.out.results,"seqtable.fasta"),
        tsv=os.path.join(dir.out.results,"sampleSeqCounts.tsv")
    params:
        samples=samples.names,
        tmpdir=lambda w, input: os.path.split(input[0])[0]
    conda:
        os.path.join(dir.env,   'pysam.yaml')
    benchmark:
        os.path.join(dir.out.bench,"merge_seq_table.txt")
    log:
        os.path.join(dir.out.stderr,'merge_seq_table.log')
    script:
        os.path.join(dir.scripts,   'mergeSeqTable.py')


rule archive_for_assembly:
    """Copy the files that will be required in the assembly steps; fastq.gz files will be generated from these"""
    input:
        os.path.join(dir.out.temp,"p01","{sample}_R1.all.fasta")
    output:
        temp(os.path.join(dir.out.assembly,"{sample}_R1.all.fasta"))
    params:
        dir.out.assembly
    shell:
        """cp {input} {params}"""