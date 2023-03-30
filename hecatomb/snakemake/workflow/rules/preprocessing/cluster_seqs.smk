rule create_host_index:
    """Create the minimap2 index for mapping to the host"""
    input:
        host = dir.dbs.host.fasta,
        cont = os.path.join(dir.dbs.contaminants, "vector_contaminants.fa")
    output:
        dir.dbs.host.index
    benchmark:
        os.path.join(dir.out.bench,"create_host_index.txt")
    log:
        os.path.join(dir.out.stderr,'create_host_index.log')
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.sml.time
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    shell:
        """
        minimap2 -t {threads} -d {output} <(zcat {input.host}; cat {input.cont}) 2> {log}
        rm {log}
        """


rule cluster_similar_sequences:
    """Preprocessing step 05: Cluster similar sequences.

     Sequences clustered at CLUSTERID in config.yaml.
    """
    input:
        fq=os.path.join(dir.out.temp,"p04","{sample}_R1.all.fastq"),
        # summ = optionalSummary[1]
    output:
        temp(os.path.join(dir.out.temp,"p05","{sample}_R1_rep_seq.fasta")),
        temp(os.path.join(dir.out.temp,"p05","{sample}_R1_cluster.tsv")),
        temp(os.path.join(dir.out.temp,"p05","{sample}_R1_all_seqs.fasta"))
    params:
        respath=lambda wildcards, output: os.path.split(output[0])[0],
        tmppath=lambda wildcards, output: os.path.join(os.path.split(output[0])[0],f"{wildcards.sample}_TMP"),
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
        os.path.join(dir.env,"mmseqs2.yaml")
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
        seqs=os.path.join(dir.out.temp,"p05","{sample}_R1_rep_seq.fasta"),
        counts=os.path.join(dir.out.temp,"p05","{sample}_R1_cluster.tsv"),
    # summ = optionalSummary[2]
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
        os.path.join(dir.env,"seqkit.yaml")
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
        seqtables=expand(os.path.join(dir.out.temp,"p06","{sample}_R1.seqtable"),sample=samples.names),
    output:
        fa=os.path.join(dir.out.results,"seqtable.fasta"),
        tsv=os.path.join(dir.out.results,"sampleSeqCounts.tsv")
    params:
        samples=samples.names,
        tmpdir=lambda wildcards, input: os.path.split(input[0])[0]
    conda:
        os.path.join(dir.env,'pysam.yaml')
    benchmark:
        os.path.join(dir.out.bench,"merge_seq_table.txt")
    log:
        os.path.join(dir.out.stderr,'merge_seq_table.log')
    script:
        os.path.join(dir.scripts,'mergeSeqTable.py')