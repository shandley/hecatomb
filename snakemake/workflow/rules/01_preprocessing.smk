"""
Snakemake rule file to preprocess Illumina sequence data for virome analysis.

What is accomplished with these rules?
    - Non-biological sequence removal (primers, adapters)
    - Host sequence removal
    - Removal of redundant sequences (clustering)
        - Creation of sequence count table
        - Calculation of sequence properties (e.g. GC content, tetramer frequencies)

Rob Edwards, Jan 2020
Updated: Scott Handley, March 2021
Updated: Michael Roach, Q2/3 2021
"""


rule remove_5prime_primer:
    """Preprocessing step 01: Remove 5' primer.
    
    Default RdA/B Primer sequences are provided in the file primerB.fa. If your lab uses other primers you will need to
    place them in CONPATH (defined in the Hecatomb.smk) and change the file name from primerB.fa to your file name below.
    """
    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + file_extension),
        r2 = os.path.join(READDIR, PATTERN_R2 + file_extension),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "p01", f"{PATTERN_R1}.s1.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "p01", f"{PATTERN_R2}.s1.out.fastq")),
        stats = os.path.join(STATS, "p01", "{sample}.s1.stats.tsv")
    benchmark:
        os.path.join(BENCH, "p01_leftmost_primerB.{sample}.txt")
    log:
        os.path.join(STDERR, "p01_leftmost_primerB.{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
            removeifeitherbad=f trimpolya=10 ordered=t rcomp=f ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        """

rule remove_3prime_contaminant:
    """Preprocessing step 02: Remove 3' read through contaminant. 
    
    This is sequence that occurs if the library fragment is shorter than 250 bases and the sequencer reads through the 
    the 3' end. We use the full length of primerB plus 6 bases of the adapter to detect this event and remove everything
    to the right of that molecule when detected.
    """
    input:
        r1 = os.path.join(TMPDIR, "p01", f"{PATTERN_R1}.s1.out.fastq"),
        r2 = os.path.join(TMPDIR, "p01", f"{PATTERN_R2}.s1.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "p02", f"{PATTERN_R1}.s2.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "p02", f"{PATTERN_R2}.s2.out.fastq")),
        stats = os.path.join(STATS, "p02", "{sample}.s2.stats.tsv")
    benchmark:
        os.path.join(BENCH, "p02_3prime_contaminant.{sample}.txt")
    log:
        os.path.join(STDERR, "p02_3prime_contaminant.{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=r removeifeitherbad=f ordered=t rcomp=f ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        """

rule remove_primer_free_adapter:
    """Preprocessing step 03: Remove primer free adapter (both orientations). 
    
    Rarely the adapter will be seen in the molecule indpendent of the primer. This removes those instances as well as 
    everything to the right of the detected primer-free adapter.
    """
    input:
        r1 = os.path.join(TMPDIR, "p02", f"{PATTERN_R1}.s2.out.fastq"),
        r2 = os.path.join(TMPDIR, "p02", f"{PATTERN_R2}.s2.out.fastq"),
        primers = os.path.join(CONPATH, "nebnext_adapters.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "p03", f"{PATTERN_R1}.s3.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "p03", f"{PATTERN_R2}.s3.out.fastq")),
        stats = os.path.join(STATS, "p03", "{sample}.s3.stats.tsv")
    benchmark:
        os.path.join(BENCH, "p03_primer_free_adapter.{sample}.txt")
    log:
        os.path.join(STDERR, "p03_primer_free_adapter.{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=10 ktrim=r removeifeitherbad=f ordered=t rcomp=t ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        """

rule remove_adapter_free_primer:
    """Preprocessing step 04: Remove adapter free primer (both orientations). 
    
    Rarely the primer is detected without the primer. This removes those instances as well as everything to the right 
    of the detected adapter-free primer. 
    """
    input:
        r1 = os.path.join(TMPDIR, "p03", f"{PATTERN_R1}.s3.out.fastq"),
        r2 = os.path.join(TMPDIR, "p03", f"{PATTERN_R2}.s3.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "p04", f"{PATTERN_R1}.s4.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "p04", f"{PATTERN_R2}.s4.out.fastq")),
        stats = os.path.join(STATS, "p04", "{sample}.s4.stats.tsv")
    benchmark:
        os.path.join(BENCH, "p04_adapter_free_primer.{sample}.txt")
    log:
        os.path.join(STDERR, "p04_adapter_free_primer.{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=0 removeifeitherbad=f ordered=t rcomp=t ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        """

rule remove_vector_contamination:
    """Preprocessing step 05: Vector contamination removal (PhiX + NCBI UniVecDB)"""
    input:
        r1 = os.path.join(TMPDIR, "p04", f"{PATTERN_R1}.s4.out.fastq"),
        r2 = os.path.join(TMPDIR, "p04", f"{PATTERN_R2}.s4.out.fastq"),
        primers = os.path.join(CONPATH, "vector_contaminants.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "p05", f"{PATTERN_R1}.s5.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "p05", f"{PATTERN_R2}.s5.out.fastq")),
        stats = os.path.join(STATS, "p05", "{sample}.s5.stats.tsv")
    benchmark:
        os.path.join(BENCH, "p05_vector_contamination_{sample}.txt")
    log:
        os.path.join(STDERR, "p05_vector_contamination_{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log};
        """
        
rule remove_low_quality:
    """Preprocessing step 06: Remove remaining low-quality bases and short reads. 
    
    Quality score can be modified in config.yaml (QSCORE).
    """
    input:
        r1 = os.path.join(TMPDIR, "p05", f"{PATTERN_R1}.s5.out.fastq"),
        r2 = os.path.join(TMPDIR, "p05", f"{PATTERN_R2}.s5.out.fastq")
    output:
        r1 = temp(os.path.join(TMPDIR, "p06", f"{PATTERN_R1}.s6.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "p06", f"{PATTERN_R2}.s6.out.fastq")),
        stats = os.path.join(STATS, "p06", "{sample}.s6.stats.tsv")
    benchmark:
        os.path.join(BENCH, "p06_low_quality_{sample}.txt")
    log:
        os.path.join(STDERR, "p06_low_quality_{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            ordered=t \
            qtrim=r maxns=2 \
            entropy={config[ENTROPY]} \
            entropywindow={config[ENTROPYWINDOW]} \
            trimq={config[QSCORE]} \
            minlength={config[READ_MINLENGTH]} \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log};
        """

rule create_host_index:
    """Create the minimap2 index for mapping to the host; this will save time."""
    input:
        HOSTFA,
    output:
        HOSTINDEX
    benchmark:
        os.path.join(BENCH, "p00_create_host_index.txt")
    log:
        os.path.join(STDERR, 'p00_create_host_index.log')
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -t {threads} -d {output} <(cat {input}) 2> {log}"

rule host_removal_mapping:
    """Preprocessing step 07a: Host removal: mapping to host. 
    
    Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    If your reference is not available you need to add it using 'Hecatomb addHost'
    """
    input:
        r1 = os.path.join(TMPDIR, "p06", f"{PATTERN_R1}.s6.out.fastq"),
        r2 = os.path.join(TMPDIR, "p06", f"{PATTERN_R2}.s6.out.fastq"),
        host = HOSTINDEX
    output:
        r1 = temp(os.path.join(TMPDIR, "p07", f"{PATTERN_R1}.unmapped.fastq")),
        r2 = temp(os.path.join(TMPDIR, "p07", f"{PATTERN_R2}.unmapped.fastq")),
        s = temp(os.path.join(TMPDIR, "p07", f"{PATTERN_R1}.unmapped.singletons.fastq")),
        o = temp(os.path.join(TMPDIR, "p07", f"{PATTERN_R1}.other.singletons.fastq"))
    benchmark:
        os.path.join(BENCH, "p07a_host_removal_mapping.{sample}.txt")
    log:
        mm = os.path.join(STDERR, "p07a_host_removal_mapping.{sample}.minimap.log"),
        sv = os.path.join(STDERR, "p07a_host_removal_mapping.{sample}.samtoolsView.log"),
        fq = os.path.join(STDERR, "p07a_host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax sr -t {threads} --secondary=no {input.host} {input.r1} {input.r2} 2> {log.mm} \
            | samtools view -f 4 -h 2> {log.sv} \
            | samtools fastq -NO -1 {output.r1} -2 {output.r2} -0 {output.o} -s {output.s} 2> {log.fq}
        """

rule nonhost_read_repair:
    """Preprocessing step 07b: Parse R1/R2 singletons (if singletons at all)"""
    input:
        s = os.path.join(TMPDIR, "p07", f"{PATTERN_R1}.unmapped.singletons.fastq"),
        o = os.path.join(TMPDIR, "p07", f"{PATTERN_R1}.other.singletons.fastq")
    output:
        sr1 = temp(os.path.join(TMPDIR, "p07", f"{PATTERN_R1}.u.singletons.fastq")),
        sr2 = temp(os.path.join(TMPDIR, "p07", f"{PATTERN_R2}.u.singletons.fastq")),
        or1 = temp(os.path.join(TMPDIR, "p07", f"{PATTERN_R1}.o.singletons.fastq")),
        or2 = temp(os.path.join(TMPDIR, "p07", f"{PATTERN_R2}.o.singletons.fastq"))
    benchmark:
        os.path.join(BENCH, "p07b_read_repair_{sample}.txt")
    log:
        os.path.join(STDERR, "p07b_read_repair_{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        {{ reformat.sh in={input.s} out={output.sr1} out2={output.sr2} \
            -Xmx{resources.mem_mb}m;
        reformat.sh in={input.o} out={output.or1} out2={output.or2} \
            -Xmx{resources.mem_mb}m; }} 2>> {log}
        """

rule nonhost_read_combine:
    """Preprocessing step 07c: Combine paired and singleton reads."""
    input:
        r1 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R1.unmapped.fastq"),
        r2 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R2.unmapped.fastq"),
        sr1 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R1.u.singletons.fastq"),
        sr2 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R2.u.singletons.fastq"),
        or1 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R1.o.singletons.fastq"),
        or2 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R2.o.singletons.fastq")
    output:
        t1 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R1.singletons.fastq"),
        t2 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R2.singletons.fastq"),
        r1 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R1.all.fastq"),
        r2 = os.path.join(TMPDIR, "p07", f"{PATTERN}_R2.all.fastq")
    benchmark:
        os.path.join(BENCH, f"p07c_read_combine_{PATTERN}.txt")
    log:
        os.path.join(STDERR, f"p07c_read_combine_{PATTERN}.log")
    shell:
        """
        {{ cat {input.sr1} {input.or1} > {output.t1};
        cat {input.sr2} {input.or2} > {output.t2};
        cat {input.r1} {output.t1} > {output.r1};
        cat {input.r2} {output.t2} > {output.r2}; }} 2> {log}
        """

rule remove_exact_dups:
    """Preprocessing step 08: Remove exact duplicates
    
    Exact duplicates are considered PCR generated and not accounted for in the count table (seqtable_all.tsv)
    """
    input:
        os.path.join(TMPDIR, "p07", f"{PATTERN_R1}.all.fastq")
    output:
        temp(os.path.join(TMPDIR, "p08", f"{PATTERN_R1}.deduped.out.fastq"))
    benchmark:
        os.path.join(BENCH, "p08_remove_exact_dups.{sample}.txt")
    log:
        os.path.join(STDERR, "p08_remove_exact_dups.{sample}.log")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        dedupe.sh in={input} out={output} \
            ac=f ow=t \
            threads={threads} -Xmx{resources.mem_mb}m 2> {log}
        """
          
rule cluster_similar_sequences: ### TODO: CHECK IF WE STILL HAVE ANY READS LEFT AT THIS POINT
    """Preprocessing step 09: Cluster similar sequences.
     
     Sequences clustered at CLUSTERID in config.yaml.
    """
    input:
        os.path.join(TMPDIR, "p08", f"{PATTERN_R1}.deduped.out.fastq")
    output:
        temp(os.path.join(TMPDIR, "p09", f"{PATTERN_R1}_rep_seq.fasta")),
        temp(os.path.join(TMPDIR, "p09", f"{PATTERN_R1}_cluster.tsv")),
        temp(os.path.join(TMPDIR, "p09", f"{PATTERN_R1}_all_seqs.fasta"))
    params:
        respath = os.path.join(TMPDIR, "p09"),
        tmppath = os.path.join(TMPDIR, "p09", "{sample}_TMP"),
        prefix = PATTERN_R1
    benchmark:
        os.path.join(BENCH, "p09_cluster_similar_sequences.{sample}.txt")
    log:
        os.path.join(STDERR, "p09_cluster_similar_sequences.{sample}.log")
    resources:
        mem_mb = MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """ 
        mmseqs easy-linclust {input} {params.respath}/{params.prefix} {params.tmppath} \
            --kmer-per-seq-scale 0.3 \
            -c {config[CLUSTERID]} --cov-mode 1 --threads {threads} &> {log};
        """
        
rule create_individual_seqtables:
    """Preprocessing step 10: Create individual seqtables. 
    
    A seqtable is a count table with each sequence as a row, each column as a sample and each cell the counts of each 
    sequence per sample.
    """
    input:
        seqs = os.path.join(TMPDIR, "p09", f"{PATTERN_R1}_rep_seq.fasta"),
        counts = os.path.join(TMPDIR, "p09", f"{PATTERN_R1}_cluster.tsv")
    output:
        seqs = temp(os.path.join(TMPDIR, "p10", f"{PATTERN_R1}.seqs")),
        counts = temp(os.path.join(TMPDIR, "p10", f"{PATTERN_R1}.counts")),
        seqtable = temp(os.path.join(TMPDIR, "p10", f"{PATTERN_R1}.seqtable"))
    benchmark:
        os.path.join(BENCH, "p10_individual_seqtables.{sample}.txt")
    log:
        os.path.join(STDERR, "p10_individual_seqtables.{sample}.txt")
    resources:
        mem_mb = BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/seqkit.yaml"
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
        """


rule merge_seq_table:
    """Preprocessing step 11: Merge seq tables
    
    Reads the sequences and counts from each samples' seqtable text file and converts to fasta format for the rest of 
    the pipline.
    """
    input:
        expand(os.path.join(TMPDIR, "p10", "{sample}_R1.seqtable"), sample=SAMPLES)
    output:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        tsv = os.path.join(RESULTS, "sampleSeqCounts.tsv")
    params:
        resultsdir = directory(RESULTS),
    benchmark:
        os.path.join(BENCH, "p11_merge_seq_table.txt")
    log:
        os.path.join(STDERR, 'p11_merge_seq_table.log')
    run:
        import logging
        logging.basicConfig(filename=log[0],filemode='w',level=logging.DEBUG)
        logging.debug('Reading in individual sample seqtables')
        outFa = open(output.fa, 'w')
        outTsv = open(output.tsv, 'w')
        for sample in SAMPLES:
            logging.debug(f'sample {sample}')
            seqId = 0
            seqCounts = 0
            st = os.path.join(TMPDIR, "p10", f"{sample}_R1.seqtable")
            counts = open(st, 'r')
            line = counts.readline() # skip header
            for line in counts:
                l = line.split()
                if len(l) == 2:
                    id = ':'.join((sample, l[1], str(seqId))) # fasta header = >sample:count:seqId
                    seqCounts += int(l[1])
                    seqId = seqId + 1
                    outFa.write(f'>{id}\n{l[0]}\n')
                else:
                    logging.warning(f'Possible malformed sample seqtable {st}')
            counts.close()
            outTsv.write(f'{sample}\t{seqCounts}\n')
        outFa.close()
        outTsv.close()
        logging.debug('Done')

