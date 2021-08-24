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
Updated: Michael Roach, Q2 2021
"""

import os
import sys
    
# NOTE: bbtools uses "threads=auto" by default that typically uses all threads, so no need to specify. 
# -Xmx is used to specify the memory allocation for bbtools operations
# Set your -Xmx specifications in your configuration file 

rule remove_5prime_primer:
    """Step 01: Remove 5' primer.
    
    Default RdA/B Primer sequences are provided in the file primerB.fa. If your lab uses other primers you will need to
    place them in CONPATH (defined in the Hecatomb.smk) and change the file name from primerB.fa to your file name below.
    """
    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + file_extension),
        r2 = os.path.join(READDIR, PATTERN_R2 + file_extension),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_01", f"{PATTERN_R1}.s1.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_01", f"{PATTERN_R2}.s1.out.fastq")),
        stats = os.path.join(STATS, "step_01", "{sample}.s1.stats.tsv")
    benchmark:
        os.path.join(BENCH, "remove_leftmost_primerB.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_leftmost_primerB.{sample}.log")
    resources:
        mem_mb=BBToolsMem
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
    """Step 02: Remove 3' read through contaminant. 
    
    This is sequence that occurs if the library fragment is shorter than 250 bases and the sequencer reads through the 
    the 3' end. We use the full length of primerB plus 6 bases of the adapter to detect this event and remove everything
    to the right of that molecule when detected.
    """
    input:
        r1 = os.path.join(TMPDIR, "step_01", f"{PATTERN_R1}.s1.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_01", f"{PATTERN_R2}.s1.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_02", f"{PATTERN_R1}.s2.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_02", f"{PATTERN_R2}.s2.out.fastq")),
        stats = os.path.join(STATS, "step_02", "{sample}.s2.stats.tsv")
    benchmark:
        os.path.join(BENCH, "remove_3prime_contaminant.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_3prime_contaminant.{sample}.log")
    resources:
        mem_mb=BBToolsMem
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
    """Step 03: Remove primer free adapter (both orientations). 
    
    Rarely the adapter will be seen in the molecule indpendent of the primer. This removes those instances as well as 
    everything to the right of the detected primer-free adapter.
    """
    input:
        r1 = os.path.join(TMPDIR, "step_02", f"{PATTERN_R1}.s2.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_02", f"{PATTERN_R2}.s2.out.fastq"),
        primers = os.path.join(CONPATH, "nebnext_adapters.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_03", f"{PATTERN_R1}.s3.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_03", f"{PATTERN_R2}.s3.out.fastq")),
        stats = os.path.join(STATS, "step_03", "{sample}.s3.stats.tsv")
    benchmark:
        os.path.join(BENCH, "remove_primer_free_adapter.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_primer_free_adapter.{sample}.log")
    resources:
        mem_mb=BBToolsMem
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
    """Step 04: Remove adapter free primer (both orientations). 
    
    Rarely the primer is detected without the primer. This removes those instances as well as everything to the right 
    of the detected adapter-free primer. 
    """
    input:
        r1 = os.path.join(TMPDIR, "step_03", f"{PATTERN_R1}.s3.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_03", f"{PATTERN_R2}.s3.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_04", f"{PATTERN_R1}.s4.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_04", f"{PATTERN_R2}.s4.out.fastq")),
        stats = os.path.join(STATS, "step_04", "{sample}.s4.stats.tsv")
    benchmark:
        os.path.join(BENCH, "remove_adapter_free_primer.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_adapter_free_primer.{sample}.log")
    resources:
        mem_mb=BBToolsMem
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
    """Step 05: Vector contamination removal (PhiX + NCBI UniVecDB)"""
    input:
        r1 = os.path.join(TMPDIR, "step_04", f"{PATTERN_R1}.s4.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_04", f"{PATTERN_R2}.s4.out.fastq"),
        primers = os.path.join(CONPATH, "vector_contaminants.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_05", f"{PATTERN_R1}.s5.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_05", f"{PATTERN_R2}.s5.out.fastq")),
        stats = os.path.join(STATS, "step_05", "{sample}.s5.stats.tsv")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s05.remove_vector_contamination_{sample}.txt")
    log:
        log = os.path.join(STDERR, "step_05", "s5_{sample}.log")
    resources:
        mem_mb=BBToolsMem
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
    """Step 06: Remove remaining low-quality bases and short reads. 
    
    Quality score can be modified in config.yaml (QSCORE).
    """
    input:
        r1 = os.path.join(TMPDIR, "step_05", f"{PATTERN_R1}.s5.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_05", f"{PATTERN_R2}.s5.out.fastq")
    output:
        r1 = temp(os.path.join(TMPDIR, f"{PATTERN_R1}.clean.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, f"{PATTERN_R2}.clean.out.fastq")),
        stats = os.path.join(STATS, "step_06", "{sample}.s6.stats.tsv")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s06.remove_low_quality_{sample}.txt")
    log:
        log = os.path.join(STDERR, "step_06", "s6_{sample}.log")
    resources:
        mem_mb=BBToolsMem
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
        # os.path.join(CONPATH, "line_sine.fasta") ########### TODO: check implementation
    output:
        HOSTINDEX
    benchmark:
        os.path.join(BENCH, "create_host_index.txt")
    log:
        os.path.join(STDERR, 'create_host_index.log')
    resources:
        mem_mb=BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/minimap2.yaml"
    shell:
        "minimap2 -d {output} <(cat {input})"

rule host_removal_mapping:
    """Step 07a: Host removal: mapping to host. 
    
    Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    If your reference is not available you need to add it using 'Hecatomb addHost'
    """
    input:
        r1 = os.path.join(TMPDIR, f"{PATTERN_R1}.clean.out.fastq"),
        r2 = os.path.join(TMPDIR, f"{PATTERN_R2}.clean.out.fastq"),
        host = HOSTINDEX
    output:
        r1=temp(os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.unmapped.fastq")),
        r2=temp(os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R2}.unmapped.fastq")),
        singletons=temp(os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.unmapped.singletons.fastq"))
    benchmark:
        os.path.join(BENCH, "host_removal_mapping.{sample}.txt")
    log:
        mm=os.path.join(STDERR, "host_removal_mapping.{sample}.minimap.log"),
        sv=os.path.join(STDERR, "host_removal_mapping.{sample}.samtoolsView.log"),
        fq=os.path.join(STDERR, "host_removal_mapping.{sample}.samtoolsFastq.log")
    resources:
        mem_mb=BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax sr -t {threads} --secondary=no {input.host} {input.r1} {input.r2} 2> {log.mm} \
            | samtools view -f 4 -h 2> {log.sv} \
            | samtools fastq -NO -1 {output.r1} -2 {output.r2} -0 /dev/null -s {output.singletons} 2> {log.fq}
        """

# rule extract_host_unmapped:
#     """Step 07b: Extract unmapped (non-host) sequences from sam files"""
#     input:
#         sam = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".sam")
#     output:
#         r1 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.fastq")),
#         r2 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".unmapped.fastq")),
#         singletons = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.singletons.fastq"))
#     benchmark:
#         os.path.join(BENCH, "PREPROCESSING", "s07b.extract_host_unmapped_{sample}.txt")
#     log:
#         log = os.path.join(LOGS, "step_07b", "s07b_{sample}.log")
#
#     resources:
#         mem_mb=100000,
#         cpus=64
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """
#         samtools fastq --threads {resources.cpus} -NO -1 {output.r1} -2 {output.r2} \
#         -0 /dev/null \
#         -s {output.singletons} \
#         {input.sam} 2> {log};
#         """

rule nonhost_read_repair:
    """Step 07b: Parse R1/R2 singletons (if singletons at all)"""
    input:
        singletons = os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.unmapped.singletons.fastq")
    output:
        r1 = temp(os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.singletons.fastq")),
        r2 = temp(os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R2}.singletons.fastq"))
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s07c.nonhost_read_repair_{sample}.txt")
    log:
        log = os.path.join(STDERR, "step_07c", "s07c_{sample}.log")
    resources:
        mem_mb=BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        reformat.sh in={input.singletons} out={output.r1} out2={output.r2} \
            -Xmx{resources.mem_mb}m 2> {log};
        """

rule nonhost_read_combine:
    """Step 07c: Combine paired and singleton reads."""
    input:
        r1 = os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.unmapped.fastq"),
        r2 = os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R2}.unmapped.fastq"),
        r1s = os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.singletons.fastq"),
        r2s = os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R2}.singletons.fastq")
    output:
        r1 = os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.all.fastq"),
        r2 = os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R2}.all.fastq")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s07d.nonhost_read_combine_{sample}.txt")
    shell:
        """
        cat {input.r1} {input.r1s} > {output.r1};
        cat {input.r2} {input.r2s} > {output.r2};
        """

rule remove_exact_dups:
    """Step 08: Remove exact duplicates
    
    Exact duplicates are considered PCR generated and not accounted for in the count table (seqtable_all.tsv)
    """
    input:
        os.path.join(QC, "HOST_REMOVED", f"{PATTERN_R1}.all.fastq")
    output:
        temp(os.path.join(QC, "CLUSTERED", f"{PATTERN_R1}.deduped.out.fastq"))
    benchmark:
        os.path.join(BENCH, "remove_exact_dups.{sample}.txt")
    log:
        os.path.join(STDERR, "remove_exact_dups.{sample}.log")
    resources:
        mem_mb=BBToolsMem
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
          
rule cluster_similar_sequences:
    """Step 09: Cluster similar sequences.
     
     Sequences clustered at CLUSTERID in config.yaml.
    """
    input:
        os.path.join(QC, "CLUSTERED", f"{PATTERN_R1}.deduped.out.fastq")
    output:
        temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_rep_seq.fasta")),
        temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_cluster.tsv")),
        temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_all_seqs.fasta"))
    params:
        respath=os.path.join(QC, "CLUSTERED", "LINCLUST"),
        tmppath=os.path.join(QC, "CLUSTERED", "LINCLUST", "{sample}_TMP"),
        prefix=PATTERN_R1
    benchmark:
        os.path.join(BENCH, "cluster_similar_sequences.{sample}.txt")
    log:
        os.path.join(STDERR, "cluster_similar_sequences.{sample}.log")
    resources:
        mem_mb=MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """ 
        mmseqs easy-linclust {input} {params.respath}/{params.prefix} {params.tmppath} \
            --kmer-per-seq-scale 0.3 \
            -c {config[CLUSTERID]} --cov-mode 1 --threads {threads} &>> {log};
        """
        
rule create_individual_seqtables:
    """Step 10: Create individual seqtables. 
    
    A seqtable is a count table with each sequence as a row, each column as a sample and each cell the counts of each 
    sequence per sample.
    """
    input:
        seqs = os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_rep_seq.fasta"),
        counts = os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}_cluster.tsv")
    output:
        seqs = temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}.seqs")),
        counts = temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}.counts")),
        seqtable = temp(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{PATTERN_R1}.seqtable"))
    benchmark:
        os.path.join(BENCH, "create_individual_seqtables.{sample}.txt")
    resources:
        mem_mb=BBToolsMem
    threads:
        BBToolsCPU
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit sort {input.seqs} --quiet -j {threads} -w 5000 -t dna \
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
        paste {output.seqs} {output.counts} > {output.seqtable};
        """

# rule merge_individual_seqtables:
#     """
#
#     Step 11: Merge individual sequence tables into combined seqtable
#
#     """
#     input:
#         files = expand(os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + ".seqtable"), sample=SAMPLES)
#     output:
#         seqtable = os.path.join(RESULTS, "seqtable_all.tsv"),
#         tab2fx = temporary(os.path.join(RESULTS, "seqtable.tab2fx"))
#     params:
#         resultsdir = directory(RESULTS),
#     benchmark:
#         os.path.join(BENCH, "PREPROCESSING", "s11.merge_seq_tables.txt")
#     log:
#         log = os.path.join(LOGS, "step_11", "s11.log")
#     resources:
#         mem_mb=100000,
#         cpus=64
#     params:
#         resultsdir = directory(RESULTS),
#     conda:
#         "../envs/R.yaml"
#     script:
#         "../scripts/seqtable_merge.R"
#
# rule convert_seqtable_tab_to_fasta:
#     """
#
#     Step 12: Convert tabular seqtable output to fasta and create index
#
#     """
#     input:
#         os.path.join(RESULTS, "seqtable.tab2fx")
#     output:
#         seqtable = os.path.join(RESULTS, "seqtable.fasta"),
#         stats = os.path.join(RESULTS, "seqtable.stats"),
#         idx = os.path.join(RESULTS, "seqtable.fasta.fai")
#     benchmark:
#         os.path.join(BENCH, "PREPROCESSING", "s12.convert_seqtable_tab_2_fasta.txt")
#     log:
#         log = os.path.join(LOGS, "step_12", "s12.log")
#     resources:
#         mem_mb=100000,
#         cpus=64
#     conda:
#         "../envs/samtools.yaml"
#     shell:
#         """
#         # Convert
#         seqkit tab2fx {input} -j {config[System][Threads]} -w 5000 -t dna -o {output.seqtable};
#
#         # Calculate seqtable statistics
#         seqkit stats {output.seqtable} -j {config[System][Threads]} -T > {output.stats};
#
#         # Create seqtable index
#         samtools faidx {output.seqtable} -o {output.idx};
#
#         """

rule merge_seq_table:
    """Step 11: Merge seq tables
    
    Reads the sequences and counts from each samples' seqtable text file and converts to fasta format for the rest of 
    the pipline.
    """
    input:
        expand(os.path.join(QC, "CLUSTERED", "LINCLUST", "{sample}_R1.seqtable"), sample=SAMPLES)
    output:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        tsv = os.path.join(RESULTS, "sampleSeqCounts.tsv")
    params:
        resultsdir = directory(RESULTS),
    benchmark:
        os.path.join(BENCH, "merge_seq_table.txt")
    run:
        outFa = open(output.fa, 'w')
        outTsv = open(output.tsv, 'w')
        for sample in SAMPLES:
            seqId = 0
            seqCounts = 0
            counts = open(os.path.join(QC, "CLUSTERED", "LINCLUST", f"{sample}_R1.seqtable"), 'r')
            line = counts.readline() # skip header
            for line in counts:
                l = line.split()
                id = ':'.join((sample, l[1], str(seqId))) # fasta header = >sample:count:seqId
                seqCounts += int(l[1])
                seqId = seqId + 1
                outFa.write(f'>{id}\n{l[0]}\n')
            counts.close()
            outTsv.write(f'{sample}\t{seqCounts}\n')
        outFa.close()
        outTsv.close()

# rule create_seqtable_index:
#     """Step 12: Index seqtable.fasta for rapid access with samtools faidx."""
#     input:
#         os.path.join(RESULTS, "seqtable.fasta")
#     output:
#         os.path.join(RESULTS, "seqtable.fasta.fai")
#     conda:
#         "../envs/samtools.yaml"
#     benchmark:
#         os.path.join(BENCH, 'create_seqtable_index.txt')
#     log:
#         os.path.join(STDERR, 'create_seqtable_index.log')
#     resources:
#         mem_mb=8000
#     shell:
#         """
#         samtools faidx {input} -o {output} 2> {log}
#         """

# rule calculate_seqtable_sequence_properties:
#     """Step 13: Calculate additional sequence properties (ie. GC-content, tetramer frequencies) per sequence"""
#     input:
#         os.path.join(RESULTS, "seqtable.fasta")
#     output:
#         gc = os.path.join(RESULTS, "seqtable_properties.gc"),
#         tetramer = os.path.join(RESULTS, "seqtable_properties.tetramer"),
#         seq_properties = os.path.join(RESULTS, "seqtable_properties.tsv")
#     benchmark:
#         os.path.join(BENCH, "s13.calculate_seqtable_sequence_properties.txt")
#     log:
#         log1 = os.path.join(STDERR, "step_13", "s13.gc.log"),
#         log2 = os.path.join(STDERR, "step_13", "s13.tetramer.log")
#     resources:
#         mem_mb=100000,
#         cpus=64
#     conda:
#         "../envs/bbmap.yaml"
#     shell:
#         """
#         # Calcualate per sequence GC content
#         countgc.sh in={input} format=2 ow=t | awk 'NF' > {output.gc};
#         sed -i '1i id\tGC' {output.gc} 2> {log.log1};
#
#         # Calculate per sequence tetramer frequency
#         tetramerfreq.sh in={input} w=0 ow=t | \
#         tail -n+2 | \
#         cut --complement -f2 > {output.tetramer} 2> {log.log2};
#
#         sed -i 's/scaffold/id/' {output.tetramer};
#
#         # Combine
#         csvtk join -f 1 {output.gc} {output.tetramer} -t -T > {output.seq_properties};
#         """


