"""
Snakemake rule file to preprocess Illumina sequence data in preperation for taxonomic assignment.

What is accomplished with this script?
    - Non-biological sequence removal (primers, adapters)

Additional Reading:
    - Hecatomb GitHub: https://github.com/shandley/hecatomb
    - Official Snakemake documentation: https://snakemake.readthedocs.io/en/stable/

Historical Notes:
- Updated Snakefile based on [contaminant_removal.sh](../base/contaminant_removal.sh)

Rob Edwards, Jan 2020
Updated: Scott Handley, Jan 2021


"""

import os
import sys

# NOTE: bbtools uses "threads=auto" by default that typically uses all threads, so no need to specify. 
# -Xmx is used to specify the memory allocation for bbtools operations
# Set your -Xmx specifications in your configuration file 

rule remove_leftmost_primerB:
    """
    
    Step 01: Remove leftmost primer.
    
    Primer sequences used in the Handley lab are included (primerB.fa). If your lab uses other primers you will need to place them in CONPATH (defined in the Snakefile) and change the file name from primerB.fa to your file name below.
    
    """
    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + file_extension),
        r2 = os.path.join(READDIR, PATTERN_R2 + file_extension),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_01", PATTERN_R1 + ".s1.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_01", PATTERN_R2 + ".s1.out.fastq")),
        stats = os.path.join(STATS, "step_01", "{sample}.s1.stats.tsv")
    benchmark:
        "BENCHMARKS/preprocessing/step_01/removeprimerB_{sample}.txt"
    log:
        "LOGS/step_01/{sample}.s1.log"
    resources:
        mem_mb=100000,
        cpus=64
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
            -Xmx{config[System][Memory]}g 2> {log}
        """

rule remove_3prime_contaminant:
    """
    
    Step 02: Remove 3' read through contaminant
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_01", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_01", PATTERN_R2 + ".s1.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_02", PATTERN_R1 + ".s2.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_02", PATTERN_R2 + ".s2.out.fastq")),
        stats = os.path.join(STATS, "step_02", "{sample}.s2.stats.tsv")
    benchmark:
        "BENCHMARKS/preprocessing/step_02/remove_3prime_contaminant_{sample}.txt"
    log:
        "LOGS/step_02/{sample}.s2.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=r removeifeitherbad=f ordered=t rcomp=f ow=t \
            -Xmx{config[System][Memory]}g 2> {log}
        """

rule remove_primer_free_adapter:
    """
    
    Step 03: Remove primer free adapter (both orientations)
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_02", PATTERN_R1 + ".s2.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_02", PATTERN_R2 + ".s2.out.fastq"),
        primers = os.path.join(CONPATH, "nebnext_adapters.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_03", PATTERN_R1 + ".s3.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_03", PATTERN_R2 + ".s3.out.fastq")),
        stats = os.path.join(STATS, "step_03", "{sample}.s3.stats.tsv")
    benchmark:
        "BENCHMARKS/preprocessing/step_03/remove_primer_free_adapter_{sample}.txt"
    log:
        "LOGS/step_03/{sample}.s3.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=10 ktrim=r removeifeitherbad=f ordered=t rcomp=t ow=t \
            -Xmx{config[System][Memory]}g 2> {log}
        """

rule remove_adapter_free_primer:
    """
    
    Step 04: Remove adapter free primer (both orientations)
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_03", PATTERN_R1 + ".s3.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_03", PATTERN_R2 + ".s3.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_04", PATTERN_R1 + ".s4.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_04", PATTERN_R2 + ".s4.out.fastq")),
        stats = os.path.join(STATS, "step_04", "{sample}.s4.stats.tsv")
    benchmark:
        "BENCHMARKS/preprocessing/step_04/remove_adapter_free_primer_{sample}.txt"
    log:
        "LOGS/step_04/{sample}.s4.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=0 removeifeitherbad=f ordered=t rcomp=t ow=t \
            -Xmx{config[System][Memory]}g 2> {log}
        """

rule remove_vector_contamination:
    """
    
    Step 05: Vector contamination removal (PhiX + NCBI UniVecDB)
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_04", PATTERN_R1 + ".s4.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_04", PATTERN_R2 + ".s4.out.fastq"),
        primers = os.path.join(CONPATH, "vector_contaminants.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_05", PATTERN_R1 + ".s5.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_05", PATTERN_R2 + ".s5.out.fastq")),
        stats = os.path.join(STATS, "step_05", "{sample}.s5.stats.tsv")
    benchmark:
        "BENCHMARKS/preprocessing/step_05/remove_vector_contamination_{sample}.txt"
    log:
        "LOGS/step_05/{sample}.s5.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t \
            -Xmx{config[System][Memory]}g 2> {log}
        """
        
rule remove_low_quality:
    """
    
    Step 06: Remove remaining low-quality bases and short reads
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_05", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_05", PATTERN_R2 + ".s5.out.fastq")
    output:
        r1 = temp(os.path.join(TMPDIR, PATTERN_R1 + ".clean.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, PATTERN_R2 + ".clean.out.fastq")),
        stats = os.path.join(STATS, "step_06", "{sample}.s6.stats.tsv")
    benchmark:
        "BENCHMARKS/preprocessing/step_06/remove_low_quality_{sample}.txt"
    log:
        "LOGS/step_06/{sample}.s6.log"
    resources:
        mem_mb=100000,
        cpus=64
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
            trimq={config[QSCORE]} \
            minlength={config[MINLENGTH]} \
            -Xmx{config[System][Memory]}g 2> {log} 
        """

rule host_removal_mapping:
    """
    
    Step 07a: Host removal. Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    
    Step 7 has several substeps (7a, 7b, 7c and 7d) to subset and fix read pairing issues
    
    If your reference is not available post an issue on GitHub requesting it to be added (https://github.com/shandley/hecatomb)
    
    """
    input:
        r1 = os.path.join(TMPDIR, PATTERN_R1 + ".clean.out.fastq"),
        r2 = os.path.join(TMPDIR, PATTERN_R2 + ".clean.out.fastq"),
        hostpath = HOSTPATH
    output:
        sam = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".sam")
    benchmark:
        "BENCHMARKS/preprocessing/step_07/host_removal_{sample}.txt"
    log:
        "LOGS/host_removal/{sample}.host_removal.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax sr -t {config[System][Threads]} {input.hostpath} {input.r1} {input.r2} 2> {log} | \
        samtools view -F 2048 -h | \
        samtools view -f 4 -h > {output.sam} 2> {log};
        """

rule extract_host_unmapped:
    """
    
    Step 07b: Extract unmapped fastq files from sam files
    
    """
    input:
        sam = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".sam")
    output:
        r1 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.fastq")),
        r2 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".unmapped.fastq")),
        singletons = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.singletons.fastq"))
    benchmark:
        "BENCHMARKS/preprocessing/step_07/extract_host_unmapped_{sample}.txt"
    log:
        "LOGS/host_removal/{sample}.extract_host_unmapped.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        samtools fastq -NO -1 {output.r1} -2 {output.r2} \
        -0 /dev/null \
        -s {output.singletons} \
        {input.sam} 2> {log}
        """

rule nonhost_read_repair:
    """
    
    Step 07c: Parse R1/R2 singletons (if singletons at all)
    
    """
    input:
        singletons = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.singletons.fastq")
    output:
        r1 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq")),
        r2 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq"))
    benchmark:
        "BENCHMARKS/preprocessing/step_07/nonhost_read_repair_{sample}.txt"
    log:
        "LOGS/host_removal/{sample}.nonhost_read_repair.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        reformat.sh in={input.singletons} out={output.r1} out2={output.r2} \
        -Xmx{config[System][Memory]}g 2> {log}
        """

rule nonhost_read_combine:
    """
    
    Step 07d: Combine R1+R1_singletons and R2+R2_singletons
    
    """
    input:
        r1 = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.fastq"),
        r2 = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".unmapped.fastq"),
        r1s = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq"),
        r2s = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq")
    output:
        r1 = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".all.fastq"),
        r2 = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".all.fastq")
    benchmark:
        "BENCHMARKS/preprocessing/step_07/nonhost_read_combine_{sample}.txt"
    resources:
        mem_mb=100000,
        cpus=64
    shell:
        """
        cat {input.r1} {input.r1s} > {output.r1};
        cat {input.r2} {input.r2s} > {output.r2}
        """

rule remove_exact_dups:
    """
    
    Step 08: Remove exact duplicates
    
    """
    input:
        os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".all.fastq")
    output:
        os.path.join(QC, "CLUSTERED", PATTERN_R1 + ".deduped.out.fastq")
    benchmark:
        "BENCHMARKS/preprocessing/step_08/remove_exact_dups_{sample}.txt"
    log:
        "LOGS/clustering/{sample}.dedupe.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        dedupe.sh in={input} out={output} \
        ac=f ow=t \
        -Xmx{config[System][Memory]}g 2> {log}
        """
          
rule cluster_similar_sequences:
    """
    
    Step 09: Cluster similar sequences at CLUSTERID in config.yaml. Default: 97% identity
    
    """
    input:
        os.path.join(QC, "CLUSTERED", PATTERN_R1 + ".deduped.out.fastq")
    output:
        os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + "_rep_seq.fasta"),
        os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + "_cluster.tsv"),
        temporary(os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + "_all_seqs.fasta"))
    params:
        respath=os.path.join(QC, "CLUSTERED", "LINCLUST"),
        tmppath=os.path.join(QC, "CLUSTERED", "LINCLUST", "TMP"),
        prefix=PATTERN_R1
    benchmark:
        "BENCHMARKS/preprocessing/step_08/cluster_similar_seqs_{sample}.txt"
    log:
        "LOGS/clustering/{sample}.linclust.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """ 
        mmseqs easy-linclust {input} {params.respath}/{params.prefix} {params.tmppath} \
        --kmer-per-seq-scale 0.3 \
        -c 0.95 --cov-mode 1 --threads {config[System][Threads]} &>> {log}
        """
        
rule create_individual_seqtables:
    """
    
    Step 10: Create individual seqtables. A seqtable is a count table with each feature (sequence) as a row, each column as a sample and each cell the counts of each sequence per sample
    
    """
    input:
        seqs=os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + "_rep_seq.fasta"),
        counts=os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + "_cluster.tsv")
    output:
        seqs=os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + ".seqs"),
        counts=os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + ".counts"),
        seqtable=os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + ".seqtable")
    params:
        prefix=PATTERN
    benchmark:
        "BENCHMARKS/preprocessing/step_10/individual_seq_tables_{sample}.txt"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit sort {input.seqs} --quiet -j {config[System][Threads]} -w 5000 -t dna | \
        seqkit fx2tab -j {config[System][Threads]} -w 5000 -t dna | \
        sed 's/\\t\\+$//' | \
        cut -f2,3 | \
        sed '1i sequence' > {output.seqs}
        
        cut -f1 {input.counts} | \
        sort | \
        uniq -c | \
        awk -F ' ' '{{print$2"\\t"$1}}' | \
        cut -f2 | \
        sed "1i {params.prefix}" > {output.counts}
        
        paste {output.seqs} {output.counts} > {output.seqtable}
        """
        
rule merge_individual_seqtables:
    """
    
    Step 11: Merge individual sequence tables into combined seqtable
    
    """
    input:
        files = expand(os.path.join(QC, "CLUSTERED", "LINCLUST", PATTERN_R1 + ".seqtable"), sample=SAMPLES)
    output:
        seqtable = os.path.join(RESULTS, "seqtable_all.tsv"),
        tab2fx = temporary(os.path.join(RESULTS, "seqtable.tab2fx"))
    params:
        resultsdir = directory(RESULTS),
    benchmark:
        "BENCHMARKS/preprocessing/step_10/merge_seq_table.txt"
    resources:
        mem_mb=100000,
        cpus=64
    params:
        resultsdir = directory(RESULTS),
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/seqtable_merge.R"

rule convert_seqtable_tab_2_fasta:
    """
    Step 12: Convert tabular seqtable output to fasta
    """
    input:
        os.path.join(RESULTS, "seqtable.tab2fx")
    output:
        os.path.join(RESULTS, "seqtable.fasta")
    benchmark:
        "BENCHMARKS/preprocessing/step_10/convert_seqtable_tab_2_fasta.txt"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit tab2fx {input} -j {config[System][Threads]} -w 5000 -t dna -o {output}
        """

rule create_seqtable_index:
    """
    
    Step 13: Index seqtable.fasta for rapid samtools access in later steps
    
    """
    input:
        os.path.join(RESULTS, "seqtable.fasta")
    output:
        os.path.join(RESULTS, "seqtable.faidx")
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools faidx {input} -o {output}
        """

rule calculate_seqtable_sequence_properties:
    """
    
    Step 14: Calculate additional sequence properties (ie. GC-content) per sequence
    
    """
    input:
        os.path.join(RESULTS, "seqtable.fasta")
    output:
        os.path.join(RESULTS, "seqtable_properties.tsv")
    benchmark:
        "BENCHMARKS/preprocessing/step_14/calculate_seqtable_sequence_properties.txt"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        seqkit fx2tab -j {config[System][Threads]} --gc -H {input} | cut -f1,4 > {output}
        """

rule assembly_kmer_normalization:
    """
    
    Step 15: Kmer normalization. Data reduction for assembly improvement
    
    """
    input:
        r1 = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.fastq"),
        r2 = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".unmapped.fastq"),
        r1s = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq"),
        r2s = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq")
    output:
        r1_norm = os.path.join(ASSEMBLY, PATTERN_R1 + ".norm.fastq"),
        r2_norm = os.path.join(ASSEMBLY, PATTERN_R2 + ".norm.fastq")
    benchmark:
        "BENCHMARKS/preprocessing/assembly/kmer_normalization_{sample}.txt"
    log:
        "LOGS/assembly/{sample}.kmer_normalization.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbnorm.sh in={input.r1} in2={input.r2} \
        extra={input.r1s},{input.r2s} \
        out={output.r1_norm} out2={output.r2_norm} \
        target=100 \
        ow=t \
        -Xmx{config[System][Memory]}g 2> {log}
    """

rule individual_sample_assembly:
    """
    
    Step 16: Individual sample assemblies
    
    """
    input:
        r1_norm = os.path.join(ASSEMBLY, PATTERN_R1 + ".norm.fastq"),
        r2_norm = os.path.join(ASSEMBLY, PATTERN_R2 + ".norm.fastq"),
        r1s = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq"),
        r2s = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq")
    output:
        contigs=os.path.join(ASSEMBLY, PATTERN, PATTERN + ".contigs.fa")
    params:
        mh_dir=directory(os.path.join(ASSEMBLY, '{sample}')),
        contig_dic=directory(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY"))
    benchmark:
        "BENCHMARKS/preprocessing/assembly/megahit_{sample}.txt"
    log:
        "LOGS/assembly/{sample}.megahit.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        rmdir {params.mh_dir};
        
        megahit -1 {input.r1_norm} -2 {input.r2_norm} -r {input.r1s},{input.r2s} \
        -o {params.mh_dir} --out-prefix {wildcards.sample} \
        --k-min 45 --k-max 225 --k-step 26 --min-count 2 &>> {log}      
        """
        
rule concatenate_contigs:
    """
    
    Step 17: Concatenate individual assembly outputs (contigs) into a single file
    
    """
    input:
        lambda wildcards: expand(os.path.join(ASSEMBLY, PATTERN, PATTERN + ".contigs.fa"), sample=SAMPLES)
    output:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.fasta")
    shell:
        """
        
        cat {input} > {output}
        
        """

rule contig_reformating_and_stats:
    """
    
    Step 18: Remove short contigs (Default: 500). Defined in config[CONTIG_SIZE_THRESH]
    
    """
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.fasta")
    output:
        rename = temporary(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.renamed.fasta")),
        size = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs_size_selected.fasta"),
        stats = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.stats"),
        sketch = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs_size_selected.sketch")
    benchmark:
        "BENCHMARKS/preprocessing/contig_reformating.txt"
    log:
        "LOGS/assembly/contig_reformating.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        rename.sh in={input} out={output.rename} \
        ow=t \
        -Xmx{config[System][Memory]}g;
        
        reformat.sh in={output.rename} out={output.size} \
        ml={config[MINLENGTH]} \
        ow=t \
        -Xmx{config[System][Memory]}g;
        
        statswrapper.sh in={input} out={output.stats} \
        format=2 \
        ow=t;
        
        sendsketch.sh in={output.size} out={output.sketch} \
        address=nt mode=sequence format=3 \
        ow=t \
        -Xmx{config[System][Memory]}g 2> {log};
        """

rule population_assembly:
    """
    
    Step 19: Create 'contig dictionary' of all unique contigs present in the study (population assembly)
    
    """
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs_size_selected.fasta")
    output:
        assembly = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta"),
        stats = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "contig_dictionary.stats"),
        sketch = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "contig_dictionary.sketch")
    params:
        flye_out=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE")
    benchmark:
        "BENCHMARKS/assembly/population_assembly.txt"
    log:
        "LOGS/assembly/population_assembly.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/metaflye.yaml"
    shell:
        """
        flye --subassemblies {input} -t {resources.cpus} --plasmids -o {params.flye_out} -g 1g;
        
        statswrapper.sh in={output.assembly} out={output.stats} \
        format=2 \
        ow=t;
        
        sendsketch.sh in={output.assembly} out={output.sketch} \
        address=nt mode=sequence format=3 \
        ow=t \
        -Xmx{config[System][Memory]}g 2> {log};
        """

rule coverage_calculations:
    """
    
    Step 20a: Calculate contig coverage and extract unmapped reads
    
    """
    input:
        r1 = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".all.fastq"),
        r2 = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".all.fastq"),
        ref = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta")
    output:
        sam=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".aln.sam.gz"),
        unmap=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".umapped.fastq"),
        covstats=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".cov_stats"),
        rpkm=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".rpkm"),
        statsfile=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".statsfile"),
        scafstats=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".scafstats")
    benchmark:
        "BENCHMARKS/assembly/coverage_calculations_{sample}.txt"
    log:
        "LOGS/assembly/coverage_calculations_{sample}.log"
    resources:
        mem_mb=100000,
        cpus=64
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbmap.sh ref={input.ref} in={input.r1} in2={input.r2} \
        nodisk \
        out={output.sam} \
        outu={output.unmap} \
        ambiguous=random \
        physcov=t \
        covstats={output.covstats} \
        rpkm={output.rpkm} \
        statsfile={output.statsfile} \
        scafstats={output.scafstats} \
        maxindel=100 minid=90 \
        ow=t 2> \
        -Xmx{config[System][Memory]}g 2> {log};
        """

rule create_contig_count_table:
    """
    
    Step 20b: Filter low coverage contigs
    
    """
    input:
        rpkm = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".rpkm"),
        covstats=os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + ".cov_stats")
    output:
        counts_tmp = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + "_counts.tmp"),
        TPM_tmp = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + "_TPM.tmp"),
        TPM = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + "_TPM"),
        TPM_final = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + "_TPM.final"),
        cov_temp = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + "_cov.tmp"),
        count_tbl = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + "_contig_counts.tsv")
    shell:
        """
        ## TPM Calculator
        # Prepare table & calculate RPK
        tail -n+6 {input.rpkm} > {output.counts_tmp} | \
        cut -f1,2,5,6,8 | \
        awk 'BEGIN{{ FS=OFS="\t" }} {{ print $0, $3/($2/1000) }}' > {output.counts_tmp}
        
        # Calculate size factor
        sizef=$(awk 'BEGIN{{ total=0 }} {{ total=total+$6 }} END{{ printf total }}' {output.counts_tmp});
        
        # Calculate TPM
        awk -v awkvar="$sizef" 'BEGIN{{ FS=OFS="\t" }} {{ print $0, $6/awkvar }}' < {output.counts_tmp} > {output.TPM_tmp};
        
        # Add sample name
        awk -v awkvar="{wildcards.sample}" 'BEGIN{{FS=OFS="\t"}} {{ print awkvar, $0 }}' < {output.TPM_tmp} > {output.TPM};
        
        # Remove RPK
        cut -f1-6,8 {output.TPM} > {output.TPM_final};
        
        ## Coverage stats modifications
        tail -n+2 {input.covstats} | cut -f2,4,5,6,9 > {output.cov_temp};
        
        ## Combine tables
        paste {output.TPM_final} {output.cov_temp} > {output.count_tbl};
        
        """

rule concatentate_contig_count_tables:
    """
    
    Rule 20c: Concatenate contig count tables
    Note: this is done as a separate rule due to how snakemake handles i/o files. It does not work well in Rule 20b as the i/o PATTERNS are different.
    
    """
    input:
        #lambda wildcards: expand(os.path.join(ASSEMBLY, PATTERN, PATTERN + ".contigs.fa"), sample=SAMPLES)
        lambda wildcards: expand(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", PATTERN + "_contig_counts.tsv"), sample=SAMPLES)
    output:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING",  "contig_count_table.tsv")
    shell:
        """
        
        cat {input} > {output};
        
        sed -i '1i sample_id\tcontig_id\tlength\treads\tRPKM\tFPKM\tTPM\tavg_fold_cov\tcontig_GC\tcov_perc\tcov_bases\tmedian_fold_cov' {output};
        
        """
    
    




