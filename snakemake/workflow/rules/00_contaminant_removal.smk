"""

Snakefile based on [contaminant_removal.sh](../base/contaminant_removal.sh)

Rob Edwards, Jan 2020
Updated: Scott Handley, Nov 2020
"""

import os

# NOTE: bbtools uses "threads=auto" by default that typically uses all threads, so no need to specify. 
# -Xmx is used to specify the memory allocation for bbtools operations
# Set your -Xmx specifications in your configuration file 

rule remove_leftmost_primerB:
    """
    Step 01: Remove leftmost primerB
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
        "benchmarks/removeprimerB_{sample}.txt"
    log:
        "LOGS/step_01/{sample}.s1.log"
    resources:
        mem_mb=20000,
        cpus=8
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
            -Xmx{resources.mem_mb}m 2> {log}
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
        "benchmarks/remove_3prime_contaminant_{sample}.txt"
    log:
        "LOGS/step_02/{sample}.s2.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=r removeifeitherbad=f ordered=t rcomp=f ow=t \
            -Xmx{resources.mem_mb}m 2> {log}
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
        "benchmarks/remove_primer_free_adapter_{sample}.txt"
    log:
        "LOGS/step_03/{sample}.s3.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=10 ktrim=r removeifeitherbad=f ordered=t rcomp=t ow=t \
            -Xmx{resources.mem_mb}m 2> {log}
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
        "benchmarks/remove_adapter_free_primer_{sample}.txt"
    log:
        "LOGS/step_04/{sample}.s4.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=0 removeifeitherbad=f ordered=t rcomp=t ow=t \
            -Xmx{resources.mem_mb}m 2> {log}
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
        "benchmarks/remove_vector_contamination_{sample}.txt"
    log:
        "LOGS/step_05/{sample}.s5.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t \
            -Xmx{resources.mem_mb}m 2> {log}
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
        "benchmarks/remove_low_quality_{sample}.txt"
    log:
        "LOGS/step_06/{sample}.s6.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            ordered=t \
            qtrim=r maxns=2 \
            trimq={config[QSCORE]} \
            minlength={config[MINLENGTH]} 2> {log} 
        """

rule host_mapping:
    """
    Step 07a: Host removal. Must define host in config file (see Host: ). Host should be masked of viral sequence
    """
    input:
        r1 = os.path.join(TMPDIR, PATTERN_R1 + ".clean.out.fastq"),
        r2 = os.path.join(TMPDIR, PATTERN_R2 + ".clean.out.fastq"),
        hostpath = HOSTPATH
    output:
        sam = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".sam")
    benchmark:
        "benchmarks/host_mapping_{sample}.txt"
    log:
        "LOGS/host_removal/{sample}.host_removal.log"
    resources:
        mem_mb=100000,
        cpus=8
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        minimap2 -ax sr -t 64 {input.hostpath} {input.r1} {input.r2} 2> {log} | \
        samtools view -F 2048 -h | \
        samtools view -f 4 -h > {output.sam};
        """

rule nonhost_read_parsing:
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
        "benchmarks/nonhost_read_parsing_{sample}.txt"
    resources:
        mem_mb=100000,
        cpus=8
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        samtools fastq -NO -1 {output.r1} -2 {output.r2} \
        -0 /dev/null \
        -s {output.singletons} \
        {input.sam}
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
        "benchmarks/nonhost_read_repair_{sample}.txt"
    resources:
        mem_mb=100000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        reformat.sh in={input.singletons} out={output.r1} out2={output.r2}
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
        "benchmarks/singleton_read_parsing_{sample}.txt"
    resources:
        mem_mb=100000,
        cpus=8
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
        "benchmarks/remove_exact_dups_{sample}.txt"
    log:
        "LOGS/clustering/{sample}.dedupe.log"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        dedupe.sh in={input} \
                out={output} \
                ac=f  ow=t \
                -Xmx{resources.mem_mb}m 2> {log}
        """
          
rule cluster_similar_sequences:
    """
    Step 09: Cluster similar sequences
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
        "benchmark/cluster_similar_seqs_{sample}.txt"
    log:
        "LOGS/clustering/{sample}.linclust.log"
    resources:
        mem_mb=64000,
        cpus=24
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """ 
        mmseqs easy-linclust {input} {params.respath}/{params.prefix} {params.tmppath} \
        --kmer-per-seq-scale 0.3 \
        -c 0.95 --cov-mode 1 --threads {resources.cpus} &>> {log}
        """
    
