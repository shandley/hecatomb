"""

Snakefile based on [contaminant_removal.sh](../base/contaminant_removal.sh)

Rob Edwards, Jan 2020
"""

import os

# NOTE: bbduk uses "threads=auto" by default that typically uses all threads, so no need to specify. -Xmx is used to
# specify the memory allocation

rule remove_leftmost_primerB:
    """
    Step 1: Remove leftmost primerB. Not the reverse complements
    """
    input:
        r1 = remove_leftmost_primerB_r1_input,
        r2 = remove_leftmost_primerB_r2_input,
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = os.path.join(QC, "step_1", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(QC, "step_1", PATTERN_R2 + ".s1.out.fastq"),
        stats = os.path.join(QC, "step_1", "{sample}.s1.stats.txt")
    benchmark:
        "benchmarks/removeprimerB_{sample}.txt"
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
            -Xmx{resources.mem_mb}m
        """

rule remove_3prime_contaminant:
    """
    Step 2: Remove 3' read through contaminant
    """
    input:
        r1 = os.path.join(QC, "step_1", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(QC, "step_1", PATTERN_R2 + ".s1.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = os.path.join(QC, "step_2", PATTERN_R1 + ".s2.out.fastq"),
        r2 = os.path.join(QC, "step_2", PATTERN_R2 + ".s2.out.fastq"),
        stats = os.path.join(QC, "step_2", "{sample}.s2.stats.txt")
    benchmark:
        "benchmarks/remove_3prime_contaminant_{sample}.txt"
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
            -Xmx{resources.mem_mb}m
        """

rule remove_primer_free_adapter:
    """
    Step 3: Remove primer free adapter (both orientations)
    """
    input:
        r1 = os.path.join(QC, "step_2", PATTERN_R1 + ".s2.out.fastq"),
        r2 = os.path.join(QC, "step_2", PATTERN_R2 + ".s2.out.fastq"),
        primers = os.path.join(CONPATH, "nebnext_adapters.fa")
    output:
        r1 = os.path.join(QC, "step_3", PATTERN_R1 + ".s3.out.fastq"),
        r2 = os.path.join(QC, "step_3", PATTERN_R2 + ".s3.out.fastq"),
        stats = os.path.join(QC, "step_3", "{sample}.s3.stats.txt")
    benchmark:
        "benchmarks/remove_primer_free_adapter_{sample}.txt"
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
            -Xmx{resources.mem_mb}m
        """

rule remove_adapter_free_primer:
    """
    Step 4: Remove adapter free primer (both orientations)
    """
    input:
        r1 = os.path.join(QC, "step_3", PATTERN_R1 + ".s3.out.fastq"),
        r2 = os.path.join(QC, "step_3", PATTERN_R2 + ".s3.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = os.path.join(QC, "step_4", PATTERN_R1 + ".s4.out.fastq"),
        r2 = os.path.join(QC, "step_4", PATTERN_R2 + ".s4.out.fastq"),
        stats = os.path.join(QC, "step_4", "{sample}.s4.stats.txt")
    benchmark:
        "benchmarks/remove_adapter_free_primer_{sample}.txt"
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
            -Xmx{resources.mem_mb}m
        """

rule remove_vector_contamination:
    """
    Step 5: Vector contamination removal (PhiX + NCBI UniVecDB)
    """
    input:
        r1 = os.path.join(QC, "step_4", PATTERN_R1 + ".s4.out.fastq"),
        r2 = os.path.join(QC, "step_4", PATTERN_R2 + ".s4.out.fastq"),
        primers = os.path.join(CONPATH, config['DatabaseFiles']['contaminants'])
    output:
        r1 = os.path.join(QC, "step_5", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(QC, "step_5", PATTERN_R2 + ".s5.out.fastq"),
        stats = os.path.join(QC, "step_5", "{sample}.s5.stats.txt")
    benchmark:
        "benchmarks/remove_vector_contamination_{sample}.txt"
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
            -Xmx{resources.mem_mb}m
        """







