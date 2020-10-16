"""
This is drop in replacement for 01_contaminant removal that uses the same initial input files and generates
the same output files. It will probably break snakemake if you include them both!

"""
import os
import sys


rule bowtie2_host_removal:
    """
    Step 6. Remove host and LINE/SINES
    These steps have been converted to bowtie2 from bbmap as it is more efficient
    """
    input:
        r1 = os.path.join(QC, "step_5", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(QC, "step_5", PATTERN_R2 + ".s5.out.fastq"),
    output:
        os.path.join(QC, "step_6", '{sample}.host.bam')
    benchmark:
        "benchmarks/bowtie2_host_removal_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    params:
        idx = HOSTBT2
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -p {resources.cpus} -x {params.idx} -1 {input.r1} -2 {input.r2} | \
        samtools view --threads {resources.cpus} -bh | \
        samtools sort --threads {resources.cpus} -o {output} -
        """

rule reads_host_mapped:
    """
    Note, in deconseq.snakefile these are all separate rules,
    but I merge them here so that they run on the same cluster node
    and it will reduce the size of the dag
    """
    input:
        os.path.join(QC, "step_6", '{sample}.host.bam')
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + '_host.mapped.fastq'),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + '_host.mapped.fastq'),
        s = os.path.join(QC, "step_6", '{sample}_singletons_host.mapped.fastq')
    benchmark:
        "benchmarks/reads_host_mapped_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -G 12 -f 65 {input} > {output.r1} &&
        samtools fastq --threads {resources.cpus} -G 12 -f 129 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus} -F 5 {input} > {output.s}
        """

rule reads_host_unmapped:
    input:
        os.path.join(QC, "step_6", '{sample}.host.bam')
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + '_host.unmapped.fastq'),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + '_host.unmapped.fastq'),
        s = os.path.join(QC, "step_6", '{sample}_singletons_host.unmapped.fastq')
    benchmark:
        "benchmarks/reads_host_unmapped_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -f 77  {input} > {output.r1} && 
        samtools fastq --threads {resources.cpus} -f 141 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus} -f 4 -F 1  {input} > {output.s}
        """

rule line_sine_bam:
    input:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + '_host.unmapped.fastq'),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + '_host.unmapped.fastq'),
        s = os.path.join(QC, "step_6", '{sample}_singletons_host.unmapped.fastq')
    output:
        os.path.join(QC, "step_6", '{sample}.linesine.bam')
    benchmark:
        "benchmarks/line_sine_bam_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    params:
        idx = os.path.join(CONPATH, "line_sine")
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -p {resources.cpus} -x {params.idx} -1 {input.r1} -2 {input.r2} \
        -U {input.s} | samtools view --threads {resources.cpus} -bh | \
        samtools sort --threads {resources.cpus} -o {output} -
        """

rule reads_linesine_mapped:
    input:
        os.path.join(QC, "step_6", '{sample}.linesine.bam')
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + '_linesine.mapped.fastq'),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + '_linesine.mapped.fastq'),
        s = os.path.join(QC, "step_6", '{sample}_singletons_linesine.mapped.fastq')
    benchmark:
        "benchmarks/reads_linesine_mapped_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -G 12 -f 65 {input} > {output.r1} &&
        samtools fastq --threads {resources.cpus} -G 12 -f 129 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus} -F 5 {input} > {output.s}
        """

rule reads_linesine_unmapped:
    input:
        os.path.join(QC, "step_6", '{sample}.linesine.bam')
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq"),
        s = os.path.join(QC, "step_6", '{sample}_singletons.s6.out.fastq')
    benchmark:
        "benchmarks/reads_linesine_unmapped_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus}  -f 77  {input} > {output.r1} && 
        samtools fastq --threads {resources.cpus}  -f 141 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus}  -f 4 -F 1  {input} > {output.s}
        """

rule remove_bacteria:
    """
    Step 8: Remove bacterial contaminants reserving viral and ambiguous sequences
    """
    input:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq"),
        s = os.path.join(QC, "step_6", '{sample}_singletons.s6.out.fastq')
    output:
        os.path.join(QC, "step_7", '{sample}.bacteria.bam')
    params:
        idx = BACBT2
    benchmark:
        "benchmarks/remove_bacteria_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bowtie2.yaml"
    shell:
        """
        bowtie2 -p {resources.cpus} -x {params.idx} -1 {input.r1} -2 {input.r2} \
        -U {input.s} | samtools view --threads {resources.cpus} -bh | \
        samtools sort --threads {resources.cpus} -o {output} -
        """

rule reads_bacteria_mapped:
    input:
        os.path.join(QC, "step_7", '{sample}.bacteria.bam')
    output:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + '.bacteria_mapped.fastq'),
        r2 = os.path.join(QC, "step_7", PATTERN_R2 + '.bacteria_mapped.fastq'),
        s1 = os.path.join(QC, "step_7", PATTERN_R1 + '.bacteria_mapped_singletons.fastq'),
        s2 = os.path.join(QC, "step_7", PATTERN_R2 + '.bacteria_mapped_singletons.fastq')
    benchmark:
        "benchmarks/reads_bacteria_mapped_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -G 12 -f 65 {input} > {output.r1} &&
        samtools fastq --threads {resources.cpus}  -G 12 -f 129 {input} > {output.r2} &&
        samtools fastq  --threads {resources.cpus} -f 72 -F 4 {input} > {output.s1} &&
        samtools fastq  --threads {resources.cpus} -f 136 -F 4 {input} > {output.s2}
        """

rule reads_bacteria_unmapped:
    input:
        os.path.join(QC, "step_7", '{sample}.bacteria.bam')
    output:
        r1 = os.path.join(QC, "step_8", PATTERN_R1 + ".bacteria_unmapped.fastq"),
        r2 = os.path.join(QC, "step_8", PATTERN_R2 + ".bacteria_unmapped.fastq"),
        s1 = os.path.join(QC, "step_8", PATTERN_R1 + ".bacteria_unmapped_singletons.fastq"),
        s2 = os.path.join(QC, "step_8", PATTERN_R2 + ".bacteria_unmapped_singletons.fastq")
    benchmark:
        "benchmarks/reads_bacteria_unmapped_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -f 77  {input} > {output.r1} && 
        samtools fastq --threads {resources.cpus} -f 141 {input} > {output.r2} &&
        samtools fastq --threads {resources.cpus} -f 68 -F 8  {input} > {output.s1} &&
        samtools fastq --threads {resources.cpus} -f 132 -F 8  {input} > {output.s2}
        """

rule concatenate_singletons:
    """
    Combine the unmapped reads and singletons into a single file for the next step
    """
    input:
        r1 = os.path.join(QC, "step_8", PATTERN_R1 + ".bacteria_unmapped.fastq"),
        s1 = os.path.join(QC, "step_8", PATTERN_R1 + ".bacteria_unmapped_singletons.fastq")
    output:
        os.path.join(QC, "step_9", PATTERN_R1 + ".viral_amb.fastq")
    shell:
        "cat {input} > {output}"