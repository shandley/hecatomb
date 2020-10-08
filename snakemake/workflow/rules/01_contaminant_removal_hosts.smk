"""
This section of contaminant removal has been abstracted as we may replace pieces and parts of it later


"""
import os
import sys

rule host_removal:
    """
    Step 6: Host removal
    """
    input:
        r1 = os.path.join(QC, "step_5", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(QC, "step_5", PATTERN_R2 + ".s5.out.fastq"),
        refpath = os.path.join(HOSTPATH, "ref")
    output:
        unmapped = os.path.join(QC, "step_6", "{sample}.host.unmapped.s6.out.fastq"),
        mapped = os.path.join(QC, "step_6", "{sample}.host.mapped.s6.out.fastq")
    params:
        hostpath = HOSTPATH
    benchmark:
        "benchmarks/remove_vector_contamination_{sample}.txt"
    resources:
        mem_mb=100000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbmap.sh in={input.r1} in2={input.r2} \
            outu={output.unmapped} outm={output.mapped} \
            path={params.hostpath} \
            semiperfectmode=t quickmatch fast ordered=t ow=t \
            -Xmx{resources.mem_mb}m
        """

rule line_sine_removal:
    """
    Step 6a. Remove any LINES and SINES in the sequences.
    """
    input:
        unmapped = os.path.join(QC, "step_6", "{sample}.host.unmapped.s6.out.fastq"),
        linesine = os.path.join(CONPATH, "line_sine.fasta")
    output:
        unmapped = os.path.join(QC, "step_6", "{sample}.linesine.unmapped.s6.out.fastq"),
        mapped   = os.path.join(QC, "step_6", "{sample}.linesine.mapped.s6.out.fastq"),
        stats    = os.path.join(QC, "step_6", "{sample}.linesine.stats")
    benchmark:
        "benchmarks/line_sine_removal_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.unmapped} out={output.unmapped} \
          outm={output.mapped} \
          ref={input.linesine} k=31 hdist=1 stats={output.stats} \
          -Xmx{resources.mem_mb}m
        """

rule repair_paired_ends:
    """
    Step 6b. Repair the paired ends
    """
    input:
        unmapped = os.path.join(QC, "step_6", "{sample}.linesine.unmapped.s6.out.fastq")
    output:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq")
    benchmark:
        "benchmarks/repair_paired_ends_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """   
        repair.sh in={input.unmapped} \
            out={output.r1} out2={output.r2} \
            ow=t \
            -Xmx{resources.mem_mb}m
        """

rule trim_low_quality:
    """
    Step 7: Trim low-quality bases
    """
    input:
        r1 = os.path.join(QC, "step_6", PATTERN_R1 + ".s6.out.fastq"),
        r2 = os.path.join(QC, "step_6", PATTERN_R2 + ".s6.out.fastq")
    output:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + ".s7.out.fastq"),
        r2 = os.path.join(QC, "step_7", PATTERN_R2 + ".s7.out.fastq"),
        singletons = os.path.join(QC, "step_7", "{sample}.singletons.s7.out.fastq"),
        stats = os.path.join(QC, "step_7", "{sample}.s7.stats.txt")
    benchmark:
        "benchmarks/trim_low_quality_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} outs={output.singletons} \
            stats={output.stats} \
            qtrim=r trimq=20 maxns=2 minlength=50 ordered=t ow=t \
            -Xmx{resources.mem_mb}m
        """

"""
In the original contaminant_removal.sh pipeline these were sequential. 
However, snakemake can easily run these in parallel as they
are independent. So I split them into a few separate rules
"""

rule get_r1_singletons:
    """
    Step 8b.i Split R1 singletons
    """
    input:
        singletons = os.path.join(QC, "step_7", "{sample}.singletons.s7.out.fastq")
    output:
        r1singletons = os.path.join(QC, "step_8", "{sample}.singletons.R1.out.fastq"),
    shell:
        """
        egrep --no-group-separator -A 3 '/1|1:N:' {input.singletons} > {output.r1singletons};
        """

rule get_r2_singletons:
    """
    Step 8b.ii Split R2 singletons
    """
    input:
        singletons = os.path.join(QC, "step_7", "{sample}.singletons.s7.out.fastq")
    output:
        r2singletons = os.path.join(QC, "step_8", "{sample}.singletons.R2.out.fastq")
    params:
        odir = os.path.join(QC, "step_8")
    shell:
        """
        egrep --no-group-separator -A 3 '/2|2:N:' {input.singletons} > {output.r2singletons};
        """

rule concat_r1:
    """
    Step 8c.i Concatenate reads
    """
    input:
        r1 = os.path.join(QC, "step_7", PATTERN_R1 + ".s7.out.fastq"),
        r1singletons = os.path.join(QC, "step_8", "{sample}.singletons.R1.out.fastq")
    output:
        r1combo = os.path.join(QC, "step_8", PATTERN_R1 + ".s8.out.fastq"),
    shell:
        """
        cat {input.r1} {input.r1singletons} > {output.r1combo};
        """

rule concat_r2:
    """
    Step 8c.ii Concatenate reads
    """
    input:
        r2 = os.path.join(QC, "step_7", PATTERN_R2 + ".s7.out.fastq"),
        r2singletons = os.path.join(QC, "step_8", "{sample}.singletons.R2.out.fastq")
    output:
        r2combo = os.path.join(QC, "step_8", PATTERN_R2 + ".s8.out.fastq")
    shell:
        """
        cat {input.r2} {input.r2singletons} > {output.r2combo};
        """

rule remove_bacteria:
    """
    Step 9: Remove bacterial contaminants reserving viral and ambiguous sequences
    """
    input:
        r1 = os.path.join(QC, "step_8", PATTERN_R1 + ".s8.out.fastq"),
        bacpath = os.path.join(BACPATH, "ref")
    output:
        mapped = os.path.join(QC, "step_9", PATTERN_R1 + ".bacterial.fastq"),
        unmapped = os.path.join(QC, "step_9", PATTERN_R1 + ".viral_amb.fastq"),
        scafstats = os.path.join(QC, "step_9", PATTERN_R1 + ".scafstats.txt")
    params:
        bacpath = BACPATH
    benchmark:
        "benchmarks/remove_bacteria_{sample}.txt"
    resources:
        mem_mb=100000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbmap.sh in={input.r1} \
            path={params.bacpath} \
            outm={output.mapped} outu={output.unmapped} \
            scafstats={output.scafstats} \
            semiperfectmode=t quickmatch fast ordered=t ow=t \
            -Xmx{resources.mem_mb}m
       """



