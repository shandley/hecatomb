"""

Snakefile based on [contaminant_removal.sh](../base/contaminant_removal.sh)

Rob Edwards, Jan 2020
"""

import os
import sys

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()

DBDIR = config['Paths']['Databases']
if not os.path.exists(DBDIR):
    os.mkdir(DBDIR)

# paths for our databases
BACPATH = os.path.join(DBDIR, "bac_giant_unique_species")
HOSTPATH = os.path.join(DBDIR, "human_masked")
CONPATH = os.path.join(DBDIR, "contaminants")

# paths for our data. This is where we will read and put things
READDIR = config['Paths']['Reads']
CLUMPED = config['Output']["Clumped"]
QC = config['Output']['QC']

SAMPLES, = glob_wildcards(os.path.join(READDIR, '{sample}_R1.fastq.gz'))
PATTERN_R1 = '{sample}_R1'
PATTERN_R2 = '{sample}_R2'

"""
Summary:
    # Step 0: Clumpify reads (https://jgi.doe.gov/data-and-tools/bbtools/bb-tools-user-guide/clumpify-guide/)
    # Step 1: Remove 5' amplification primer
    # Step 2: Remove 3' read through contaminant (Reverse complement of amplification primer + 6 bases of the adapter)
    # Step 3: Remove primer free adapter (both orientations)
    # Step 4: Remove adapter free primer (both orientations)
    # Step 5: PhiX Removal and vector contamination removal
    # Step 6: Host-removal
    # Step 7: Trim low-quality bases
    # Step 8: Remove bacterial contaminants reserving viral and aambiguous sequences
"""

rule all:
    input:
        # the database directories
        os.path.join(BACPATH, "ref"),
        os.path.join(HOSTPATH, "ref"),
        # expand(os.path.join(CLUMPED, "{sample}_R1.clumped.fastq.gz"), sample=SAMPLES)
        # expand(os.path.join(QC, "step_1", "{sample}.s1.stats.txt"), sample=SAMPLES)

        # step 9 output
        expand(os.path.join(QC, "step_9", "{sample}.viral_amb.fastq"), sample=SAMPLES)

# NOTE: bbduk uses "threads=auto" by default that typically uses all threads, so no need to specify. -Xmx is used to
# specify the memory allocation

rule remove_leftmost_primerB:
    """
    Step 1: Remove leftmost primerB. Not the reverse complements
    """
    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + ".fastq.gz"),
        r2 = os.path.join(READDIR, PATTERN_R2 + ".fastq.gz"),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = os.path.join(QC, "step_1", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(QC, "step_1", PATTERN_R2 + ".s1.out.fastq"),
        stats = os.path.join(QC, "step_1", "{sample}.s1.stats.txt")
    benchmark:
        "benchmarks/removeprimerB_{sample}.txt"
    resources:
        time_min = 240,
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
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
        time_min = 240,
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
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
        time_min = 240,
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
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
        time_min = 240,
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
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
        time_min = 240,
        mem_mb=20000,
        cpus=8
    conda:
        "envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t \
            -Xmx{resources.mem_mb}m
        """







