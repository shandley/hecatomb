"""

Snakefile based on [cluster_count.sh](../base/cluster_count.sh)

Rob Edwards, Jan 2020
"""

import os
import sys

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()

DBDIR = config['Paths']['Databases']


CLUMPED = config['Output']["Clumped"]
QC = config['Output']['QC']

DATADIR = os.path.join(QC, "step_9")
if not os.path.exists(DATADIR):
    sys.stderr.write("ERROR: Please run contaminant_removal.snakefile for the first data steps")
    # we could just add this as a rule ...
    sys.exit(1)


# Step 1: Remove exact duplicates (consider to be PCR artifacts)
# Step 2: Deduplicate
# Step 3: Reformat and prepare for sequence table generation

SAMPLES, = glob_wildcards(os.path.join(DATADIR, '{sample}_viral_amb.fastq'))

rule all:
    input:
        expand(os.path.join(QC, "counts", "{sample}_seqtable.txt"), sample=SAMPLES)


rule remove_exact_dups:
    """
    Step 1: Remove exact duplicates
    """
    input:
        os.path.join(DATADIR, "{sample}_viral_amb.fastq")
    output:
        os.path.join(QC, "step_10", "{sample}_R1.s9.deduped.out.fastq")
    shell:
        """
        dedupe.sh in={input} \
                out={output} \
                ac=f  ow=t -Xmx128g;
        """

rule deduplicate:
    """
    Step 2: Dereplicate
    """
    input:
        os.path.join(QC, "step_10", "{sample}_R1.s9.deduped.out.fastq")
    output:
        fa = os.path.join(QC, "step_11", "{sample}_R1.best.fasta")
        stats = os.path.join(QC, "step_11", "{sample}_stats.txt")
    shell:
        """
        dedupe.sh in={input} \
            csf={output.stats} out={output.fa} \
            ow=t s=4 rnc=t pbr=t -Xmx128g;
        """

rule extract_seq_counts:
    """
    Step 3: Extract sequences and counts for seqtable (count table)
    """
    input:
        os.path.join(QC, "step_11", "{sample}_R1.best.fasta")
    output:
        os.path.join(QC, "step_12", "{sample}_reformated.fasta")
    shell:
        """
        reformat.sh in={input} out={output} \
            deleteinput=t fastawrap=0 \
            ow=t \
            -Xmx128g;
        """

rule extract_counts:
    """
    Parse and combine stats and contig files
    """
    input:
        os.path.join(QC, "step_12", "{sample}_reformated.fasta") 
    output:
        os.path.join(QC, "counts", "{sample}_seqs.txt")
    shell:
        """
        grep -v '>' {input} | sed '1i sequence' > {output}
        """

rule extract_counts_ids:
    """
    Extract sequence IDs
    """
    input:
        os.path.join(QC, "step_12", "{sample}_reformated.fasta")
    output:
        os.path.join(QC, "counts", "{sample}_contig_ids.txt")
    shell:
        """
        grep '>' {input} | sed 's|>Cluster_||' | awk -F "," '{{ print$1 }}' | sort -n | sed '1i contig_ids' > {output}
        """

rule exract_count_stats:
    """
    Extract counts
    """
    input:
        os.path.join(QC, "step_11", "{sample}_stats.txt")
    output:
        os.path.join(QC, "counts", "{sample}_counts.txt")
    shell:
        """
        cut -f 2 {input} | sed "1s/size/$F/" > {output}
        """

rule create_seq_table:
    """
    Create sequence table
    """
    input:
        seq = os.path.join(QC, "counts", "{sample}_seqs.txt")
        cnt = os.path.join(QC, "counts", "{sample}_counts.txt")
    output:
        os.path.join(QC, "counts", "{sample}_seqtable.txt")
    shell:
        """
        paste {input.seq} {input.cnt} > {output}
        """



