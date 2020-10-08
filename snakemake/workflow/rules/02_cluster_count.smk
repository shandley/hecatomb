"""

Snakefile based on [cluster_count.sh](../base/cluster_count.sh)

Rob Edwards, Jan 2020
"""

import os
import sys

rule remove_exact_dups:
    """
    Step 1: Remove exact duplicates
    """
    input:
        os.path.join(QC, "step_9", PATTERN_R1 + ".viral_amb.fastq")
    output:
        os.path.join(QC, "step_10", PATTERN_R1 + ".s9.deduped.out.fastq")
    benchmark:
        "benchmarks/remove_exact_dups_{sample}.txt"
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
                -Xmx{resources.mem_mb}m
        """

rule deduplicate:
    """
    Step 2: Dereplicate
    """
    input:
        os.path.join(QC, "step_10", PATTERN_R1 + ".s9.deduped.out.fastq")
    output:
        fa = os.path.join(QC, "step_11", PATTERN_R1 + ".best.fasta"),
        stats = os.path.join(QC, "step_11", "{sample}_stats.txt")
    benchmark:
        "benchmarks/deduplicate_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        dedupe.sh in={input} \
            csf={output.stats} out={output.fa} \
            ow=t s=4 rnc=t pbr=t \
            -Xmx{resources.mem_mb}m
        """

rule extract_seq_counts:
    """
    Step 3: Extract sequences and counts for seqtable (count table)
    """
    input:
        os.path.join(QC, "step_11", PATTERN_R1 + ".best.fasta")
    output:
        os.path.join(QC, "step_12", PATTERN_R1 + ".reformated.fasta")
    benchmark:
        "benchmarks/extract_seq_counts_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        reformat.sh in={input} out={output} \
            deleteinput=t fastawrap=0 \
            ow=t \
            -Xmx{resources.mem_mb}m
        """

rule extract_counts:
    """
    Parse and combine stats and contig files
    """
    input:
        os.path.join(QC, "step_12", PATTERN_R1 + ".reformated.fasta")
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
        os.path.join(QC, "step_12", PATTERN_R1 + ".reformated.fasta")
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
    params:
        s = "{sample}"
    shell:
        """
        cut -f 2 {input} | sed "1s/size/{params.s}/" > {output}
        """

rule create_seq_table:
    """
    Create sequence table
    """
    input:
        seq = os.path.join(QC, "counts", "{sample}_seqs.txt"),
        cnt = os.path.join(QC, "counts", "{sample}_counts.txt")
    output:
        os.path.join(QC, "counts", "{sample}_seqtable.txt")
    shell:
        """
        paste {input.seq} {input.cnt} > {output}
        """



