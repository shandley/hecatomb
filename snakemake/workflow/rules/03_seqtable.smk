"""

Snakefile specifically to run the R code.

Rob Edwards, Feb 2020
"""

import os
import sys


rule seqtable_first:
    input:
        os.path.join(RESULTS, "seqtable_all.tsv"),
        os.path.join(RESULTS, "seqtable.tab2fx")

rule merge_seq_table:
    """
    Merge seq counts
    """
    input:
        files = expand(os.path.join(QC, "counts", "{sample}_seqtable.txt"), sample=SAMPLES)
    output:
        seqtable = os.path.join(RESULTS, "seqtable_all.tsv"),
        tab2fx = os.path.join(RESULTS, "seqtable.tab2fx")
    params:
        resultsdir = directory(RESULTS),
    benchmark:
        "benchmarks/merge_seq_table.txt"
    resources:
        mem_mb=20000,
        cpus=8
    params:
        resultsdir = directory(RESULTS),
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/seqtable_merge.R"



