###################################################################
#                                                                 #
# Run prinseq if we want to clean the files before starting       #
#                                                                 #
###################################################################

import os
import sys

rule prinseq_low_quality:
    """
    Step 00: Trim low-quality bases.

    Here we use prinseq++ to trim and dereplicate the sequences by quaity
    """


    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + file_extension),
        r2 = os.path.join(READDIR, PATTERN_R2 + file_extension)
    output:
        r1 = os.path.join(QC, "step_0", PATTERN_R1 + ".good_out.s0.fastq"),
        r2 = os.path.join(QC, "step_0", PATTERN_R2 + ".good_out.s0.fastq"),
        s1 = os.path.join(QC, "step_0", PATTERN_R1 + ".single_out.s0.fastq"),
        s2 = os.path.join(QC, "step_0", PATTERN_R2 + ".single_out.s0.fastq"),
        b1 = temporary(os.path.join(QC, "step_0", PATTERN_R1 + ".bad_out_R1.fastq")),
        b2 = temporary(os.path.join(QC, "step_0", PATTERN_R2 + ".bad_out_R2.fastq"))
    benchmark:
        "benchmarks/trim_low_quality_{sample}.txt"
    resources:
        mem_mb=20000,
        cpus=8
    params:
        o = os.path.join(QC, "step_0")
    conda:
        "../envs/prinseq.yaml"
    shell:
        """
            prinseq++ -min_len 60 -min_qual_mean 25 -ns_max_n 1 -derep 1 \
                    -out_format 0 -trim_tail_left 5 -trim_tail_right 5 \
                    -ns_max_n 5  -trim_qual_type min -trim_qual_left 30 \
                    -trim_qual_right 30 -trim_qual_window 10 \
                    -threads {resources.cpus} \
                    -out_good {output.r1} -out_single {output.s1} -out_bad {output.b1} \
                    -out_good2 {output.r2} -out_single2 {output.s2} -out_bad2 {output.b2} \
                    -fastq {input.r1} \
                    -fastq2 {input.r2};
        """
