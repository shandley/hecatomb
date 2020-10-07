"""
Based on mmseqs_pviral_nt_check.sh

Rob Edwards August, 2020

"""


import os
import sys


rule mmseqs_pviral_nt_check_first:
    input:
        os.path.join(NT_CHECKED_OUT, "phage_nt_seqs.fasta"),
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_seqs.fasta"),
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage.tsv")



rule find_nt_phages:
    input:
        # this input comes from mmseqs_pviral_nt.snakefile
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_lineage.tsv")
    output:
        os.path.join(NT_CHECKED_OUT, "phage_nt_table.tsv")
    shell:
        """
        tail -n+2 {input} | grep -f {PHAGE_LINEAGES} |  sort -n -k1 \
                > {output}
        """

rule list_nt_phages:
    input:
        os.path.join(NT_CHECKED_OUT, "phage_nt_table.tsv")
    output:
         os.path.join(NT_CHECKED_OUT, "phage_nt_table.list")
    shell:
        """
        cut -f1 {input} > {output}
        """

rule pull_nt_phage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(NT_CHECKED_OUT, "phage_nt_table.list")
    output:
        os.path.join(NT_CHECKED_OUT, "phage_nt_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule non_nt_phage_viruses:
    input:
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_lineage.tsv")
    output:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_table.tsv")
    shell:
        """
        tail -n+2 {input} | grep -vf {PHAGE_LINEAGES} |  sort -n -k1 \
                > {output}
        """

rule list_non_nt_viruses:
    input:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_table.tsv")
    output:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_table.list")
    shell:
        """
        cut -f1 {input} > {output}
        """

rule pull_nt_non_phage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_table.list")
    output:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule create_nt_db:
    input:
        os.path.join(NT_CHECKED_OUT, "pviral_virus_nt_seqs.fasta")
    output:
         idx = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.index"),
         dbt = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.dbtype")
    params:
         st = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB")
    benchmark: "benchmarks/create_nt_db.txt"
    resources:
        mem_mb=100000,
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createdb {input} {params.st} --dbtype 2
        """

rule nt_search_checked:
    input:
        idx = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.dbtype")
    output:
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.dbtype")
    params:
        st = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB"),
        rdb = os.path.join(NT_CHECKED_OUT, "resultDB")
    benchmark: "benchmarks/nt_search_checked.txt"
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs search {params.st} {NTDB} {params.rdb} $(mktemp -d -p {TMPDIR}) \
        -a -e 0.000001 --search-type 3 --cov-mode 2 -c 0.95 --threads {resources.cpus}
        """

rule filter_nt_db:
    input:
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.dbtype")
    output:
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.dbtype")
    params:
        rdb = os.path.join(NT_CHECKED_OUT, "resultDB"),
        rfh = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit")
    benchmark: "benchmarks/filter_nt_db.txt"
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs filterdb {params.rdb} {params.rfh}  --extract-lines 1 --threads {resources.cpus}
        """

rule convert_nt_alias:
    input:
        sqi = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.index"),
        sqd = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB.dbtype"),
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.dbtype")
    output:
        os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.m8")
    params:
        fh = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit"),
        sq = os.path.join(NT_CHECKED_OUT, "seqtable_queryDB"),
    benchmark: "benchmarks/convert_nt_alias.txt"
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis {params.sq} {BVMDB} {params.fh} {output} --threads {resources.cpus}
        """

rule annotate_checked_nt:
    input:
        fhtbl = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.m8")
    output:
        linout = os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage.tsv")
    params:
        taxtax = TAXTAX
    benchmark: "benchmarks/annotate_checked_nt.txt"
    resources:
        mem_mb=100000,
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/mmseqs_pviral_nt_check_annotate.R"


