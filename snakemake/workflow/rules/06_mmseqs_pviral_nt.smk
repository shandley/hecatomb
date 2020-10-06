"""
Based on mmseqs_pviral_nt.sh

Rob Edwards April, 2020

"""

import os
import sys

rule mmseqs_pviral_nt_first:
    input:
        os.path.join(NT_OUT, "resultDB.firsthit.m8"),
        os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_lineage.tsv")

rule create_nt_querydb:
    input:
        os.path.join(AA_OUT, "pviral_aa_unclassified_seqs.fasta")
    output:
        idx = os.path.join(NT_OUT, "seqtable_queryDB.index"),
        dbt = os.path.join(NT_OUT, "seqtable_queryDB.dbtype")
    params:
        st = os.path.join(NT_OUT, "seqtable_queryDB")
    benchmark: "benchmarks/create_nt_querydb.txt"
    resources:
        mem_mb=100000,
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createdb {input} {params.st}  --dbtype 2
        """

rule nt_search:
    input:
        idx = os.path.join(NT_OUT, "seqtable_queryDB.index"),
        dbt = os.path.join(NT_OUT, "seqtable_queryDB.dbtype")
    output:
        idx = os.path.join(NT_OUT, "resultDB.index"),
        dbt = os.path.join(NT_OUT, "resultDB.dbtype")
    params:
        st = os.path.join(NT_OUT, "seqtable_queryDB"),
        rdb = os.path.join(NT_OUT, "resultDB")
    benchmark: "benchmarks/nt_search.txt"
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs search {params.st} {NTDB} {params.rdb} $(mktemp -d -p {TMPDIR}) \
        -a -e 0.000001 --search-type 3 --cov-mode 2 -c 0.95  --threads {resources.cpus}
        """

rule nt_top_hit:
    input:
        idx = os.path.join(NT_OUT, "resultDB.index"),
        dbt = os.path.join(NT_OUT, "resultDB.dbtype")
    output:
        os.path.join(NT_OUT, "resultDB.firsthit.dbtype"),
        os.path.join(NT_OUT, "resultDB.firsthit.index")
    params:
        rdb = os.path.join(NT_OUT, "resultDB"),
        rdbfh = os.path.join(NT_OUT, "resultDB.firsthit")
    benchmark: "benchmarks/nt_top_hit.txt"
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs filterdb {params.rdb} {params.rdbfh} --extract-lines 1  --threads {resources.cpus}
        """

rule nt_to_m8:
    input:
        sidx = os.path.join(NT_OUT, "seqtable_queryDB.index"),
        sdbt = os.path.join(NT_OUT, "seqtable_queryDB.dbtype"),
        ridx = os.path.join(NT_OUT, "resultDB.firsthit.dbtype"),
        rdbt = os.path.join(NT_OUT, "resultDB.firsthit.index")
    output:
        os.path.join(NT_OUT, "resultDB.firsthit.m8")
    params:
        st = os.path.join(NT_OUT, "seqtable_queryDB"),
        rdbfh = os.path.join(NT_OUT, "resultDB.firsthit")
    benchmark: "benchmarks/nt_to_m8.txt"
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis {params.st} {NTDB} {params.rdbfh} {output}  --threads {resources.cpus}
        """

rule nt_annotate:
    input:
        fhtbl = os.path.join(NT_OUT, "resultDB.firsthit.m8")
    output:
        linout = os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_lineage.tsv")
    params:
        taxtax = TAXTAX
    benchmark: "benchmarks/nt_annotate.txt"
    resources:
        mem_mb=100000,
    conda:
        "../envs/R.yaml"
    script:
        "../scripts/mmseqs_pviral_nt_annotate.R"






