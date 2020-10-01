"""
Based on mmseqs_pviral_nt.sh

Rob Edwards April, 2020

"""


import os
import sys

### These imports are from hecatomb.snakefile

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


DBDIR = config['Paths']['Databases']
TMPDIR = config['Paths']['Temp']
RESULTS = config['Output']['Results']
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR, exist_ok=True)

AA_OUT  = os.path.join(RESULTS, "mmseqs_aa_out")
if not os.path.exists(AA_OUT):
    os.makedirs(AA_OUT, exist_ok=True)

AA_OUT_CHECKED  = os.path.join(RESULTS, "mmseqs_aa_checked_out")
if not os.path.exists(AA_OUT_CHECKED):
    os.makedirs(AA_OUT_CHECKED, exist_ok=True)



### These are new imports/definitions

NUCLPATH = os.path.join(DBDIR, "nucleotides")
NTDB = os.path.join(NUCLPATH, "refseq_virus_nt_UniVec_masked", "nt.fnaDB")
if not os.path.exists(NTDB):
    sys.stderr.write(f"FATAL: You appear not to have the nucleotide ")
    sys.stderr.write(f"database {NTDB} installed.\n")
    sys.stderr.write(f"Please download the databases using the download_databases.snakefile\n")
    sys.exit()


NT_OUT = os.path.join(RESULTS, "mmseqs_nt_out")
if not os.path.exists(NT_OUT):
    os.makedirs(NT_OUT)

NT_CHECKED_OUT = os.path.join(RESULTS, "mmseqs_nt_checked_out")
if not os.path.exists(NT_CHECKED_OUT):
    os.makedirs(NT_CHECKED_OUT)

# taxonomizr taxa
TAXPATH  = os.path.join(DBDIR, "taxonomy")
TAXTAX = os.path.join(TAXPATH, "taxonomizr_accessionTaxa.sql")
if not os.path.exists(TAXTAX):
    sys.stderr.write(f"FATAL: You appear not to have the taxonomizr ")
    sys.stderr.write(f"database {TAXTAX} installed.\n")
    sys.stderr.write(f"Please download the databases using the download_databases.snakefile\n")
    sys.exit()




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
    shell:
        """
        mmseqs search {params.st} {NTDB} {params.rdb} $(mktemp -d -p {TMPDIR}) \
        -a -e 0.000001 --search-type 3 --cov-mode 2 -c 0.95
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
    shell:
        """
        mmseqs filterdb {params.rdb} {params.rdbfh} --extract-lines 1
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
    shell:
        """
        mmseqs convertalis {params.st} {NTDB} {params.rdbfh} {output}
        """

rule nt_annotate:
    input:
        fhtbl = os.path.join(NT_OUT, "resultDB.firsthit.m8")
    output:
        linout = os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_lineage.tsv")
    params:
        taxtax = TAXTAX
    script:
        "scripts/mmseqs_pviral_nt_annotate.R"






