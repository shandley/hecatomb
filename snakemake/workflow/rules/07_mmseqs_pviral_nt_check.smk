"""
Based on mmseqs_pviral_nt_check.sh

Rob Edwards August, 2020

"""


import os
import sys

### These imports are from hecatomb.snakefile

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


RESULTS = config['Output']['Results']
DBDIR = config['Paths']['Databases']
TMPDIR = config['Paths']['Temp']
if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR, exist_ok=True)

AA_OUT  = os.path.join(RESULTS, "mmseqs_aa_out")
if not os.path.exists(AA_OUT):
    os.makedirs(AA_OUT, exist_ok=True)

AA_OUT_CHECKED  = os.path.join(RESULTS, "mmseqs_aa_checked_out")
if not os.path.exists(AA_OUT_CHECKED):
    os.makedirs(AA_OUT_CHECKED, exist_ok=True)


PHAGE_LINEAGES = os.path.join(DBDIR, "phages", "phage_taxonomic_lineages.txt")
if not os.path.exists(PHAGE_LINEAGES):
    sys.stderr.write("FATAL: phages/phage_taxonomic_lineages.txt not ")
    sys.stderr.write("found in the databases directory. Please check ")
    sys.stderr.write("you have the latest version of the databases\n")
    sys.exit()



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




### These are new imports/definitions

BVMDB = os.path.join(NUCLPATH, "bac_virus_masked", "nt.fnaDB")
if not os.path.exists(BVMDB):
    sys.stderr.write(f"FATAL: You appear not to have the nucleotide ")
    sys.stderr.write(f"database {BVMDB} installed.\n")
    sys.stderr.write(f"Please download the databases using the download_databases.snakefile\n")
    sys.exit()

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
    shell:
        """
        mmseqs search {params.st} {NTDB} {params.rdb} $(mktemp -d -p {TMPDIR}) \
        -a -e 0.000001 --search-type 3 --cov-mode 2 -c 0.95
        """

rule filter_nt_db:
    input:
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.dbtype")
    output:
        idx = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.index"),
        dbt = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.dbtype")
    params:
        rdb = os.path.join(NT_CHECKED_OUT, "resultDB")
    shell:
        """
        mmseqs filterdb {params.rdb} {output}  --extract-lines 1
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
    shell:
        """
        mmseqs convertalis {params.sq} {BVMDB} {params.fh} {output}
        """

rule annotate_checked_nt:
    input:
        fhtbl = os.path.join(NT_CHECKED_OUT, "resultDB.firsthit.m8")
    output:
        linout = os.path.join(NT_CHECKED_OUT, "mmseqs_pviral_nt_checked_lineage.tsv")
    params:
        taxtax = TAXTAX
    script:
        "../scripts/mmseqs_pviral_nt_check_annotate.R"


