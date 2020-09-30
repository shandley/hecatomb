"""
Snakefile to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

This is based on [mmseqs_pviral_aa.sh](../base/mmseqs_pviral_aa.sh)


References:

Heavy reliance on:
        # mmseqs2: https://github.com/soedinglab/MMseqs2
        # pullseq: https://github.com/bcthomas/pullseq
        # SeqKit: https://bioinf.shenwei.me/seqkit/

REQUIRES that targetDB has already been indexed
If it has not been index then run the following script in the directory of your choice: uniprot_viral_DB_build.sh (found in /accessory)                                                                  
Note: mmseqs2 taxonomy is currently most useful if you have UniProt formatted fasta databases
more details about database building can be found at: https://github.com/soedinglab/mmseqs2/wiki#taxonomy-assignment-using-mmseqs-taxonomy                                                               


Rob Edwards, March 2020

"""

"""
Before we begin, we need to make sure the databases have been downloaded
and installed.

If not, we throw an error here rather than proceeding to download them.
That way, the downloading is asynchronous.
"""

import os
import sys

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


DBDIR = config['Paths']['Databases']

# paths for our databases
PROTPATH = os.path.join(DBDIR, "proteins")

if not os.path.exists(PROTPATH):
    sys.stderr.write("FATAL: You appear not to have the protein databases. Please download the databases using the download_databases.snakefile\n")
    sys.exit()

# The results directory we create with hecatomb
RESULTS = config['Output']['Results']
AA_OUT  = os.path.join(RESULTS, "mmseqs_aa_out")
if not os.path.exists(AA_OUT):
    os.mkdir(AA_OUT)

# how much memory we have
XMX = config['System']['Memory']


# The virus database, clustered at 99% with cd-hit and then compiled
# with mmseqs

VIRDB = os.path.join(PROTPATH, "uniprot_virus_c99.db")
if not os.path.exists(VIRDB):
    sys.stderr.write(f"FATAL: {VIRDB} does not exist. Please ensure you")
    sys.stderr.write(" have installed the databases\n")
    sys.exit()

# Phage lineages. We normally read this from a text file
#TODO do we want it as a file or a set?
def phage_lineages:
    pl = set()
    with open("../base/phage_taxonomic_lineages.txt", 'r') as f:
        for l in f:
            pl.add(l.strip())
    return pl

















rule convert_seqtable_to_fasta:
    input:
        os.path.join(RESULTS, "seqtable.tab2fx")
    output:
        os.path.join(RESULTS, "seqtable.fasta")
    shell:
        "seqkit tab2fx {input} -w 5000 -o {output}"

rule create_seqtable_db:
    input:
        os.path.join(RESULTS, "seqtable.fasta")
    output:
        os.path.join(AA_OUT, "seqtable_query.db")
    shell:
        "mmseqs createdb {input} {output} --dont-shuffle 0 --dbtype 0"

rule seqtable_taxsearch:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
        db = os.path.join(PROTPATH, "uniprot_virus_c99.db")
    output:
        tr = os.path.join(AA_OUT, "taxonomyResult"),
        aa = os.path.join(AA_OUT, "tmp_aa")
    shell:
        """
        mmseqs taxonomy {input.sq} {input.db} {output.tr} {output.aa} \
        -a --start-sens 1 --sens-steps 3 -s 7 \
        --search-type 2 --tax-output-mode 1
        """

rule seqtable_convert_alignments:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
        db = os.path.join(PROTPATH, "uniprot_virus_c99.db"),
        tr = os.path.join(AA_OUT, "taxonomyResult")
    output:
        os.path.join(AA_OUT, "aln.m8")
    shell:
        """
        mmseqs convertalis {input.sq} {input.db} {input.tr} {output} \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
        """

rule seqtable_lca:
    input:
        db = os.path.join(PROTPATH, "uniprot_virus_c99.db"),
        tr = os.path.join(AA_OUT, "taxonomyResult")
    output:
        os.path.join(AA_OUT, "lca.db")
    shell:
        """
        mmseqs lca {input.db} {input.tr} {output} --tax-lineage true \
        --lca-ranks "superkingdom:phylum:class:order:family:genus:species"
        """
    








