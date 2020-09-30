"""
Snakefile to check the amino acid searches

Based upon ../base/mmseqs_pviral_aa_check.sh

Rob Edwards, April 2020

"""

import os
import sys

if not config:
    sys.stderr.write("FATAL: Please define a config file using the --configfile command line option.\n")
    sys.stderr.write("examples are provided in the Git repo\n")
    sys.exit()


DBDIR = config['Paths']['Databases']
TMPDIR = config['Paths']['Temp']
if not os.path.exists(TMPDIR):
    os.mkdir(TMPDIR)

# paths for our databases
PROTPATH = os.path.join(DBDIR, "proteins")

if not os.path.exists(PROTPATH):
    sys.stderr.write("FATAL: You appear not to have the protein databases. Please download the databases using the download_databases.snakefile\n")
    sys.exit()

# The results directory we create with hecatomb
RESULTS = config['Output']['Results']
AA_OUT  = os.path.join(RESULTS, "mmseqs_aa_out")
if not os.path.exists(AA_OUT):
    sys.stderr.write(f"FATAL: You do not appear to have directory {AA_OUT}\n")
    sys.stderr.write("Did you run mmseqs_pviral.snakefile?")
    sys.exit()

AA_OUT_CHECKED  = os.path.join(RESULTS, "mmseqs_aa_checked_out")
if not os.path.exists(AA_OUT_CHECKED):
    os.mkdir(AA_OUT_CHECKED)

# how much memory we have
XMX = config['System']['Memory']


# The virus database, clustered at 99% with cd-hit and then compiled
# with mmseqs

VIRDB = os.path.join(PROTPATH, "uniprot_virus_c99.db")
if not os.path.exists(VIRDB):
    sys.stderr.write(f"FATAL: {VIRDB} does not exist. Please ensure you")
    sys.stderr.write(" have installed the databases\n")
    sys.exit()

PHAGE_LINEAGES = os.path.join(DBDIR, "phages", "phage_taxonomic_lineages.txt")
if not os.path.exists(PHAGE_LINEAGES):
    sys.stderr.write("FATAL: phages/phage_taxonomic_lineages.txt not ")
    sys.stderr.write("found in the databases directory. Please check ")
    sys.stderr.write("you have the latest version of the databases\n")
    sys.exit()

# uniref50 + viruses
URVPATH = os.path.join(PROTPATH, "uniref_plus_virus")
URVDB = os.path.join(URVPATH, "uniref50_virus.db") # uniref50 + viruses database
if not os.path.exists(URVDB):
    sys.stderr.write("FATAL: {URVDB} not found.\n")
    sys.stderr.write("Please make sure that you have run ")
    sys.stderr.write("download_databases.snakefile before commencing\n")
    sys.exit()

rule all:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.report"),
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv"),
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_tax_table.tsv"),
        os.path.join(AA_OUT_CHECKED, "unclassified_checked_aa_seqs.fasta")


rule create_viral_seqs_db:
    input:
        os.path.join(AA_OUT, "viruses_seqs.fasta")
    output:
        os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB")
    shell:
        """
        mmseqs createdb {input} {output} --shuffle 0 --dbtype 0
        """

rule viral_seqs_tax_search:
    """
    The taxnomy result maybe several disjoint files
    """
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.index")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult")
    shell:
        """
        mmseqs taxonomy {input.vqdb} {URVDB} {params.tr} \
            $(mktemp -d -p {TMPDIR}) \
            -a -s 7 --search-type 2 --tax-output-mode 1
        """

rule viral_seqs_convertalis:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult")
    output:
        os.path.join(AA_OUT_CHECKED, "aln.m8")
    shell:
        """
        mmseqs convertalis {input.vqdb} {URVDB} {params.tr} {output} \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
        """

rule viral_seqs_lca:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult")
    output:
        os.path.join(AA_OUT_CHECKED, "lca.db.dbtype"),
        os.path.join(AA_OUT_CHECKED, "lca.db.index")
    shell:
        """
        mmseqs lca {URVDB} {params.tr} {output} \
        --tax-lineage true \
        --lca-ranks "superkingdom,phylum,class,order,family,genus,species"
        """

rule extract_top_hit:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.dbtype"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.index")
    shell:
        """
        mmseqs filterdb {params.tr} {output} --extract-lines 1
        """

rule convertalis_vsqd:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        trfhdb = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.dbtype")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8")
    params:
        trfh = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit")
    shell:
        """
        mmseqs convertalis {input.vqdb} {URVDB} {params.trfh} {output} 
        """

rule create_taxtable_vsqd:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db.dbtype")
    params:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.tsv")
    shell:
        """
        mmseqs createtsv {input.vqdb} {params.lcadb} {output}
        """

rule create_kraken_vsqd:
    input:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db.dbtype")
    params:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.report")
    shell:
        """
        mmseqs taxonomyreport {URVDB} {params.lcadb} {output}
        """

rule nonphages_to_pyloseq_table:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv")
    shell:
        """
        grep -v 'Bacteria;' {input} | \
            grep 'Viruses;' | \
            grep -v -f  {PHAGE_LINEAGES} | cut -f1,5 | \
            sed 's/:/\t/g' | \
                sort -n -k1 > {output}
        """

rule nonphages_to_pyloseq_list:
    input:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.list")
    shell:
        """
        cut -f1 {input} > {output}
        """

rule pull_nonphage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.list")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule nonphage_seqs_to_tab:
    input:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.fasta")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.tab")
    shell:
        """
        perl -pe 'if (s/^>//) {{chomp; s/$/\t/}}' {input} > {output}
        """

rule nonphage_to_tax_table:
    input:
        tab = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_seqs.tab"),
        tsv = os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_table.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "viruses_checked_aa_tax_table.tsv")
    shell:
        """
        join {input.tab} {input.tsv} | \
        cut -d ' ' --output-delimiter=$'\t' -f 2-9 | \
        sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' \
        > {output}
        """

rule non_viral_lineages:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.tsv")
    output:
        os.path.join(AA_OUT_CHECKED, "unclassified_checked_aa_seqs.list")
    shell:
        """
        grep -v 'Viruses;' {input} | cut -f1 | \
           sort -n -k1 > {output}
        """
 
rule pull_non_viral_lineages:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT_CHECKED, "unclassified_checked_aa_seqs.list")
    output:
        os.path.join(AA_OUT_CHECKED, "unclassified_checked_aa_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

