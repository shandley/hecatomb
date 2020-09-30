"""
Snakefile to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

This is based on [mmseqs_pviral_aa.sh](../base/mmseqs_pviral_aa.sh)


References:

Heavy reliance on:
        # mmseqs2: https://github.com/soedinglab/MMseqs2
        # pullseq: https://github.com/bcthomas/pullseq

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

PHAGE_LINEAGES = os.path.join(DBDIR, "phages", "phage_taxonomic_lineages.txt")
if not os.path.exists(PHAGE_LINEAGES):
    sys.stderr.write("FATAL: phages/phage_taxonomic_lineages.txt not ")
    sys.stderr.write("found in the databases directory. Please check ")
    sys.stderr.write("you have the latest version of the databases\n")
    sys.exit()



rule all:
    input:
        os.path.join(AA_OUT, "phage_tax_table.tsv"),
        os.path.join(AA_OUT, "viruses_tax_table.tsv"),
        os.path.join(AA_OUT, "unclassified_seqs.fasta")

rule convert_seqtable_to_fasta:
    input:
        os.path.join(RESULTS, "seqtable.tab2fx")
    output:
        os.path.join(RESULTS, "seqtable.fasta")
    shell:
        "sed -e 's/^/>/; s/\\t/\\n/' {input} > {output}"

rule create_seqtable_db:
    input:
        os.path.join(RESULTS, "seqtable.fasta")
    output:
        os.path.join(AA_OUT, "seqtable_query.db")
    shell:
        "mmseqs createdb --shuffle 0 --dbtype 0 {input} {output}"

rule seqtable_taxsearch:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
        db = os.path.join(PROTPATH, "uniprot_virus_c99.db")
    output:
        tr = os.path.join(AA_OUT, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT, "taxonomyResult")
    shell:
        """
        mmseqs taxonomy {input.sq} {input.db} {params.tr} $(mktemp -d -p {TMPDIR}) \
        -a --start-sens 1 --sens-steps 3 -s 7 \
        --search-type 2 --tax-output-mode 1
        """

rule seqtable_convert_alignments:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
        db = os.path.join(PROTPATH, "uniprot_virus_c99.db"),
        tr = os.path.join(AA_OUT, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT, "taxonomyResult")
    output:
        os.path.join(AA_OUT, "aln.m8")
    shell:
        """
        mmseqs convertalis {input.sq} {input.db} {params.tr} {output} \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
        """

rule seqtable_lca:
    input:
        db = os.path.join(PROTPATH, "uniprot_virus_c99.db"),
        tr = os.path.join(AA_OUT, "taxonomyResult.dbtype")
    output:
        os.path.join(AA_OUT, "lca.db.dbtype")
    params:
        lc = os.path.join(AA_OUT, "lca.db"),
        tr = os.path.join(AA_OUT, "taxonomyResult")
    shell:
        """
        mmseqs lca {input.db} {params.tr} {params.lc} --tax-lineage true \
        --lca-ranks "superkingdom,phylum,class,order,family,genus,species";
        """

rule seqtable_taxtable_tsv:
    input:
        sq = os.path.join(AA_OUT, "seqtable_query.db"),
        lc = os.path.join(AA_OUT, "lca.db.dbtype")
    params:
        lc = os.path.join(AA_OUT, "lca.db")
    output:
        os.path.join(AA_OUT, "taxonomyResult.tsv")
    shell:
        """
        mmseqs createtsv {input.sq} {params.lc} {output}
        """

rule seqtable_create_kraken:
    input:
        db = os.path.join(PROTPATH, "uniprot_virus_c99.db"),
        lc = os.path.join(AA_OUT, "lca.db")
    output:
        os.path.join(AA_OUT, "taxonomyResult.report")
    shell:
        """
        mmseqs taxonomyreport {input.db} {input.lc} {output}
        """

## Adjust taxonomy table and extract viral lineages
# Extract all (virus + phage) potential viral sequences

rule find_viruses:
    input:
        os.path.join(AA_OUT, "taxonomyResult.tsv")
    output:
        os.path.join(AA_OUT, "all_viruses_table.tsv")
    shell:
        """
        grep 'Viruses;' {input} | cut -f1,5 | sed 's/phi14:2/phi14_2/g' | \
                sed 's/;/\\t/g' | \
                sort -n -k1 > {output}
        """

# Extract phage viral lineages and generate taxonomy table for import into R as PhyloSeq object
rule find_phages:
    input:
        av = os.path.join(AA_OUT, "all_viruses_table.tsv")
    output:
        os.path.join(AA_OUT, "phage_table.tsv")
    shell:
        "grep -f {PHAGE_LINEAGES} {input.av} > {output}"

rule find_phage_seqs:
    input:
        os.path.join(AA_OUT, "phage_table.tsv")
    output:
        os.path.join(AA_OUT, "phage_seqs.list")
    shell:
        "cut -f1 {input} > {output}"

rule pull_phage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT, "phage_seqs.list")
    output:
        os.path.join(AA_OUT, "phage_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule phage_seqs_to_tab:
    input:
        os.path.join(AA_OUT, "phage_seqs.fasta")
    output:
        os.path.join(AA_OUT, "phage_seqs.tab")
    shell:
        """
        perl -pe 'if (s/^>//) {{chomp; s/$/\t/}}' {input} > {output}
        """

rule phage_to_tax_table:
    input:
        tab = os.path.join(AA_OUT, "phage_seqs.tab"),
        tsv = os.path.join(AA_OUT, "phage_table.tsv")
    output:
        os.path.join(AA_OUT, "phage_tax_table.tsv")
    shell:
        """
        join {input.tab} {input.tsv} | \
        cut -d ' ' --output-delimiter=$'\t' -f 2-9 | \
        sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' \
        > {output}
        """

# Extract non-phage viral lineages and generate taxonomy table for import into R as PhyloSeq object
rule find_non_phages:
    input:
        av = os.path.join(AA_OUT, "all_viruses_table.tsv")
    output:
        os.path.join(AA_OUT, "viruses_table.tsv")
    shell:
        "grep -vf {PHAGE_LINEAGES} {input.av} > {output}"

rule find_non_phage_seqs:
    input:
        os.path.join(AA_OUT, "viruses_table.tsv")
    output:
        os.path.join(AA_OUT, "viruses_seqs.list")
    shell:
        "cut -f1 {input} > {output}"

rule pull_non_phage_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT, "viruses_seqs.list")
    output:
        os.path.join(AA_OUT, "viruses_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """

rule non_phage_seqs_to_tab:
    input:
        os.path.join(AA_OUT, "viruses_seqs.fasta")
    output:
        os.path.join(AA_OUT, "viruses_seqs.tab")
    shell:
        """
        perl -pe 'if (s/^>//) {{chomp; s/$/\t/}}' {input} > {output}
        """

rule non_phage_to_tax_table:
    input:
        tab = os.path.join(AA_OUT, "viruses_seqs.tab"),
        tsv = os.path.join(AA_OUT, "viruses_table.tsv")
    output:
        os.path.join(AA_OUT, "viruses_tax_table.tsv")
    shell:
        """
        join {input.tab} {input.tsv} | \
        cut -d ' ' --output-delimiter=$'\t' -f 2-9 | \
        sed '1isequence\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies' \
        > {output}
        """


# Extract unclassified lineages
rule unclassified_lineages:
    input:
        os.path.join(AA_OUT, "taxonomyResult.tsv")
    output:
        os.path.join(AA_OUT, "pviral_unclassified_seqs.list")
    shell:
        """
        grep -v 'Viruses;' {input} | cut -f1 | \
                sort -n -k1 > {output}
        """

rule pull_unclassified_seqs:
    input:
        fa = os.path.join(RESULTS, "seqtable.fasta"),
        ls = os.path.join(AA_OUT, "pviral_unclassified_seqs.list")
    output:
        os.path.join(AA_OUT, "unclassified_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """


