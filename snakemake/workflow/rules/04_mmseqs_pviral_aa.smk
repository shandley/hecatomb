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

import os
import sys

if not os.path.exists(TMPDIR):
    os.makedirs(TMPDIR, exist_ok=True)

if not os.path.exists(AA_OUT):
    os.makedirs(AA_OUT, exist_ok=True)

rule mmseqs_pviral_aa_first:
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
    benchmark:
        "benchmarks/create_seqtable_db.txt"
    resources:
        mem_mb=20000,
    conda:
        "../envs/mmseqs2.yaml"
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
    benchmark:
        "benchmarks/seqtable_taxsearch.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs taxonomy {input.sq} {input.db} {params.tr} $(mktemp -d -p {TMPDIR}) \
        -a --start-sens 1 --sens-steps 3 -s 7 --threads {resources.cpus} \
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
    benchmark:
        "benchmarks/seqtable_convert_alignments.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis {input.sq} {input.db} {params.tr} {output} --threads {resources.cpus} \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
        """

rule seqtable_lca:
    input:
        tr = os.path.join(AA_OUT, "taxonomyResult.dbtype")
    output:
        os.path.join(AA_OUT, "lca.db.dbtype")
    params:
        lc = os.path.join(AA_OUT, "lca.db"),
        tr = os.path.join(AA_OUT, "taxonomyResult")
    benchmark:
        "benchmarks/seqtable_lca.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs lca  {VIRDB} {params.tr} {params.lc}  --tax-lineage 1 --threads {resources.cpus} \
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
    benchmark:
        "benchmarks/seqtable_taxtable_tsv.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createtsv --threads {resources.cpus} {input.sq} {params.lc} {output}
        """

rule seqtable_create_kraken:
    input:
        db = os.path.join(PROTPATH, "uniprot_virus_c99.db"),
        lc = os.path.join(AA_OUT, "lca.db")
    output:
        os.path.join(AA_OUT, "taxonomyResult.report")
    benchmark:
        "benchmarks/seqtable_create_kraken.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs taxonomyreport --threads {resources.cpus} {VIRDB} {input.lc} {output}
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
        os.path.join(AA_OUT, "pviral_aa_unclassified_seqs.fasta")
        #os.path.join(AA_OUT, "unclassified_seqs.fasta")
    shell:
        """
        grep --no-group-separator -A 1 -Fwf {input.ls} {input.fa} > {output}
        """


