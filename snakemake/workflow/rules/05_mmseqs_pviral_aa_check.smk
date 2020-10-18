"""
Snakefile to check the amino acid searches

Based upon ../base/mmseqs_pviral_aa_check.sh

Rob Edwards, April 2020

"""

import os
import sys




if not os.path.exists(AA_OUT_CHECKED):
    os.makedirs(AA_OUT_CHECKED, exist_ok=True)


rule mmseqs_pviral_aa_check_first:
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
    benchmark: "benchmarks/create_viral_seqs_db.txt"
    resources:
        mem_mb=20000
    conda:
        "../envs/mmseqs2.yaml"
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
    benchmark: "benchmarks/viral_seqs_tax_search.txt"
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult")
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs taxonomy {input.vqdb} {URVDB} {params.tr} \
            $(mktemp -d -p {TMPDIR}) --threads {resources.cpus} \
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
    benchmark: "benchmarks/viral_seqs_convertalis.txt"
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis {input.vqdb} {URVDB} {params.tr} {output} --threads {resources.cpus} \
        --format-output "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits,qaln,taln"
        """

rule viral_seqs_lca:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult"),
        ot = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "lca.db.dbtype"),
        os.path.join(AA_OUT_CHECKED, "lca.db.index")
    benchmark: "benchmarks/viral_seqs_lca.txt"
    resources:
        mem_mb=100000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs lca {URVDB} {params.tr} {params.ot} \
        --tax-lineage 1 --threads {resources.cpus} \
        --lca-ranks "superkingdom,phylum,class,order,family,genus,species"
        """

rule extract_top_hit:
    input:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.dbtype")
    params:
        tr = os.path.join(AA_OUT_CHECKED, "taxonomyResult"),
        fh = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.dbtype"),
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.index")
    benchmark: "benchmarks/extract_top_hit.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs filterdb {params.tr} {params.fh} --extract-lines 1  --threads {resources.cpus}
        """

rule convertalis_vsqd:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        trfhdb = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.dbtype")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit.m8")
    params:
        trfh = os.path.join(AA_OUT_CHECKED, "taxonomyResult.firsthit")
    benchmark: "benchmarks/convertalis_vsqd.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs convertalis {input.vqdb} {URVDB} {params.trfh} {output} --threads {resources.cpus}
        """

rule create_taxtable_vsqd:
    input:
        vqdb = os.path.join(AA_OUT_CHECKED, "viral_seqs_queryDB"),
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db.dbtype")
    params:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.tsv")
    benchmark: "benchmarks/create_taxtable_vsqd.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs createtsv {input.vqdb} {params.lcadb} {output} --threads {resources.cpus}
        """

rule create_kraken_vsqd:
    input:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db.dbtype")
    params:
        lcadb = os.path.join(AA_OUT_CHECKED, "lca.db")
    output:
        os.path.join(AA_OUT_CHECKED, "taxonomyResult.report")
    benchmark: "benchmarks/create_kraken_vsqd.txt"
    resources:
        mem_mb=20000,
        cpus=16
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        mmseqs taxonomyreport {URVDB} {params.lcadb} {output} --threads {resources.cpus}
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
            sed 's/;/\t/g' | \
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

