"""
Hecatomb.smk to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

History: This is based on [mmseqs_pviral_aa.sh](../base/mmseqs_pviral_aa.sh)

Rob Edwards, March 2020
Overhauled - Michael Roach, Q2/Q3 2021
"""


rule primary_aa_search:
    """Taxon step 01: Assign taxonomy to dir["out"]["results"]/seqtable.fasta sequences using mmseqs2 with primary AA database"""
    input:
        seqs = os.path.join(dir["out"]["results"], "seqtable.fasta"),
        db = dir["dbs"]["primaryAA"]
    output:
        os.path.join(dir["out"]["primaryAA"], "mmseqs.primary.aa.alignments.tsv")
    params:
        alnRes = os.path.join(dir["out"]["primaryAA"], "mmseqs_aa_primary"),
        filtaa = config["mmseqs"]["filtAA"],
        sensaa = config["mmseqs"]["sens"],
        memsplit = str(int(0.75 * int(resources["big"]["mem"]))) + "M",
    benchmark:
        os.path.join(dir["out"]["bench"], "primary_aa_taxonomy_assignment.txt")
    log:
        os.path.join(dir["out"]["stderr"], "primary_aa_taxonomy_assignment.log")
    resources:
        mem_mb = resources["big"]["mem"],
        mem = str(resources["big"]["mem"]) + "MB",
        time = resources["big"]["time"]
    threads:
        resources["big"]["cpu"]
    conda:
        os.path.join(dir["env"], "mmseqs2.yaml")
    shell:
        "mmseqs easy-search {input.seqs} {input.db} {output} {params.alnRes} "
            " {params.filtaa} {params.sensaa} "
            " --threads {threads} --split-memory-limit {params.memsplit} &> {log} "



rule primary_aa_parsing:
    """Taxon step 02: Parse primary AA search results for classified (potentially viral) and unclassified sequences"""
    input:
        alnsort = os.path.join(dir["out"]["primaryAA"], "mmseqs.primary.aa.alignments.tsv"),
        seqs = os.path.join(dir["out"]["results"], "seqtable.fasta"),
    output:
        class_seqs = os.path.join(dir["out"]["primaryAA"], "primary.aa.classified.fasta"),
    resources:
        time = resources["sml"]["time"]
    benchmark:
        os.path.join(dir["out"]["bench"], "primary_aa_parsing.txt")
    log:
        os.path.join(dir["out"]["stderr"], "primary_aaA_parsing.log")
    script:
        os.path.join(dir["scripts"],  "aaPrimaryParse.py")


rule secondary_aa_taxonomy_assignment:
    """Taxon step 03: Check taxonomic assignments in MMSEQS_AA_PRIMARY_classified.fasta using mmseqs2"""
    input:
        seqs = os.path.join(dir["out"]["primaryAA"], "primary.aa.classified.fasta"),
        db = os.path.join(dir["dbs"]["secondaryAA"], "sequenceDB")
    output:
        lca = os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary_lca.tsv"),
        report = os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary_report"),
        tophit_report = os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary_tophit_report"),
        aln = os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary_tophit_aln")
    params:
        alnRes=os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary"),
        tmppath=os.path.join(dir["out"]["secondaryAA"], "tmp"),
        filtaa = config["mmseqs"]["filtAA"],
        formataa = config["immutable"]["reqAA"],
        sensaa = config["mmseqs"]["sens"],
        memsplit = str(int(0.75 * int(resources["big"]["mem"]))) + "M",
        aaHeader = config["immutable"]["mmseqsHeaderAA"]
    benchmark:
        os.path.join(dir["out"]["bench"], "secondary_aa_taxonomy_assignment.txt")
    log:
        os.path.join(dir["out"]["stderr"], "secondary_aa_taxonomy_assignment.log")
    resources:
        mem_mb = resources["big"]["mem"],
        mem = str(resources["big"]["mem"]) + "MB",
        time = resources["big"]["time"]
    threads:
        resources["big"]["cpu"]
    conda:
        os.path.join(dir["env"], "mmseqs2.yaml")
    shell:
        "{{ "
        "mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} "
            "{params.filtaa} {params.sensaa} {params.formataa} "
            "--lca-mode 2 --threads {threads} --split-memory-limit {params.memsplit}; "
        "}} &> {log} "


rule secondary_aa_tophit_lineage:
    """Taxon step 04: Add/reformat tophit viral lineages with up-to-date* NCBI taxonomy"""
    input:
        db = dir["dbs"]["taxonomy"],
        tophit = os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary_tophit_aln")
    output:
        tophit_lineage_refomated = os.path.join(dir["out"]["secondaryAA"], "tophit.lineage.reformated"),
    conda:
        os.path.join(dir["env"], "seqkit.yaml")
    resources:
        time = resources["sml"]["time"],
        mem = resources["ram"]["mem"]
    params:
        taxonFormat = lambda wildcards: config["immutable"]["taxonkitReformat"]
    benchmark:
        os.path.join(dir["out"]["bench"], "secondary_aa_tophit_lineage.txt")
    log:
        os.path.join(dir["out"]["stderr"], "secondary_aa_tophit_lineage.log")
    group:
        "secondary_aa_parsing"
    shell:
        "{{ "
        "cut -f1,20 {input.tophit} "
            "| taxonkit lineage --data-dir {input.db} -i 2 "
            "| taxonkit reformat --data-dir {input.db} -i 3 {params.taxonFormat} "
            "| cut --complement -f3 "
            "> {output.tophit_lineage_refomated}; }} &> {log} "


rule secondary_aa_refactor_finalize:
    """Taxon step 05: Remove sequences to be refactored from LCA table and recombine with updated taxonomies."""
    input:
        db = dir["dbs"]["taxonomy"],
        lca = os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary_lca.tsv"),
    output:
        lca_reformated = os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary_lca.reformatted"),
    conda:
        os.path.join(dir["env"], "seqkit.yaml")
    resources:
        time=resources["sml"]["time"],
        mem=resources["ram"]["mem"]
    params:
        taxonFormat = lambda wildcards: config["immutable"]["taxonkitReformat"]
    benchmark:
        os.path.join(dir["out"]["bench"], "secondary_aa_refactor_finalize.txt")
    log:
        os.path.join(dir["out"]["stderr"], "secondary_aa_refactor_finalize.log")
    group:
        "secondary_aa_parsing"
    shell:
        "{{ cut -f1,2 {input.lca} "
            "| taxonkit lineage --data-dir {input.db} -i 2 "
            "| taxonkit reformat --data-dir {input.db} -i 3 {params.taxonFormat} "
            "| cut --complement -f3 "
            "> {output.lca_reformated}; }} &> {log} "


rule secondary_aa_output_table:
    """Taxon step 06: Join sequence info, tophit align info, and LCA or tophit lineage info into the output format table"""
    input:
        aln = os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary_tophit_aln"),
        lca = os.path.join(dir["out"]["secondaryAA"], "mmseqs.aa.secondary_lca.reformatted"),
        top = os.path.join(dir["out"]["secondaryAA"], "tophit.lineage.reformated"),
        counts = os.path.join(dir["out"]["results"], "sampleSeqCounts.tsv"),
        balt = os.path.join(dir["dbs"]["tables"], "2020_07_27_Viral_classification_table_ICTV2019.txt")
    output:
        vir = os.path.join(dir["out"]["secondaryAA"], "AA_bigtable.tsv"),
        nonvir = os.path.join(dir["out"]["secondaryAA"], "AA_bigtable.nonviral.tsv")
    resources:
        time=resources["sml"]["time"],
        mem=resources["ram"]["mem"]
    benchmark:
        os.path.join(dir["out"]["bench"], "secondary_aa_generate_output_table.txt")
    log:
        os.path.join(dir["out"]["stderr"], "secondary_aa_generate_output_table.log")
    group:
        "secondary_aa_parsing"
    params:
        taxIdIgnore = config["mmseqs"]["taxIdIgnore"].split(),
        bigtableHeader = config["immutable"]["bigtableHeader"]
    script:
        os.path.join(dir["scripts"],  "aaBigtable.py")


rule secondary_aa_parsing:
    """Taxon step 07: Parse out all sequences that remain unclassified following the Secondary AA search and refactoring."""
    input:
        bigtable = os.path.join(dir["out"]["secondaryAA"],"AA_bigtable.tsv"),
        seqs = os.path.join(dir["out"]["results"], "seqtable.fasta")
    output:
        unclass_seqs = os.path.join(dir["out"]["primaryAA"], "primary.aa.unclassified.fasta")
    resources:
        time=resources["sml"]["time"],
        mem=resources["ram"]["mem"]
    benchmark:
        os.path.join(dir["out"]["bench"], "secondary_aa_parsing.txt")
    log:
        os.path.join(dir["out"]["stderr"], "secondary_aa_parsing.log")
    group:
        "secondary_aa_parsing"
    script:
        os.path.join(dir["scripts"],  "aaSecondaryParse.py")

        
rule primary_nt_search:
    """Taxon step 08: Primary nucleotide search of unclassified viral-like sequences from aa search"""
    input:
        seqs = os.path.join(dir["out"]["primaryAA"], "primary.aa.unclassified.fasta"),
        db = dir["dbs"]["primaryNT"]
    output:
        os.path.join(dir["out"]["primaryNT"], "mmseqs.primary.nt.alignments.tsv"),
    params:
        tmppath = os.path.join(dir["out"]["primaryNT"], "mmseqs_aa_tmp"),
        filtnt = config["mmseqs"]["filtNT"],
        ntsens = config["mmseqs"]["sens"],
        memsplit = str(int(0.75 * int(resources["big"]["mem"]))) + "M"
    benchmark:
        os.path.join(dir["out"]["bench"], "primary_nt_taxonomic_assignment.txt")
    log:
        os.path.join(dir["out"]["stderr"], "primary_nt_taxonomic_assignment.log")
    resources:
        mem_mb = resources["big"]["mem"],
        mem = str(resources["big"]["mem"]) + "MB",
        time = resources["big"]["time"]
    threads:
        resources["big"]["cpu"]
    conda:
        os.path.join(dir["env"], "mmseqs2.yaml")
    shell:
        "mmseqs easy-search {input.seqs} {input.db} {output} {params.tmppath} "
            "{params.ntsens} {params.filtnt} "
            "--search-type 3 --threads {threads} --split-memory-limit {params.memsplit} &> {log} "


rule primary_nt_parsing:
    """Taxon step 10: Extract unclassified sequences from secondary translated search
    
    xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
    """
    input:
        seqs = os.path.join(dir["out"]["primaryAA"], "primary.aa.unclassified.fasta"),
        align = os.path.join(dir["out"]["primaryNT"], "mmseqs.primary.nt.alignments.tsv")
    output:
        class_seqs = os.path.join(dir["out"]["primaryNT"], "primary.nt.classified.fasta"),
        unclass_seqs = os.path.join(dir["out"]["primaryNT"], "primary.nt.unclassified.fasta")
    resources:
        time = resources["sml"]["time"]
    benchmark:
        os.path.join(dir["out"]["bench"], "primary_nt_parsing.txt")
    log:
        os.path.join(dir["out"]["stderr"], "primary_nt_parsing.log")
    script:
        os.path.join(dir["scripts"],  "ntPrimaryParse.py")


rule secondary_nt_search:
    """Taxon step 11: Secondary nucleotide search of viral hits from primary nucleotide search"""
    input:
        seqs = os.path.join(dir["out"]["primaryNT"], "primary.nt.classified.fasta"),
        db = os.path.join(dir["dbs"]["secondaryNT"], "sequenceDB")
    output:
        aln = os.path.join(dir["out"]["secondaryNT"], "mmseqs.secondary.nt.alignments.tsv"),
        tax = temp(os.path.join(dir["out"]["secondaryNT"], "all.taxid"))
    params:
        tmp = os.path.join(dir["out"]["secondaryNT"], "tmp"),
        ntfilt = config["mmseqs"]["filtNT"],
        sensnt = config["mmseqs"]["sens"],
        format = config["immutable"]["secondaryNtFormat"],
        memsplit = str(int(0.75 * int(resources["big"]["mem"]))) + "M"
    benchmark:
        os.path.join(dir["out"]["bench"], "secondary_nt_taxonomic_assignment.txt")
    log:
        os.path.join(dir["out"]["stderr"], "secondary_nt_taxonomic_assignment.log")
    resources:
        mem_mb=resources["big"]["mem"],
        mem=str(resources["big"]["mem"]) + "MB",
        time = resources["big"]["time"]
    threads:
        resources["big"]["cpu"]
    conda:
        os.path.join(dir["env"], "mmseqs2.yaml")
    shell:
        "{{ "
        "if [[ -d {params.tmp} ]]; then rm -r {params.tmp}; fi; "
        "mmseqs easy-search {input.seqs} {input.db} {output.aln} {params.tmp} "
            "{params.sensnt} {params.ntfilt} {params.format} "
            "--search-type 3 --threads {threads} --split-memory-limit {params.memsplit}; "
        "cut -f1,2 {output.aln} "
            "| sed 's/tid|//' "
            "| sed 's/|.*//' "
            "> {output.tax}; "
        "}} &> {log} "


rule secondary_nt_lca_table:
    """Taxon step 14: Create table for taxonkit lineage for secondary NT search"""
    input:
        align = os.path.join(dir["out"]["secondaryNT"], "all.taxid")
    output:
        lin = temp(os.path.join(dir["out"]["secondaryNT"], "all.lin")),
        top = temp(os.path.join(dir["out"]["secondaryNT"], "top.lin"))
    benchmark:
        os.path.join(dir["out"]["bench"], "secondary_nt_lca_table.txt")
    log:
        os.path.join(dir["out"]["stderr"], "secondary_nt_lca_table.log")
    resources:
        time = resources["sml"]["time"],
        mem = resources["ram"]["mem"]
    group:
        "secondary_nt_parsing"
    script:
        os.path.join(dir["scripts"],  "ntSecondaryLca.py")


rule secondary_nt_calc_lca:
    """Taxon step 15: Calculate the lca for the secondary NT search"""
    input:
        lin = os.path.join(dir["out"]["secondaryNT"], "all.lin"),
        top = os.path.join(dir["out"]["secondaryNT"], "top.lin"),
        taxdb = dir["dbs"]["taxonomy"]
    output:
        lca_lineage = os.path.join(dir["out"]["secondaryNT"], "lca_lineage.tsv"),
        top_lineage = os.path.join(dir["out"]["secondaryNT"], "top_lineage.tsv"),
    resources:
        time=resources["sml"]["time"],
        mem=resources["ram"]["mem"]
    params:
        taxonFormat = lambda wildcards: config["immutable"]["taxonkitReformat"],
    conda:
        os.path.join(dir["env"], "mmseqs2.yaml")
    benchmark:
        os.path.join(dir["out"]["bench"], "secondary_nt_calc_lca.txt")
    log:
        os.path.join(dir["out"]["stderr"], "secondary_nt_calc_lca.log")
    group:
        "secondary_nt_parsing"
    shell:
        "{{ "
        "taxonkit lca -i 2 -s ';' --data-dir {input.taxdb} {input.lin} | "
            "awk -F '\t' '$3 != 0' | "
            "taxonkit lineage -i 3 --data-dir {input.taxdb} | "
            "taxonkit reformat --data-dir {input.taxdb} -i 4 {params.taxonFormat} | "
            "cut --complement -f 3,4 "
            "> {output.lca_lineage}; "
        "taxonkit lineage -i 2 --data-dir {input.taxdb} {input.top} | "
            "taxonkit reformat --data-dir {input.taxdb} -i 3 {params.taxonFormat} | "
            "cut --complement -f 3 > {output.top_lineage}; "
        "}} &> {log}; "


rule secondary_nt_generate_output_table:
    """Taxon step 16: Create the NT portion of the bigtable"""
    input:
        aln = os.path.join(dir["out"]["secondaryNT"], "mmseqs.secondary.nt.alignments.tsv"),
        lca = os.path.join(dir["out"]["secondaryNT"], "lca_lineage.tsv"),
        top = os.path.join(dir["out"]["secondaryNT"], "top_lineage.tsv"),
        counts = os.path.join(dir["out"]["results"], "sampleSeqCounts.tsv"),
        balt = os.path.join(dir["dbs"]["tables"], "2020_07_27_Viral_classification_table_ICTV2019.txt")
    output:
        vir = os.path.join(dir["out"]["secondaryNT"], "NT_bigtable.tsv"),
        nonvir = os.path.join(dir["out"]["secondaryNT"], "NT_bigtable.nonviral.tsv")
    resources:
        time=resources["sml"]["time"],
        mem=resources["ram"]["mem"]
    params:
        taxIdIgnore = config["mmseqs"]["taxIdIgnore"].split(),
        bigtableHeader = config["immutable"]["bigtableHeader"]
    benchmark:
        os.path.join(dir["out"]["bench"], "secondary_nt_generate_output_table.txt")
    log:
        os.path.join(dir["out"]["stderr"], "secondary_nt_generate_output_table.log")
    group:
        "secondary_nt_parsing"
    script:
        os.path.join(dir["scripts"],  "ntBigtable.py")


rule combine_aa_nt:
    """Taxon step 17: Combine the AA and NT bigtables"""
    input:
        aa = os.path.join(dir["out"]["secondaryAA"],"AA_bigtable.tsv"),
        nt = os.path.join(dir["out"]["secondaryNT"], "NT_bigtable.tsv")
    output:
        os.path.join(dir["out"]["results"], "bigtable.tsv")
    benchmark:
        os.path.join(dir["out"]["bench"], "combine_AA_NT.txt")
    log:
        os.path.join(dir["out"]["stderr"], "combine_AA_NT.log")
    resources:
        time=resources["sml"]["time"],
        mem=resources["ram"]["mem"]
    group:
        "secondary_nt_parsing"
    shell:
        "{{ cat {input.aa} > {output}; "
        "tail -n+2 {input.nt} >> {output}; }} &> {log}; "
