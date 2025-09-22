"""
Hecatomb.smk to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

History: This is based on [mmseqs_pviral_aa.sh](../base/mmseqs_pviral_aa.sh)

Rob Edwards, March 2020
Overhauled - Michael Roach, Q2/Q3 2021
"""


rule read_annotation_primary_aa_search:
    """Taxon step 01: Assign taxonomy to config["hecatomb"]["args"]["output_paths"]["results"]/seqtable.fasta sequences using mmseqs2 with primary AA database"""
    input:
        seqs = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "seqtable.fasta"),
        db = os.path.join(config["hecatomb"]["args"]["database_paths"]["primaryAA"], "sequenceDB")
    output:
        os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryAA"], "mmseqs.primary.aa.alignments.tsv")
    params:
        alnRes = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryAA"], "mmseqs_aa_primary"),
        filtaa = config["hecatomb"]["mmseqs"]["filtAA"],
        sensaa = config["hecatomb"]["mmseqs"][config["hecatomb"]["args"]["search"]],
        memsplit = str(int(0.75 * int(config["resources"]["big"]["mem"]))) + "M",
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "primary_aa_taxonomy_assignment.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "primary_aa_taxonomy_assignment.log")
    resources:
        **config["resources"]["big"]
    threads:
        config["resources"]["big"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    shell:
        "mmseqs easy-search {input.seqs} {input.db} {output} {params.alnRes} "
            "{params.filtaa} {params.sensaa} "
            "--threads {threads} --split-memory-limit {params.memsplit} &> {log}; "


rule read_annotation_primary_aa_parsing:
    """Taxon step 02: Parse primary AA search results for classified (potentially viral) and unclassified sequences"""
    input:
        alnsort = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryAA"], "mmseqs.primary.aa.alignments.tsv"),
        seqs = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "seqtable.fasta"),
    output:
        class_seqs = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryAA"], "primary.aa.classified.fasta"),
    resources:
        **config["resources"]["sml"]
    threads:
        config["resources"]["sml"]["cpu"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "primary_aa_parsing.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "primary_aaA_parsing.log")
    script:
        os.path.join("..", "..", "scripts",  "aaPrimaryParse.py")


rule read_annotation_secondary_aa_taxonomy_assignment:
    """Taxon step 03: Check taxonomic assignments in MMSEQS_AA_PRIMARY_classified.fasta using mmseqs2"""
    input:
        seqs = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryAA"], "primary.aa.classified.fasta"),
        db = os.path.join(config["hecatomb"]["args"]["database_paths"]["secondaryAA"], "sequenceDB")
    output:
        lca = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary_lca.tsv"),
        report = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary_report"),
        tophit_report = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary_tophit_report"),
        aln = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary_tophit_aln")
    params:
        alnRes=os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary"),
        tmppath=os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "tmp"),
        filtaa = config["hecatomb"]["mmseqs"]["filtAA"],
        formataa = config["hecatomb"]["immutable"]["reqAA"],
        sensaa = config["hecatomb"]["mmseqs"][config["hecatomb"]["args"]["search"]],
        memsplit = str(int(0.75 * int(config["resources"]["big"]["mem"]))) + "M",
        aaHeader = config["hecatomb"]["immutable"]["mmseqsHeaderAA"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "secondary_aa_taxonomy_assignment.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "secondary_aa_taxonomy_assignment.log")
    resources:
        **config["resources"]["big"]
    threads:
        config["resources"]["big"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    shell:
        "{{ mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} "
            "{params.filtaa} {params.sensaa} {params.formataa} "
            "--lca-mode 2 --threads {threads} --split-memory-limit {params.memsplit}; "
        "}} &> {log} "


rule read_annotation_secondary_aa_tophit_lineage:
    """Taxon step 04: Add/reformat tophit viral lineages with up-to-date* NCBI taxonomy"""
    input:
        db = config["hecatomb"]["args"]["database_paths"]["taxonomy"],
        tophit = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary_tophit_aln")
    output:
        tophit_lineage_refomated = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "tophit.lineage.reformated"),
    conda:
        os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "secondary_aa_tophit_lineage.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "secondary_aa_tophit_lineage.log")
    group:
        "secondary_aa_parsing"
    shell:
        "{{ cut -f1,20 {input.tophit} "
            "| taxonkit reformat2 --data-dir {input.db} -I 2 --miss-rank-repl NA "
                r"-f '{{domain|acellular root|superkingdom}}\t{{phylum}}\t{{class}}\t{{order}}\t{{family}}\t{{genus}}\t{{species}}' "
            "> {output.tophit_lineage_refomated}; }} 2> {log} "


rule read_annotation_secondary_aa_refactor_finalize:
    """Taxon step 05: Remove sequences to be refactored from LCA table and recombine with updated taxonomies."""
    input:
        db = config["hecatomb"]["args"]["database_paths"]["taxonomy"],
        lca = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary_lca.tsv"),
    output:
        lca_reformated = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary_lca.reformatted"),
    conda:
        os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "secondary_aa_refactor_finalize.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "secondary_aa_refactor_finalize.log")
    group:
        "secondary_aa_parsing"
    shell:
        "{{ cut -f1,2 {input.lca} "
            "| taxonkit reformat2 --data-dir {input.db} -I 2 --miss-rank-repl NA "
                r"-f '{{domain|acellular root|superkingdom}}\t{{phylum}}\t{{class}}\t{{order}}\t{{family}}\t{{genus}}\t{{species}}' "
            "> {output.lca_reformated}; }} &> {log} "


rule read_annotation_secondary_aa_output_table:
    """Taxon step 06: Join sequence info, tophit align info, and LCA or tophit lineage info into the output format table"""
    input:
        aln = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary_tophit_aln"),
        lca = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "mmseqs.aa.secondary_lca.reformatted"),
        top = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "tophit.lineage.reformated"),
        counts = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "sampleSeqCounts.tsv"),
        balt = os.path.join(config["hecatomb"]["args"]["database_paths"]["tables"], "2020_07_27_Viral_classification_table_ICTV2019.txt")
    output:
        vir = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "AA_bigtable.tsv"),
        nonvir = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"], "AA_bigtable.nonviral.tsv")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "secondary_aa_generate_output_table.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "secondary_aa_generate_output_table.log")
    group:
        "secondary_aa_parsing"
    params:
        taxIdIgnore = config["hecatomb"]["mmseqs"]["taxIdIgnore"].split(),
        bigtableHeader = config["hecatomb"]["immutable"]["bigtableHeader"]
    script:
        os.path.join("..", "..", "scripts",  "aaBigtable.py")


rule read_annotation_secondary_aa_parsing:
    """Taxon step 07: Parse out all sequences that remain unclassified following the Secondary AA search and refactoring."""
    input:
        bigtable = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"],"AA_bigtable.tsv"),
        seqs = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "seqtable.fasta")
    output:
        unclass_seqs = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryAA"], "primary.aa.unclassified.fasta")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "secondary_aa_parsing.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "secondary_aa_parsing.log")
    group:
        "secondary_aa_parsing"
    script:
        os.path.join("..", "..", "scripts",  "aaSecondaryParse.py")

        
rule read_annotation_primary_nt_search:
    """Taxon step 08: Primary nucleotide search of unclassified viral-like sequences from aa search"""
    input:
        seqs = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryAA"], "primary.aa.unclassified.fasta"),
        db = os.path.join(config["hecatomb"]["args"]["database_paths"]["primaryNT"], "sequenceDB")
    output:
        os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryNT"], "mmseqs.primary.nt.alignments.tsv"),
    params:
        tmppath = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryNT"], "mmseqs_aa_tmp"),
        filtnt = config["hecatomb"]["mmseqs"]["filtNT"],
        ntsens = config["hecatomb"]["mmseqs"][config["hecatomb"]["args"]["search"]],
        memsplit = str(int(0.75 * int(config["resources"]["big"]["mem"]))) + "M"
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "primary_nt_taxonomic_assignment.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "primary_nt_taxonomic_assignment.log")
    resources:
        **config["resources"]["big"]
    threads:
        config["resources"]["big"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    shell:
        "mmseqs easy-search {input.seqs} {input.db} {output} {params.tmppath} "
            "{params.ntsens} {params.filtnt} "
            "--search-type 3 --threads {threads} --split-memory-limit {params.memsplit} &> {log} "


rule read_annotation_primary_nt_parsing:
    """Taxon step 10: Extract unclassified sequences from secondary translated search
    
    xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
    """
    input:
        seqs = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryAA"], "primary.aa.unclassified.fasta"),
        align = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryNT"], "mmseqs.primary.nt.alignments.tsv")
    output:
        class_seqs = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryNT"], "primary.nt.classified.fasta"),
        unclass_seqs = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryNT"], "primary.nt.unclassified.fasta")
    resources:
        **config["resources"]["sml"]
    threads:
        config["resources"]["sml"]["cpu"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "primary_nt_parsing.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "primary_nt_parsing.log")
    script:
        os.path.join("..", "..", "scripts",  "ntPrimaryParse.py")


rule read_annotation_secondary_nt_search:
    """Taxon step 11: Secondary nucleotide search of viral hits from primary nucleotide search"""
    input:
        seqs = os.path.join(config["hecatomb"]["args"]["temp_paths"]["primaryNT"], "primary.nt.classified.fasta"),
        db = os.path.join(config["hecatomb"]["args"]["database_paths"]["secondaryNT"], "sequenceDB")
    output:
        aln = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "mmseqs.secondary.nt.alignments.tsv"),
        tax = temp(os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "all.taxid"))
    params:
        tmp = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "tmp"),
        ntfilt = config["hecatomb"]["mmseqs"]["filtNT"],
        sensnt = config["hecatomb"]["mmseqs"][config["hecatomb"]["args"]["search"]],
        format = config["hecatomb"]["immutable"]["secondaryNtFormat"],
        memsplit = str(int(0.75 * int(config["resources"]["big"]["mem"]))) + "M"
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "secondary_nt_taxonomic_assignment.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "secondary_nt_taxonomic_assignment.log")
    resources:
        **config["resources"]["big"]
    threads:
        config["resources"]["big"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    shell:
        "{{ if [[ -d {params.tmp} ]]; then rm -r {params.tmp}; fi; "
        "mmseqs easy-search {input.seqs} {input.db} {output.aln} {params.tmp} "
            "{params.sensnt} {params.ntfilt} {params.format} "
            "--search-type 3 --threads {threads} --split-memory-limit {params.memsplit}; "
        "cut -f1,2 {output.aln} "
            "| sed 's/tid|//' "
            "| sed 's/|.*//' "
            "> {output.tax}; "
        "}} &> {log} "


rule read_annotation_secondary_nt_lca_table:
    """Taxon step 14: Create table for taxonkit lineage for secondary NT search"""
    input:
        align = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "all.taxid")
    output:
        lin = temp(os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "all.lin")),
        top = temp(os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "top.lin"))
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "secondary_nt_lca_table.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "secondary_nt_lca_table.log")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    group:
        "secondary_nt_parsing"
    script:
        os.path.join("..", "..", "scripts",  "ntSecondaryLca.py")


rule read_annotation_secondary_nt_calc_lca:
    """Taxon step 15: Calculate the lca for the secondary NT search"""
    input:
        lin = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "all.lin"),
        top = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "top.lin"),
        taxdb = config["hecatomb"]["args"]["database_paths"]["taxonomy"]
    output:
        lca_lineage = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "lca_lineage.tsv"),
        top_lineage = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "top_lineage.tsv"),
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "secondary_nt_calc_lca.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "secondary_nt_calc_lca.log")
    group:
        "secondary_nt_parsing"
    shell:
        "{{ taxonkit lca -i 2 -s ';' --data-dir {input.taxdb} {input.lin} "
            r"| awk -F '\t' '$3 != 0' "
            "| taxonkit reformat2 --data-dir {input.taxdb} -I 3 --miss-rank-repl NA "
                r"-f '{{domain|acellular root|superkingdom}}\t{{phylum}}\t{{class}}\t{{order}}\t{{family}}\t{{genus}}\t{{species}}' "
            "| cut --complement -f2 "
            "> {output.lca_lineage}; "
        "taxonkit reformat2 -I 2 --data-dir {input.taxdb} {input.top} --miss-rank-repl NA "
            r"-f '{{domain|acellular root|superkingdom}}\t{{phylum}}\t{{class}}\t{{order}}\t{{family}}\t{{genus}}\t{{species}}' "
            "> {output.top_lineage}; "
        "}} &> {log}; "


rule read_annotation_secondary_nt_generate_output_table:
    """Taxon step 16: Create the NT portion of the bigtable"""
    input:
        aln = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "mmseqs.secondary.nt.alignments.tsv"),
        lca = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "lca_lineage.tsv"),
        top = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "top_lineage.tsv"),
        counts = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "sampleSeqCounts.tsv"),
        balt = os.path.join(config["hecatomb"]["args"]["database_paths"]["tables"], "2020_07_27_Viral_classification_table_ICTV2019.txt")
    output:
        vir = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "NT_bigtable.tsv"),
        nonvir = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "NT_bigtable.nonviral.tsv")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    params:
        taxIdIgnore = config["hecatomb"]["mmseqs"]["taxIdIgnore"].split(),
        bigtableHeader = config["hecatomb"]["immutable"]["bigtableHeader"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "secondary_nt_generate_output_table.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "secondary_nt_generate_output_table.log")
    group:
        "secondary_nt_parsing"
    script:
        os.path.join("..", "..", "scripts",  "ntBigtable.py")


rule read_annotation_combine_aa_nt:
    """Taxon step 17: Combine the AA and NT bigtables"""
    input:
        aa = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryAA"],"AA_bigtable.tsv"),
        nt = os.path.join(config["hecatomb"]["args"]["temp_paths"]["secondaryNT"], "NT_bigtable.tsv")
    output:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "bigtable.tsv")
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "combine_AA_NT.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "combine_AA_NT.log")
    resources:
        **config["resources"]["ram"]
    threads:
        config["resources"]["ram"]["cpu"]
    group:
        "secondary_nt_parsing"
    shell:
        "{{ cat {input.aa} > {output}; "
        "tail -n+2 {input.nt} >> {output}; }} 2> {log}; "
