"""
Hecatomb.smk to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

History: This is based on [mmseqs_pviral_aa.sh](../base/mmseqs_pviral_aa.sh)

Rob Edwards, March 2020
Overhauled - Michael Roach, Q2/Q3 2021
"""


rule PRIMARY_AA_taxonomy_assignment:
    """Taxon step 01: Assign taxonomy to dir.out.results/seqtable.fasta sequences using mmseqs2 with primary AA database"""
    input:
        seqs = os.path.join(dir.out.results, "seqtable.fasta"),
        db = os.path.join(dir.dbs.primaryAA, "sequenceDB")
    output:
        lca = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_lca.tsv"),
        report = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_report"),
        tophit_report = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_tophit_report"),
        aln = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_tophit_aln"),
        alnsort = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_tophit_aln_sorted")
    params:
        alnRes = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY"),
        tmppath = os.path.join(dir.out.primaryAA, "mmseqs_aa_tmp"),
        filtaa = config.mmseqs.filtAAprimary,
        formataa = config.immutable.reqAA,
        sensaa = config.mmseqs.sensAA,
        memsplit = str(int(0.75 * int(config.resources.big.mem))) + 'M',
        aaHeader = config.immutable.mmseqsHeaderAA
    benchmark:
        os.path.join(dir.out.bench, "PRIMARY_AA_taxonomy_assignment.txt")
    log:
        os.path.join(dir.out.stderr, "PRIMARY_AA_taxonomy_assignment.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    group:
        "primaryaa"
    shell:
        """
        {{ # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
            {params.filtaa} {params.sensaa} {params.formataa} \
            --threads {threads} --split-memory-limit {params.memsplit};
            
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i {params.aaHeader}' > {output.alnsort};
        }} &> {log}
        rm {log}
        """


rule PRIMARY_AA_parsing:
    """Taxon step 02: Parse primary AA search results for classified (potentially viral) and unclassified sequences"""
    input:
        alnsort = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_tophit_aln_sorted"),
        seqs = os.path.join(dir.out.results, "seqtable.fasta"),
    output:
        class_seqs = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_classified.fasta"),
    resources:
        time = config.resources.sml.time
    benchmark:
        os.path.join(dir.out.bench, 'PRIMARY_AA_parsing.txt')
    log:
        os.path.join(dir.out.stderr, 'PRIMARY_AA_parsing.log')
    group:
        "primaryaa"
    script:
        os.path.join(dir.scripts,  'aaPrimaryParse.py')


rule SECONDARY_AA_taxonomy_assignment:
    """Taxon step 03: Check taxonomic assignments in MMSEQS_AA_PRIMARY_classified.fasta using mmseqs2"""
    input:
        seqs = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_classified.fasta"),
        db = os.path.join(dir.dbs.secondaryAA, "sequenceDB")
    output:
        lca = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_lca.tsv"),
        report = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_report"),
        tophit_report = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_tophit_report"),
        aln = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_tophit_aln"),
        alnsort = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_tophit_aln_sorted")
    params:
        alnRes=os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY"),
        tmppath=os.path.join(dir.out.secondaryAA, "mmseqs_aa_tmp"),
        filtaa = config.mmseqs.filtAAsecondary,
        formataa = config.immutable.reqAA,
        sensaa = config.mmseqs.sensAA,
        memsplit = str(int(0.75 * int(config.resources.big.mem))) + 'M',
        aaHeader = config.immutable.mmseqsHeaderAA
    benchmark:
        os.path.join(dir.out.bench, "SECONDARY_AA_taxonomy_assignment.txt")
    log:
        os.path.join(dir.out.stderr, "SECONDARY_AA_taxonomy_assignment.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    group:
        "secondaryaa"
    shell:
        """
        {{ # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
            {params.filtaa} {params.sensaa} {params.formataa} \
            --threads {threads} --split-memory-limit {params.memsplit};
            
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i {params.aaHeader}' > {output.alnsort};
        }} &> {log}
        rm {log}
        """


rule SECONDARY_AA_tophit_lineage:
    """Taxon step 04: Add/reformat tophit viral lineages with up-to-date* NCBI taxonomy"""
    input:
        db = dir.dbs.taxonomy,
        tophit = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_tophit_aln")
    output:
        tophit_lineage_refomated = os.path.join(dir.out.secondaryAA, "tophit.lineage.reformated"),
    conda:
        os.path.join(dir.env, "seqkit.yaml")
    resources:
        time = config.resources.sml.time
    params:
        taxonFormat = lambda wildcards: config.immutable.taxonkitReformat
    benchmark:
        os.path.join(dir.out.bench, "SECONDARY_AA_tophit_lineage.txt")
    log:
        os.path.join(dir.out.stderr, "SECONDARY_AA_tophit_lineage.log")
    group:
        "secondaryaa"
    shell:
        """
        {{ # Make a table: SeqID <tab> taxID
        cut -f1,20 {input.tophit} \
            | taxonkit lineage --data-dir {input.db} -i 2 \
            | taxonkit reformat --data-dir {input.db} -i 3 {params.taxonFormat} \
            | cut --complement -f3 \
            > {output.tophit_lineage_refomated}; }} &> {log}
        rm {log}
        """


rule SECONDARY_AA_refactor_finalize:
    """Taxon step 05: Remove sequences to be refactored from LCA table and recombine with updated taxonomies."""
    input:
        db = dir.dbs.taxonomy,
        lca = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_lca.tsv"),
    output:
        lca_reformated = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_lca.reformated"),
    conda:
        os.path.join(dir.env, "seqkit.yaml")
    resources:
        time = config.resources.sml.time
    params:
        taxonFormat = lambda wildcards: config.immutable.taxonkitReformat
    benchmark:
        os.path.join(dir.out.bench, "SECONDARY_AA_refactor_finalize.txt")
    log:
        os.path.join(dir.out.stderr, "SECONDARY_AA_refactor_finalize.log")
    group:
        "secondaryaa"
    shell:
        """
        {{ cut -f1,2 {input.lca} \
            | taxonkit lineage --data-dir {input.db} -i 2 \
            | taxonkit reformat --data-dir {input.db} -i 3 {params.taxonFormat} \
            | cut --complement -f3 \
            > {output.lca_reformated}; }} &> {log}
        rm {log}
        """


rule SECONDARY_AA_generate_output_table:
    """Taxon step 06: Join sequence info, tophit align info, and LCA or tophit lineage info into the output format table"""
    input:
        aln = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_tophit_aln"),
        lca = os.path.join(dir.out.secondaryAA, "MMSEQS_AA_SECONDARY_lca.reformated"),
        top = os.path.join(dir.out.secondaryAA, "tophit.lineage.reformated"),
        counts = os.path.join(dir.out.results, "sampleSeqCounts.tsv"),
        balt = os.path.join(dir.dbs.tables, '2020_07_27_Viral_classification_table_ICTV2019.txt')
    output:
        vir = os.path.join(dir.out.secondaryAA, "AA_bigtable.tsv"),
        nonvir = os.path.join(dir.out.secondaryAA, "AA_bigtable.nonviral.tsv")
    resources:
        time = config.resources.sml.time
    benchmark:
        os.path.join(dir.out.bench, "SECONDARY_AA_generate_output_table.txt")
    log:
        os.path.join(dir.out.stderr, "SECONDARY_AA_generate_output_table.log")
    group:
        "secondaryaa"
    params:
        taxIdIgnore = config.mmseqs.taxIdIgnore.split(),
        bigtableHeader = config.immutable.bigtableHeader
    script:
        os.path.join(dir.scripts,  'aaBigtable.py')


rule SECONDARY_AA_parsing:
    """Taxon step 07: Parse out all sequences that remain unclassified following the Secondary AA search and refactoring."""
    input:
        bigtable = os.path.join(dir.out.secondaryAA,"AA_bigtable.tsv"),
        seqs = os.path.join(dir.out.results, "seqtable.fasta")
    output:
        unclass_seqs = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_unclassified.fasta")
    resources:
        time = config.resources.sml.time
    benchmark:
        os.path.join(dir.out.bench, "SECONDARY_AA_parsing.txt")
    log:
        os.path.join(dir.out.stderr, 'SECONDARY_AA_parsing.log')
    group:
        "secondaryaa"
    script:
        os.path.join(dir.scripts,  'aaSecondaryParse.py')

        
rule PRIMARY_NT_taxonomic_assignment:
    """Taxon step 08: Primary nucleotide search of unclassified viral-like sequences from aa search"""
    input:
        seqs = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_unclassified.fasta"),
        db = os.path.join(dir.dbs.primaryNT, "sequenceDB")
    output:
        queryDB = os.path.join(dir.out.primaryNT, "queryDB"),
        result = os.path.join(dir.out.primaryNT, "results", "result.index"),
        type = os.path.join(dir.out.primaryNT, "results", "result.dbtype")
    params:
        respath = os.path.join(dir.out.primaryNT, "results", "result"),
        tmppath = os.path.join(dir.out.primaryNT, "mmseqs_aa_tmp"),
        filtnt = config.mmseqs.filtNTprimary,
        ntsens = config.mmseqs.sensNT,
        memsplit = str(int(0.75 * int(config.resources.big.mem))) + 'M'
    benchmark:
        os.path.join(dir.out.bench, "PRIMARY_NT_taxonomic_assignment.txt")
    log:
        os.path.join(dir.out.stderr, "PRIMARY_NT_taxonomic_assignment.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    group:
        "primarynt"
    shell:
        """
        {{ # Create query database
        rm -f {output.type}
        mmseqs createdb {input.seqs} {output.queryDB} --dbtype 2;
        # mmseqs search
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
            {params.ntsens} {params.filtnt} \
            --search-type 3 --threads {threads} --split-memory-limit {params.memsplit};
        }} &> {log}
        rm {log}
        """


rule PRIMARY_NT_reformat:
    """Taxon step 09: Collect some summary statistics on primay nucleotide search"""
    input:
        queryDB = os.path.join(dir.out.primaryNT, "queryDB"),
        db = os.path.join(dir.dbs.primaryNT, "sequenceDB"),
        taxdb = dir.dbs.taxonomy
    output:
        result = os.path.join(dir.out.primaryNT, "results", "firsthit.index"),
        align = os.path.join(dir.out.primaryNT, "results", "result.m8"),
        lineage = os.path.join(dir.out.primaryNT, "primary_nt.lineage"),
        reformated = os.path.join(dir.out.primaryNT, "primary_nt.tsv"),
    params:
        inputpath = os.path.join(dir.out.primaryNT, "results", "result"),
        respath = os.path.join(dir.out.primaryNT, "results", "firsthit"),
        taxonFormat = lambda wildcards: config.immutable.taxonkitReformat
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    resources:
        time = config.resources.sml.time
    benchmark:
        os.path.join(dir.out.bench, "PRIMARY_NT_reformat.txt")
    log:
        os.path.join(dir.out.stderr, 'PRIMARY_NT_reformat.log')
    group:
        "primarynt"
    shell:
        """
        {{ # Filter TopHit results
        mmseqs filterdb {params.inputpath} {params.respath} --extract-lines 1;
        # Convert to alignments
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align};
        # Assign taxonomy
        cut -f2 {output.align} | \
            awk -F '|' '{{ print$2 }}' | \
            taxonkit lineage --data-dir {input.taxdb} > {output.lineage};
        # Reformat TopHit viral lineage information
        taxonkit reformat --data-dir {input.taxdb} {output.lineage} -i 2 {params.taxonFormat} | \
            cut --complement -f2 > {output.reformated};
        }} &> {log}
        rm {log}
        """

        
rule PRIMARY_NT_parsing:
    """Taxon step 10: Extract unclassified sequences from secondary translated search
    
    xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
    """
    input:
        seqs = os.path.join(dir.out.primaryAA, "MMSEQS_AA_PRIMARY_unclassified.fasta"),
        align = os.path.join(dir.out.primaryNT, "results", "result.m8")
    output:
        class_seqs = os.path.join(dir.out.primaryNT, "classified_seqs.fasta"),
        unclass_seqs = os.path.join(dir.out.primaryNT, "unclassified_seqs.fasta")
    resources:
        time = config.resources.sml.time
    benchmark:
        os.path.join(dir.out.bench, "PRIMARY_NT_parsing.txt")
    log:
        os.path.join(dir.out.stderr, 'PRIMARY_NT_parsing.log')
    group:
        "primarynt"
    script:
        os.path.join(dir.scripts,  'ntPrimaryParse.py')


rule SECONDARY_NT_taxonomic_assignment:
    """Taxon step 11: Secondary nucleotide search of viral hits from primary nucleotide search"""
    input:
        seqs = os.path.join(dir.out.primaryNT, "classified_seqs.fasta"),
        db = os.path.join(dir.dbs.secondaryNT, "sequenceDB")
    output:
        queryDB = os.path.join(dir.out.secondaryNT, "queryDB"),
        result = os.path.join(dir.out.secondaryNT, "results", "result.index"),
        type = os.path.join(dir.out.secondaryNT, "results", "result.dbtype")
    params:
        respath = os.path.join(dir.out.secondaryNT, "results", "result"),
        tmppath = os.path.join(dir.out.secondaryNT, "mmseqs_aa_tmp"),
        ntfilt = config.mmseqs.filtNTsecondary,
        sensnt = config.mmseqs.sensNT,
        memsplit = str(int(0.75 * int(config.resources.big.mem))) + 'M'
    benchmark:
        os.path.join(dir.out.bench, "SECONDARY_NT_taxonomic_assignment.txt")
    log:
        log = os.path.join(dir.out.stderr, "SECONDARY_NT_taxonomic_assignment.log")
    resources:
        mem_mb=config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    group:
        "secondarynt"
    shell:
        """
        {{ # Create query database
        rm -f {output.type}
        mmseqs createdb {input.seqs} {output.queryDB} --dbtype 2
        # mmseqs search
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
            {params.sensnt} {params.ntfilt} \
            --search-type 3 --threads {threads} --split-memory-limit {params.memsplit};
        }} &> {log}
        rm {log}
        """


rule SECONDARY_NT_summary:
    """Taxon step 12: Summary statistics for secondary nucleotide search"""
    input:
        queryDB = os.path.join(dir.out.secondaryNT, "queryDB"),
        db = os.path.join(dir.dbs.secondaryNT, "sequenceDB"),
        taxdb = dir.dbs.taxonomy
    output:
        result = os.path.join(dir.out.secondaryNT, "results", "tophit.index"),
        align = os.path.join(dir.out.secondaryNT, "results", "tophit.m8"),
        lineage = os.path.join(dir.out.secondaryNT, "SECONDARY_nt.lineage"),
        reformated = os.path.join(dir.out.secondaryNT, "SECONDARY_nt.tsv"),
    params:
        inputpath = os.path.join(dir.out.secondaryNT, "results", "result"),
        respath = os.path.join(dir.out.secondaryNT, "results", "tophit"),
        taxonFormat = lambda wildcards: config.immutable.taxonkitReformat,
        convertAli = config.immutable.mmseqConvertAliFormat
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    resources:
        time = config.resources.sml.time
    benchmark:
        os.path.join(dir.out.bench, "SECONDARY_NT_summary.txt")
    log:
        os.path.join(dir.out.stderr, 'SECONDARY_NT_summary.log')
    group:
        "secondarynt"
    shell:
        """
        {{ # Filter TopHit results
        mmseqs filterdb {params.inputpath} {params.respath} --extract-lines 1;
        # Convert to alignments
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} \
            --format-output {params.convertAli};
        # Assign taxonomy
        cut -f1,2 {output.align} \
            | sed 's/tid|//' \
            | sed 's/|.*//' \
            | taxonkit lineage -i 2 --data-dir {input.taxdb} > {output.lineage};
        # Reformat TopHit viral lineage information
        taxonkit reformat --data-dir {input.taxdb} {output.lineage} -i 3 {params.taxonFormat} | \
            cut --complement -f3 \
            > {output.reformated};
        }} &> {log}
        rm {log}
        """

        
rule SECONDARY_NT_convert:
    """Taxon step 13: Reformat secondary search LCA taxon assignment"""
    input:
        queryDB = os.path.join(dir.out.secondaryNT, "queryDB"),
        db = os.path.join(dir.dbs.secondaryNT, "sequenceDB")
    output:
        align = temp(os.path.join(dir.out.secondaryNT, "results", "all.m8")),
    params:
        respath = os.path.join(dir.out.secondaryNT, "results", "result")
    resources:
        time = config.resources.sml.time
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    benchmark:
        os.path.join(dir.out.bench, "SECONDARY_NT_convert.txt")
    log:
        os.path.join(dir.out.stderr, 'SECONDARY_NT_convert.log')
    group:
        "secondarynt"
    shell:
        """
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} \
            --format-output "query,target" \
            2> {log}
        rm {log}
        """


rule secondary_nt_lca_table:
    """Taxon step 14: Create table for taxonkit lineage for secondary NT search"""
    input:
        align = os.path.join(dir.out.secondaryNT, "results", "all.m8")
    output:
        lin = os.path.join(dir.out.secondaryNT, "results", "all.lin")
    benchmark:
        os.path.join(dir.out.bench, "secondary_nt_lca_table.txt")
    log:
        os.path.join(dir.out.stderr, 'secondary_nt_lca_table.log')
    resources:
        time = config.resources.sml.time
    group:
        "secondarynt"
    script:
        os.path.join(dir.scripts,  'ntSecondaryLca.py')


rule secondary_nt_calc_lca:
    """Taxon step 15: Calculate the lca for the secondary NT search"""
    input:
        lin = os.path.join(dir.out.secondaryNT, "results", "all.lin"),
        taxdb = dir.dbs.taxonomy
    output:
        lca_lineage = os.path.join(dir.out.secondaryNT, "results", "lca.lineage"),
        reformated = os.path.join(dir.out.secondaryNT, "results", "secondary_nt_lca.tsv")
    resources:
        time = config.resources.sml.time
    params:
        taxonFormat = lambda wildcards: config.immutable.taxonkitReformat,
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    benchmark:
        os.path.join(dir.out.bench, "secondary_nt_calc_lca.txt")
    log:
        os.path.join(dir.out.stderr, 'secondary_nt_calc_lca.log')
    group:
        "secondarynt"
    shell:
        """
        {{
        # calculate lca and lineage
        taxonkit lca -i 2 -s ';' --data-dir {input.taxdb} {input.lin} | \
            taxonkit lineage -i 3 --data-dir {input.taxdb} | \
            cut --complement -f 2 \
            > {output.lca_lineage} 2> {log}
        
        # Reformat lineages
        awk -F '\t' '$2 != 0' {output.lca_lineage} | \
            taxonkit reformat --data-dir {input.taxdb} -i 3 {params.taxonFormat} 2>> {log} |
            cut --complement -f3 > {output.reformated}
        }} &> {log}
        rm {log}
        """


rule SECONDARY_NT_generate_output_table:
    """Taxon step 16: Create the NT portion of the bigtable"""
    input:
        aln = os.path.join(dir.out.secondaryNT, "results", "tophit.m8"),
        top = os.path.join(dir.out.secondaryNT, "SECONDARY_nt.tsv"),
        lca = os.path.join(dir.out.secondaryNT, "results", "secondary_nt_lca.tsv"),
        counts = os.path.join(dir.out.results,"sampleSeqCounts.tsv"),
        balt = os.path.join(dir.dbs.tables,'2020_07_27_Viral_classification_table_ICTV2019.txt')
    output:
        vir = os.path.join(dir.out.secondaryNT, "NT_bigtable.tsv"),
        nonvir = os.path.join(dir.out.secondaryNT, "NT_bigtable.nonviral.tsv")
    resources:
        time = config.resources.sml.time
    params:
        taxIdIgnore = config.mmseqs.taxIdIgnore.split(),
        bigtableHeader = config.immutable.bigtableHeader
    benchmark:
        os.path.join(dir.out.bench, "SECONDARY_NT_generate_output_table.txt")
    log:
        os.path.join(dir.out.stderr, 'SECONDARY_NT_generate_output_table.log')
    group:
        "secondarynt"
    script:
        os.path.join(dir.scripts,  'ntBigtable.py')


rule combine_AA_NT:
    """Taxon step 17: Combine the AA and NT bigtables"""
    input:
        aa = os.path.join(dir.out.secondaryAA,"AA_bigtable.tsv"),
        nt = os.path.join(dir.out.secondaryNT, "NT_bigtable.tsv")
    output:
        os.path.join(dir.out.results, "bigtable.tsv")
    benchmark:
        os.path.join(dir.out.bench, "combine_AA_NT.txt")
    log:
        os.path.join(dir.out.stderr, 'combine_AA_NT.log')
    resources:
        time = config.resources.sml.time
    group:
        "secondarynt"
    shell:
        """
        {{ cat {input.aa} > {output};
        tail -n+2 {input.nt} >> {output}; }} &> {log}
        rm {log}
        """
