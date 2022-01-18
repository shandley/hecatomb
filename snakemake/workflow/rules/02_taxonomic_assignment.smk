"""
Hecatomb.smk to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

History: This is based on [mmseqs_pviral_aa.sh](../base/mmseqs_pviral_aa.sh)

Rob Edwards, March 2020
Overhauled - Michael Roach, Q2/Q3 2021
"""


rule PRIMARY_AA_taxonomy_assignment:
    """Taxon step 01: Assign taxonomy to RESULTS/seqtable.fasta sequences using mmseqs2 with primary AA database"""
    input:
        seqs = os.path.join(RESULTS, "seqtable.fasta"),
        db = os.path.join(UNIVIRDB, "sequenceDB")
    output:
        lca = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.tsv"),
        report = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_report"),
        tophit_report = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_report"),
        aln = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln"),
        alnsort = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted")
    params:
        alnRes = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY"),
        tmppath = os.path.join(PRIMARY_AA_OUT, "mmseqs_aa_tmp")
    benchmark:
        os.path.join(BENCH, "PRIMARY_AA_taxonomy_assignment.txt")
    log:
        os.path.join(STDERR, "PRIMARY_AA_taxonomy_assignment.log")
    resources:
        mem_mb = MMSeqsMem,
        time = MMSeqsTimeMin
    threads:
        MMSeqsCPU
    conda:
        os.path.join("..", "envs", "mmseqs2.yaml")
    shell:
        """
        {{ # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
            {config[filtAAprimary]} {MMSeqsSensAA} {config[reqAA]} \
            --threads {threads} --split-memory-limit {MMSeqsMemSplit};
            
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        }} &> {log}
        rm {log}
        """


rule PRIMARY_AA_parsing:
    """Taxon step 02: Parse primary AA search results for classified (potentially viral) and unclassified sequences"""
    input:
        alnsort = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted"),
        seqs = os.path.join(RESULTS, "seqtable.fasta"),
    output:
        class_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
    resources:
        mem_mb = MiscMem
    threads:
        MiscCPU
    benchmark:
        os.path.join(BENCH, 'PRIMARY_AA_parsing.txt')
    log:
        os.path.join(STDERR, 'PRIMARY_AA_parsing.log')
    script:
        os.path.join('..', 'scripts', 'aaPrimaryParse.py')


rule SECONDARY_AA_taxonomy_assignment:
    """Taxon step 03: Check taxonomic assignments in MMSEQS_AA_PRIMARY_classified.fasta using mmseqs2"""
    input:
        seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
        db = os.path.join(UNIREF50VIR, "sequenceDB")
    output:
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv"),
        report = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_report"),
        tophit_report = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_report"),
        aln = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln"),
        alnsort = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln_sorted")
    params:
        alnRes=os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY"),
        tmppath=os.path.join(SECONDARY_AA_OUT, "mmseqs_aa_tmp")
    benchmark:
        os.path.join(BENCH, "SECONDARY_AA_taxonomy_assignment.txt")
    log:
        os.path.join(STDERR, "SECONDARY_AA_taxonomy_assignment.log")
    resources:
        mem_mb = MMSeqsMem,
        time = MMSeqsTimeMin
    threads:
        MMSeqsCPU
    conda:
        os.path.join("..", "envs", "mmseqs2.yaml")
    shell:
        """
        {{ # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
            {config[filtAAsecondary]} {MMSeqsSensAA} {config[reqAA]} \
            --threads {threads} --split-memory-limit {MMSeqsMemSplit};
            
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        }} &> {log}
        rm {log}
        """


rule SECONDARY_AA_tophit_lineage:
    """Taxon step 04: Add/reformat tophit viral lineages with up-to-date* NCBI taxonomy"""
    input:
        db = TAX,
        tophit = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln")
    output:
        tophit_lineage_refomated = os.path.join(SECONDARY_AA_OUT, "tophit.lineage.reformated"),
    conda:
        os.path.join("..", "envs", "seqkit.yaml")
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    benchmark:
        os.path.join(BENCH, "SECONDARY_AA_tophit_lineage.txt")
    log:
        os.path.join(STDERR, "SECONDARY_AA_tophit_lineage.log")
    shell:
        """
        {{ # Make a table: SeqID <tab> taxID
        cut -f1,20 {input.tophit} \
            | taxonkit lineage --data-dir {input.db} -i 2 \
            | taxonkit reformat --data-dir {input.db} -i 3 \
                -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F --fill-miss-rank \
            | cut --complement -f3 \
            > {output.tophit_lineage_refomated}; }} &> {log}
        rm {log}
        """


rule SECONDARY_AA_refactor_finalize:
    """Taxon step 05: Remove sequences to be refactored from LCA table and recombine with updated taxonomies."""
    input:
        db = TAX,
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv"),
    output:
        lca_reformated = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.reformated"),
    conda:
        os.path.join("..", "envs", "seqkit.yaml")
    resources:
        mem_mb = MiscMem
    threads:
        MiscCPU
    benchmark:
        os.path.join(BENCH, "SECONDARY_AA_refactor_finalize.txt")
    log:
        os.path.join(STDERR, "SECONDARY_AA_refactor_finalize.log")
    shell:
        """
        {{ cut -f1,2 {input.lca} \
            | taxonkit lineage --data-dir {input.db} -i 2 \
            | taxonkit reformat --data-dir {input.db} -i 3 \
                -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F --fill-miss-rank \
            | cut --complement -f3 \
            > {output.lca_reformated}; }} &> {log}
        rm {log}
        """


rule SECONDARY_AA_generate_output_table:
    """Taxon step 06: Join sequence info, tophit align info, and LCA or tophit lineage info into the output format table"""
    input:
        aln = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln"),
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.reformated"),
        top = os.path.join(SECONDARY_AA_OUT, "tophit.lineage.reformated"),
        counts = os.path.join(RESULTS, "sampleSeqCounts.tsv"),
        balt = os.path.join(TABLES, '2020_07_27_Viral_classification_table_ICTV2019.txt')
    output:
        os.path.join(SECONDARY_AA_OUT, "AA_bigtable.tsv")
    resources:
        mem_mb = MiscMem
    threads:
        MiscCPU
    benchmark:
        os.path.join(BENCH, "SECONDARY_AA_generate_output_table.txt")
    log:
        os.path.join(STDERR, "SECONDARY_AA_generate_output_table.log")
    params:
        taxIdIgnore = config['taxIdIgnore'].split()
    script:
        os.path.join('..', 'scripts', 'aaBigtable.py')


rule SECONDARY_AA_parsing:
    """Taxon step 07: Parse out all sequences that remain unclassified following the Secondary AA search and refactoring."""
    input:
        bigtable = os.path.join(SECONDARY_AA_OUT,"AA_bigtable.tsv"),
        seqs = os.path.join(RESULTS, "seqtable.fasta")
    output:
        unclass_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta")
    resources:
        mem_mb = MiscMem
    threads:
        MiscCPU
    benchmark:
        os.path.join(BENCH, "SECONDARY_AA_parsing.txt")
    log:
        os.path.join(STDERR, 'SECONDARY_AA_parsing.log')
    script:
        os.path.join('..', 'scripts', 'aaSecondaryParse.py')

        
rule PRIMARY_NT_taxonomic_assignment:
    """Taxon step 08: Primary nucleotide search of unclassified viral-like sequences from aa search"""
    input:
        seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta"),
        db = os.path.join(NCBIVIRDB, "sequenceDB")
    output:
        queryDB = os.path.join(PRIMARY_NT_OUT, "queryDB"),
        result = os.path.join(PRIMARY_NT_OUT, "results", "result.index"),
        type = os.path.join(PRIMARY_NT_OUT, "results", "result.dbtype")
    params:
        respath = os.path.join(PRIMARY_NT_OUT, "results", "result"),
        tmppath = os.path.join(PRIMARY_NT_OUT, "mmseqs_aa_tmp")
    benchmark:
        os.path.join(BENCH, "PRIMARY_NT_taxonomic_assignment.txt")
    log:
        os.path.join(STDERR, "PRIMARY_NT_taxonomic_assignment.log")
    resources:
        mem_mb = MMSeqsMem,
        time = MMSeqsTimeMin
    threads:
        MMSeqsCPU
    conda:
        os.path.join("..", "envs", "mmseqs2.yaml")
    shell:
        """
        {{ # Create query database
        rm -f {output.type}
        mmseqs createdb {input.seqs} {output.queryDB} --dbtype 2;
        # mmseqs search
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
            {MMSeqsSensNT} {config[filtNTprimary]} \
            --search-type 3 --threads {threads} --split-memory-limit {MMSeqsMemSplit};
        }} &> {log}
        rm {log}
        """


rule PRIMARY_NT_reformat:
    """Taxon step 09: Collect some summary statistics on primay nucleotide search"""
    input:
        queryDB = os.path.join(PRIMARY_NT_OUT, "queryDB"),
        db = os.path.join(NCBIVIRDB, "sequenceDB"),
        taxdb = TAX
    output:
        result = os.path.join(PRIMARY_NT_OUT, "results", "firsthit.index"),
        align = os.path.join(PRIMARY_NT_OUT, "results", "result.m8"),
        lineage = os.path.join(PRIMARY_NT_OUT, "primary_nt.lineage"),
        reformated = os.path.join(PRIMARY_NT_OUT, "primary_nt.tsv"),
    params:
        inputpath = os.path.join(PRIMARY_NT_OUT, "results", "result"),
        respath = os.path.join(PRIMARY_NT_OUT, "results", "firsthit")
    conda:
        os.path.join("..", "envs", "mmseqs2.yaml")
    resources:
        mem_mb = MiscMem
    threads:
        MiscCPU
    benchmark:
        os.path.join(BENCH, "PRIMARY_NT_reformat.txt")
    log:
        os.path.join(STDERR, 'PRIMARY_NT_reformat.log')
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
        taxonkit reformat --data-dir {input.taxdb} {output.lineage} -i 2 \
            -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F --fill-miss-rank |
            cut --complement -f2 > {output.reformated};
        }} &> {log}
        rm {log}
        """

        
rule PRIMARY_NT_parsing:
    """Taxon step 10: Extract unclassified sequences from secondary translated search
    
    xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
    """
    input:
        seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta"),
        align = os.path.join(PRIMARY_NT_OUT, "results", "result.m8")
    output:
        class_seqs = os.path.join(PRIMARY_NT_OUT, "classified_seqs.fasta"),
        unclass_seqs = os.path.join(PRIMARY_NT_OUT, "unclassified_seqs.fasta")
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    benchmark:
        os.path.join(BENCH, "PRIMARY_NT_parsing.txt")
    log:
        os.path.join(STDERR, 'PRIMARY_NT_parsing.log')
    script:
        os.path.join('..', 'scripts', 'ntPrimaryParse.py')


rule SECONDARY_NT_taxonomic_assignment:
    """Taxon step 11: Secondary nucleotide search of viral hits from primary nucleotide search"""
    input:
        seqs = os.path.join(PRIMARY_NT_OUT, "classified_seqs.fasta"),
        db = os.path.join(POLYMICRODB, "sequenceDB")
    output:
        queryDB = os.path.join(SECONDARY_NT_OUT, "queryDB"),
        result = os.path.join(SECONDARY_NT_OUT, "results", "result.index"),
        type = os.path.join(SECONDARY_NT_OUT, "results", "result.dbtype")
    params:
        respath = os.path.join(SECONDARY_NT_OUT, "results", "result"),
        tmppath = os.path.join(SECONDARY_NT_OUT, "mmseqs_aa_tmp")
    benchmark:
        os.path.join(BENCH, "SECONDARY_NT_taxonomic_assignment.txt")
    log:
        log = os.path.join(STDERR, "SECONDARY_NT_taxonomic_assignment.log")
    resources:
        mem_mb=MMSeqsMem,
        time=MMSeqsTimeMin
    threads:
        MMSeqsCPU
    conda:
        os.path.join("..", "envs", "mmseqs2.yaml")
    shell:
        """
        {{ # Create query database
        rm -f {output.type}
        mmseqs createdb {input.seqs} {output.queryDB} --dbtype 2
        # mmseqs search
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
            {MMSeqsSensNT} {config[filtNTsecondary]} \
            --search-type 3 --threads {threads} --split-memory-limit {MMSeqsMemSplit};
        }} &> {log}
        rm {log}
        """


rule SECONDARY_NT_summary:
    """Taxon step 12: Summary statistics for secondary nucleotide search"""
    input:
        queryDB = os.path.join(SECONDARY_NT_OUT, "queryDB"),
        db = os.path.join(POLYMICRODB, "sequenceDB"),
        taxdb = TAX
    output:
        result = os.path.join(SECONDARY_NT_OUT, "results", "tophit.index"),
        align = os.path.join(SECONDARY_NT_OUT, "results", "tophit.m8"),
        lineage = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt.lineage"),
        reformated = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt.tsv"),
    params:
        inputpath = os.path.join(SECONDARY_NT_OUT, "results", "result"),
        respath = os.path.join(SECONDARY_NT_OUT, "results", "tophit")
    conda:
        os.path.join("..", "envs", "mmseqs2.yaml")
    resources:
        mem_mb = MiscMem
    threads:
        MiscCPU
    benchmark:
        os.path.join(BENCH, "SECONDARY_NT_summary.txt")
    log:
        os.path.join(STDERR, 'SECONDARY_NT_summary.log')
    shell:
        """
        {{ # Filter TopHit results
        mmseqs filterdb {params.inputpath} {params.respath} --extract-lines 1;
        # Convert to alignments
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} \
            --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader";
        # Assign taxonomy
        cut -f1,2 {output.align} \
            | sed 's/tid|//' \
            | sed 's/|.*//' \
            | taxonkit lineage -i 2 --data-dir {input.taxdb} > {output.lineage};
        # Reformat TopHit viral lineage information
        taxonkit reformat --data-dir {input.taxdb} {output.lineage} -i 3 \
            -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F --fill-miss-rank |
            cut --complement -f3 \
            > {output.reformated};
        }} &> {log}
        rm {log}
        """

        
rule SECONDARY_NT_convert:
    """Taxon step 13: Reformat secondary search LCA taxon assignment"""
    input:
        queryDB = os.path.join(SECONDARY_NT_OUT, "queryDB"),
        db = os.path.join(POLYMICRODB, "sequenceDB")
    output:
        align = os.path.join(SECONDARY_NT_OUT, "results", "all.m8"),
    params:
        respath = os.path.join(SECONDARY_NT_OUT, "results", "result")
    resources:
        mem_mb = MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        os.path.join("..", "envs", "mmseqs2.yaml")
    benchmark:
        os.path.join(BENCH, "SECONDARY_NT_convert.txt")
    log:
        os.path.join(STDERR, 'SECONDARY_NT_convert.log')
    shell:
        """
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} \
            --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader" \
            2> {log}
        rm {log}
        """


rule secondary_nt_lca_table:
    """Taxon step 14: Create table for taxonkit lineage for secondary NT search"""
    input:
        align = os.path.join(SECONDARY_NT_OUT, "results", "all.m8")
    output:
        lin = os.path.join(SECONDARY_NT_OUT, "results", "all.lin")
    benchmark:
        os.path.join(BENCH, "secondary_nt_lca_table.txt")
    log:
        os.path.join(STDERR, 'secondary_nt_lca_table.log')
    resources:
        mem_mb = MiscMem
    script:
        os.path.join('..', 'scripts', 'ntSecondaryLca.py')


rule secondary_nt_calc_lca:
    """Taxon step 15: Calculate the lca for the secondary NT search"""
    input:
        lin = os.path.join(SECONDARY_NT_OUT, "results", "all.lin"),
        taxdb = TAX
    output:
        lca_lineage = os.path.join(SECONDARY_NT_OUT, "results", "lca.lineage"),
        reformated = os.path.join(SECONDARY_NT_OUT, "results", "secondary_nt_lca.tsv")
    resources:
        mem_mb = MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        os.path.join("..", "envs", "mmseqs2.yaml")
    benchmark:
        os.path.join(BENCH, "secondary_nt_calc_lca.txt")
    log:
        os.path.join(STDERR, 'secondary_nt_calc_lca.log')
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
            taxonkit reformat --data-dir {input.taxdb} -i 3 \
                -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F --fill-miss-rank 2>> {log} |
            cut --complement -f3 > {output.reformated}
        }} &> {log}
        rm {log}
        """


rule SECONDARY_NT_generate_output_table:
    """Taxon step 16: Create the NT portion of the bigtable"""
    input:
        aln = os.path.join(SECONDARY_NT_OUT, "results", "tophit.m8"),
        top = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt.tsv"),
        lca = os.path.join(SECONDARY_NT_OUT, "results", "secondary_nt_lca.tsv"),
        counts = os.path.join(RESULTS,"sampleSeqCounts.tsv"),
        balt = os.path.join(TABLES,'2020_07_27_Viral_classification_table_ICTV2019.txt')
    output:
        os.path.join(SECONDARY_NT_OUT, "NT_bigtable.tsv")
    resources:
        mem_mb = MiscMem
    threads:
        MiscCPU
    params:
        taxIdIgnore = config['taxIdIgnore'].split()
    benchmark:
        os.path.join(BENCH, "SECONDARY_NT_generate_output_table.txt")
    log:
        os.path.join(STDERR, 'SECONDARY_NT_generate_output_table.log')
    script:
        os.path.join('..', 'scripts', 'ntBigtable.py')


rule combine_AA_NT:
    """Taxon step 17: Combine the AA and NT bigtables"""
    input:
        aa = os.path.join(SECONDARY_AA_OUT,"AA_bigtable.tsv"),
        nt = os.path.join(SECONDARY_NT_OUT, "NT_bigtable.tsv")
    output:
        os.path.join(RESULTS, "bigtable.tsv")
    benchmark:
        os.path.join(BENCH, "combine_AA_NT.txt")
    log:
        os.path.join(STDERR, 'combine_AA_NT.log')
    shell:
        """
        {{ cat {input.aa} > {output};
        tail -n+2 {input.nt} >> {output}; }} &> {log}
        rm {log}
        """
