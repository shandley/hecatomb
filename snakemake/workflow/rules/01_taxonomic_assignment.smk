"""
Snakefile to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

This is based on [mmseqs_pviral_aa.sh](../base/mmseqs_pviral_aa.sh)

REQUIRES that targetDB has already been indexed
If it has not been index then run the following script in the directory of your choice: uniprot_viral_DB_build.sh (found in /accessory)                                                                  
Note: mmseqs2 taxonomy is currently most useful if you have UniProt formatted fasta databases
more details about database building can be found at: https://github.com/soedinglab/mmseqs2/wiki#taxonomy-assignment-using-mmseqs-taxonomy                                                               


Rob Edwards, March 2020

"""

import os
import sys


rule mmseqs_pviral_aa:
    input:
        os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.tsv"),
        os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_report"),
        os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_report"),
        os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln"),
        os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted")
        #os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.ids"),
        #os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted.ids"),
        #os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.ids"),
        #os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
        #os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta"),
        #os.path.join(PRIMARY_AA_OUT, "mmseqs_primary_aa_search_summary.txt")

rule seqtable_translated_taxonomy_PRIMARY_search:
    """
    Assign taxonomy to RESULTS/seqtable.fasta sequences using mmseqs2 to query a viral protein reference databases.
    
    This is a nt-to-aa translated search, similar to blastX.
    
    All sequences assigned to a viral lineage will be reserved for the secondary query against a trans-kingdom database to confirm their viral lineage.
    
    Reference database: All Uniprot viral proteins clustered at 99% identity.
    
    Note: The sequences checked in the next step are taken from the tophit_aln file. This means they will not have LCA assigned taxonomy, but it is the most inclusive group to send to the next step where any false-positives will be sorted out. For now, we want to capture anything that potentially looks viral.

    """
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
        alnRes=os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY"),
        tmppath=os.path.join(PRIMARY_AA_OUT, "mmseqs_aa_tmp")
    benchmark:
        "BENCHMARKS/mmseqs_primary_aa.txt"
    log:
        "LOGS/mmseqs/mmseqs_primary_aa.log"
    resources:
        mem_mb=64000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
        -a --start-sens 1 --sens-steps 3 -s 7 --search-type 2 --tax-output-mode 2 --lca-mode 4 --shuffle 0 \
        --lca-ranks "superkingdom,phylum,class,order,family,genus,species" \
        --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
        --tax-lineage 1 \
        -e 0.1;
        
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        """
        
rule seqtable_translated_taxonomy_PRIMARY_tabulation:
    """
    1) Parse potential viral and unclassified sequences
    
    2) Tabulate some statistics
        - Number of input sequences: seqtable.fasta
        - Number of sequences classified as having a viral lineage: MMSEQS_AA_PRIMARY_classified.fasta
        - Number of sequences having no inidication of having a viral lineage: MMSEQS_AA_PRIMARY_unclassified.fasta
        - Number of sequences with a viral LCA lineage
        - Number of sequences with a viral tophit lineage
        - etc.
        
    MMSEQS_AA_PRIMARY_classified.fasta will be subjected to a secondary translated search against a trans-kingdom database to confirm true virality.
    
    MMSEQS_AA_PRIMARY_unclassified.fasta will be queried against a viral nucleotide database to detect similarity to non-coding regions of viral genomes or to sequences not represented in the UniProt protein databases.
        
    """
    input:
        lca = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.tsv"),
        alnsort = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted"),
        seqs = os.path.join(RESULTS, "seqtable.fasta")
    output:
        allseq_ids = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.ids"),
        tophit_ids = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted.ids"),
        unclass_ids = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.ids"),
        class_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
        unclass_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta"),
        summary = os.path.join(RESULTS, "SEARCH_SUMMARY", "MMSEQS_AA_PRIMARY_summary.txt"),
        virfam = os.path.join(RESULTS, "SEARCH_SUMMARY", "MMSEQS_AA_PRIMARY_virus_family_summary.txt"),
        incomplete_lineage = os.path.join(RESULTS, "SEARCH_SUMMARY", "MMSEQS_AA_PRIMARY_incomplete_lineage.txt")
    conda:
        "../envs/samtools.yaml"
    log:
        "LOGS/mmseqs/mmseqs_primary_aa_parsing.log"
    shell:
        """
        # Extract full entry ID list (all sequences)
        # Note: the lca output table has one entry per input sequence
        sort -k1 -n {input.lca} | cut -f1 > {output.allseq_ids};
        
        # Extract tophit entry ids
        tail -n+2 {input.alnsort} | cut -f1 > {output.tophit_ids};
        
        # Compare tophit ids to all ids to extract ids with no hits
        comm -23 <(sort {output.allseq_ids}) <(sort {output.tophit_ids}) | sort -n > {output.unclass_ids};
        
        # Extract classified sequencess (Antyhing with a tophit. Just any signal it may be viral at this point)
        xargs samtools faidx {input.seqs} -n 5000 < {output.tophit_ids} > {output.class_seqs};
        
        # Extract unclassified sequencess
        xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
        
        # Viral family summary
        grep "d_Viruses" {input.lca} | cut -f5 | awk -F ';' '{{ print$5 }}' | sort | uniq -c | sort -k1 -nr > {output.virfam};
        
        # Taxonomy IDs without full lineage information in the reference database
        tail -n+2 {input.alnsort} | grep -v "d_Viruses" | awk -F '\t' '{{ print$20"\\t"$19 }}' | sort -k1 -nr | uniq -c | sort -k1 -nr | sed '1iCount\tTaxID\tUniProt_ID' > {output.incomplete_lineage};
        
        ## Summary table variables
        # Calculate number of input sequences (INPUT) from seqtable
        INPUT=$(grep -c ">" {input.seqs})
        
        # Calculate number of classified sequences
        CLASS=$(grep -c ">" {output.class_seqs})
        
        # Calculate number of unclassified sequences
        UNCLASS=$(grep -c ">" {output.unclass_seqs})
        
        # Calculate % of classified (pCLASS) and 9(pUNCLASS) unclassified sequences
        pCLASS=$(echo "scale=4 ; 100*($CLASS / $INPUT)" | bc)
        pUNCLASS=$(echo "scale=4 ; 100*($UNCLASS / $INPUT)" | bc)
        
        # Calculate the number of LCA calls with viral lineage
        LCA_ALL=$(wc -l {input.lca})
        LCA_VIRAL=$(grep -c "d_Viruses" {input.lca})
        LCA_NONVIRAL=$(grep -cv "d_Viruses" {input.lca})
        pLCA_VIRAL=$(echo "scale=4 ; 100*($LCA_VIRAL / $INPUT)" | bc)
        pLCA_NONVIRAL=$(echo "scale=4 ; 100*($LCA_NONVIRAL / $INPUT)" | bc)
        
        # Calculate the number of tophits with complete viral lineages
        TH_ALL=$(tail -n+2 {input.alnsort} | wc -l)
        pTH_ALL=$(echo "scale=4 ; 100*($TH_ALL / $INPUT)" | bc)
        TH_COMPLETE=$(tail -n+2 {input.alnsort} | grep -c "d_Viruses")
        TH_INCOMPLETE=$(tail -n+2 {input.alnsort} | grep -cv "d_Viruses")
        pTH_COMPLETE=$(echo "scale=4 ; 100*($TH_COMPLETE / $INPUT)" | bc)
        pTH_INCOMPLETE=$(echo "scale=6 ; 100*($TH_INCOMPLETE / $INPUT)" | bc)
        
        # Create summary table
        # NOTE: Here docs are explicit and will be formated as shown. Do not add tabs in front of doc as you normally would in snakemake code! Format it exactly how you want to see it!
        cat << EOF > {output.summary}
PARAMATER    VALUE    PERCENT    NOTES
Input    $INPUT    NA    Number of sequences input into primary mmseqs amino acid search
Classified    $CLASS    $pCLASS    Number and % of sequences classified with a viral (== d_Viruses) lineage
Unclassified    $UNCLASS    $pUNCLASS    Number and % of sequences classified without a viral (!= d_Viruses) lineage
Viral LCA    $LCA_VIRAL    $pLCA_VIRAL    Number and % of sequences with a viral lineage (== d_Viruses) after LCA calculation
Nonviral LCA    $LCA_NONVIRAL    $pLCA_NONVIRAL    Number and % of sequences without a viral lineage (!= d_Viruses) after LCA calculation
Viral TopHit    $TH_ALL    $pTH_ALL    Number and % of sequences with a tophit viral (== d_Viruses) lineage
Viral TopHit with complete lineage    $TH_COMPLETE    $pTH_COMPLETE    Number and % of sequences with complete viral lineage information in the reference taxonomy
Viral TopHit without complete lineage    $TH_INCOMPLETE    $pTH_INCOMPLETE    Number and % of sequences without complete viral lineage information in the reference taxonomy
        """

rule seqtable_translated_taxonomy_SECONDARY_search:
    """
    Check taxonomic assignments made in rule seqtable_translated_taxonomy_primary_search using mmseqs2 to a transkingdom reference amino acid database.
    
    Sequences assigned to a viral lineage ("d_Viruses") from prior rule are checked.
    
    This is a nt-to-aa translated search, similar to blastX.
    
    Reference database: UniRef50 combined with all Uniprot viral proteins clustered at 99% identity (UniRef50 + virus).
        - The origianl viral protein sequence are added back to ensure no viral proteins are lost in the constrction of the UniRef50 database
    
    """
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
        "BENCHMARKS/mmseqs_secondary_aa.txt"
    log:
        "LOGS/mmseqs/mmseqs_secondary_aa.log"
    resources:
        mem_mb=64000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
        -a --start-sens 1 --sens-steps 3 -s 7 --search-type 2 --tax-output-mode 2 --lca-mode 2 --shuffle 0 \
        --lca-ranks "superkingdom,phylum,class,order,family,genus,species" \
        --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
        --tax-lineage 1 \
        -e 0.01;
        
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        """
        
rule seqtable_translated_taxonomy_SECONDARY_tabulation:
    """

        
    """
    input:
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv"),
        alnsort = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln_sorted"),
        seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta")
    output:
        allseq_ids = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.ids"),
        tophit_ids = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln_sorted.ids"),
        unclass_ids = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_unclassified.ids"),
        class_seqs = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_classified.fasta"),
        unclass_seqs = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_unclassified.fasta"),
        summary = os.path.join(RESULTS, "SEARCH_SUMMARY", "MMSEQS_AA_SECONDARY_summary.txt"),
        virfam = os.path.join(RESULTS, "SEARCH_SUMMARY", "MMSEQS_AA_SECONDARY_virus_family_summary.txt"),
        incomplete_lineage = os.path.join(RESULTS, "SEARCH_SUMMARY", "MMSEQS_AA_SECONDARY_incomplete_lineage.txt")
    conda:
        "../envs/samtools.yaml"
    log:
        "LOGS/mmseqs/mmseqs_SECONDARY_aa_parsing.log"
    shell:
        """
        # Extract full entry ID list (all sequences)
        # Note: the lca output table has one entry per input sequence
        sort -k1 -n {input.lca} | cut -f1 > {output.allseq_ids};
        
        # Extract tophit entry ids
        tail -n+2 {input.alnsort} | cut -f1 > {output.tophit_ids};
        
        # Compare tophit ids to all ids to extract ids with no hits
        comm -23 <(sort {output.allseq_ids}) <(sort {output.tophit_ids}) | sort -n > {output.unclass_ids};
        
        # Extract classified sequencess
        xargs samtools faidx {input.seqs} -n 5000 < {output.tophit_ids} > {output.class_seqs};
        
        # Extract unclassified sequencess
        xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
        
        # Viral family summary
        grep "d_Viruses" {input.lca} | cut -f5 | awk -F ';' '{{ print$5 }}' | sort | uniq -c | sort -k1 -nr > {output.virfam};
        
        # Taxonomy IDs without full lineage information in the reference database
        tail -n+2 {input.alnsort} | grep -v "d_Viruses" | awk -F '\t' '{{ print$20"\\t"$19 }}' | sort -k1 -nr | uniq -c | sort -k1 -nr | sed '1iCount\tTaxID\tUniProt_ID' > {output.incomplete_lineage};
        
        ## Summary table variables
        # Calculate number of input sequences (INPUT) from seqtable
        INPUT=$(grep -c ">" {input.seqs})
        
        # Calculate number of classified sequences
        CLASS=$(grep -c ">" {output.class_seqs})
        
        # Calculate number of unclassified sequences
        UNCLASS=$(grep -c ">" {output.unclass_seqs})
        
        # Calculate % of classified (pCLASS) and 9(pUNCLASS) unclassified sequences
        pCLASS=$(echo "scale=4 ; 100*($CLASS / $INPUT)" | bc)
        pUNCLASS=$(echo "scale=4 ; 100*($UNCLASS / $INPUT)" | bc)
        
        # Calculate the number of LCA calls with viral lineage
        LCA_ALL=$(wc -l {input.lca})
        LCA_VIRAL=$(grep -c "d_Viruses" {input.lca})
        LCA_NONVIRAL=$(grep -cv "d_Viruses" {input.lca})
        pLCA_VIRAL=$(echo "scale=4 ; 100*($LCA_VIRAL / $INPUT)" | bc)
        pLCA_NONVIRAL=$(echo "scale=4 ; 100*($LCA_NONVIRAL / $INPUT)" | bc)
        
        # Calculate the number of tophits with complete viral lineages
        TH_ALL=$(tail -n+2 {input.alnsort} | wc -l)
        pTH_ALL=$(echo "scale=4 ; 100*($TH_ALL / $INPUT)" | bc)
        TH_COMPLETE=$(tail -n+2 {input.alnsort} | grep -c "d_Viruses")
        TH_INCOMPLETE=$(tail -n+2 {input.alnsort} | grep -cv "d_Viruses")
        pTH_COMPLETE=$(echo "scale=4 ; 100*($TH_COMPLETE / $INPUT)" | bc)
        pTH_INCOMPLETE=$(echo "scale=6 ; 100*($TH_INCOMPLETE / $INPUT)" | bc)
        
        # Create summary table
        # NOTE: Here docs are explicit and will be formated as shown. Do not add tabs in front of doc as you normally would in snakemake code! Format it exactly how you want to see it!
        cat << EOF > {output.summary}
PARAMATER    VALUE    PERCENT    NOTES
Input    $INPUT    NA    Number of sequences input into SECONDARY mmseqs amino acid search
Classified    $CLASS    $pCLASS    Number and % of sequences classified with a viral (== d_Viruses) lineage
Unclassified    $UNCLASS    $pUNCLASS    Number and % of sequences classified without a viral (!= d_Viruses) lineage
Viral LCA    $LCA_VIRAL    $pLCA_VIRAL    Number and % of sequences with a viral lineage (== d_Viruses) after LCA calculation
Nonviral LCA    $LCA_NONVIRAL    $pLCA_NONVIRAL    Number and % of sequences without a viral lineage (!= d_Viruses) after LCA calculation
Viral TopHit    $TH_ALL    $pTH_ALL    Number and % of sequences with a tophit viral (== d_Viruses) lineage
Viral TopHit with complete lineage    $TH_COMPLETE    $pTH_COMPLETE    Number and % of sequences with complete viral lineage information in the reference taxonomy
Viral TopHit without complete lineage    $TH_INCOMPLETE    $pTH_INCOMPLETE    Number and % of sequences without complete viral lineage information in the reference taxonomy
        """


