"""
Snakefile to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2                                                             

History: This is based on [mmseqs_pviral_aa.sh](../base/mmseqs_pviral_aa.sh)

Rob Edwards, March 2020

"""

import os
import sys


rule PRIMARY_AA_taxonomy_assignment:
    """
    Assign taxonomy to RESULTS/seqtable.fasta sequences using mmseqs2 
    
    - Reference database: all UniProt viral (UNIVIRDB) protein sequences clustered at 99% ID
        
    - This is a nt-to-aa translated search, similar to blastX
    
    - All sequences assigned to a viral lineage will be reserved for the secondary query against a trans-kingdom database (UniRef50 + UNIVIRDB) to confirm their viral lineage
    
    - The sequences checked in the next step are selected if their tophit is in a viral lineage (tophit_aln file). This means they may not have LCA assigned taxonomy, but it is the most inclusive group to send to the next step (secondary search) where any false-positives will be sorted out. For now, we want to capture anything that potentially looks viral

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
        os.path.join(BENCH, "MMSEQS", "mmseqs_primary_AA.txt")
    log:
        log = os.path.join(LOGS, "MMSEQS", "mmseqs_primary_AA.log")
    resources:
        mem_mb=64000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
        -a --start-sens 1 --sens-steps 3 -s 7 \
        --tax-output-mode 2 --search-type 2 --lca-mode 2 --shuffle 0 \
        --lca-ranks "superkingdom,phylum,class,order,family,genus,species" \
        --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
        --tax-lineage 1 \
        -e {config[PRIMAAE]} &>> {log};
        
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        """
        
rule PRIMARY_AA_parsing:
    """
    
    Parse primary AA search results for classified (potentially viral) and unclassified sequences
        
    - MMSEQS_AA_PRIMARY_classified.fasta will be subjected to a secondary translated search (rule SECONDARY_AA_taxonomy_assignment:) against a trans-kingdom database (UniRef50 + UNIVIRDB) to confirm true viral lineage
    
    - MMSEQS_AA_PRIMARY_unclassified.fasta will be queried against a viral nucleotide database (nt-to-nt) to detect similarity to non-coding regions of viral genomes or to sequences not represented in the UniProt protein databases
        
    """
    input:
        lca = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.tsv"),
        alnsort = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted"),
        seqs = os.path.join(RESULTS, "seqtable.fasta")
    output:
        allseq_ids = temporary(os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.ids")),
        tophit_ids = temporary(os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted.ids")),
        unclass_ids = temporary(os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.ids")),
        class_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
        unclass_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta")
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        # Extract full entry ID list (all sequences)
        # Note: the lca output table has one entry per input sequence
        sort -k1 -n {input.lca} | cut -f1 > {output.allseq_ids};
        
        # Extract tophit entry ids
        tail -n+2 {input.alnsort} | cut -f1 > {output.tophit_ids};
        
        # For the Primary search we want to take anything that looks potentially viral to the Secondary search
        # There are several reasons in which a sequence may not be classified due to LCA issues (dealt with in subsequent rules
        # Therefore, we will not rely on LCA to call potential viral hits
        # And will just extract anything that looks potentially viral as a tophit
        # LCA will be used more rigorously for the final calls after the Secondary search
        
        # Compare tophit ids to all ids to extract ids with no tophits to the virusDB
        comm -23 <(sort {output.allseq_ids}) <(sort {output.tophit_ids}) | sort -n > {output.unclass_ids};
        
        # Extract classified sequences (Antyhing with a tophit. Just any signal it may be viral at this point)
        xargs samtools faidx {input.seqs} -n 5000 < {output.tophit_ids} > {output.class_seqs};
        
        # Extract unclassified sequencess
        xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
    
        """

rule PRIMARY_AA_summary:
    """
    Tabulate some summary statistics from the primary AA mmseqs search against virusDB
    
        - Number of input sequences: seqtable.fasta
        - Number and % of sequences classified as having a viral lineage: MMSEQS_AA_PRIMARY_classified.fasta
        - Number and % of sequences having no inidication of having a viral lineage: MMSEQS_AA_PRIMARY_unclassified.fasta
        - Number and % of sequences with a viral LCA lineage
        - Number and % of sequences with a viral tophit lineage
        - Number and % of seqeunces with a TopHit lineage but lacking an LCA lineage
        - Potentially more to be added
    """
    input:
        lca = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.tsv"),
        tophit = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln"),
        seqs = os.path.join(RESULTS, "seqtable.fasta"),
        class_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
        unclass_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta")
    output:
        virord = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_virus_order_summary.tsv"),
        virfam = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_virus_family_summary.tsv"),
        virgen = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_virus_genus_summary.tsv"),
        virspe = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_virus_species_summary.tsv"),
        tophit_only_ids = temporary(os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_only.ids")),
        tophit_only_tsv = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_only.tsv"),
        tmpsummary = temporary(os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_summary.tmp")),
        tophit_broken = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_broken_taxID.tsv"),
        summary = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_summary.tsv")
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        # Viral order summary
        grep "d_Viruses" {input.lca} | cut -f5 | awk -F ';' '{{ print$4 }}' | csvtk freq -H -n -r -T -t | sed '1i Family\tPrimary_AA_Order_Frequency' > {output.virord};
        
        # Viral family summary
        grep "d_Viruses" {input.lca} | cut -f5 | awk -F ';' '{{ print$5 }}' | csvtk freq -H -n -r -T -t | sed '1i Family\tPrimary_AA_Family_Frequency' > {output.virfam};
        
        # Viral genus summary
        grep "d_Viruses" {input.lca} | cut -f5 | awk -F ';' '{{ print$6 }}' | csvtk freq -H -n -r -T -t | sed '1i Family\tPrimary_AA_Order_Frequency' > {output.virgen};
        
        # Viral species summary
        grep "d_Viruses" {input.lca} | cut -f5 | awk -F ';' '{{ print$7 }}' | csvtk freq -H -n -r -T -t | sed '1i Family\tPrimary_AA_Order_Frequency' > {output.virspe};
        
        ## Summary table variables
        # Calculate number of input sequences (INPUT) from seqtable
        INPUT=$(grep -c ">" {input.seqs})
        
        # Calculate the number of LCA calls with viral lineage
        LCA_ALL=$(wc -l < {input.lca})
        pLCA_ALL=$(echo "scale=4 ; 100*($LCA_ALL / $INPUT)" | bc)
        LCA_VIRAL=$(grep -c "d_Viruses" {input.lca})
        pLCA_VIRAL=$(echo "scale=4 ; 100*($LCA_VIRAL / $INPUT)" | bc)
        LCA_NONVIRAL=$(grep -cv "d_Viruses" {input.lca})
        pLCA_NONVIRAL=$(echo "scale=4 ; 100*($LCA_NONVIRAL / $INPUT)" | bc)
        
        # Calculate the number of LCA calls that went to the virus root
        LCA_VIRAL_ROOT=$(awk '$2 == 10239' {input.lca} | wc -l)
        pLCA_VIRAL_ROOT=$(echo "scale=4 ; 100*($LCA_VIRAL_ROOT / $LCA_VIRAL)" | bc)
        
        # Calculate the number of tophits with complete viral lineages
        TH_ALL=$(wc -l < {input.tophit})
        pTH_ALL=$(echo "scale=4 ; 100*($TH_ALL / $INPUT)" | bc)
        TH_VIRAL=$(grep -c "d_Viruses" {input.tophit})
        pTH_VIRAL=$(echo "scale=4 ; 100*($TH_VIRAL / $INPUT)" | bc)
        TH_NONVIRAL=$(grep -cv "d_Viruses" {input.tophit})
        pTH_NONVIRAL=$(echo "scale=6 ; 100*($TH_NONVIRAL / $INPUT)" | bc)
        
        # Identify sequences with a tophit, but lost in the lca
        comm -13 <(grep "d_Viruses" {input.lca} | cut -f1 | sort) <(cut -f1 {input.tophit} | sort) > {output.tophit_only_ids};
        csvtk join -f1 {output.tophit_only_ids} {input.tophit} -H -t -T > {output.tophit_only_tsv};
            
        TH_ONLY=$(wc -l < {output.tophit_only_tsv})
        pTH_ONLY=$(echo "scale=10 ; 100*($TH_ONLY / $INPUT)" | bc)
        
        # Create list of 'broken" NCBI taxon IDs in the tophit table
        # Broken = lineage information unavailable for this taxon ID (deleted, moved, etc.)
        grep -v "d_Viruses" {input.tophit} | \
        cut -f20 | \
        csvtk freq -H -n -r -T -t | \
        sed '1i taxID\tFrequency' > {output.tophit_broken};
        
        # Calculate number of classified sequences
        CLASS=$(grep -c ">" {input.class_seqs})
        
        # Calculate number of unclassified sequences
        UNCLASS=$(grep -c ">" {input.unclass_seqs})
        
        # Calculate % of classified (pCLASS) and (pUNCLASS) unclassified sequences
        pCLASS=$(echo "scale=4 ; 100*($CLASS / $INPUT)" | bc)
        pUNCLASS=$(echo "scale=4 ; 100*($UNCLASS / $INPUT)" | bc)
        
        ## Create summary table
        rm -f {output.tmpsummary};
        touch {output.tmpsummary};
        
        printf "PARAMATER\tCOUNT\tPERCENT\tNOTES\n" >> {output.tmpsummary};
        printf "Input\t%d\t%s\t%s\n" "$INPUT" "NA" "Sequences input into primary mmseqs2 translated search (# of sequences in seqtable.fasta)" >> {output.tmpsummary};
        printf "LCA entries \t%d\t%.4f\t%s\n" "$LCA_ALL" "$pLCA_ALL" "# of LCA entries" >> {output.tmpsummary};
        printf "LCA Viral\t%d\t%.4f\t%s\n" "$LCA_VIRAL" "$pLCA_VIRAL" "Sequences with a viral LCA lineage" >> {output.tmpsummary};
        printf "LCA Nonviral\t%d\t%.4f\t%s\n" "$LCA_NONVIRAL" "$pLCA_NONVIRAL" "Sequences without a viral LCA lineage" >> {output.tmpsummary};
        printf "LCA Virus-Root LCA\t%d\t%.4f\t%s\n" "$LCA_VIRAL_ROOT" "$pLCA_VIRAL_ROOT" "Sequences classified to LCA virus-root taxonomy" >> {output.tmpsummary};
        printf "LCA Nonviral with a TopHit\t%d\t%.4f\t%s\n" "$TH_ONLY" "$pTH_ONLY" "Nonviral LCA but with a viral TopHit (listed in: MMSEQS_AA_PRIMARY_tophit_only.tsv)" >> {output.tmpsummary};
        printf "TopHit entries \t%d\t%.4f\t%s\n" "$TH_ALL" "$pTH_ALL" "# of TopHit entries" >> {output.tmpsummary};
        printf "TopHit Viral\t%d\t%.4f\t%s\n" "$TH_VIRAL" "$pTH_VIRAL" "Sequences with a tophit viral lineage" >> {output.tmpsummary};
        printf "TopHit with broken taxID\t%d\t%.4f\t%s\n" "$TH_NONVIRAL" "$pTH_NONVIRAL" "Sequences with broken lineage information in the reference database" >> {output.tmpsummary};
        printf "Classified\t%d\t%.4f\t%s\n" "$CLASS" "$pCLASS" "Sequences classified with a viral lineage to be checked in a secondary mmseqs search against a trans-kingdom database" >> {output.tmpsummary};
        printf "Unclassified\t%d\t%.4f\t%s\n" "$UNCLASS" "$pUNCLASS" "Sequences not classified with a viral lineage. Will be checked in a nt-vs-nt search" >> {output.tmpsummary};
 
        # Prettify
        csvtk pretty -t {output.tmpsummary} > {output.summary};
        
        
        """

rule SECONDARY_AA_taxonomy_assignment:
    """
    Check taxonomic assignments in MMSEQS_AA_PRIMARY_classified.fasta using mmseqs2
    
    - Reference database: UniRef50 + UNIVIRDB. All UniProtKB protein entries (all domains of life) clustered at 50% ID (https://www.uniprot.org/help/uniref) concatenated to UNIVIRDB
    
    - This is a nt-to-aa translated search, similar to blastX.
    
    - All sequences assiged to a viral lineage will be reserved as our final translated taxonomic calls
    
    - All sequences not classified as viral will be subjected to an untranslated search (nt-vs-nt) in later rules
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
        os.path.join(BENCH, "MMSEQS", "mmseqs_secondary_AA.txt")
    log:
        log = os.path.join(LOGS, "MMSEQS", "mmseqs_secondary_AA.log")
    resources:
        mem_mb=64000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
        -a --start-sens 1 --sens-steps 3 -s 7 \
        --tax-output-mode 2 --search-type 2 --lca-mode 2 --shuffle 0 \
        --lca-ranks "superkingdom,phylum,class,order,family,genus,species" \
        --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
        --tax-lineage 1 \
        -e {config[SECAAE]} &>> {log};
        
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        """
        
rule SECONDARY_AA_tophit_refactor:
    """
    
    Add/reformat tophit viral lineages with up-to-date* NCBI taxonomy
    
    - UniRef50 (https://www.uniprot.org/help/uniref) and other UniRef databases are a mixture of protein entries from UniProtKB (https://www.uniprot.org/help/uniprotkb) and UniParc (https://www.uniprot.org/help/uniparc)
    
    - The taxonomic lineage assignment algorithm of mmseqs2 will not provide full annotation of UniParc entries
    
    - We can updated those here using taxonkit lineage
    
    - This is also to ensure that the taxa are annotated based on the NCBI taxonomy which you can rapidly update in the databases/tax/taxonomy directory
    
    """
    input:
        db = TAX,
        tophit = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln")
    output:
        tophit_seqids = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.seq.ids")),
        tophit_taxids = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.tax.ids")),
        tophit_seq_taxids = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.seq_tax.ids")),
        tophit_lineage = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.lineage")),
        tophit_lineage_refomated = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.lineage.reformated")),
        tophit_tmp = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.lineage.vir.reformated")),
        tophit_tmp_updated = os.path.join(SECONDARY_AA_OUT, "tophit.tax_tmp_updated.tsv"),
        tophit_kingdom_freq = os.path.join(SECONDARY_AA_OUT, "tophit.kingdom.freq"),
        uncl_superkingdom_freq = os.path.join(SECONDARY_AA_OUT, "tophit.unclass_superkingdom.freq"),
        tophit_keyword_nonviral_ids = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.keyword_nonviral.ids")),
        tophit_keyword_nonviral_freq = os.path.join(SECONDARY_AA_OUT, "tophit.keyword_nonviral.freq"),
        tophit_keyword_nonviral_tmp = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.keyword_nonviral.tmp")),
        tophit_keyword_nonviral_list = os.path.join(SECONDARY_AA_OUT, "tophit.keyword_nonviral.list"),
        tophit_keyword_bac_freq = os.path.join(SECONDARY_AA_OUT, "tophit.keyword_bac.freq"),
        tophit_keyword_vir_freq_tr_ids = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.keyword_vir_tr.ids")),
        tophit_keyword_vir_freq_tr_tmp = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.keyword_vir_tr.tmp")),
        tophit_keyword_vir_freq_tr = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.keyword_vir_tr.freq")),
        tophit_keyword_vir_freq_sp_ids = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.keyword_vir_sp.ids")),
        tophit_keyword_vir_freq_sp_tmp = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.keyword_vir_sp.tmp")),
        tophit_keyword_vir_freq_sp = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.keyword_vir_sp.freq")),
        tophit_keyword_vir_freq = os.path.join(SECONDARY_AA_OUT, "tophit.keyword_vir.freq"),
        tophit_keyword_vir_list = os.path.join(SECONDARY_AA_OUT, "tophit.keyword_vir.list"),
        tophit_keyword_all_list = os.path.join(SECONDARY_AA_OUT, "tophit.keyword_all.list")
    conda:
        "../envs/seqkit.yaml"
    log:
        log = os.path.join(LOGS, "MMSEQS", "mmseqs_secondary_tophit_refactor.log")
    shell:
        """
        # Make a list of all sequence IDs represented in TopHit table
        cut -f1 {input.tophit} > {output.tophit_seqids};
        
        # Pull all UniProt TaxIDs from TopHit table
        cut -f20 {input.tophit} > {output.tophit_taxids};
        
        # Combine sequence IDs with taxIDs
        paste {output.tophit_seqids} {output.tophit_taxids} > {output.tophit_seq_taxids};
        
        # Add NCBI lineage information and collect all sequences with a viral lineage
        taxonkit lineage --data-dir {input.db} {output.tophit_seq_taxids} -i 2 | \
        cut --complement -f2 > {output.tophit_lineage};
        
        # Reformat TopHit viral lineage information
        taxonkit reformat --data-dir {input.db} {output.tophit_lineage} -i 2 \
        -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F --fill-miss-rank |
        cut --complement -f2 > {output.tophit_lineage_refomated};
        
        # Remove redundant lineage information
        cut --complement -f21,22 {input.tophit} > {output.tophit_tmp};
        
        # Create final TopHit viral sequence table with alignment information
        csvtk join -f1 {output.tophit_tmp} {output.tophit_lineage_refomated} -H -t -T > {output.tophit_tmp_updated} 2> {log};
        
        ## Summarize kingdom frequency information
            
        # Create table of kingdom frequencies
        cut -f21 {output.tophit_tmp_updated} | csvtk freq -H -n -r -T -t > {output.tophit_kingdom_freq};
        
        sed -i '1i kingdom\tfrequency' {output.tophit_kingdom_freq};
        
        # Create unclassified superkingdom frequency table
        grep 'unclassified unclassified entries superkingdom' {output.tophit_tmp_updated} | \
            cut -f27 | \
            csvtk freq -H -n -r -T -t > {output.uncl_superkingdom_freq};
            
        sed -i '1i description\tfrequency' {output.uncl_superkingdom_freq};
        
        ## Identify and tag provirus sequences based on UniRef protein name annotations
        # Separate all TopHit non-viral and viral sequences
        # Add column indicating there classification
        
        # Create frequency table of reference database key words for all nonviral TopHit entries
        awk -F '\t' '$21 != "Viruses"' {output.tophit_tmp_updated} | \
            cut -f19 | \
            sed 's/ n=[0-9]//g' | \
            awk -F "|" '{{ print $3 }}' | \
            awk -F "OS=" '{{ print$1 }}' > {output.tophit_keyword_nonviral_freq};
        
        sed -i '1i keyword\tfrequency' {output.tophit_keyword_nonviral_freq};
        
        # Create TopHit keyword list for all TopHit nonviral entries
        # Create id list
        awk -F '\t' '$21 != "Viruses"' {output.tophit_tmp_updated} | \
        cut -f1 > {output.tophit_keyword_nonviral_ids};
        
        awk -F '\t' '$21 != "Viruses"' {output.tophit_tmp_updated} | \
            cut -f19 | \
            sed 's/ n=[0-9]//g' | \
            awk -F "|" '{{ print $3 }}' | \
            awk -F "OS=" '{{ print$1 }}' > {output.tophit_keyword_nonviral_tmp};
            
        paste {output.tophit_keyword_nonviral_ids} {output.tophit_keyword_nonviral_tmp} > {output.tophit_keyword_nonviral_list};
        
        # Create frequency table of reference database key words for all TopHit bacteria
        awk -F '\t' '$21 == "Bacteria"' {output.tophit_tmp_updated} | \
            sed 's/ n=[0-9]//g' | \
            cut -f19 | \
            awk -F "|" '{{ print $3 }}' | \
            awk -F "OS=" '{{ print$1 }}' | \
            csvtk freq -H -n -r -T -t > {output.tophit_keyword_bac_freq};
            
        sed -i '1i keyword\tfrequency' {output.tophit_keyword_bac_freq};
            
        # Create frequency table of reference database key words for all TopHit viruses
        # First for the 'tr|' annotated viruses
        awk -F '\t' '$21 == "Viruses"' {output.tophit_tmp_updated} | \
            grep 'tr|' | \
            cut -f1 > {output.tophit_keyword_vir_freq_tr_ids};
            
        awk -F '\t' '$21 == "Viruses"' {output.tophit_tmp_updated} | \
            sed 's/ n=[0-9]//g' | \
            grep 'tr|' | \
            cut -f19 | \
            awk -F "OS=" '{{ print$1 }}' | \
            sed 's/ /:/' | \
            awk -F ":" '{{ print$2 }}' > {output.tophit_keyword_vir_freq_tr_tmp};
            
        paste {output.tophit_keyword_vir_freq_tr_ids} {output.tophit_keyword_vir_freq_tr_tmp} > {output.tophit_keyword_vir_freq_tr};
        
        # Next for the 'sp|' annotated viruses
        awk -F '\t' '$21 == "Viruses"' {output.tophit_tmp_updated} | \
            grep 'sp|' | \
            cut -f1 > {output.tophit_keyword_vir_freq_sp_ids};
        
        awk -F '\t' '$21 == "Viruses"' {output.tophit_tmp_updated} | \
            sed 's/ n=[0-9]//g' | \
            grep 'sp|' | \
            cut -f19 | \
            awk -F "|" '{{ print$3 }}' | \
            awk -F "OS=" '{{ print$1 }}' > {output.tophit_keyword_vir_freq_sp_tmp};
        
        paste {output.tophit_keyword_vir_freq_sp_ids} {output.tophit_keyword_vir_freq_sp_tmp} > {output.tophit_keyword_vir_freq_sp};
            
        cat {output.tophit_keyword_vir_freq_tr} {output.tophit_keyword_vir_freq_sp} > {output.tophit_keyword_vir_list};
            
        # Summarize viral keyword frequency table
        cat {output.tophit_keyword_vir_freq_tr} {output.tophit_keyword_vir_freq_sp} | \
        csvtk freq -H -n -r -T -t > {output.tophit_keyword_vir_freq};
            
        sed -i '1i keyword\tfrequency' {output.tophit_keyword_vir_freq};
        
        # Create all entry keyword list
        cat {output.tophit_keyword_vir_list} {output.tophit_keyword_nonviral_list} > {output.tophit_keyword_all_list};
        
        """
    
rule SECONDARY_AA_LCA_virus_root_refactor:
    """
    Update LCA virus root taxonomy to tophit taxonomy.
    
    - Some LCA virus taxonomies will default to: Viruses;uc_Viruses;uc_Viruses;uc_Viruses;uc_Viruses;uc_Viruses;uc_Viruses
    
    - We call this issue 'virus root taxonomy'
    
    - These sequences will be defaulted to their tophit taxonomy and marked accordingly
    """
    input:
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv"),
        tophit_updated = os.path.join(SECONDARY_AA_OUT, "tophit.tax_tmp_updated.tsv")
    output:
        lca_virus_root_seqids = temporary(os.path.join(SECONDARY_AA_OUT, "lca_virus_root.seq.ids")),
        lca_virus_root_vir_tsv = os.path.join(SECONDARY_AA_OUT, "lca_virus_root.vir.tsv")
    conda:
        "../envs/seqkit.yaml"
    log:
        log = os.path.join(LOGS, "MMSEQS", "mmseqs_secondary_lca_virus_root_refactor.log")
    shell:
        """
        # Isolate sequences pushed to virus root by LCA (Viruses;uc_Viruses;uc_Viruses;uc_Viruses;uc_Viruses;uc_Viruses;uc_Viruses)
        awk '$2 == 10239' {input.lca} | cut -f1 > {output.lca_virus_root_seqids};
        
        # Join LCA-root sequnece IDs to updated virus tophit table
        csvtk join -f1 {output.lca_virus_root_seqids} {input.tophit_updated} -H -t -T > {output.lca_virus_root_vir_tsv};
        
        """

rule SECONDARY_AA_LCA_virus_unclassified_refactor:
    """
    Update LCA unclassified taxonomy to tophit taxonomy.
    
    - There are a variety of reasons that a sequence may be assigned a tophit viral taxonomic lineage and not an LCA taxonomy (e.g. the LCA of a sequence assigned to a bacteria and a virus is 'root')
    
    - These sequences will be defaulted to their tophit taxonomy and marked accordingly.
    """
    input:
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv"),
        tophit_updated = os.path.join(SECONDARY_AA_OUT, "tophit.tax_tmp_updated.tsv")
    output:
        lca_unclass_seqids_0 = temporary(os.path.join(SECONDARY_AA_OUT, "lca_unclass.seq.ids_0")),
        lca_unclass_seqids_1 = temporary(os.path.join(SECONDARY_AA_OUT, "lca_unclass.seq.ids_1")),
        lca_unclass_seqids = temporary(os.path.join(SECONDARY_AA_OUT, "lca_unclass.seq.ids")),
        lca_unclass_vir_tsv = os.path.join(SECONDARY_AA_OUT, "lca_unclass.vir.tsv")
    conda:
        "../envs/seqkit.yaml"
    log:
        log = os.path.join(LOGS, "MMSEQS", "mmseqs_secondary_lca_unclassified_refactor.log")
    shell:
        """
        # Isolate unclassified sequences from LCA
        awk '$2 == 0' {input.lca} | cut -f1 > {output.lca_unclass_seqids_0};
        awk '$2 == 1' {input.lca} | cut -f1 > {output.lca_unclass_seqids_1};
        
        cat {output.lca_unclass_seqids_0} {output.lca_unclass_seqids_1} | cut -f1 > {output.lca_unclass_seqids};
        
        # Join unclassified LCA sequence IDs to virus tophit viral table
        # Save only the potential viral hits
        csvtk join -f1 {output.lca_unclass_seqids} {input.tophit_updated} -H -t -T > {output.lca_unclass_vir_tsv};
          
        """

rule SECONDARY_AA_refactor_finalize:
    """
    
    Remove sequences to be refactored from LCA table and recombine with updated taxonomies.
    
    """
    input:
        db = TAX,
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv"),
        tophit_updated = os.path.join(SECONDARY_AA_OUT, "tophit.tax_tmp_updated.tsv"),
        tophit_keyword_all_list = os.path.join(SECONDARY_AA_OUT, "tophit.keyword_all.list"),
        lca_virus_root_vir_tsv = os.path.join(SECONDARY_AA_OUT, "lca_virus_root.vir.tsv"),
        lca_unclass_vir_tsv = os.path.join(SECONDARY_AA_OUT, "lca_unclass.vir.tsv")
    output:
        lca_filt = temporary(os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca_filt.tmp")),
        lca_filt_lineage = temporary(os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca_filt.lineage")),
        lca_filt_linegage_reformated = temporary(os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca_filt.reformated")),
        tophit_updated_filt = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.tax_tmp_updated.filt")),
        tophit_lca_tmp = temporary(os.path.join(SECONDARY_AA_OUT, "tophit_lca.tmp")),
        tophit_lca_key = temporary(os.path.join(SECONDARY_AA_OUT, "tophit_lca.key")),
        tophit_kingdom = temporary(os.path.join(SECONDARY_AA_OUT, "tophit.kingdom")),
        tophit_lca_tsv = temporary(os.path.join(SECONDARY_AA_OUT, "tophit_lca.tsv")),
        lca_root_key = temporary(os.path.join(SECONDARY_AA_OUT, "lca_root.key")),
        lca_root_king = temporary(os.path.join(SECONDARY_AA_OUT, "lca_root.kingdom")),
        lca_unclass_key = temporary(os.path.join(SECONDARY_AA_OUT, "lca_unclass.key")),
        lca_unclass_king = temporary(os.path.join(SECONDARY_AA_OUT, "lca_unclass.kingdom")),
        translated_tmp = temporary(os.path.join(SECONDARY_AA_OUT, "translated_final.tmp")),
        translated_final = os.path.join(SECONDARY_AA_OUT, "translated_final.tsv")
    conda:
        "../envs/seqkit.yaml"
    log:
        log = os.path.join(LOGS, "MMSEQS", "mmseqs_secondary_lca_refactor_final.log")
    shell:
        """
        ## Create new base LCA table. Which will be recombined with the updated taxonomy tables from the previous refactoring rules
        
        # Remove virus root and unclassified entries from virus LCA table
        awk '$2 != 0' {input.lca} | awk '$2 != 10239' | awk '$2 != 1' > {output.lca_filt};
        
        # Pull all sequence and tax ids from filtered LCA table and reformat to updated NCBI taxonomy
        cut -f1,2 {output.lca_filt} | \
            taxonkit lineage --data-dir {input.db} -i 2 > {output.lca_filt_lineage};
            
        cut --complement -f2 {output.lca_filt_lineage} | taxonkit reformat --data-dir {input.db} -i 2 \
        -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F --fill-miss-rank | \
        cut --complement -f2 > {output.lca_filt_linegage_reformated};
        
        # Combine TopHit alignment information with updated LCA lineage information
        cut --complement -f21- {input.tophit_updated} > {output.tophit_updated_filt}; 
        
        csvtk join -f1 {output.tophit_updated_filt} {output.lca_filt_linegage_reformated} -T -t -H > {output.tophit_lca_tmp};
        
        # Add keyword to tophit_lca.tsv
        csvtk join -f1 {output.tophit_lca_tmp} {input.tophit_keyword_all_list} -T -t -H > {output.tophit_lca_key};
        
        # Add tophit kingdom category to tophit_lca.tsv
        cut -f1,21 {input.tophit_updated} | awk -F"\t" '{{ print$0"\t""tophit_"$2 }}' | \
            cut -f1,3 > {output.tophit_kingdom};
        
        csvtk join -f1 {output.tophit_lca_key} {output.tophit_kingdom} -T -t -H > {output.tophit_lca_tsv};
        
        sed -i 's/$/\\tlca_aa_classified/' {output.tophit_lca_tsv};
        
        ## Add keyword & tophit kingdom category to update lca files
        # First for LCA virus-root
        # Add keyword
        csvtk join -f1 {input.lca_virus_root_vir_tsv} {input.tophit_keyword_all_list} -T -t -H > {output.lca_root_key};
        
        # Add TopHit kingdom
        csvtk join -f1 {output.lca_root_key} {output.tophit_kingdom} -T -t -H > {output.lca_root_king};
        
        sed -i 's/$/\\tlca_aa_virus_root_updated/' {output.lca_root_king};
        
        # Second for lca unclassified
        # Add keyword
        csvtk join -f1 {input.lca_unclass_vir_tsv} {input.tophit_keyword_all_list} -T -t -H > {output.lca_unclass_key};
        
        csvtk join -f1 {output.lca_unclass_key} {output.tophit_kingdom} -T -t -H > {output.lca_unclass_king};
        
        sed -i 's/$/\\tlca_aa_unclassified_updated/' {output.lca_unclass_king};
        
        # Combine into final AA table
        cat {output.tophit_lca_tsv} {output.lca_root_king} {output.lca_unclass_king} > {output.translated_tmp};
        
        sed -i 's/$/\\ttranslated/' {output.translated_tmp};
        
        sort -k1 -n {output.translated_tmp} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatches\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tquery_header\ttarget_header\ttaxID\tKingdom\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies\tprotein_name\ttophit_kingdom\tclassification\tquery_type' > {output.translated_final};
        """

rule SECONDARY_AA_parsing:
    """
    
    Parse out all sequences that remain unclassified following the Secondary AA search and refactoring.
        
    """
    input:
        lca = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.tsv"),
        translated_final = os.path.join(SECONDARY_AA_OUT, "translated_final.tsv"),
        seqs = os.path.join(RESULTS, "seqtable.fasta")
    output:
        allseq_ids = temporary(os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_PRIMARY_lca.ids")),
        class_ids = temporary(os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca_classified.ids")),
        unclass_ids = temporary(os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca_unclassified.ids")),
        unclass_seqs = os.path.join(SECONDARY_AA_OUT, "translated_unclassified.fasta")
    conda:
        "../envs/samtools.yaml"
    log:
        "LOGS/mmseqs/mmseqs_SECONDARY_aa_parsing.log"
    shell:
        """
        # Extract full entry ID list (all sequences)
        # Note: the lca output table has one entry per input sequence
        # This is a list of all sequence id's that went into the primary search
        cut -f1 {input.lca} > {output.allseq_ids};
        
        # Extract sequences IDs for everything that was classified as virus following the secondary AA search and refactoring
        tail -n+2 {input.translated_final} | cut -f1 > {output.class_ids};
        
        # Compare all classified AA search sequences to input classified sequences to create a list of unclassified sequences
        comm -23 <(sort {output.allseq_ids}) <(sort {output.class_ids}) | sort -n > {output.unclass_ids};
        
        # Extract unclassified sequences from seocndary translated search
        xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
             
        """
        
rule PRIMARY_NT_taxonomic_assignment:
    """
    
    For sequences not assigned a translated (aa-to-nt) taxonomy attempt to assing an untranslated (nt-to-nt) taxonomy to a viral nucleic acid database.
    
    """
    input:
        seqs = os.path.join(SECONDARY_AA_OUT, "translated_unclassified.fasta"),
        db = os.path.join(NCBIVIRDB, "sequenceDB")
    output:
        queryDB = os.path.join(PRIMARY_NT_OUT, "queryDB"),
        result = os.path.join(PRIMARY_NT_OUT, "results", "result.index")
    params:
        respath = os.path.join(PRIMARY_NT_OUT, "results", "result"),
        tmppath = os.path.join(PRIMARY_NT_OUT, "mmseqs_aa_tmp")
    benchmark:
            os.path.join(BENCH, "MMSEQS", "mmseqs_primary_NT.txt")
    log:
        log = os.path.join(LOGS, "MMSEQS", "mmseqs_primary_NT.log")
    resources:
        mem_mb=64000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        # Create query database
        mmseqs createdb {input.seqs} {output.queryDB} --dbtype 2;
        
        # mmseqs search
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
        --start-sens 2 -s 7 --sens-steps 3 \
        --search-type 3 \
        -e {config[PRIMNTE]} &>> {log};
    
        """
        
rule PRIMARY_NT_summary:
    """
    
    Summarize primary untranslated (nt-to-nt) search results.

    """
    input:
        queryDB = os.path.join(PRIMARY_NT_OUT, "queryDB"),
        db = os.path.join(NCBIVIRDB, "sequenceDB"),
        taxdb = TAX
    output:
        result = os.path.join(PRIMARY_NT_OUT, "results", "firsthit.index"),
        align = os.path.join(PRIMARY_NT_OUT, "results", "result.m8"),
        lineage = temporary(os.path.join(PRIMARY_NT_OUT, "primary_nt.lineage")),
        reformated = os.path.join(PRIMARY_NT_OUT, "PRIMARY_nt.tsv"),
        phyl_sum = os.path.join(PRIMARY_NT_OUT, "PRIMARY_nt_phylum_summary.tsv"),
        class_sum = os.path.join(PRIMARY_NT_OUT, "PRIMARY_nt_class_summary.tsv"),
        ord_sum = os.path.join(PRIMARY_NT_OUT, "PRIMARY_nt_order_summary.tsv"),
        fam_sum = os.path.join(PRIMARY_NT_OUT, "PRIMARY_nt_family_summary.tsv"),
        gen_sum = os.path.join(PRIMARY_NT_OUT, "PRIMARY_nt_genus_summary.tsv"),
        spe_sum = os.path.join(PRIMARY_NT_OUT, "PRIMARY_nt_species_summary.tsv")
    params:
        inputpath = os.path.join(PRIMARY_NT_OUT, "results", "result"),
        respath = os.path.join(PRIMARY_NT_OUT, "results", "firsthit")
    conda:
        "../envs/mmseqs2.yaml"
    log:
        "LOGS/mmseqs/mmseqs_PRIMARY_nt_summary.log"
    shell:
        """
        # Filter TopHit results
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
        
        ## Generate summary tables
        # Phylum
        cut -f3 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Phylum\tPrimary_NT_Phylum_Frequency'> {output.phyl_sum};
            
        # Class
        cut -f4 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Class\tPrimary_NT_Class_Frequency'> {output.class_sum};
        
        # Order
        cut -f5 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Order\tPrimary_NT_Order_Frequency'> {output.ord_sum};
            
        # Family
        cut -f6 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Family\tPrimary_NT_Family_Frequency'> {output.fam_sum};
        
        # Genus
        cut -f7 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Genus\tPrimary_NT_Genus_Frequency'> {output.gen_sum};
            
        # Species
        cut -f3 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Species\tPrimary_NT_Species_Frequency'> {output.spe_sum};
        """
        
rule PRIMARY_NT_parsing:
    """
    
    Parse and summarize primary untranlsated (nt-to-nt) results.
    
    """
    input:
        seqs = os.path.join(SECONDARY_AA_OUT, "translated_unclassified.fasta"),
        align = os.path.join(PRIMARY_NT_OUT, "results", "result.m8")
    output:
        allseq_ids = temporary(os.path.join(PRIMARY_NT_OUT, "all_input.ids")),
        class_ids = temporary(os.path.join(PRIMARY_NT_OUT, "primary_classified.ids")),
        class_seqs = os.path.join(PRIMARY_NT_OUT, "classified_seqs.fasta"),
        unclass_ids = temporary(os.path.join(PRIMARY_NT_OUT, "unclassified.ids")),
        unclass_seqs = os.path.join(PRIMARY_NT_OUT, "unclassified_seqs.fasta")
    conda:
        "../envs/samtools.yaml"
    log:
        "LOGS/mmseqs/mmseqs_PRIMARY_nt_parsing.log"
    resources:
        mem_mb=64000,
        cpus=64
    shell:
        """
        # Extract full entry ID list (all sequences)
        seqkit seq -n {input.seqs} > {output.allseq_ids};
        
        # Extract primary nt hit ids
        cut -f1 {input.align} > {output.class_ids};
        
        # Extract unclassified ids
        comm -23 <(sort {output.allseq_ids}) <(sort {output.class_ids}) > {output.unclass_ids};
        
        # Extract classified sequences from primary untranslated search
        xargs samtools faidx {input.seqs} -n 5000 < {output.class_ids} > {output.class_seqs};
        
        # Extract unclassified sequences from seocndary translated search
        xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
        
        """

rule SECONDARY_NT_taxonomic_assignment:
    """
    
    Confrim primary untranslated (nt-to-nt) taxonomic assignments by checking against a more comprehensive, transkingdom nucleotide database.
    
    """
    input:
        seqs = os.path.join(PRIMARY_NT_OUT, "classified_seqs.fasta"),
        db = os.path.join(POLYMICRODB, "sequenceDB")
    output:
        queryDB = os.path.join(SECONDARY_NT_OUT, "queryDB"),
        result = os.path.join(SECONDARY_NT_OUT, "results", "result.index")
    params:
        respath = os.path.join(SECONDARY_NT_OUT, "results", "result"),
        tmppath = os.path.join(SECONDARY_NT_OUT, "mmseqs_aa_tmp")
    benchmark:
            os.path.join(BENCH, "MMSEQS", "mmseqs_secondary_NT.txt")
    log:
        log = os.path.join(LOGS, "MMSEQS", "mmseqs_secondary_NT.log")
    resources:
        mem_mb=64000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        # Create query database
        mmseqs createdb {input.seqs} {output.queryDB} --dbtype 2;
        
        # mmseqs search
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
        -a 1 --search-type 3 \
        --start-sens 2 -s 7 --sens-steps 3 \
        -e {config[SECNTE]} &>> {log};
    
        """

rule SECONDARY_NT_summary:
    """
    
    Summarize secondary untranslated (nt-to-nt) confirmatory search results.
    
    """
    input:
        queryDB = os.path.join(SECONDARY_NT_OUT, "queryDB"),
        db = os.path.join(POLYMICRODB, "sequenceDB"),
        taxdb = TAX
    output:
        result = os.path.join(SECONDARY_NT_OUT, "results", "tophit.index"),
        align = os.path.join(SECONDARY_NT_OUT, "results", "tophit.m8"),
        lineage = temporary(os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt.lineage")),
        reformated = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt.tsv"),
        king_sum = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt_kingdom_summary.tsv"),
        phyl_sum = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt_phylum_summary.tsv"),
        class_sum = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt_class_summary.tsv"),
        ord_sum = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt_order_summary.tsv"),
        fam_sum = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt_family_summary.tsv"),
        gen_sum = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt_genus_summary.tsv"),
        spe_sum = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt_species_summary.tsv")
    params:
        inputpath = os.path.join(SECONDARY_NT_OUT, "results", "result"),
        respath = os.path.join(SECONDARY_NT_OUT, "results", "tophit")
    resources:
        mem_mb=64000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    log:
        "LOGS/mmseqs/mmseqs_SECONDARY_nt_summary.log"
    shell:
        """
        # Filter TopHit results
        mmseqs filterdb {params.inputpath} {params.respath} --extract-lines 1;
        
        # Convert to alignments
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} \
        --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader";
        
        # Assign taxonomy
        cut -f1,2 {output.align} | \
            awk -F '|' '{{ print$1$2 }}' | \
            sed 's/tid//g' | \
            taxonkit lineage --data-dir {input.taxdb} -i 2 > {output.lineage};
            
        # Reformat TopHit viral lineage information
        taxonkit reformat --data-dir {input.taxdb} {output.lineage} -i 3 \
        -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F --fill-miss-rank |
        cut --complement -f3 > {output.reformated};
        
        ## Generate summary tables
        # Kingdom
        cut -f3 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Kingdom\tSECONDARY_NT_Kingdom_Frequency'> {output.king_sum};
        
        # Phylum
        cut -f4 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Phylum\tSECONDARY_NT_Phylum_Frequency'> {output.phyl_sum};
            
        # Class
        cut -f5 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Class\tSECONDARY_NT_Class_Frequency'> {output.class_sum};
        
        # Order
        cut -f6 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Order\tSECONDARY_NT_Order_Frequency'> {output.ord_sum};
            
        # Family
        cut -f7 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Family\tSECONDARY_NT_Family_Frequency'> {output.fam_sum};
            
        # Genus
        cut -f8 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Genus\tSECONDARY_NT_Genus_Frequency'> {output.gen_sum};
            
        # Species
        cut -f9 {output.reformated} | \
            csvtk freq -H -n -r -T -t | \
            sed '1i Species\tSECONDARY_NT_Species_Frequency'> {output.spe_sum};
        """
        
rule SECONDARY_NT_calculate_LCA:
    """
    
    Calculate LCA from secondary untranslated (nt-to-nt) search results.
    
    """
    input:
        queryDB = os.path.join(SECONDARY_NT_OUT, "queryDB"),
        db = os.path.join(POLYMICRODB, "sequenceDB"),
        tophit = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt.tsv"),
        align = os.path.join(SECONDARY_NT_OUT, "results", "tophit.m8"),
        taxdb = TAX
    output:
        align = os.path.join(SECONDARY_NT_OUT, "results", "all.m8"),
        lca_lineage = temporary(os.path.join(SECONDARY_NT_OUT, "results", "lca.lineage")),
        reformated = os.path.join(SECONDARY_NT_OUT, "results", "secondary_nt_lca.tsv"),
        unclass_ids = os.path.join(SECONDARY_NT_OUT, "results", "lca_unclassified.ids"),
        tophit_viral = os.path.join(SECONDARY_NT_OUT, "results", "secondary_tophit_viral.tsv"),
        unclass_annot = os.path.join(SECONDARY_NT_OUT, "results", "lca_unclassified_annotated.tsv"),
        unclass_removed = os.path.join(SECONDARY_NT_OUT, "results", "secondary_nt_lca_unclassified_removed.tsv"),
        remaining_unclass_ids = os.path.join(SECONDARY_NT_OUT, "results", "lca_remaining_unclassified_removed.tsv"),
        remaining_lineage = os.path.join(SECONDARY_NT_OUT, "results", "lca_remaining_unclassified_lineage.tsv"),
        lca_final_lin = os.path.join(SECONDARY_NT_OUT, "results", "lca_final_lineages.tsv"),
        lca_final_aln = os.path.join(SECONDARY_NT_OUT, "results", "secondary_lca_final.tsv")
    params:
        respath = os.path.join(SECONDARY_NT_OUT, "results", "result")
    resources:
        mem_mb=64000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    log:
        "LOGS/mmseqs/mmseqs_SECONDARY_nt_lca.log"
    shell:
        """
        # Convert to alignments
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} \
        --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader";
        
        # Calculate LCA for each sequence
        cut -f1,2 {output.align} | \
            sed 's/tid|//' | \
            awk -F '|' '{{ print$1 }}' | \
            sort -k1,2 | \
            uniq | \
            csvtk fold -f 1 -v 2 -H -t -s ';' | \
            taxonkit lca -i 2 -s ';' --data-dir {input.taxdb} | \
            cut -f1,3 | \
            taxonkit lineage -i 2 --data-dir {input.taxdb} > {output.lca_lineage}
            
        # Reformat lineages
        awk -F '\t' '$2 != 0' {output.lca_lineage} | \
        taxonkit reformat --data-dir {input.taxdb} -i 3 \
        -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F --fill-miss-rank |
        cut --complement -f2,3 > {output.reformated};
        
        # Create list of unclassified ids
        grep "unclassified root" {output.reformated} | cut -f1 > {output.unclass_ids};
        
        # Create table of tophit virus sequences
        awk -F '\t' '$3 == "Viruses"' {input.tophit} | cut --complement -f2 > {output.tophit_viral};
        
        # Annotate LCA unclassifieds with tophit viral lineages
        csvtk join -f1 {output.unclass_ids} {output.tophit_viral} -t -T -H > {output.unclass_annot};
        
        # Remove unclassified sequences from lca table
        grep -v "unclassified root" {output.reformated} > {output.unclass_removed};
        
        # Create list of ids that remain unclassified
        comm -23 <(sort {output.unclass_ids}) <(cut -f1 {output.unclass_annot} | sort) > {output.remaining_unclass_ids};
        
        # Extract remaining unclassified lineages
        csvtk join -f1 {output.remaining_unclass_ids} {output.reformated} -H -t -T > {output.remaining_lineage};
        
        # Combine LCA + LCA-tophit assigned + remaining lineages
        cat {output.unclass_removed} {output.unclass_annot} {output.remaining_lineage} > {output.lca_final_lin};
        
        # Combint alignment information with final lineages
        csvtk join -f1 {input.align} {output.lca_final_lin} -H -T -t > {output.lca_final_aln};
        
        """
