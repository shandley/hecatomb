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

rule PRIMARY_AA_taxonomy_assignment:
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
        
rule PRIMARY_AA_taxonomy_tabulation:
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
        
        # Calculate % of classified (pCLASS) and (pUNCLASS) unclassified sequences
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

rule SECONDARY_AA_taxonomy_assignment:
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

rule SECONDARY_AA_LCA_refactor:
    """
    
    TBD
    
    """
    input:
        db = TAX,
        tophit = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln"),
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv")
    output:
        tophit_seqids = os.path.join(SECONDARY_AA_OUT, "tophit.seq.ids"),
        tophit_taxids = os.path.join(SECONDARY_AA_OUT, "tophit.tax.ids"),
        tophit_seq_taxids = os.path.join(SECONDARY_AA_OUT, "tophit.seq_tax.ids"),
        tophit_lineage_vir = os.path.join(SECONDARY_AA_OUT, "tophit.lineage.vir"),
        tophit_lineage_vir_refomated = os.path.join(SECONDARY_AA_OUT, "tophit.lineage.vir.reformated"),
        tophit_vir_tsv = os.path.join(SECONDARY_AA_OUT, "tophit.vir.tsv"),
        lca_virus_root_seqids = os.path.join(SECONDARY_AA_OUT, "lca_virus_root.seq.ids"),
        lca_virus_root_vir_tsv = os.path.join(SECONDARY_AA_OUT, "lca_virus_root.vir.tsv"),
        lca_virus_root_nonvir_seqids = os.path.join(SECONDARY_AA_OUT, "lca_virus_root.nonvir.seq.ids"),
        lca_unclass_seqids = os.path.join(SECONDARY_AA_OUT, "lca_unclass.seq.ids"),
        lca_unclass_vir_tsv = os.path.join(SECONDARY_AA_OUT, "lca_unclass.vir.tsv"),
        lca_unclass_nonvir_seqids = os.path.join(SECONDARY_AA_OUT, "lca_unclass.nonvir.seq.ids"),
        mmseqs_novir_root_tsv = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca_no_virus_root.tsv"),
        mmseqs_lca_filt_tsv = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.filtered.tsv"),
        lca_filt_seqids = os.path.join(SECONDARY_AA_OUT, "lca.filtered.seq.ids"),
        lca_filt_taxids = os.path.join(SECONDARY_AA_OUT, "lca.filtered.tax.ids"),
        lca_filt_seq_taxids = os.path.join(SECONDARY_AA_OUT, "lca.filtered.seq_tax.ids"),
        lca_filt_lineage = os.path.join(SECONDARY_AA_OUT, "lca.filtered.lineage"),
        lca_filt_linegage_reformated = os.path.join(SECONDARY_AA_OUT, "lca.filtered.lineage.reformated"),
        lca_aa_tsv = os.path.join(SECONDARY_AA_OUT, "lca_aa.tsv"),
        lca_concat_tsv = os.path.join(SECONDARY_AA_OUT, "lca_concat_aa.tsv"),
        lca_final_tsv = os.path.join(SECONDARY_AA_OUT, "lca_final_aa.tsv"),
        aa_unannotated_seqids = os.path.join(SECONDARY_AA_OUT, "aa_unannotated.ids"),
        viroDB_headers = os.path.join(SECONDARY_AA_OUT, "viroDB.headers"),
        viroDB_ids = os.path.join(SECONDARY_AA_OUT, "viroDB.ids"),
        uniref_headers = os.path.join(SECONDARY_AA_OUT, "uniref.headers"),
        uniref_ids = os.path.join(SECONDARY_AA_OUT, "uniref.ids"),
        viroDB_proteins = os.path.join(SECONDARY_AA_OUT, "viroDB.proteins"),
        viroDB_protein_tsv = os.path.join(SECONDARY_AA_OUT, "viroDB.proteins.tsv"),
        uniref_proteins = os.path.join(SECONDARY_AA_OUT, "uniref.proteins"),
        uniref_proteins_tsv = os.path.join(SECONDARY_AA_OUT, "uniref.proteins.tsv"),
        protein_annotations = os.path.join(SECONDARY_AA_OUT, "protein.annotations.tsv"),
        lca_aa_final_protein_annotated = os.path.join(SECONDARY_AA_OUT, "lca_aa_final_protein_annotated.tsv")
    conda:
        "../envs/seqkit.yaml"
    log:
        "LOGS/mmseqs/mmseqs_lca_refactor.log"
    shell:
        """
        ##### Create updated tophit table of viral sequenes
        # Make a list of all sequence IDs represented in tophit table
        cut -f1 {input.tophit} > {output.tophit_seqids};
        
        # Pull all UniPort TaxIDs from tophit table
        cut -f19 {input.tophit} | awk -F 'TaxID=' '{{ print$2 }}' | awk -F ' ' '{{ print$1 }}' > {output.tophit_taxids};
        
        # Combine sequence IDs with taxIDs
        paste {output.tophit_seqids} {output.tophit_taxids} > {output.tophit_seq_taxids};
        
        # Add NCBI lineage information and collect all sequences with a viral lineage
        taxonkit lineage --data-dir {input.db} {output.tophit_seq_taxids} -i 2 | \
        grep 'Viruses;' | \
        cut --complement -f2 > {output.tophit_lineage_vir};
        
        # Reformat tophit viral lineage information
        taxonkit reformat --data-dir {input.db} {output.tophit_lineage_vir} -i 2 \
        -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F --fill-miss-rank |
        cut --complement -f2 > {output.tophit_lineage_vir_refomated};
        
        # Create final tophit viral sequence table with alignment information
        csvtk join -f1 {output.tophit_lineage_vir_refomated} {input.tophit} -H -t -T --left-join > {output.tophit_vir_tsv};
        
        ##### Update LCA virus root taxonomy to tophit taxonomy
        # Isolate sequences pushed to virus root by LCA (Viruses;uc_Viruses;uc_Viruses;uc_Viruses;uc_Viruses;uc_Viruses;uc_Viruses)
        awk '$2 == 10239' {input.lca} | cut -f1 > {output.lca_virus_root_seqids}
        
        # Join LCA-root sequnece IDs to updated tophit viral table (tophit.vir.tsv)
        csvtk join -f1 {output.lca_virus_root_seqids} {output.tophit_vir_tsv} -H -t -T > {output.lca_virus_root_vir_tsv};
        
        # Add LCA category to unpdated LCA-root table
        sed -i 's/$/\tlca_aa_virus_root_updated/' {output.lca_virus_root_vir_tsv};
        
        # Extract Seq IDs for sequences not updated for downsteram nt-to-nt search
        comm -23 <(sort {output.lca_virus_root_seqids}) <(cut -f1 {output.lca_virus_root_vir_tsv} | sort) > {output.lca_virus_root_nonvir_seqids};
        
        ##### Update LCA unclassified taxonomy to tophit taxonomy
        # Isolate unclassified sequences from LCA
        awk '$2 == 0' {input.lca} | cut -f1 > {output.lca_unclass_seqids};
        
        # Join unclassified LCA sequence IDs to updated tophit viral table (tophit.vir.tsv)
        csvtk join -f1 {output.lca_unclass_seqids} {output.tophit_vir_tsv} -H -t -T > {output.lca_unclass_vir_tsv};
        
        # Add LCA category to unpdated LCA-unclassified table
        sed -i 's/$/\tlca_aa_unclassified_updated/' {output.lca_unclass_vir_tsv};
        
        # Extract Seq IDs for sequences not updated for downsteram nt-to-nt search
        comm -23 <(sort {output.lca_unclass_seqids}) <(cut -f1 {output.lca_unclass_vir_tsv} | sort) > {output.lca_unclass_nonvir_seqids};   
        
        ##### Seperate sequences that will be updated/parsed from lca table. Parsed sequences will either be updated here or sent to the nt-to-nt search
        # Separate virus root taxonomies from LCA table
        awk '$2 != 10239' {input.lca} > {output.mmseqs_novir_root_tsv};
        
        # Separate LCA unlcassified sequences from LCA table
        awk '$2 != 0' {output.mmseqs_novir_root_tsv} > {output.mmseqs_lca_filt_tsv};
        
        ## This is our new base LCA table. We will not superced any tophit classifications with these LCA classificiations

        # Add LCA category to unpdated base LCA table
        sed -i 's/$/\tlca_aa/' {output.mmseqs_lca_filt_tsv};
        
        ##### Update taxonomic lineages and formats of LCA assignments
        # Collect all filtered LCA sequence IDs
        cut -f1 {output.mmseqs_lca_filt_tsv} > {output.lca_filt_seqids};
        
        # Collect all filtered LCA NCBI taxIDs
        cut -f2 {output.mmseqs_lca_filt_tsv} > {output.lca_filt_taxids};
        
        # Combine sequence IDs and taxIDs
        paste {output.lca_filt_seqids} {output.lca_filt_taxids} > {output.lca_filt_seq_taxids};
        
        # Add NCBI taxonomy
        taxonkit lineage --data-dir {input.db} {output.lca_filt_seq_taxids} -i 2 > {output.lca_filt_lineage};
        
        # Reformat lca lineages
        cut --complement -f2 {output.lca_filt_lineage} | taxonkit reformat --data-dir {input.db} -i 2 \
                -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F --fill-miss-rank | \
                cut --complement -f2 > {output.lca_filt_linegage_reformated};
                
        # Combine with alignment information
        csvtk join -f1 {output.lca_filt_linegage_reformated} {input.tophit} -H -t -T > {output.lca_aa_tsv};
        
        # Add taxonomic assignment category
        sed -i 's/$/\tlca_aa/' {output.lca_aa_tsv};
        
        # Concatenate lca_root, unclassified and lca_aa tsv tables
        cat {output.lca_aa_tsv} {output.lca_virus_root_vir_tsv} {output.lca_unclass_vir_tsv} > {output.lca_concat_tsv};
        
        # Add headers
        sort -k1 -n {output.lca_concat_tsv} | \
                sed '1i query\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\ttarget\tevalue\tpident\tfident\tnident\tmismatches\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage\tquery_type' > {output.lca_final_tsv};
        
        # Concatentate all sequence IDs to go to nt-to-nt search
        cat {output.lca_virus_root_nonvir_seqids} {output.lca_unclass_nonvir_seqids} > {output.aa_unannotated_seqids};
        
        ##### PARSE PROTEIN ANNOTATIONS
        # UniRef50 and the Virus Protein DB have different annotation schemes / headers so they need to be treated separately

        # Seperate virusDB headers
        tail -n+2 {output.lca_final_tsv} | cut -f26 | grep -v "UniRef50" > {output.viroDB_headers};
        tail -n+2 {output.lca_final_tsv} | grep -v "UniRef50" | cut -f1 > {output.viroDB_ids};
        
        # Seperate UniRef50 headers
        tail -n+2 {output.lca_final_tsv} | cut -f26 | grep "UniRef50" > {output.uniref_headers};
        tail -n+2 {output.lca_final_tsv} | grep "UniRef50" | cut -f1 > {output.uniref_ids};
        
        # Isolate protein annotations for virusDB
        awk -F '|' '{{print$3}}' {output.viroDB_headers} | awk -F "OS=" '{{print$1}}' | cut --complement -f1 -d " " > {output.viroDB_proteins};
        paste {output.viroDB_ids} {output.viroDB_proteins} > {output.viroDB_protein_tsv};
        
        # Isolate protein annotations for UniRef50
        awk -F 'n=' '{{print$1}}' {output.uniref_headers} | cut -f1 --complement -d " " > {output.uniref_proteins};
        paste {output.uniref_ids} {output.uniref_proteins} > {output.uniref_proteins_tsv};
        
        # Recombine
        cat {output.viroDB_protein_tsv} {output.uniref_proteins_tsv} | \
                sort -k1 -n | \
                sed '1i query\tprotein' > {output.protein_annotations};
                
        # Combine
        csvtk join -f1 {output.lca_final_tsv} {output.protein_annotations} -t -T > {output.lca_aa_final_protein_annotated};
        """
    
rule SECONDARY_AA_taxonomy_tabulation:
    """
    
    TBD
        
    """
    input:
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv"),
        alnsort = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln_sorted"),
        seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta")
    output:
        allseq_ids = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.ids"),
        class_ids = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca_classified.ids"),
        unclass_ids = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca_unclassified.ids"),
        class_seqs = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_classified.fasta"),
        unclass_seqs = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_unclassified.fasta"),
        summary = os.path.join(RESULTS, "SEARCH_SUMMARY", "MMSEQS_AA_SECONDARY_summary.txt"),
        virfam = os.path.join(RESULTS, "SEARCH_SUMMARY", "MMSEQS_AA_SECONDARY_virus_family_summary.txt")
    conda:
        "../envs/samtools.yaml"
    log:
        "LOGS/mmseqs/mmseqs_SECONDARY_aa_parsing.log"
    shell:
        """
        # Extract full entry ID list (all sequences)
        # Note: the lca output table has one entry per input sequence
        sort -k1 -n {input.lca} | cut -f1 > {output.allseq_ids};
        
        # Extract all sequence IDs assigned a viral LCA taxonomy
        grep "d_Viruses" {input.lca} | cut -f1 | sort -n > {output.class_ids};
        
        # Extract all sequence IDs not assigned a viral LCA taxonomy
        grep -v "d_Viruses" {input.lca} | cut -f1 | sort -n > {output.unclass_ids};
        
        # Extract classified sequences
        xargs samtools faidx {input.seqs} -n 5000 < {output.class_ids} > {output.class_seqs};
        
        # Extract unclassified sequences
        xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
        
        # Viral family summary
        grep "d_Viruses" {input.lca} | cut -f5 | awk -F ';' '{{ print$5 }}' | sort | uniq -c | sort -k1 -nr > {output.virfam};
        
        ## Summary table variables
        # Calculate number of input sequences (INPUT) from seqtable
        INPUT=$(grep -c ">" {input.seqs})
        
        # Calculate the number of LCA calls with viral lineage
        LCA_ALL=$(wc -l {input.lca})
        LCA_VIRAL=$(grep -c "d_Viruses" {input.lca})
        LCA_NONVIRAL=$(grep -cv "d_Viruses" {input.lca})
        pLCA_VIRAL=$(echo "scale=4 ; 100*($LCA_VIRAL / $INPUT)" | bc)
        pLCA_NONVIRAL=$(echo "scale=4 ; 100*($LCA_NONVIRAL / $INPUT)" | bc)
        
        # Create summary table
        # NOTE: Here docs are explicit and will be formated as shown. Do not add tabs in front of doc as you normally would in snakemake code! Format it exactly how you want to see it!
        cat << EOF > {output.summary}
PARAMATER    VALUE    PERCENT    NOTES
Input    $INPUT    NA    Number of sequences input into SECONDARY mmseqs amino acid search (taken from MMSEQS_AA_PRIMARY_classified.fasta)
Viral LCA    $LCA_VIRAL    $pLCA_VIRAL    Number and % of sequences with a viral lineage (== d_Viruses) after LCA calculation
Nonviral LCA    $LCA_NONVIRAL    $pLCA_NONVIRAL    Number and % of sequences without a viral lineage (!= d_Viruses) after LCA calculation        
        """
        
rule mmseqs_AA_diagnostics:
    """
    
    TBD
    
    """
    input:
        prim_class_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
        tophit_aln = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln"),
        secondary_class_seqs = os.path.join(SECONDARY_AA_OUT, "lca_aa_final_protein_annotated.tsv"),
        db = TAX,
        baltimore = os.path.join(TABLES, "2020_07_27_Viral_classification_table_ICTV2019.txt")
    output:
        prim_class_ids = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta.ids"),
        primary_AA_classified_tophit_aln = os.path.join(PRIMARY_AA_OUT, "primary_AA_classified.tophit.aln"),
        primary_AA_classified_tophit_aln_taxids = os.path.join(PRIMARY_AA_OUT, "primary_AA_classified.tophit.aln.taxIDs"),
        primary_AA_classified_tophit_aln_seq_taxids = os.path.join(PRIMARY_AA_OUT, "primary_AA_classified.tophit.aln.seq_taxIDs"),
        primary_AA_classified_lineage = os.path.join(PRIMARY_AA_OUT, "primary_AA_classified.lineage"),
        primary_AA_classified_lineage_reformated = os.path.join(PRIMARY_AA_OUT, "primary_AA_classified.lineage.reformated"),
        primary_AA_order_summary = os.path.join(PRIMARY_AA_OUT, "primary_AA_order.summary"),
        primary_AA_family_summary = os.path.join(PRIMARY_AA_OUT, "primary_AA_family.summary"),
        primary_AA_genus_summary = os.path.join(PRIMARY_AA_OUT, "primary_AA_genus.summary"),
        primary_AA_species_summary = os.path.join(PRIMARY_AA_OUT, "primary_AA_species.summary"),
        secondary_AA_order_summary = os.path.join(PRIMARY_AA_OUT, "secondary_AA_order.summary"),
        secondary_AA_family_summary = os.path.join(PRIMARY_AA_OUT, "secondary_AA_family.summary"),
        secondary_AA_genus_summary = os.path.join(PRIMARY_AA_OUT, "secondary_AA_genus.summary"),
        secondary_AA_species_summary = os.path.join(PRIMARY_AA_OUT, "secondary_AA_species.summary"),
        order_compare = os.path.join(PRIMARY_AA_OUT, "order.compare"),
        family_compare = os.path.join(PRIMARY_AA_OUT, "family.compare"),
        genus_compare = os.path.join(PRIMARY_AA_OUT, "genus.compare"),
        species_compare = os.path.join(PRIMARY_AA_OUT, "species.compare"),
        family_compare_baltimore = os.path.join(PRIMARY_AA_OUT, "family.compare.baltimore")
    benchmark:
        "BENCHMARKS/mmseqs_diagnostics.txt"
    log:
        "../LOGS/mmseqs/mmseqs_diagnostics.log"
    conda:
        "../envs/seqkit.yaml"
    shell:
        """
        # Pull all Primary AA search classified seqIDs
        seqkit seq -n {input.prim_class_seqs} > {output.prim_class_ids};
        
        # Join with alignment statistics
        csvtk join -f 1 {output.prim_class_ids} {input.tophit_aln} -H -t -T > {output.primary_AA_classified_tophit_aln};
        
        # Pull taxIDs
        cut -f20 {output.primary_AA_classified_tophit_aln} > {output.primary_AA_classified_tophit_aln_taxids};
        
        # Combine
        paste {output.prim_class_ids} {output.primary_AA_classified_tophit_aln_taxids} > {output.primary_AA_classified_tophit_aln_seq_taxids}
        
        # Add lineage information
        taxonkit lineage {output.primary_AA_classified_tophit_aln_seq_taxids} --data-dir {input.db} -i 2 > {output.primary_AA_classified_lineage};
        
        # Reformat lineage information
        taxonkit reformat --data-dir {input.db} {output.primary_AA_classified_lineage} -i 3 \
                -f "{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}" -F --fill-miss-rank | \
                cut --complement -f2,3 > {output.primary_AA_classified_lineage_reformated};
                
        ### Generate summaries at different taxonomic levels
        ### Primary AA Summaries
        # Order summary
        cut -f5 {output.primary_AA_classified_lineage_reformated} | sort | uniq -c | sed 's/^\s*//' | sed 's/ /\t/' | sort -k1 -nr | awk -F '\t' '{{ print$2"\t"$1 }}' > {output.primary_AA_order_summary};
        
        # Family summary
        cut -f6 {output.primary_AA_classified_lineage_reformated} | sort | uniq -c | sed 's/^\s*//' | sed 's/ /\t/' | sort -k1 -nr | awk -F '\t' '{{ print$2"\t"$1 }}' > {output.primary_AA_family_summary};
        
        # Genus summary
        cut -f7 {output.primary_AA_classified_lineage_reformated} | sort | uniq -c | sed 's/^\s*//' | sed 's/ /\t/' | sort -k1 -nr | awk -F '\t' '{{ print$2"\t"$1 }}' > {output.primary_AA_genus_summary};
        
        # Species summary
        cut -f8 {output.primary_AA_classified_lineage_reformated} | sort | uniq -c | sed 's/^\s*//' | sed 's/ /\t/' | sort -k1 -nr | awk -F '\t' '{{ print$2"\t"$1 }}' > {output.primary_AA_species_summary};
        
        ### Secondary AA Summaries
        # Order summary
        tail -n+2 {input.secondary_class_seqs} | \
            cut -f5 | sort | uniq -c | \
            sed 's/^\\s*//' | sed 's/ /\\t/' | sort -k1 -nr | \
            awk -F '\t' '{{ print$2"\t"$1 }}' > {output.secondary_AA_order_summary};
        
        # Family summary
        tail -n+2 {input.secondary_class_seqs} | \
            cut -f6 | sort | uniq -c | \
            sed 's/^\\s*//' | sed 's/ /\\t/' | sort -k1 -nr | \
            awk -F '\t' '{{ print$2"\t"$1 }}' > {output.secondary_AA_family_summary};
            
        # Genus summary
        tail -n+2 {input.secondary_class_seqs} | \
            cut -f7 | sort | uniq -c | \
            sed 's/^\\s*//' | sed 's/ /\\t/' | sort -k1 -nr | \
            awk -F '\t' '{{ print$2"\t"$1 }}' > {output.secondary_AA_genus_summary};
            
        # Species summary
        tail -n+2 {input.secondary_class_seqs} | \
            cut -f8 | sort | uniq -c | \
            sed 's/^\\s*//' | sed 's/ /\\t/' | sort -k1 -nr | \
            awk -F '\t' '{{ print$2"\t"$1 }}' > {output.secondary_AA_species_summary};
            
        # Join to compare primary to secondary
        # Order
        csvtk join -f1 {output.primary_AA_order_summary} {output.secondary_AA_order_summary} -t -T > {output.order_compare};
        
        # Family
        csvtk join -f1 {output.primary_AA_family_summary} {output.secondary_AA_family_summary} -t -T > {output.family_compare};
        
        # Genus
        csvtk join -f1 {output.primary_AA_genus_summary} {output.secondary_AA_genus_summary} -t -T > {output.genus_compare};
        
        # Species
        csvtk join -f1 {output.primary_AA_species_summary} {output.secondary_AA_species_summary} -t -T > {output.species_compare};
        
        # Add Baltimore classification to family comparison
        sed -i '1i Family\tPrimary\tSecondary' {output.family_compare};
        csvtk join -f 1 {output.family_compare} {input.baltimore} -t -T --left-join --na "NA" > {output.family_compare_baltimore};
        """

rule seqtable_untranslated_taxonomy_PRIMARY_search:
    """
    Combine sequences not assigned a virus taxonomy in the primary or secondary translated searches and query a viral nuecleotide database.
    
    This is a nt-to-nt untranslated search, similar to blastn.
    
    Reference database: All viral NCBI GenBank entries clustered at 100% identity
    
    """
    input:
        primary_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta"),
        secondary_seqs = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_unclassified.fasta"),
        db = os.path.join(NCBIVIRDB, "sequenceDB")
    output:
        aa_unclassified_seqs = os.path.join(PRIMARY_NT_OUT, "aa_unclassified.fasta")
        #inputDB = os.path.join(PRIMARY_NT_OUT, "inputDB")
        #tophit = os.path.join(PRIMARY_NT_OUT, "resultDB.firsthit.m8")
    params:
        alnRes=os.path.join(PRIMARY_NT_OUT, "MMSEQS_NT_PRIMARY"),
        tmppath=os.path.join(PRIMARY_NT_OUT, "mmseqs_nt_tmp")
    benchmark:
        "BENCHMARKS/mmseqs_nt.txt"
    log:
        "LOGS/mmseqs/mmseqs_nt.log"
    resources:
        mem_mb=64000,
        cpus=64
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        # Concatenate unclassified sequences from primary and secondary translated searches
        cat {input.primary_seqs} {input.secondary_seqs} > {output.aa_unclassified_seqs}
        
        # Create query database
        mmseqs createdb {output.aa_unclassified_seqs} {params.alnRes};
        """


