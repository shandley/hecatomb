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
        os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_lca.tsv"),
        os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_report"),
        os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_tophit_report"),
        os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln"),
        os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted")
        #os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_lca.ids"),
        #os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted.ids"),
        #os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.ids"),
        #os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
        #os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta"),
        #os.path.join(AA_OUT, "mmseqs_primary_aa_search_summary.txt")

rule seqtable_translated_taxonomy_primary_search:
    """
    Assign taxonomy to RESULTS/seqtable.fasta sequences using mmseqs2 to query a viral protein reference databases.
    
    This is a nt-to-aa translated search, similar to blastX.
    
    All sequences assigned to a viral lineage will be reserved for the secondary query against a trans-kingdom database to confirm their viral lineage.
    
    Reference database: All Uniprot viral proteins clustered at 99% identity.

    """
    input:
        seqs = os.path.join(RESULTS, "seqtable.fasta"),
        db = os.path.join(UNIVIRDB, "sequenceDB")
    output:
        lca = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_lca.tsv"),
        report = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_report"),
        tophit_report = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_tophit_report"),
        aln = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln"),
        alnsort = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted")
    params:
        alnRes=os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY"),
        tmppath=os.path.join(AA_OUT, "mmseqs_aa_tmp")
    benchmark:
        "BENCHMARKS/mmseqs_pviral_aa.txt"
    log:
        "LOGS/mmseqs/mmseqs_pviral_aa.log"
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
        --tax-lineage 1
        
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        """
        
rule seqtable_translated_taxonomy_primary_tabulation:
    """
    1) Parse potential viral and unclassified sequences
    
    2) Tabulate some statistics
        - Number of input sequences: seqtable.fasta
        - Number of sequences classified as having a viral lineage: MMSEQS_AA_PRIMARY_classified.fasta
        - Number of sequences having no inidication of having a viral lineage: MMSEQS_AA_PRIMARY_unclassified.fasta
        
    MMSEQS_AA_PRIMARY_classified.fasta will be subjected to a secondary translated search against a trans-kingdom database to confirm true virality
    
    MMSEQS_AA_PRIMARY_unclassified.fasta will be queried against a viral nucleotide database to detect similarity to non-coding regions of viral genomes or to sequences not represented in the UniProt protein databases.
        
    """
    input:
        lca = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_lca.tsv"),
        alnsort = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted"),
        seqs = os.path.join(RESULTS, "seqtable.fasta")
    output:
        allseq_ids = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_lca.ids"),
        tophit_ids = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted.ids"),
        unclass_ids = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.ids"),
        class_seqs = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
        unclass_seqs = os.path.join(AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta"),
        summary = os.path.join(AA_OUT, "mmseqs_primary_aa_search_summary.txt")
    conda:
        "../envs/samtools.yaml"
    log:
        "LOGS/mmseqs/classified_sequence_parsing.log"
    shell:
        """
        # Extract full entry ID list (all sequences)
        sort -k1 -n {input.lca} | cut -f1 > {output.allseq_ids};
        
        # Extract tophit entry ids
        tail -n+2 {input.alnsort} | cut -f1 > {output.tophit_ids};
        
        # Extract ids that failed LCA
        awk '$2 =='
        
        # Compare tophit ids to all ids to extract ids with no hits
        comm -23 <(sort {output.allseq_ids}) <(sort {output.tophit_ids}) | sort -n > {output.unclass_ids};
        
        # Extract classified sequencess
        xargs samtools faidx {input.seqs} -n 5000 < {output.tophit_ids} > {output.class_seqs};
        
        # Extract unclassified sequencess
        xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
        
        # Summary table
        INPUT=$(grep -c ">" {input.seqs})
        CLASS=$(grep -c ">" {output.class_seqs})
        UNCLASS=$(grep -c ">" {output.unclass_seqs})
        LCA=$(wc -l {input.lca})
        pCLASS=$(echo "scale=4 ; 100*($CLASS / $INPUT)" | bc)
        pUNCLASS=$(echo "scale=4 ; 100*($UNCLASS / $INPUT)" | bc)
        pLCA=$(echo "scale=4 ; 100*($LCA / $INPUT)" | bc)
        
        cat <<- EOF | column -t -s $'\t' > {output.summary}
        FILE    VALUE    PERCENT    NOTES
        seqtable.fasta    $INPUT    1    Number of input sequenes
        MMSEQS_AA_PRIMARY_classified.fasta    $CLASS    $pCLASS    Number of classifiable sequences
        MMSEQS_AA_PRIMARY_unclassified.fasta $UNCLASS    $pUNCLASS    Number of unclassifiable sequences
        MMSEQS_AA_PRIMARY_lca.tsv    $LCA    $pLCA    Number of sequences with classifiable LCA
        
        """

rule seqtable_translated_taxonomy_primary_search:
    """
    Check taxonomic assignments made in rule seqtable_translated_taxonomy_primary_search using mmseqs2 to a transkingdom reference amino acid database.
    
    Sequences assigned to a viral lineage ("d_Viruses") from prior rule are checked.
    
    This is a nt-to-aa translated search, similar to blastX.
    
    Reference database: UniRef50 combined with all Uniprot viral proteins clustered at 99% identity (UniRef50 + virus).
        - The origianl viral protein sequence are added back to ensure no viral proteins are lost in the constrction of the UniRef50 database
    
    """
    input:
        seqs = os.path.join(AA_OUT, "MMSEQS_AA_SECONDARY_classified.fasta"),
        db = os.path.join(UNIREF50VIR, "sequenceDB")
    output:
        lca = os.path.join(AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv"),
        report = os.path.join(AA_OUT, "MMSEQS_AA_SECONDARY_report"),
        tophit_report = os.path.join(AA_OUT, "MMSEQS_AA_SECONDARY_tophit_report"),
        aln = os.path.join(AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln"),
        alnsort = os.path.join(AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln_sorted")
    params:
        alnRes=os.path.join(AA_OUT, "MMSEQS_AA_SECONDARY"),
        tmppath=os.path.join(AA_OUT, "mmseqs_aa_tmp")
    benchmark:
        "BENCHMARKS/mmseqs_pviral_aa.txt"
    log:
        "LOGS/mmseqs/mmseqs_pviral_aa.log"
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
        --tax-lineage 1
        
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        """
        
    
