"""
Hecatomb.smk to query target amino acid sequence database with reduced (seqtab) clustered sequences from merge_seqtable.sh using mmseqs2

History: This is based on [mmseqs_pviral_aa.sh](../base/mmseqs_pviral_aa.sh)

Rob Edwards, March 2020
Overhauled - Michael Roach, Q2/Q3 2021

"""

import os
import sys


rule PRIMARY_AA_taxonomy_assignment:
    """Assign taxonomy to RESULTS/seqtable.fasta sequences using mmseqs2 
    
    - Reference database: all UniProt viral (UNIVIRDB) protein sequences clustered at 99% ID
    - This is a nt-to-aa translated search, similar to blastX
    - All sequences assigned to a viral lineage will be reserved for the secondary query against a trans-kingdom database 
      (UniRef50 + UNIVIRDB) to confirm their viral lineage
    - The sequences checked in the next step are selected if their tophit is in a viral lineage (tophit_aln file). This 
      means they may not have LCA assigned taxonomy, but it is the most inclusive group to send to the next step 
      (secondary search) where any false-positives will be sorted out. For now, we want to capture anything that 
      potentially looks viral.
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
        log = os.path.join(STDERR, "MMSEQS", "mmseqs_primary_AA.log")
    resources:
        mem_mb=MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/mmseqs2.yaml"
    shell: # run easy taxonomy, add header
        """
        # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
            -a --start-sens 1 --sens-steps 3 -s 7 \
            --tax-output-mode 2 --search-type 2 --lca-mode 2 --shuffle 0 \
            --lca-ranks "superkingdom,phylum,class,order,family,genus,species" \
            --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
            --tax-lineage 1 --min-length {config[AAMINLEN]} \
            --threads {threads} --split-memory-limit {MMSeqsMemSplit} \
            -e {config[PRIMAAE]} &>> {log};
        
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        """
        
rule PRIMARY_AA_parsing:
    """Parse primary AA search results for classified (potentially viral) and unclassified sequences
        
    - MMSEQS_AA_PRIMARY_classified.fasta will be subjected to a secondary translated search (rule 
      SECONDARY_AA_taxonomy_assignment:) against a trans-kingdom database (UniRef50 + UNIVIRDB) to confirm true 
      viral lineage
    - MMSEQS_AA_PRIMARY_unclassified.fasta will be queried against a viral nucleotide database (nt-to-nt) to detect 
      similarity to non-coding regions of viral genomes or to sequences not represented in the UniProt protein databases
    """
    input:
        alnsort = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_tophit_aln_sorted"),
        seqs = os.path.join(RESULTS, "seqtable.fasta"),
    output:
        class_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_classified.fasta"),
        unclass_seqs = os.path.join(PRIMARY_AA_OUT, "MMSEQS_AA_PRIMARY_unclassified.fasta")
    # conda:
    #     "../envs/samtools.yaml"
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    run:
        topHit = {}
        for l in stream_tsv(input.alnsort):
            topHit[l[0]] = 1
        outClass = open(output.class_seqs, 'w')
        outUnclass = open(output.unclass_seqs, 'w')
        inFa = open(input.seqs,'r')
        for line in inFa:
            if line.startswith('>'):
                id = line.strip().replace('>','')
                seq = inFa.readline().strip()
                try:
                    topHit[id]
                    outClass.write(f'>{id}\n{seq}\n')
                except KeyError:
                    outUnclass.write(f'>{id}\n{seq}\n')
            else:
                sys.stderr.write(f'malformed {input.seqs} file, complain to Mike.')
                exit(1)
        inFa.close()
        outClass.close()
        outUnclass.close()


rule SECONDARY_AA_taxonomy_assignment:
    """Check taxonomic assignments in MMSEQS_AA_PRIMARY_classified.fasta using mmseqs2
    
    - Reference database: UniRef50 + UNIVIRDB. All UniProtKB protein entries (all domains of life) clustered at 50% ID 
      (https://www.uniprot.org/help/uniref) concatenated to UNIVIRDB
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
        log = os.path.join(STDERR, "MMSEQS", "mmseqs_secondary_AA.log")
    resources:
        mem_mb=MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/mmseqs2.yaml"
    shell: # secondary easy-tax search, add header
        """
        # Run mmseqs taxonomy module
        mmseqs easy-taxonomy {input.seqs} {input.db} {params.alnRes} {params.tmppath} \
            -a --start-sens 1 --sens-steps 3 -s 7 \
            --tax-output-mode 2 --search-type 2 --lca-mode 2 --shuffle 0 \
            --lca-ranks "superkingdom,phylum,class,order,family,genus,species" \
            --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader,taxid,taxname,taxlineage" \
            --tax-lineage 1 --split-memory-limit {MMSeqsMemSplit} --min-length {config[AAMINLEN]} \
            --threads {threads} \
            -e {config[SECAAE]} &>> {log};
        
        # Add headers
        sort -k1 -n {output.aln} | \
            sed '1i query\ttarget\tevalue\tpident\tfident\tnident\tmismatch\tqcov\ttcov\tqstart\tqend\tqlen\ttstart\ttend\ttlen\talnlen\tbits\tqheader\ttheader\ttaxid\ttaxname\tlineage' > {output.alnsort};
        """
        
rule SECONDARY_AA_tophit_lineage:
    """Add/reformat tophit viral lineages with up-to-date* NCBI taxonomy
    
    - UniRef50 (https://www.uniprot.org/help/uniref) and other UniRef databases are a mixture of protein entries from 
      UniProtKB (https://www.uniprot.org/help/uniprotkb) and UniParc (https://www.uniprot.org/help/uniparc)
    - The taxonomic lineage assignment algorithm of mmseqs2 will not provide full annotation of UniParc entries
    - We can updated those here using taxonkit lineage
    - This is also to ensure that the taxa are annotated based on the NCBI taxonomy which you can rapidly update in the 
      databases/tax/taxonomy directory
    """
    input:
        db = TAX,
        tophit = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln")
    output:
        tophit_lineage_refomated = os.path.join(SECONDARY_AA_OUT, "tophit.lineage.reformated"),
    conda:
        "../envs/seqkit.yaml"
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    log:
        log = os.path.join(STDERR, "MMSEQS", "mmseqs_secondary_tophit_refactor.log")
    shell:
        """
        # Make a table: SeqID <tab> taxID
        cut -f1,20 {input.tophit} \
            | taxonkit lineage --data-dir {input.db} -i 2 \
            | taxonkit reformat --data-dir {input.db} -i 3 \
                -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F --fill-miss-rank \
            | cut --complement -f3 \
            > {output.tophit_lineage_refomated};
        """


rule SECONDARY_AA_refactor_finalize:
    """Remove sequences to be refactored from LCA table and recombine with updated taxonomies."""
    input:
        db = TAX,
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.tsv"),
    output:
        lca_reformated = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.reformated"),
    conda:
        "../envs/seqkit.yaml"
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    log:
        log = os.path.join(STDERR, "MMSEQS", "mmseqs_secondary_lca_refactor_final.log")
    shell:
        """
        cut -f1,2 {input.lca} \
            | taxonkit lineage --data-dir {input.db} -i 2 \
            | taxonkit reformat --data-dir {input.db} -i 3 \
            -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F --fill-miss-rank \
            | cut --complement -f3 \
            > {output.lca_reformated};
        """


rule SECONDARY_AA_generate_output_table:
    """Join sequence info, tophit align info, and LCA or tophit lineage info into the output format table"""
    input:
        aln = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_tophit_aln"),
        lca = os.path.join(SECONDARY_AA_OUT, "MMSEQS_AA_SECONDARY_lca.reformated"),
        top = os.path.join(SECONDARY_AA_OUT, "tophit.lineage.reformated"),
    output:
        os.path.join(SECONDARY_AA_OUT,"AA_bigtable.tsv"),
    run:
        lcaLin = {}
        for l in stream_tsv(input.lca):
            # dont use lca lineage for taxids of 0, 1, or 10239
            if l[1] != '0' and l[1] != '1' and l[1] != '10239':
                lcaLin[l[0]] = '\t'.join((l[2:]))
        topLin = {}
        for l in stream_tsv(input.top):
            # skip if using lca lineage
            try:
                lcaLin[l[0]]
            except KeyError:
                topLin[l[0]] = '\t'.join((l[2:]))
        # iterate the tophit alignments and print the table on the fly
        out = open(output[0],'w')
        out.write('\t'.join(('seqID',
                             'sampleID',
                             'count',
                             'alnType',     # aa or nt
                             'target',
                             'evalue',
                             'pident',
                             'fident',
                             'nident',
                             'mismatches',
                             'qcov',
                             'tcov',
                             'qstart',
                             'qend',
                             'qlen',
                             'tstart',
                             'tend',
                             'tlen',
                             'alnlen',
                             'bits',
                             'taxMethod',
                             'kingdom',
                             'phylum',
                             'class',
                             'order',
                             'family',
                             'genus',
                             'species')))
        out.write('\n')
        # parse the alignments and attach appropriate info
        for l in stream_tsv(input.aln):
            try:
                taxOut = 'LCA\t' + lcaLin[l[0]]
            except KeyError:
                try:
                    taxOut = 'TopHit\t' + topLin[l[0]]
                except KeyError:
                    taxOut = '\t'.join((['NA'] * 8))
            seqInf = l[0].split(':')        # seq ID = sample:count:seqNum
            seqOut = '\t'.join((l[0], seqInf[0], seqInf[1]))
            alnOut = 'aa\t' + '\t'.join((l[1:17]))
            out.write('\t'.join((seqOut, alnOut, taxOut)))
            out.write('\n')
        out.close()

rule SECONDARY_AA_parsing:
    """Parse out all sequences that remain unclassified following the Secondary AA search and refactoring."""
    input:
        bigtable = os.path.join(SECONDARY_AA_OUT,"AA_bigtable.tsv"),
        seqs = os.path.join(RESULTS, "seqtable.fasta")
    output:
        unclass_seqs = os.path.join(SECONDARY_AA_OUT, "translated_unclassified.fasta")
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    log:
        os.path.join(STDERR, 'mmseqs', 'mmseqs_SECONDARY_aa_parsing.log')
    run:
        virSeqs = {}
        for l in stream_tsv(input.bigtable):
            virSeqs[l[0]] = 1
        inFa = open(input.seqs, 'r')
        outFa = open(output.unclass_seqs, 'w')
        for line in inFa:
            if line.startswith('>'):
                id = line.strip().replace('>','')
                seq = inFa.readline().strip()
                try:
                    virSeqs[id]
                except KeyError:
                    outFa.write(f'>{id}\n{seq}\n')
            else:
                sys.stderr.write(f'malformed {input.seqs} file, complain to Mike')
                exit(1)
        outFa.close()
        inFa.close()

        
rule PRIMARY_NT_taxonomic_assignment:
    """Primary nucleotide search of unclassified viral-like sequences from aa search"""
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
        log = os.path.join(STDERR, "MMSEQS", "mmseqs_primary_NT.log")
    resources:
        mem_mb=MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        # Create query database
        mmseqs createdb {input.seqs} {output.queryDB} --dbtype 2;
        
        # mmseqs search
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
            --start-sens 2 -s 7 --sens-steps 3 \
            --search-type 3 --min-length {config[NTMINLEN]} \
            --threads {threads} --split-memory-limit {MMSeqsMemSplit} \
            -e {config[PRIMNTE]} &>> {log};
        """
        
rule PRIMARY_NT_reformat:
    """Collect some summary statistics on primay nucleotide search"""
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
        "../envs/mmseqs2.yaml"
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    log:
        os.path.join(STDERR, 'MMSEQS', 'mmseqs_PRIMARY_nt_summary.log')
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
        """

        
rule PRIMARY_NT_parsing:
    """Extract unclassified sequences from secondary translated search
    
    xargs samtools faidx {input.seqs} -n 5000 < {output.unclass_ids} > {output.unclass_seqs};
    """
    input:
        seqs = os.path.join(SECONDARY_AA_OUT, "translated_unclassified.fasta"),
        align = os.path.join(PRIMARY_NT_OUT, "results", "result.m8")
    output:
        class_seqs = os.path.join(PRIMARY_NT_OUT, "classified_seqs.fasta"),
        unclass_seqs = os.path.join(PRIMARY_NT_OUT, "unclassified_seqs.fasta")
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    run:
        hit = {}
        inAln = open(input.align,'r')
        for l in stream_tsv(input.align):
            hit[l[0]] = 1
        inFa = open(input.seqs, 'r')
        outClass = open(output.class_seqs, 'w')
        outUnclass = open(output.unclass_seqs, 'w')
        for line in inFa:
            if line.startswith('>'):
                id = line.strip().replace('>','')
                seq = inFa.readline().strip()
                try:
                    hit[id]
                    outClass.write(f'>{id}\n{seq}\n')
                except KeyError:
                    outUnclass.write(f'>{id}\n{seq}\n')
            else:
                sys.stderr.write(f'malformed {input.seqs} file, complain to Mike')
                exit(1)
        inFa.close()
        outClass.close()
        outUnclass.close()


rule SECONDARY_NT_taxonomic_assignment:
    """Secondary nucleotide search of viral hits from primary nucleotide search"""
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
        log = os.path.join(STDERR, "MMSEQS", "mmseqs_secondary_NT.log")
    resources:
        mem_mb=MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/mmseqs2.yaml"
    shell:
        """
        # Create query database
        mmseqs createdb {input.seqs} {output.queryDB} --dbtype 2;
        
        # mmseqs search
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
            -a 1 --search-type 3 \
            --start-sens 2 -s 7 --sens-steps 3 --min-length {config[NTMINLEN]} \
            --threads {threads} --split-memory-limit {MMSeqsMemSplit} \
            -e {config[SECNTE]} &>> {log};
    
        """

rule SECONDARY_NT_summary:
    """Summary statistics for secondary nucleotide search"""
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
        "../envs/mmseqs2.yaml"
    resources:
        mem_mb=MiscMem
    threads:
        MiscCPU
    log:
        os.path.join(STDERR, 'MMSEQS', 'mmseqs_SECONDARY_nt_summary.log')
    shell:
        """
        # Filter TopHit results
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
        """

        
rule SECONDARY_NT_convert:
    """Reformat secondary search LCA taxon assignment"""
    input:
        queryDB = os.path.join(SECONDARY_NT_OUT, "queryDB"),
        db = os.path.join(POLYMICRODB, "sequenceDB")
    output:
        align = os.path.join(SECONDARY_NT_OUT, "results", "all.m8"),
    params:
        respath = os.path.join(SECONDARY_NT_OUT, "results", "result")
    resources:
        mem_mb=MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/mmseqs2.yaml"
    log:
        os.path.join(STDERR, 'MMSEQS', 'mmseqs_SECONDARY_nt_convert.log')
    shell:
        """
        # Convert to alignments
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} \
            --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,qheader,theader" \
            2> {log}
        """


rule secondary_nt_lca_table:
    """Create table for taxonkit lineage for secondary NT search"""
    input:
        align = os.path.join(SECONDARY_NT_OUT, "results", "all.m8")
    output:
        lin = temp(os.path.join(SECONDARY_NT_OUT, "results", "all.lin"))
    run:
        lin = {}
        for l in stream_tsv(input.align):
            t = l[1].split('|')     # e.g. tid|2293023|NZ_QUCO01000049.1
            try:
                if not t[1] in lin[l[0]]:
                    lin[l[0]].append(t[1])
            except KeyError:
                lin[l[0]] = [t[1]]
        out = open(output.lin, 'w')
        for s in lin.keys():
            tOut = ':'.join(lin[s])
            out.write(f'{s}\t{tOut}\n')
        out.close()


rule secondary_nt_calc_lca:
    """Calculate the lca for the secondary NT search"""
    input:
        lin = os.path.join(SECONDARY_NT_OUT, "results", "all.lin"),
        taxdb = TAX
    output:
        lca_lineage = os.path.join(SECONDARY_NT_OUT, "results", "lca.lineage"),
        reformated = os.path.join(SECONDARY_NT_OUT, "results", "secondary_nt_lca.tsv")
    resources:
        mem_mb=MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        "../envs/mmseqs2.yaml"
    log:
        os.path.join(STDERR, 'MMSEQS', 'mmseqs_SECONDARY_nt_calc.log')
    shell:
        """
        # calculate lineages
        taxonkit lineage -i 2 --data-dir {input.taxdb} {input.lin} > {output.lca_lineage} 2> {log}
        
        # Reformat lineages
        awk -F '\t' '$2 != 0' {output.lca_lineage} | \
            taxonkit reformat --data-dir {input.taxdb} -i 3 \
                -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F --fill-miss-rank 2>> {log} |
            cut --complement -f3 > {output.reformated}
        """


rule SECONDARY_NT_generate_output_table:
    input:
        aln = os.path.join(SECONDARY_NT_OUT, "results", "tophit.m8"),
        top = os.path.join(SECONDARY_NT_OUT, "SECONDARY_nt.tsv"),
        lca = os.path.join(SECONDARY_NT_OUT, "results", "secondary_nt_lca.tsv")
    output:
        os.path.join(SECONDARY_NT_OUT, "NT_bigtable.tsv")
    run:
        lcaLin = {}
        for l in stream_tsv(input.lca):
            # dont use lca lineage for taxids of 0, 1, or 10239
            if l[1] != '0' and l[1] != '1' and l[1] != '10239':
                lcaLin[l[0]] = '\t'.join((l[2:]))
        topLin = {}
        for l in stream_tsv(input.top):
            # skip if using lca lineage
            try:
                lcaLin[l[0]]
            except KeyError:
                topLin[l[0]] = '\t'.join((l[2:]))
        # output
        out = open(output[0], 'w')
        out.write('\t'.join(('seqID',
                             'sampleID',
                             'count',
                             'alnType',     # aa or nt
                             'target',
                             'evalue',
                             'pident',
                             'fident',
                             'nident',
                             'mismatches',
                             'qcov',
                             'tcov',
                             'qstart',
                             'qend',
                             'qlen',
                             'tstart',
                             'tend',
                             'tlen',
                             'alnlen',
                             'bits',
                             'taxMethod',
                             'kingdom',
                             'phylum',
                             'class',
                             'order',
                             'family',
                             'genus',
                             'species')))
        out.write('\n')
        for l in stream_tsv(input.aln):
            try:
                taxOut = 'LCA\t' + lcaLin[l[0]]
            except KeyError:
                try:
                    taxOut = 'TopHit\t' + topLin[l[0]]
                except KeyError:
                    taxOut = '\t'.join((['NA'] * 8))
            seqInf = l[0].split(':')        # seq ID = sample:count:seqNum
            seqOut = '\t'.join((l[0], seqInf[0], seqInf[1]))
            alnOut = 'nt\t' + '\t'.join((l[1:17]))
            out.write('\t'.join((seqOut, alnOut, taxOut)))
            out.write('\n')
        out.close()

rule combine_AA_NT:
    input:
        aa = os.path.join(SECONDARY_AA_OUT,"AA_bigtable.tsv"),
        nt = os.path.join(SECONDARY_NT_OUT, "NT_bigtable.tsv")
    output:
        os.path.join(RESULTS, "bigtable.tsv")
    shell:
        """
        cat {input.aa} > {output};
        tail -n+2 {input.nt} >> {output};
        """

rule krona_text_format:
    input:
        os.path.join(RESULTS, "bigtable.tsv")
    output:
        os.path.join(RESULTS, "kronaText.tsv")
    run:
        counts = {}
        for l in stream_tsv(input[0]):
            if l[0]=="seqID":
                continue
            t = '\t'.join(l[21:])
            try:
                counts[t] += int(l[2])
            except KeyError:
                counts[t] = int(l[2])
        outFH = open(output[0],'w')
        for k in sorted(counts.keys()):
            outFH.write(f'{counts[k]}\t{k}\n')
        outFH.close()

rule krona_plot:
    input:
        os.path.join(RESULTS,"kronaText.tsv")
    output:
        os.path.join(RESULTS, "kronaPlot.html")
    conda:
        os.path.join('../', 'envs', 'krona.yaml')
    shell:
        """
        ktImportText {input} -o {output}
        """
