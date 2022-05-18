rule mmseqs_contig_annotation:
    """Contig annotation step 01: Assign taxonomy to contigs in contig_dictionary using mmseqs
    
    Database: NCBI virus assembly with taxID added
    """
    input:
        contigs=os.path.join(RESULTS,"assembly.fasta"),
        db=os.path.join(POLYMICRODB, "sequenceDB")
    output:
        queryDB=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","queryDB"),
        result=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","results","result.index")
    params:
        respath=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","results","result"),
        tmppath=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","mmseqs_nt_tmp")
    benchmark:
        os.path.join(BENCH, "mmseqs_contig_annotation.txt")
    log:
        os.path.join(STDERR, "mmseqs_contig_annotation.log")
    resources:
        mem_mb=MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        os.path.join("../", "envs", "mmseqs2.yaml")
    shell:
        """
        {{
        mmseqs createdb {input.contigs} {output.queryDB} --dbtype 2;
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
            {MMSeqsSensNT} --split-memory-limit {MMSeqsMemSplit} {config[filtNTsecondary]} \
            --search-type 3 ; }} &> {log}
        rm {log}
        """


rule mmseqs_contig_annotation_summary:
    """Contig annotation step 02: Summarize mmseqs contig annotation results"""
    input:
        queryDB=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","queryDB"),
        db=os.path.join(POLYMICRODB, "sequenceDB"),
        taxdb=TAX
    output:
        result=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","results","tophit.index"),
        align=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","results","tophit.m8"),
        tsv = os.path.join(RESULTS, "contigAnnotations.tsv")
        # lineage=temporary(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt.lineage")),
        # reformated=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt.tsv"),
        # phyl_sum=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_phylum_summary.tsv"),
        # class_sum=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_class_summary.tsv"),
        # ord_sum=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_order_summary.tsv"),
        # fam_sum=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_family_summary.tsv"),
        # gen_sum=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_genus_summary.tsv"),
        # spe_sum=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","SECONDARY_nt_species_summary.tsv")
    params:
        inputpath=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","results","result"),
        respath=os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","FLYE","results","tophit"),
        header='\\\t'.join(['contigID',
                           'evalue',
                           'pident',
                           'fident',
                           'nident',
                           'mismatch',
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
                           'target',
                           'kingdom',
                           'phylum',
                           'class',
                           'order',
                           'family',
                           'genus',
                           'species\\\n'])
    benchmark:
        os.path.join(BENCH, "mmseqs_contig_annotation_summary.txt")
    log:
        os.path.join(STDERR, "mmseqs_contig_annotation_summary.log")
    resources:
        mem_mb=MMSeqsMem
    threads:
        MMSeqsCPU
    conda:
        os.path.join("../", "envs", "mmseqs2.yaml")
    shell:
        """
        {{ 
        # Filter TopHit results
        mmseqs filterdb {params.inputpath} {params.respath} --extract-lines 1;
        
        # Convert to alignments
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} \
            --format-output "query,target,evalue,pident,fident,nident,mismatch,qcov,tcov,qstart,qend,qlen,tstart,tend,tlen,alnlen,bits,target";
        
        # Header for output table
        printf {params.header} > {output.tsv};
        
        # Assign taxonomy
        sed 's/tid|//' {output.align} | \
            sed 's/|\S*//' | \
            taxonkit lineage --data-dir {input.taxdb} -i 2 | \
            taxonkit reformat --data-dir {input.taxdb} -i 19 -f "{{k}}\\t{{p}}\\t{{c}}\\t{{o}}\\t{{f}}\\t{{g}}\\t{{s}}" -F --fill-miss-rank | \
            cut --complement -f2,19 >> {output.tsv};
        }} 2> {log}
        rm {log}
        """
        # ## Generate summary tables
        # # Phylum
        # cut -f3 {output.reformated} | \
        #     csvtk freq -H -n -r -T -t | \
        #     sed '1i Phylum\tSECONDARY_NT_Phylum_Frequency'> {output.phyl_sum};
        # # Class
        # cut -f4 {output.reformated} | \
        #     csvtk freq -H -n -r -T -t | \
        #     sed '1i Class\tSECONDARY_NT_Class_Frequency'> {output.class_sum};
        # # Order
        # cut -f5 {output.reformated} | \
        #     csvtk freq -H -n -r -T -t | \
        #     sed '1i Order\tSECONDARY_NT_Order_Frequency'> {output.ord_sum};
        # # Family
        # cut -f6 {output.reformated} | \
        #     csvtk freq -H -n -r -T -t | \
        #     sed '1i Family\tSECONDARY_NT_Family_Frequency'> {output.fam_sum};
        # # Genus
        # cut -f7 {output.reformated} | \
        #     csvtk freq -H -n -r -T -t | \
        #     sed '1i Genus\tSECONDARY_NT_Genus_Frequency'> {output.gen_sum};
        # # Species
        # cut -f3 {output.reformated} | \
        #     csvtk freq -H -n -r -T -t | \
        #     sed '1i Species\tSECONDARY_NT_Species_Frequency'> {output.spe_sum}; }} &> {log}
        # rm {log}
        # """


