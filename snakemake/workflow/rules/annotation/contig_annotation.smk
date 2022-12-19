rule mmseqs_contig_annotation:
    """Contig annotation step 01: Assign taxonomy to contigs in contig_dictionary using mmseqs
    
    Database: NCBI virus assembly with taxID added
    """
    input:
        contigs=os.path.join(dir.out.results,"assembly.fasta"),
        db=os.path.join(dir.dbs.secondaryNT, "sequenceDB")
    output:
        queryDB=os.path.join(dir.out.assembly,"FLYE","queryDB"),
        result=os.path.join(dir.out.assembly,"FLYE","results","result.index")
    params:
        respath=os.path.join(dir.out.assembly,"FLYE","results","result"),
        tmppath=os.path.join(dir.out.assembly,"FLYE","mmseqs_nt_tmp"),
        sensnt = config.mmseqs.sensNT,
        memsplit = str(int(0.75 * int(config.resources.big.mem))) + 'M',
        filtnt = config.mmseqs.filtNTsecondary
    benchmark:
        os.path.join(dir.out.bench, "mmseqs_contig_annotation.txt")
    log:
        os.path.join(dir.out.stderr, "mmseqs_contig_annotation.log")
    resources:
        mem_mb=config.resources.big.mem
    threads:
        config.resources.big.cpu
    conda:
        os.path.join("../", "envs", "mmseqs2.yaml")
    shell:
        """
        {{
        mmseqs createdb {input.contigs} {output.queryDB} --dbtype 2;
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
            {params.sensnt} --split-memory-limit {params.memsplit} {params.filtnt} \
            --search-type 3 ; }} &> {log}
        rm {log}
        """


rule mmseqs_contig_annotation_summary:
    """Contig annotation step 02: Summarize mmseqs contig annotation results"""
    input:
        queryDB=os.path.join(dir.out.assembly,"FLYE","queryDB"),
        db=os.path.join(dir.dbs.secondaryNT, "sequenceDB"),
        taxdb=dir.dbs.taxonomy
    output:
        result=os.path.join(dir.out.assembly,"FLYE","results","tophit.index"),
        align=os.path.join(dir.out.assembly,"FLYE","results","tophit.m8"),
        tsv = os.path.join(dir.out.results, "contigAnnotations.tsv")
    params:
        inputpath=os.path.join(dir.out.assembly,"FLYE","results","result"),
        respath=os.path.join(dir.out.assembly,"FLYE","results","tophit"),
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
                           'species\\\\n'])
    benchmark:
        os.path.join(dir.out.bench, "mmseqs_contig_annotation_summary.txt")
    log:
        os.path.join(dir.out.stderr, "mmseqs_contig_annotation_summary.log")
    resources:
        mem_mb=config.resources.big.mem
    threads:
        config.resources.big.cpu
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
