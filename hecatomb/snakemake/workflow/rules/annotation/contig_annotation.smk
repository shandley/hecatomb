rule mmseqs_contig_annotation:
    """Contig annotation step 01: Assign taxonomy to contigs in contig_dictionary using mmseqs
    
    Database: NCBI virus assembly with taxID added
    """
    input:
        contigs=os.path.join(dir.out.results,f"{config.args.assembly}_assembly.fasta"),
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
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    group:
        "contigannot"
    shell:
        """
        {{
        mmseqs createdb {input.contigs} {output.queryDB} --dbtype 2;
        mmseqs search {output.queryDB} {input.db} {params.respath} {params.tmppath} \
            {params.sensnt} --split-memory-limit {params.memsplit} {params.filtnt} \
            --search-type 3 --threads {threads} ; }} &> {log}
        rm {log}
        """


rule mmseqs_contig_annotation_summary:
    """Contig annotation step 02: Summarize mmseqs contig annotation results"""
    input:
        queryDB=os.path.join(dir.out.assembly,"FLYE","queryDB"),
        db=os.path.join(dir.dbs.secondaryNT,"sequenceDB"),
        taxdb=dir.dbs.taxonomy
    output:
        result=os.path.join(dir.out.assembly,"FLYE","results","tophit.index"),
        align=os.path.join(dir.out.assembly,"FLYE","results","tophit.m8"),
        tsv = os.path.join(dir.out.results,"contigAnnotations.tsv")
    params:
        inputpath=os.path.join(dir.out.assembly,"FLYE","results","result"),
        respath=os.path.join(dir.out.assembly,"FLYE","results","tophit"),
        header=config.immutable.contigAnnotHeader,
        taxonFormat=lambda wildcards: config.immutable.taxonkitReformat,
        convertAliSummary=config.immutable.convertAliSummary
    benchmark:
        os.path.join(dir.out.bench, "mmseqs_contig_annotation_summary.txt")
    log:
        os.path.join(dir.out.stderr, "mmseqs_contig_annotation_summary.log")
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env, "mmseqs2.yaml")
    group:
        "contigannot"
    shell:
        """
        {{ 
        # Filter TopHit results
        mmseqs filterdb {params.inputpath} {params.respath} --extract-lines 1;
        
        # Convert to alignments
        mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} {params.convertAliSummary};
        
        # Header for output table
        printf "{params.header}\n" > {output.tsv};
        
        # Assign taxonomy
        sed 's/tid|//' {output.align} | \
            sed 's/|\S*//' | \
            taxonkit lineage --data-dir {input.taxdb} -i 2 | \
            taxonkit reformat --data-dir {input.taxdb} -i 19 {params.taxonFormat} | \
            cut --complement -f2,19 >> {output.tsv};
        }} &> {log}
        rm {log}
        """
