rule mmseqs_contig_annotation:
    """Contig annotation step 01: Assign taxonomy to contigs in contig_dictionary using mmseqs
    
    Database: NCBI virus assembly with taxID added
    """
    input:
        contigs=os.path.join(dir["out"]["results"], config["args"]["assembly"] + '_assembly.fasta'),
        db=os.path.join(dir["dbs"]["secondaryNT"], "sequenceDB")
    output:
        queryDB=os.path.join(dir["out"]["assembly"],"FLYE","queryDB"),
        result=os.path.join(dir["out"]["assembly"],"FLYE","results","result.index")
    params:
        resdir=os.path.join(dir["out"]["assembly"],"FLYE","results"),
        prefix=os.path.join(dir["out"]["assembly"],"FLYE","results", "result"),
        tmppath=os.path.join(dir["out"]["assembly"],"FLYE","mmseqs_nt_tmp"),
        sensnt = config["mmseqs"]["sens"],
        memsplit = str(int(0.75 * int(resources["big"]["mem"]))) + "M",
        filtnt = config["mmseqs"]["filtNT"]
    benchmark:
        os.path.join(dir["out"]["bench"], "mmseqs_contig_annotation.txt")
    log:
        os.path.join(dir["out"]["stderr"], "mmseqs_contig_annotation.log")
    resources:
        mem_mb = resources["big"]["mem"],
        mem = str(resources["big"]["mem"]) + "MB",
        time = resources["big"]["time"]
    threads:
        resources["big"]["cpu"]
    conda:
        os.path.join(dir["env"], "mmseqs2.yaml")
    group:
        "contigannot"
    shell:
        "{{ "
        "if [[ -d {params.resdir} ]]; then rm -r {params.resdir}; fi; "
        "if [[ -d {params.tmppath} ]]; then rm -r {params.tmppath}; fi; "
        "mkdir -p {params.resdir}; "
        "mmseqs createdb {input.contigs} {output.queryDB} --dbtype 2; "
        "mmseqs search {output.queryDB} {input.db} {params.prefix} {params.tmppath} "
            "{params.sensnt} --split-memory-limit {params.memsplit} {params.filtnt} "
            "--search-type 3 --threads {threads} ; }} &> {log}"


rule mmseqs_contig_annotation_summary:
    """Contig annotation step 02: Summarize mmseqs contig annotation results"""
    input:
        queryDB=os.path.join(dir["out"]["assembly"],"FLYE","queryDB"),
        db=os.path.join(dir["dbs"]["secondaryNT"],"sequenceDB"),
        taxdb=dir["dbs"]["taxonomy"]
    output:
        result=os.path.join(dir["out"]["assembly"],"FLYE","results","tophit.index"),
        align=os.path.join(dir["out"]["assembly"],"FLYE","results","tophit.m8"),
        tsv = os.path.join(dir["out"]["results"],"contigAnnotations.tsv")
    params:
        inputpath=os.path.join(dir["out"]["assembly"],"FLYE","results","result"),
        respath=os.path.join(dir["out"]["assembly"],"FLYE","results","tophit"),
        header=config["immutable"]["contigAnnotHeader"],
        taxonFormat=lambda wildcards: config["immutable"]["taxonkitReformat"],
        secondaryNtFormat=config["immutable"]["secondaryNtFormat"]
    benchmark:
        os.path.join(dir["out"]["bench"], "mmseqs_contig_annotation_summary.txt")
    log:
        os.path.join(dir["out"]["stderr"], "mmseqs_contig_annotation_summary.log")
    resources:
        mem_mb = resources["big"]["mem"],
        mem = str(resources["big"]["mem"]) + "MB",
        time = resources["big"]["time"]
    threads:
        resources["big"]["cpu"]
    conda:
        os.path.join(dir["env"], "mmseqs2.yaml")
    group:
        "contigannot"
    shell:
        "{{ "
        "mmseqs filterdb {params.inputpath} {params.respath} --extract-lines 1; "
        "mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} {params.secondaryNtFormat}; "
        "printf '{params.header}\n' > {output.tsv}; "
        "sed 's/tid|//' {output.align} | "
            r"sed 's/|\S*//' | "
            "taxonkit lineage --data-dir {input.taxdb} -i 2 | "
            "taxonkit reformat --data-dir {input.taxdb} -i 18 {params.taxonFormat} | "
            "cut --complement -f2,19 >> {output.tsv}; "
        "}} &> {log}; "
