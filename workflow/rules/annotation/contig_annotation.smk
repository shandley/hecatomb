rule contig_annotation_mmseqs_search:
    """Contig annotation step 01: Assign taxonomy to contigs in contig_dictionary using mmseqs
    
    Database: NCBI virus assembly with taxID added
    """
    input:
        contigs=os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], config["hecatomb"]["args"]["assembly"] + '_assembly.fasta'),
        db=os.path.join(config["hecatomb"]["args"]["database_paths"]["secondaryNT"], "sequenceDB")
    output:
        queryDB=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","queryDB"),
        result=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","results","result.index")
    params:
        resdir=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","results"),
        prefix=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","results", "result"),
        tmppath=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","mmseqs_nt_tmp"),
        sensnt = config["hecatomb"]["mmseqs"][config["hecatomb"]["args"]["search"]],
        memsplit = str(int(0.75 * int(config["resources"]["big"]["mem"]))) + "M",
        filtnt = config["hecatomb"]["mmseqs"]["filtNT"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "mmseqs_contig_annotation.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "mmseqs_contig_annotation.log")
    resources:
        **config["resources"]["big"]
    threads:
        config["resources"]["big"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    group:
        "contigannot"
    shell:
        "{{ if [[ -d {params.resdir} ]]; then rm -r {params.resdir}; fi; "
        "if [[ -d {params.tmppath} ]]; then rm -r {params.tmppath}; fi; "
        "mkdir -p {params.resdir}; "
        "mmseqs createdb {input.contigs} {output.queryDB} --dbtype 2; "
        "mmseqs search {output.queryDB} {input.db} {params.prefix} {params.tmppath} "
            "{params.sensnt} --split-memory-limit {params.memsplit} {params.filtnt} "
            "--search-type 3 --threads {threads} ; }} &> {log}"


rule contig_annotation_mmseqs_summary:
    """Contig annotation step 02: Summarize mmseqs contig annotation results"""
    input:
        queryDB=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","queryDB"),
        db=os.path.join(config["hecatomb"]["args"]["database_paths"]["secondaryNT"],"sequenceDB"),
        taxdb=config["hecatomb"]["args"]["database_paths"]["taxonomy"]
    output:
        result=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","results","tophit.index"),
        align=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","results","tophit.m8"),
        tsv = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"],"contigAnnotations.tsv")
    params:
        inputpath=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","results","result"),
        respath=os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"FLYE","results","tophit"),
        header=config["hecatomb"]["immutable"]["contigAnnotHeader"],
        secondaryNtFormat=config["hecatomb"]["immutable"]["secondaryNtFormat"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "mmseqs_contig_annotation_summary.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "mmseqs_contig_annotation_summary.log")
    resources:
        **config["resources"]["big"]
    threads:
        config["resources"]["big"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    group:
        "contigannot"
    shell:
        "{{ mmseqs filterdb {params.inputpath} {params.respath} --extract-lines 1; "
        "mmseqs convertalis {input.queryDB} {input.db} {params.respath} {output.align} {params.secondaryNtFormat}; "
        "printf '{params.header}\n' > {output.tsv}; "
        "sed 's/tid|//' {output.align} | "
            r"sed 's/|/\t/' | "
            "taxonkit lineage --data-dir {input.taxdb} -i 2 | "
            "taxonkit reformat --data-dir {input.taxdb} -i 19 "
                r"-f '{{k}}\t{{p}}\t{{c}}\t{{o}}\t{{f}}\t{{g}}\t{{s}}' -F --fill-miss-rank | "
            "cut --complement -f2,19 >> {output.tsv}; "
        "}} &> {log}; "
