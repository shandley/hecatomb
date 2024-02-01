rule buildEnv:
    output:
        os.path.join(dir["out"]["temp"], "{env}.done")
    conda:
        lambda wildcards: os.path.join(dir["env"], wildcards.env)
    shell:
        "touch {output}"


rule trimnami_config:
    output:
        os.path.join(dir["out"]["temp"], "trimnami.config.yaml")
    params:
        config = trimnami,
    localrule:
        True
    run:
        import yaml
        with open(output[0],"w") as f:
            yaml.dump(params.config, f)


rule run_trimnami:
    """Get coverage statistics with Koverage"""
    input:
        tsv = os.path.join(dir["out"]["results"], "hecatomb.samples.tsv"),
        config = os.path.join(dir["out"]["temp"], "trimnami.config.yaml")
    output:
        targets["trimnami"]
    params:
        out_dir = os.path.join(dir["out"]["base"], "trimnami"),
        host = lambda w: "--host " + dir["dbs"]["hostFasta"] if not config["args"]["host"].lower() == "none" else "",
        trim = config["args"]["trim"],
        minimap_mode = lambda w: "map-ont " if config["args"]["trim"] == "nanopore" else "sr ",
        profile= lambda wildcards: "--profile " + config["args"]["profile"] if config["args"]["profile"] else "",
        fastqc = lambda wildcards: "--fastqc " if config["args"]["fastqc"] else "",
    threads:
        lambda wildcards: resources["sml"]["cpu"] if config["args"]["profile"] else resources["big"]["cpu"]
    resources:
        mem_mb = lambda wildcards: resources["sml"]["mem"] if config["args"]["profile"] else resources["big"]["mem"],
        mem = lambda wildcards: str(resources["sml"]["mem"]) + "MB" if config["args"]["profile"] else str(resources["big"]["mem"]) + "MB",
        time = resources["big"]["time"]
    conda:
        os.path.join(dir["env"], "trimnami.yaml")
    shell:
        "trimnami run {params.trim} "
            "--reads {input.tsv} "
            "--configfile {input.config} "
            "{params.host} "
            "--output {params.out_dir} "
            "--threads {threads} "
            "--minimap {params.minimap_mode} "
            "{params.fastqc} "
            "{params.profile}; "


rule cluster_sequences:
    input:
        fq= lambda wildcards: samples["trimmed"][wildcards.sample]["R1"],
    output:
        temp(os.path.join(dir["out"]["temp"],"{sample}_R1_rep_seq.fasta")),
        temp(os.path.join(dir["out"]["temp"],"{sample}_R1_cluster.tsv")),
        temp(os.path.join(dir["out"]["temp"],"{sample}_R1_all_seqs.fasta"))
    params:
        respath=lambda wildcards, output: os.path.split(output[0])[0],
        tmppath=lambda wildcards, output: os.path.join(os.path.split(output[0])[0],f"{wildcards.sample}_TMP"),
        prefix="{sample}_R1",
        config=config["mmseqs"]["linclustParams"]
    benchmark:
        os.path.join(dir["out"]["bench"],"cluster_similar_sequences.{sample}.txt")
    log:
        os.path.join(dir["out"]["stderr"],"cluster_similar_sequences.{sample}.log")
    resources:
        mem_mb=resources["big"]["mem"],
        mem=str(resources["big"]["mem"]) + "MB",
    threads:
        resources["big"]["cpu"]
    conda:
        os.path.join(dir["env"],"mmseqs2.yaml")
    shell:
        "mmseqs easy-linclust {input.fq} {params.respath}/{params.prefix} {params.tmppath} "
            "{params.config} "
            "--threads {threads} &> {log}; "


rule create_individual_seqtables:
    """Preprocessing step 06: Create individual seqtables. 

    A seqtable is a count table with each sequence as a row, each column as a sample and each cell the counts of each 
    sequence per sample.
    """
    input:
        seqs=os.path.join(dir["out"]["temp"],"{sample}_R1_rep_seq.fasta"),
        counts=os.path.join(dir["out"]["temp"],"{sample}_R1_cluster.tsv"),
    output:
        seqs=temp(os.path.join(dir["out"]["temp"],"{sample}_R1.seqs")),
        counts=temp(os.path.join(dir["out"]["temp"],"{sample}_R1.counts")),
        seqtable=temp(os.path.join(dir["out"]["temp"],"{sample}_R1.seqtable"))
    benchmark:
        os.path.join(dir["out"]["bench"],"individual_seqtables.{sample}.txt")
    log:
        os.path.join(dir["out"]["stderr"],"individual_seqtables.{sample}.txt")
    resources:
        mem_mb=resources["big"]["mem"],
        mem=str(resources["big"]["mem"]) + "MB",
    threads:
        resources["big"]["cpu"]
    conda:
        os.path.join(dir["env"],"seqkit.yaml")
    shell:
        "{{ seqkit sort {input.seqs} --quiet -j {threads} -w 5000 -t dna "
            "| seqkit fx2tab -w 5000 -t dna "
            "| sed 's/\\t\\+$//' "
            "| cut -f2,3 "
            "| sed '1i sequence' > {output.seqs}; "
        "cut -f1 {input.counts} "
            "| sort "
            "| uniq -c "
            "| awk -F ' ' '{{print$2\"\\t\"$1}}' "
            "| cut -f2 "
            "| sed '1i {wildcards.sample}' > {output.counts}; "
        "paste {output.seqs} {output.counts} > {output.seqtable}; }} 2> {log}; "


rule merge_seq_tables:
    """Preprocessing step 07: Merge seq tables

    Reads the sequences and counts from each samples' seqtable text file and converts to fasta format for the rest of 
    the pipline.
    """
    input:
        seqtables=expand(os.path.join(dir["out"]["temp"],"{sample}_R1.seqtable"),sample=samples["names"]),
    output:
        fa=os.path.join(dir["out"]["results"],"seqtable.fasta"),
        tsv=os.path.join(dir["out"]["results"],"sampleSeqCounts.tsv")
    params:
        samples=samples["names"],
        tmpdir=lambda wildcards, input: os.path.split(input[0])[0]
    conda:
        os.path.join(dir["env"],"pysam.yaml")
    benchmark:
        os.path.join(dir["out"]["bench"],"merge_seq_table.txt")
    log:
        os.path.join(dir["out"]["stderr"],"merge_seq_table.log")
    script:
        os.path.join(dir["scripts"],"mergeSeqTable.py")