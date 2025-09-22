rule preprocessing_build_env:
    output:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], "{env}.done")
    conda:
        lambda wildcards: os.path.join("..", "..", "envs", wildcards.env)
    shell:
        "touch {output}"


rule preprocessing_cluster_sequences:
    input:
        fq= lambda wildcards: config["trimnami"]["trimmed"][wildcards.sample]["R1"],
    output:
        temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"],"{sample}_R1_rep_seq.fasta")),
        temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"],"{sample}_R1_cluster.tsv")),
        temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"],"{sample}_R1_all_seqs.fasta"))
    params:
        respath=lambda wildcards, output: os.path.split(output[0])[0],
        tmppath=lambda wildcards, output: os.path.join(os.path.split(output[0])[0], wildcards.sample + "_TMP"),
        prefix="{sample}_R1",
        config=config["hecatomb"]["mmseqs"]["linclustParams"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"],"cluster_similar_sequences.{sample}.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"],"cluster_similar_sequences.{sample}.log")
    resources:
        **config["resources"]["lrg"]
    threads:
        config["resources"]["lrg"]["cpu"]
    conda:
        os.path.join("..", "..", "envs","mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
    shell:
        "mmseqs easy-linclust {input.fq} "
            "{params.respath}/{params.prefix} "
            "{params.tmppath} "
            "{params.config} "
            "--threads {threads} "
            "&> {log}; "


rule preprocessing_create_individual_seqtables:
    """Preprocessing step 06: Create individual seqtables. 

    A seqtable is a count table with each sequence as a row, each column as a sample and each cell the counts of each 
    sequence per sample.
    """
    input:
        seqs=os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"],"{sample}_R1_rep_seq.fasta"),
        counts=os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"],"{sample}_R1_cluster.tsv"),
    output:
        seqs=temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"],"{sample}_R1.seqs")),
        counts=temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"],"{sample}_R1.counts")),
        seqtable=temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"],"{sample}_R1.seqtable"))
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"],"individual_seqtables.{sample}.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"],"individual_seqtables.{sample}.txt")
    resources:
        **config["resources"]["lrg"]
    threads:
        config["resources"]["lrg"]["cpu"]
    conda:
        os.path.join("..", "..", "envs","mmseqs2_seqkit_taxonkit_csvtk.yaml")
    container:
        config["hecatomb"]["container"]["mmseqs2_seqkit_taxonkit_csvtk"]
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


rule preprocessing_merge_seq_tables:
    """Preprocessing step 07: Merge seq tables

    Reads the sequences and counts from each samples' seqtable text file and converts to fasta format for the rest of 
    the pipline.
    """
    input:
        seqtables=expand(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"],"{sample}_R1.seqtable"),sample=samples["names"]),
    output:
        fa=os.path.join(config["hecatomb"]["args"]["output_paths"]["results"],"seqtable.fasta"),
        tsv=os.path.join(config["hecatomb"]["args"]["output_paths"]["results"],"sampleSeqCounts.tsv")
    params:
        samples=samples["names"],
        tmpdir=lambda wildcards, input: os.path.split(input[0])[0]
    conda:
        os.path.join("..", "..", "envs","krona_curl_zstd_pysam.yaml")
    container:
        config["hecatomb"]["container"]["krona_curl_zstd_pysam"]
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"],"merge_seq_table.txt")
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"],"merge_seq_table.log")
    script:
        os.path.join("..", "..", "scripts","mergeSeqTable.py")