
rule assembly_lr_cross_assembly:
    """Alternative to merged assembly; assemble everything together in one hit."""
    input:
        [ x for x in config["trimnami"]["targets"]["reads"] if not x.endswith("multiqc.html") ]
    output:
        assembly = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "cross_assembly.fasta"),
        graph = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "cross_assembly.gfa"),
        tar = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"crossAssembly.tar.zst")
    params:
        params = config["hecatomb"]["assembly"]["metaflye"],
        dir = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"], "crossAssembly"),
        assembly = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"], "crossAssembly", "assembly.fasta"),
        graph = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"], "crossAssembly", "assembly_graph.gfa"),
    resources:
        **config["resources"]["big"]
    threads:
        config["resources"]["big"]["cpu"]
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "canu_cross_assembly.log")
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "canu_cross_assembly.txt")
    conda:
        os.path.join("..", "..", "envs", "bbmap_bedtools_pigz_flye.yaml")
    container:
        config["hecatomb"]["container"]["bbmap_bedtools_pigz_flye"]
    shell:
        "flye -o {params.dir} -t {threads} {params.params} {input} 2> {log}; "
        "mv {params.assembly} {output.assembly}; "
        "mv {params.graph} {output.graph}; "
        "tar cf - {params.dir} | zstd -T{threads} -9 > {output.tar} 2> {log}; "


rule assembly_canu_sample_assembly:
    """Per-sample assembly with canu; also works for unmapped rescue reads"""
    input:
        lambda wildcards: config["trimnami"]["trimmed"][wildcards.sample]["S"] if \
            "S" in config["trimnami"]["trimmed"][wildcards.sample] else config["trimnami"]["trimmed"][wildcards.sample]["R1"]
    output:
        ctg = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"{sample}","{sample}.contigs.fasta"),
        ctgq = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"{sample}","{sample}.contigs.uniq.fasta"),
        un = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"{sample}","{sample}.unassembled.fasta"),
        unq = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"{sample}","{sample}.unassembled.uniq.fasta"),
        tar = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"{sample}.tar.zst")
    params:
        settings = config["hecatomb"]["assembly"]["canu"],
        canu_dir = lambda w, output: os.path.split(output.ctg)[0]
    resources:
        **config["resources"]["lrg"]
    threads:
        config["resources"]["lrg"]["cpu"]
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "canu_sample_assembly.{sample}.log")
    conda:
        os.path.join("..", "..", "envs", "canu.yaml")
    container:
        config["hecatomb"]["container"]["krona_curl_zstd_pysam"]
    shell:
        "{{ canu {params.settings} {input} "
            "batThreads={threads} "
            "batMemory={resources.mem_mb}M "
            "-p {wildcards.sample} "
            "-d {params.canu_dir}; "
        "sed 's/>tig/>{wildcards.sample}./' {output.ctg} > {output.ctgq}; "
        "sed 's/>tig/>{wildcards.sample}./' {output.un} > {output.unq}; "
        "tar cf - {params.canu_dir} | zstd -T{threads} -9 > {output.tar}; }} 2> {log}; "


rule assembly_combine_canu_unassembled:
    """Combine the unassembled reads from all canu assemblies"""
    input:
        expand(os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"{sample}","{sample}.unassembled.uniq.fasta"), sample=samples["names"])
    output:
        temp(os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"unmappedRescue_R1.all.fasta.gz"))
    resources:
        **config["resources"]["lrg"]
    threads:
        config["resources"]["lrg"]["cpu"]
    group:
        "assembly"
    shell:
        "cat {input} > {output}"


rule assembly_combine_canu_contigs:
    """Combine contigs from all samples plus unmapped rescue assembly"""
    input:
        expand(os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"{sample}","{sample}.contigs.uniq.fasta"), sample=samples["names"])
    output:
        os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"],"all_sample_contigs.fasta.gz")
    resources:
        **config["resources"]["lrg"]
    threads:
        config["resources"]["lrg"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "bbmap_bedtools_pigz_flye.yaml")
    container:
        config["hecatomb"]["container"]["bbmap_bedtools_pigz_flye"]
    group:
        "assembly"
    shell:
        "cat {input} | pigz -p {threads} -1 -c > {output}"
