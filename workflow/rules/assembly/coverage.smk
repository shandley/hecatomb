rule assembly_population_subassembly:
    """Assembly step 05: Create "contig dictionary" of all unique contigs present in the study (aka: population assembly)"""
    input:
        os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"], "all_sample_contigs.fasta.gz")
    output:
        assembly = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "merged_assembly.fasta"),
        graph = os.path.join(config["hecatomb"]["args"]["output_paths"]["results"], "merged_assembly.gfa"),
        stats = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"], "FLYE", "contig_dictionary.stats")
    params:
        flye_out = lambda w, output: os.path.split(output.stats)[0],
        flye_params = config["hecatomb"]["assembly"]["flye"],
        assembly = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"], "FLYE", "assembly.fasta"),
        graph = os.path.join(config["hecatomb"]["args"]["temp_paths"]["assembly"], "FLYE", "assembly_graph.gfa")
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "population_assembly.txt")
    log:
        log1 = os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "population_assembly.flye.log"),
        log2 = os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "population_assembly.stats.log")
    resources:
        **config["resources"]["lrg"]
    threads:
        config["resources"]["lrg"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "bbmap_bedtools_pigz_flye.yaml")
    container:
        config["hecatomb"]["container"]["bbmap_bedtools_pigz_flye"]
    shell:
        "flye --subassemblies {input} "
            "-t {threads} "
            "--plasmids "
            "-o {params.flye_out} "
            "{params.flye_params} "
            "&> {log.log1}; "
        "mv {params.assembly} {output.assembly}; "
        "mv {params.graph} {output.graph}; "
        "statswrapper.sh in={output.assembly} "
            "out={output.stats} "
            "format=2 "
            "ow=t &> {log.log2}; "

