rule population_assembly:
    """Assembly step 05: Create "contig dictionary" of all unique contigs present in the study (aka: population assembly)"""
    input:
        os.path.join(dir["out"]["assembly"], "all_sample_contigs.fasta.gz")
    output:
        assembly = os.path.join(dir["out"]["results"], "merged_assembly.fasta"),
        graph = os.path.join(dir["out"]["results"], "merged_assembly.gfa"),
        stats = os.path.join(dir["out"]["assembly"], "FLYE", "contig_dictionary.stats")
    params:
        flye_out = lambda w, output: os.path.split(output.stats)[0],
        flye_params = config["assembly"]["flye"],
        assembly = os.path.join(dir["out"]["assembly"], "FLYE", "assembly.fasta"),
        graph = os.path.join(dir["out"]["assembly"], "FLYE", "assembly_graph.gfa")
    benchmark:
        os.path.join(dir["out"]["bench"], "population_assembly.txt")
    log:
        log1 = os.path.join(dir["out"]["stderr"], "population_assembly.flye.log"),
        log2 = os.path.join(dir["out"]["stderr"], "population_assembly.stats.log")
    resources:
        mem_mb = resources["lrg"]["mem"],
        mem = str(resources["lrg"]["mem"]) + "MB",
        time = resources["lrg"]["time"]
    threads:
        resources["lrg"]["cpu"]
    conda:
        os.path.join(dir["env"], "metaflye.yaml")
    shell:
        " flye --subassemblies {input} "
            "-t {threads} --plasmids -o {params.flye_out} {params.flye_params} "
            "&>> {log.log1}; "
        "mv {params.assembly} {output.assembly}; "
        "mv {params.graph} {output.graph}; "
        "statswrapper.sh in={output.assembly} out={output.stats} "
            "format=2 "
            "ow=t 2> {log.log2}; "


rule koverage_samples:
    """Generate samples TSV for Koverage"""
    input:
        targets["trimnami"]
    output:
        temp(os.path.join(dir["out"]["temp"], "samples_trimmed.tsv"))
    params:
        samples["trimmed"]
    run:
        from metasnek import fastq_finder
        fastq_finder.write_samples_tsv(params[0],output[0])


rule koverage_calculations:
    """Get coverage statistics with Koverage"""
    input:
        tsv = os.path.join(dir["out"]["temp"], "samples_trimmed.tsv"),
        ref = os.path.join(dir["out"]["results"], config["args"]["assembly"] + '_assembly.fasta'),
        req = targets["preprocessing"]
    output:
        os.path.join(dir["out"]["results"], "sample_coverage.tsv"),
        os.path.join(dir["out"]["results"], "all_coverage.tsv")
    params:
        out_dir = dir["out"]["base"],
        minimap_mode = lambda wildcards: "map-ont" if config["args"]["trim"] == "nanopore" else "sr",
        profile= lambda wildcards: "--profile " + config["args"]["profile"] if config["args"]["profile"] else "",
    threads:
        lambda wildcards: resources["sml"]["cpu"] if config["args"]["profile"] else resources["big"]["cpu"]
    resources:
        mem_mb = lambda wildcards: resources["sml"]["mem"] if config["args"]["profile"] else resources["big"]["mem"],
        mem = lambda wildcards: str(resources["sml"]["mem"]) + "MB" if config["args"]["profile"] else str(resources["big"]["mem"]) + "MB",
        time = resources["big"]["time"]
    conda:
        os.path.join(dir["env"], "koverage.yaml")
    shell:
        "koverage run "
            "--reads {input.tsv} "
            "--ref {input.ref} "
            "--output {params.out_dir} "
            "--threads {threads} "
            "--minimap {params.minimap_mode} "
            "{params.profile}; "
