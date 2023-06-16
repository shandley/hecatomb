rule population_assembly:
    """Assembly step 05: Create 'contig dictionary' of all unique contigs present in the study (aka: population assembly)"""
    input:
        os.path.join(dir.out.assembly, "all_sample_contigs.fasta.gz")
    output:
        assembly = os.path.join(dir.out.results, "merged_assembly.fasta"),
        graph = os.path.join(dir.out.results, "merged_assembly_graph.gfa"),
        stats = os.path.join(dir.out.assembly, "FLYE", "contig_dictionary.stats")
    params:
        flye_out = lambda w, output: os.path.split(output.stats)[0],
        flye_params = config.assembly.flye,
        assembly = os.path.join(dir.out.assembly, "FLYE", "assembly.fasta"),
        graph = os.path.join(dir.out.assembly, "FLYE", "assembly_graph.gfa")
    benchmark:
        os.path.join(dir.out.bench, "population_assembly.txt")
    log:
        log1 = os.path.join(dir.out.stderr, "population_assembly.flye.log"),
        log2 = os.path.join(dir.out.stderr, "population_assembly.stats.log")
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "metaflye.yaml")
    shell:
        """
        flye --subassemblies {input} -t {threads} --plasmids -o {params.flye_out} {params.flye_params} &>> {log.log1}
        rm {log.log1}
        mv {params.assembly} {output.assembly}
        mv {params.graph} {output.graph}
        statswrapper.sh in={output.assembly} out={output.stats} \
            format=2 \
            ow=t 2> {log.log2}
        rm {log.log2}
        """


rule koverage_samples:
    """Generate samples TSV for Koverage"""
    output:
        os.path.join(dir.out.temp, "samples_trimmed.tsv")
    params:
        names = samples.names,
        dir = dir.out.assembly,
        R12 = lambda w: ["_R1", "_R2"] if config.args.library in ["paired", "roundAB"] else ["_R1"]
    run:
        with open(output[0], 'w') as out_fh:
            for sample in params.names:
                out_fh.write(sample)
                for R in params.R12:
                    out_fh.write("\t" + os.path.join(params.dir, sample + R + ".all.fastq.gz"))
                out_fh.write("\n")


rule koverage_calculations:
    """Get coverage statistics with Koverage"""
    input:
        tsv = os.path.join(dir.out.temp, "samples_trimmed.tsv"),
        ref = os.path.join(dir.out.results, f"{config.args.assembly}_assembly.fasta"),
        req = targets.preprocessing
    output:
        os.path.join(dir.out.results, "sample_coverage.tsv"),
        os.path.join(dir.out.results, "all_coverage.tsv")
    params:
        out_dir = dir.out.base,
        minimap_mode = lambda w: "map-ont" if config.args.library == "longread" else "sr"
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    conda:
        os.path.join(dir.env, "koverage.yaml")
    shell:
        """
        koverage run \
            --reads {input.tsv} \
            --ref {input.ref} \
            --output {params.out_dir} \
            --threads {threads} \
            --minimap {params.minimap_mode}
        """
