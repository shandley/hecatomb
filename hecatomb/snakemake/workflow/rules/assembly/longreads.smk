
rule lr_cross_assembly:
    """Alternative to merged assembly; assemble everything together in one hit."""
    input:
        targets["trimnami"]
    output:
        assembly = os.path.join(dir["out"]["results"], "cross_assembly.fasta"),
        graph = os.path.join(dir["out"]["results"], "cross_assembly.gfa"),
        tar = os.path.join(dir["out"]["assembly"],"crossAssembly.tar.zst")
    params:
        params = config["assembly"]["metaflye"],
        dir = os.path.join(dir["out"]["assembly"], "crossAssembly"),
        assembly = os.path.join(dir["out"]["assembly"], "crossAssembly", "assembly.fasta"),
        graph = os.path.join(dir["out"]["assembly"], "crossAssembly", "assembly_graph.gfa"),
    resources:
        mem_mb = resources["big"]["mem"],
        mem = str(resources["big"]["mem"]) + "MB",
        time = resources["big"]["time"]
    threads:
        resources["big"]["cpu"]
    log:
        os.path.join(dir["out"]["stderr"], "canu_cross_assembly.log")
    benchmark:
        os.path.join(dir["out"]["bench"], "canu_cross_assembly.txt")
    conda:
        os.path.join(dir["env"], "metaflye.yaml")
    shell:
        "flye -o {params.dir} -t {threads} {params.params} {input} 2> {log} "
        "mv {params.assembly} {output.assembly} "
        "mv {params.graph} {output.graph} "
        "tar cf - {params.dir} | zstd -T{threads} -9 > {output.tar} 2> {log} "


rule canu_sample_assembly:
    """Per-sample assembly with canu; also works for unmapped rescue reads"""
    input:
        os.path.join(dir["out"]["trim"], "{sample}" + "_S" + config["args"]["hostStr"] + ".fastq.gz")
    output:
        ctg = os.path.join(dir["out"]["assembly"],"{sample}","{sample}.contigs.fasta"),
        ctgq = os.path.join(dir["out"]["assembly"],"{sample}","{sample}.contigs.uniq.fasta"),
        un = os.path.join(dir["out"]["assembly"],"{sample}","{sample}.unassembled.fasta"),
        unq = os.path.join(dir["out"]["assembly"],"{sample}","{sample}.unassembled.uniq.fasta"),
        tar = os.path.join(dir["out"]["assembly"],"{sample}.tar.zst")
    params:
        settings = config["assembly"]["canu"],
        canu_dir = lambda w, output: os.path.split(output.ctg)[0]
    resources:
        mem_mb = resources["lrg"]["mem"],
        mem = str(resources["lrg"]["mem"]) + "MB",
        time = resources["lrg"]["time"]
    threads:
        resources["lrg"]["cpu"]
    log:
        os.path.join(dir["out"]["stderr"], "canu_sample_assembly.{sample}.log")
    conda:
        os.path.join(dir["env"], "canu.yaml")
    shell:
        "canu {params.settings} {input} "
            "batThreads={threads} "
            "batMemory={resources.mem_mb}M "
            "-p {wildcards.sample} "
            "-d {params.canu_dir} "
            "&>> {log} "
        "sed 's/>tig/>{wildcards.sample}./' {output.ctg} > {output.ctgq} "
        "sed 's/>tig/>{wildcards.sample}./' {output.un} > {output.unq} "
        "tar cf - {params.canu_dir} | zstd -T{threads} -9 > {output.tar} 2> {log} "


rule combine_canu_unassembled:
    """Combine the unassembled reads from all canu assemblies"""
    input:
        expand(os.path.join(dir["out"]["assembly"],"{sample}","{sample}.unassembled.uniq.fasta"), sample=samples["names"])
    output:
        temp(os.path.join(dir["out"]["assembly"],"unmappedRescue_R1.all.fasta.gz"))
    threads:
        resources["lrg"]["cpu"]
    resources:
        time = resources["sml"]["time"]
    group:
        "assembly"
    shell:
        "cat {input} > {output}"


rule combine_canu_contigs:
    """Combine contigs from all samples plus unmapped rescue assembly"""
    input:
        expand(os.path.join(dir["out"]["assembly"],"{sample}","{sample}.contigs.uniq.fasta"), sample=samples["names"])
    output:
        os.path.join(dir["out"]["assembly"],"all_sample_contigs.fasta.gz")
    threads:
        resources["lrg"]["cpu"]
    resources:
        time = resources["sml"]["time"]
    conda:
        os.path.join(dir["env"], "pigz.yaml")
    group:
        "assembly"
    shell:
        "cat {input} | pigz -p {threads} -1 -c > {output}"
