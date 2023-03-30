"""
Per-sample assemblies for longreads
    Take all host-removed reads, create pooled contigs "all_samples_contigs_size_selected.fasta" for use
    in combine_sample_assemblies.smk
"""

rule lr_co_assembly:
    """Alternative to cross assembly; assemble everything together in one hit."""
    input:
        expand(os.path.join(dir.out.assembly,"{sample}_R1.all.fasta.gz"), sample=samples.names)
    output:
        assembly = os.path.join(dir.out.results, "co_assembly.fasta"),
        graph = os.path.join(dir.out.results, "co_assembly_graph.gfa"),
        tar = os.path.join(dir.out.assembly,"coAssembly.tar.zst")
    params:
        params = config.assembly.metaflye,
        dir = os.path.join(dir.out.assembly, "coAssembly"),
        assembly = os.path.join(dir.out.assembly, "coAssembly", "assembly.fasta"),
        graph = os.path.join(dir.out.assembly, "coAssembly", "assembly_graph.gfa"),
    resources:
        mem_mb = config.resources.big.mem,
        time = config.resources.big.time
    threads:
        config.resources.big.cpu
    log:
        os.path.join(dir.out.stderr, "canu_co_assembly.log")
    benchmark:
        os.path.join(dir.out.bench, "canu_co_assembly.txt")
    conda:
        os.path.join(dir.env, "metaflye.yaml")
    shell:
        """
        flye -o {params.dir} -t {threads} {params.params} {input} 2> {log}
        mv {params.assembly} {output.assembly}
        mv {params.graph} {output.graph}
        tar cf - {params.dir} | zstd -T{threads} -9 > {output.tar} 2> {log}
        rm {log}
        """


rule canu_sample_assembly:
    """Per-sample assembly with canu; also works for unmapped rescue reads"""
    input:
        os.path.join(dir.out.assembly,"{sample}_R1.all.fasta.gz")
    output:
        ctg = os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.fasta"),
        ctgq = os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.uniq.fasta"),
        un = os.path.join(dir.out.assembly,"{sample}","{sample}.unassembled.fasta"),
        unq = os.path.join(dir.out.assembly,"{sample}","{sample}.unassembled.uniq.fasta"),
        tar = os.path.join(dir.out.assembly,"{sample}.tar.zst")
    params:
        settings = config.assembly.canu,
        canu_dir = lambda w, output: os.path.split(output.ctg)[0]
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    log:
        os.path.join(dir.out.stderr, "canu_sample_assembly.{sample}.log")
    conda:
        os.path.join(dir.env, "canu.yaml")
    shell:
        """
        canu {params.settings} {input} \
            batThreads={threads} \
            batMemory={resources.mem_mb}M \
            -p {wildcards.sample} \
            -d {params.canu_dir} \
            &>> {log}
        sed 's/>tig/>{wildcards.sample}./' {output.ctg} > {output.ctgq}
        sed 's/>tig/>{wildcards.sample}./' {output.un} > {output.unq}
        tar cf - {params.canu_dir} | zstd -T{threads} -9 > {output.tar} 2> {log}
        rm {log}
        """


rule combine_canu_unassembled:
    """Combine the unassembled reads from all canu assemblies"""
    input:
        expand(os.path.join(dir.out.assembly,"{sample}","{sample}.unassembled.uniq.fasta"), sample=samples.names)
    output:
        temp(os.path.join(dir.out.assembly,"unmappedRescue_R1.all.fasta.gz"))
    resources:
        time = config.resources.sml.time
    group:
        "assembly"
    shell:
        """cat {input} > {output}"""


rule combine_canu_contigs:
    """Combine contigs from all samples plus unmapped rescue assembly"""
    input:
        expand(os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.uniq.fasta"), sample=samples.names + ["unmappedRescue"])
    output:
        os.path.join(dir.out.assembly,"all_sample_contigs.fasta.gz")
    threads:
        config.resources.med.cpu
    resources:
        time = config.resources.sml.time
    params:
        compression = '-' + str(config.qc.compression)
    conda:
        os.path.join(dir.env, "pigz.yaml")
    group:
        "assembly"
    shell:
        """cat {input} | pigz -p {threads} {params.compression} -c > {output}"""


rule coverage_calculations:
    """Assembly step 07: Calculate per sample contig coverage and extract unmapped reads"""
    input:
        r1 = os.path.join(dir.out.assembly,"{sample}_R1.all.fasta.gz"),
        ref = os.path.join(dir.out.results, f"{config.args.assembly}_assembly.fasta")
    output:
        sam = temp(os.path.join(dir.out.mapping, "{sample}.aln.sam.gz")),
        unmap = temp(os.path.join(dir.out.mapping, "{sample}.unmapped.fastq")),
        covstats = temp(os.path.join(dir.out.mapping, "{sample}.cov_stats")),
        rpkm = temp(os.path.join(dir.out.mapping, "{sample}.rpkm")),
        statsfile = temp(os.path.join(dir.out.mapping, "{sample}.statsfile")),
        scafstats = temp(os.path.join(dir.out.mapping, "{sample}.scafstats"))
    benchmark:
        os.path.join(dir.out.bench, "coverage_calculations.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "coverage_calculations.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        javaAlloc = int(0.9 * config.resources.med.mem),
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
    group:
        "assembly"
    shell:
        """
        bbmap.sh ref={input.ref} in={input.r1} \
            nodisk \
            out={output.sam} \
            outu={output.unmap} \
            ambiguous=random \
            slow=t \
            physcov=t \
            covstats={output.covstats} \
            rpkm={output.rpkm} \
            statsfile={output.statsfile} \
            scafstats={output.scafstats} \
            maxindel=100 minid=90 \
            ow=t \
            threads={threads} -Xmx{resources.javaAlloc}m 2> {log}
        rm {log}
        """