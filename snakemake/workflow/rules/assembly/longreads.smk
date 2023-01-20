"""
Per-sample assemblies for longreads
    Take all host-removed reads, create pooled contigs "all_samples_contigs_size_selected.fasta" for use
    in combine_sample_assemblies.smk
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
        tar = os.path.join(dir.out.assembly,'{sample}.tar.zst')
    params:
        settings = config.assembly.canu,
        canu_dir = lambda w, output: os.path.split(output.ctg)[0]
    resources:
        mem_mb = config.resources.med.mem
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
        tar cf - {params.canu_dir} | zstd -T8 -9 > {output.tar} 2> {log}
        rm {log}
        """


rule combine_canu_unassembled:
    """Combine the unassembled reads from all canu assemblies"""
    input:
        expand(os.path.join(dir.out.assembly,"{sample}","{sample}.unassembled.uniq.fasta"), sample=samples.names)
    output:
        temp(os.path.join(dir.out.assembly,"unmappedRescue_R1.all.fasta.gz"))
    shell:
        """cat {input} > {output}"""


rule combine_canu_contigs:
    """Combine contigs from all samples plus unmapped rescue assembly"""
    input:
        expand(os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.uniq.fasta"), sample=samples.names + ['unmappedRescue'])
    output:
        temp(os.path.join(dir.out.assembly,"all_sample_contigs.fasta"))
    shell:
        """cat {input} > {output}"""


rule coverage_calculations:
    """Assembly step 07: Calculate per sample contig coverage and extract unmapped reads"""
    input:
        r1 = os.path.join(dir.out.assembly,"{sample}_R1.all.fasta.gz"),
        ref = os.path.join(dir.out.results, "assembly.fasta")
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
        javaAlloc = int(0.9 * config.resources.med.mem)
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "bbmap.yaml")
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