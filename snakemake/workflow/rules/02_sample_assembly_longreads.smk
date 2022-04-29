"""
Per-sample assemblies for longreads
    Take all host-removed reads, create pooled contigs "all_samples_contigs_size_selected.fasta" for use
    in 03_population_assembly.smk
"""

rule canu_sample_assembly:
    """Per-sample assembly with canu; also works for unmapped rescue reads"""
    input:
        os.path.join(TMPDIR,"p01","{sample}_R1.all.fasta")
    output:
        ctg = os.path.join(ASSEMBLY,"{sample}","{sample}.contigs.fasta"),
        ctgq = os.path.join(ASSEMBLY,"{sample}","{sample}.contigs.uniq.fasta"),
        un = os.path.join(ASSEMBLY,"{sample}","{sample}.unassembled.fasta"),
        unq = os.path.join(ASSEMBLY,"{sample}","{sample}.unassembled.uniq.fasta"),
        tar = os.path.join(ASSEMBLY,'{sample}.tar.zst')
    params:
        settings = config['canuSettings'],
        canu_dir = os.path.join(ASSEMBLY,"{sample}")
    resources:
        mem_mb = MhitMem
    threads:
        MhitCPU
    log:
        os.path.join(STDERR, "canu_sample_assembly.{sample}.log")
    conda:
        "../envs/canu.yaml"
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
        expand(os.path.join(ASSEMBLY,"{sample}","{sample}.unassembled.uniq.fasta"), sample=SAMPLES)
    output:
        temp(os.path.join(TMPDIR,"p01","unmappedRescue_R1.all.fasta"))
    shell:
        """cat {input} > {output}"""

rule combine_canu_contigs:
    """Combine contigs from all samples plus unmapped rescue assembly"""
    input:
        expand(os.path.join(ASSEMBLY,"{sample}","{sample}.contigs.uniq.fasta"), sample=list(SAMPLES) + ['unmappedRescue'])
    output:
        temp(os.path.join(ASSEMBLY,"CONTIG_DICTIONARY","all_samples_contigs_size_selected.fasta"))
    shell:
        """cat {input} > {output}"""

rule coverage_calculations:
    """Assembly step 07: Calculate per sample contig coverage and extract unmapped reads"""
    input:
        r1 = os.path.join(TMPDIR,"p01","{sample}_R1.all.fasta"),
        ref = os.path.join(RESULTS, "assembly.fasta")
    output:
        sam = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.aln.sam.gz")),
        unmap = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.unmapped.fastq")),
        covstats = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.cov_stats")),
        rpkm = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.rpkm")),
        statsfile = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.statsfile")),
        scafstats = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "MAPPING", "{sample}.scafstats"))
    benchmark:
        os.path.join(BENCH, "coverage_calculations.{sample}.txt")
    log:
        os.path.join(STDERR, "coverage_calculations.{sample}.log")
    resources:
        mem_mb = BBToolsMem,
        javaAlloc = int(0.95 * BBToolsMem)
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
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