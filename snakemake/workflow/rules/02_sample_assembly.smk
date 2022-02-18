"""
Per-sample assemblies for short paired reads
    Take all trimmed reads, create pooled contigs "all_samples_contigs_size_selected.fasta" for use
    in 03_population_assembly.smk
"""

rule assembly_kmer_normalization:
    """Assembly step 01: Kmer normalization. Data reduction for assembly improvement"""
    input:
        r1 = os.path.join(TMPDIR, "p07", "{sample}_R1.unmapped.fastq"),
        r2 = os.path.join(TMPDIR, "p07", "{sample}_R2.unmapped.fastq"),
        r1s = os.path.join(TMPDIR, "p07", "{sample}_R1.singletons.fastq"),
        r2s = os.path.join(TMPDIR, "p07", "{sample}_R2.singletons.fastq")
    output:
        r1_norm = os.path.join(ASSEMBLY, "{sample}_R1.norm.fastq"),
        r2_norm = os.path.join(ASSEMBLY, "{sample}_R2.norm.fastq")
    benchmark:
        os.path.join(BENCH, "kmer_normalization_{sample}.txt")
    log:
        os.path.join(STDERR, "kmer_norm_{sample}.log")
    resources:
        mem_mb = BBToolsMem,
        javaAlloc = int(0.95 * BBToolsMem)
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbnorm.sh in={input.r1} in2={input.r2} \
            extra={input.r1s},{input.r2s} \
            out={output.r1_norm} out2={output.r2_norm} \
            target=100 \
            ow=t \
            threads={threads} -Xmx{resources.javaAlloc}m 2> {log}
        rm {log}
        """

rule individual_sample_assembly:
    """Assembly step 02: Individual sample assemblies

    Megahit: https://github.com/voutcn/megahit
    """
    input:
        r1_norm = os.path.join(ASSEMBLY, "{sample}_R1.norm.fastq"),
        r2_norm = os.path.join(ASSEMBLY, "{sample}_R2.norm.fastq"),
        r1s = os.path.join(TMPDIR, "p07", "{sample}_R1.singletons.fastq"),
        r2s = os.path.join(TMPDIR, "p07", "{sample}_R2.singletons.fastq")
    output:
        contigs = os.path.join(ASSEMBLY, "{sample}", "{sample}.contigs.fa"),
        tar = os.path.join(ASSEMBLY,'{sample}.tar.zst')
    params:
        mh_dir = directory(os.path.join(ASSEMBLY,'{sample}'))
    benchmark:
        os.path.join(BENCH, "megahit_{sample}.txt")
    log:
        os.path.join(STDERR, "megahit_{sample}.log")
    resources:
        mem_mb = MhitMem
    threads:
        MhitCPU
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        if [ -d {params.mh_dir} ]; then
            rm -rf {params.mh_dir}
        fi
        megahit -1 {input.r1_norm} -2 {input.r2_norm} -r {input.r1s},{input.r2s} \
            -o {params.mh_dir} --out-prefix {wildcards.sample} \
            --k-min 45 --k-max 225 --k-step 26 --min-count 2 -t {threads} &>> {log}
        tar cf - {params.mh_dir} | zstd -T8 -9 > {output.tar} 2> {log}
        rm {log}
        """

rule mapSampleAssemblyPairedReads:
    """Map the sample paired reads to the sample assembly"""
    input:
        r1 = os.path.join(TMPDIR,"p07","{sample}_R1.unmapped.fastq"),
        r2 = os.path.join(TMPDIR,"p07","{sample}_R2.unmapped.fastq"),
        contigs = os.path.join(ASSEMBLY, "{sample}", "{sample}.contigs.fa"),
    output:
        temp(os.path.join(ASSEMBLY, '{sample}', '{sample}.pe.bam'))
    conda:
        os.path.join('..', 'envs', 'minimap2.yaml')
    threads:
        BBToolsCPU
    resources:
        mem_mb = BBToolsMem
    log:
        os.path.join(STDERR, 'sampleAssemblyMapPe.{sample}.log')
    benchmark:
        os.path.join(BENCH, 'sampleAssemblyMapPe.{sample}.txt')
    shell:
        """
        minimap2 -t {threads} -ax sr {input.contigs} {input.r1} {input.r2} | \
            samtools sort -n -o {output[0]}
        """

rule mapSampleAssemblyUnpairedReads:
    """Map the sample unpaired reads to the sample assembly"""
    input:
        r1s = os.path.join(TMPDIR,"p07","{sample}_R1.singletons.fastq"),
        r2s = os.path.join(TMPDIR,"p07","{sample}_R2.singletons.fastq"),
        contigs = os.path.join(ASSEMBLY,"{sample}","{sample}.contigs.fa")
    output:
        temp(os.path.join(ASSEMBLY,'{sample}','{sample}.assemblyUnmapped.s.fastq'))
    conda:
        os.path.join('..', 'envs', 'minimap2.yaml')
    threads:
        BBToolsCPU
    resources:
        mem_mb = BBToolsMem
    log:
        os.path.join(STDERR, 'sampleAssemblyMapS.{sample}.log')
    benchmark:
        os.path.join(BENCH, 'sampleAssemblyMapS.{sample}.txt')
    shell:
        """
        minimap2 -t {threads} -ax sr {input.contigs} {input.r1s} {input.r2s} | \
            samtools sort -n | \
            samtools fastq -f 4 > {output[0]}
        """

rule pullPairedUnmappedReads:
    """Grab the paired unmapped reads (neither pair mapped)"""
    input:
        os.path.join(ASSEMBLY,'{sample}','{sample}.pe.bam')
    output:
        r1 = temp(os.path.join(ASSEMBLY, '{sample}', '{sample}.assemblyUnmapped_R1.fastq')),
        r2 = temp(os.path.join(ASSEMBLY, '{sample}', '{sample}.assemblyUnmapped_R2.fastq')),
    conda:
        os.path.join('..','envs','samtools.yaml')
    threads:
        BBToolsCPU
    resources:
        mem_mb = BBToolsMem
    log:
        os.path.join(STDERR, 'pullPairedUnmappedReads.{sample}.log')
    benchmark:
        os.path.join(BENCH, 'pullPairedUnmappedReads.{sample}.txt')
    shell:
        """
        samtools fastq -f 77 {input} > {output.r1}
        samtools fastq -f 141 {input} > {output.r2}
        """

rule pullPairedUnmappedReadsMateMapped:
    """Grab the paired unmapped reads (mate is mapped)"""
    input:
        os.path.join(ASSEMBLY,'{sample}','{sample}.pe.bam')
    output:
        temp(os.path.join(ASSEMBLY, '{sample}', '{sample}.assemblyUnmapped.pe.s.fastq'))
    conda:
        os.path.join('..','envs','samtools.yaml')
    threads:
        BBToolsCPU
    resources:
        mem_mb = BBToolsMem
    log:
        os.path.join(STDERR, 'pullPairedUnmappedReads.{sample}.log')
    benchmark:
        os.path.join(BENCH, 'pullPairedUnmappedReads.{sample}.txt')
    shell:
        """
        samtools fastq -f5 -F8 {input[0]} > {output[0]}
        """

rule poolR1Unmapped:
    """Concatenate the unmapped, paired R1 reads for all samples"""
    input:
        expand(os.path.join(ASSEMBLY, '{sample}', '{sample}.assemblyUnmapped_R1.fastq'), sample=SAMPLES)
    output:
        temp(os.path.join(ASSEMBLY, 'unmapRescue_R1.fastq'))
    shell:
        """cat {input} > {output}"""

rule poolR2Unmapped:
    """Concatenate the unmapped, paired R2 reads for all samples"""
    input:
        expand(os.path.join(ASSEMBLY, '{sample}', '{sample}.assemblyUnmapped_R2.fastq'), sample=SAMPLES)
    output:
        temp(os.path.join(ASSEMBLY, 'unmapRescue_R2.fastq'))
    shell:
        """cat {input} > {output}"""

rule poolUnpairedUnmapped:
    """Concatenate the unmapped, unpaired reads for all samples"""
    input:
        fq = expand(os.path.join(ASSEMBLY,'{sample}','{sample}.assemblyUnmapped.{sPe}.fastq'), sample=SAMPLES, sPe = ['s','pe.s']),
    output:
        temp(os.path.join(ASSEMBLY, 'unmapRescue.s.fastq'))
    shell:
        """
        cat {input.fq} > {output}
        """

# rule archive_mhitDir:
#     """tar and zip the megahit assembly directories.
#
#     We add the count table as a requirement to make sure the directory is only archive once it's no longer needed.
#     """
#     input:
#         dir = os.path.join(ASSEMBLY,'{sample}'),
#         req = os.path.join(RESULTS, "contig_count_table.tsv")
#     output:
#         os.path.join(ASSEMBLY,'{sample}.tar.zst')
#     threads:
#         BBToolsCPU
#     resources:
#         mem_mb = BBToolsMem
#     shell:
#         """
#         tar cf - {input.dir} | zstd -T8 -9 > {output}
#         rm -rf {input.dir}
#         """

rule rescue_read_kmer_normalization:
    """Assembly step 01: Kmer normalization. Data reduction for assembly improvement"""
    input:
        r1 = os.path.join(ASSEMBLY, 'unmapRescue_R1.fastq'),
        r2 = os.path.join(ASSEMBLY, 'unmapRescue_R2.fastq'),
        s = os.path.join(ASSEMBLY, 'unmapRescue.s.fastq')
    output:
        r1_norm = temp(os.path.join(ASSEMBLY, 'unmapRescueNorm_R1.fastq')),
        r2_norm = temp(os.path.join(ASSEMBLY, 'unmapRescueNorm_R2.fastq'))
    benchmark:
        os.path.join(BENCH, "rescue_read_kmer_normalization.txt")
    log:
        os.path.join(STDERR, "rescue_read_kmer_normalization.log")
    resources:
        mem_mb = BBToolsMem,
        javaAlloc = int(0.95 * BBToolsMem)
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbnorm.sh in={input.r1} in2={input.r2} \
            extra={input.s} \
            out={output.r1_norm} out2={output.r2_norm} \
            target=100 \
            ow=t \
            threads={threads} -Xmx{resources.javaAlloc}m 2> {log}
        rm {log}
        """

rule unmapped_read_rescue_assembly:
    """Assemble the unmapped reads from all samples

    Megahit: https://github.com/voutcn/megahit
    """
    input:
        r1_norm = os.path.join(ASSEMBLY, 'unmapRescueNorm_R1.fastq'),
        r2_norm = os.path.join(ASSEMBLY, 'unmapRescueNorm_R2.fastq'),
        s = os.path.join(ASSEMBLY, 'unmapRescue.s.fastq')
    output:
        contigs = os.path.join(ASSEMBLY, 'rescue', "rescue.contigs.fa")
    params:
        mh_dir = directory(os.path.join(ASSEMBLY, 'rescue'))
    benchmark:
        os.path.join(BENCH, "unmapped_read_rescue_assembly.txt")
    log:
        os.path.join(STDERR, "unmapped_read_rescue_assembly.log")
    resources:
        mem_mb = MhitMem
    threads:
        MhitCPU
    conda:
        "../envs/megahit.yaml"
    shell:
        """
        if [ -d {params.mh_dir} ]; then
            rm -rf {params.mh_dir}
        fi
        megahit -1 {input.r1_norm} -2 {input.r2_norm} -r {input.s} \
            -o {params.mh_dir} --out-prefix rescue \
            --k-min 45 --k-max 225 --k-step 26 --min-count 2 -t {threads} &>> {log}
        rm {log}
        """

rule concatenate_contigs:
    """Assembly step 03: Concatenate individual assembly outputs (contigs) into a single file"""
    input:
        expand(os.path.join(ASSEMBLY, "{sample}", "{sample}.contigs.fa"), sample=SAMPLES),
        os.path.join(ASSEMBLY,'rescue',"rescue.contigs.fa")
    output:
        temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.fasta"))
    params:
        expand(os.path.join(ASSEMBLY,'{sample}'), sample=SAMPLES),
    benchmark:
        os.path.join(BENCH, "concatenate_assemblies.txt")
    log:
        os.path.join(STDERR, "concatenate_assemblies.log")
    shell:
        """
        cat {input} > {output} 2> {log} && rm {log}
        rm -rf {params}
        """

rule contig_reformating_and_stats:
    """Assembly step 04: Remove short contigs (Default: 1000). Defined in config[CONTIG_MINLENGTH]"""
    input:
        os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.fasta")
    output:
        rename = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.renamed.fasta")),
        size = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_samples_contigs_size_selected.fasta")),
        stats = temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_megahit_contigs.stats"))
    benchmark:
        os.path.join(BENCH, "contig_reformating.txt")
    log:
        log1 = os.path.join(STDERR, "contig_reformating_and_stats.rename.log"),
        log2 = os.path.join(STDERR, "contig_reformating_and_stats.reformat.log"),
        log3 = os.path.join(STDERR, "contig_reformating_and_stats.stats.log")
    resources:
        mem_mb = BBToolsMem,
        javaAlloc = int(0.95 * BBToolsMem)
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        rename.sh in={input} out={output.rename} \
            prefix=contig_ \
            ow=t \
            -Xmx{resources.javaAlloc}m 2> {log.log1}
        rm {log.log1}
        reformat.sh in={output.rename} out={output.size} \
            ml={config[CONTIG_MINLENGTH]} \
            ow=t \
            -Xmx{resources.javaAlloc}m 2> {log.log2}
        rm {log.log2}
        statswrapper.sh in={input} out={output.stats} \
            format=2 \
            ow=t 2> {log.log3}
        rm {log.log3}
        """


rule coverage_calculations:
    """Assembly step 07: Calculate per sample contig coverage and extract unmapped reads"""
    input:
        r1 = os.path.join(TMPDIR, "p07", "{sample}_R1.all.fastq"),
        r2 = os.path.join(TMPDIR, "p07", "{sample}_R2.all.fastq"),
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
        bbmap.sh ref={input.ref} in={input.r1} in2={input.r2} \
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