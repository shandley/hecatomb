"""
Per-sample assemblies for short paired reads
    Take all trimmed reads, create pooled contigs "all_samples_contigs_size_selected.fasta" for use
    in 03_population_assembly.smk
"""

rule assembly_kmer_normalization:
    """Assembly step 01: Kmer normalization. Data reduction for assembly improvement"""
    input:
        r1 = os.path.join(ASSEMBLY, "{sample}_R1.unmapped.fastq.gz"),
        r1s = os.path.join(ASSEMBLY, "{sample}_R1.singletons.fastq.gz")
    output:
        r1_norm = temp(os.path.join(ASSEMBLY, "{sample}_R1.norm.fastq"))
    benchmark:
        os.path.join(BENCH, "kmer_normalization_{sample}.txt")
    log:
        os.path.join(STDERR, "kmer_norm_{sample}.log")
    resources:
        mem_mb = BBToolsMem,
        javaAlloc = int(0.9 * BBToolsMem)
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbnorm.sh in={input.r1} \
            extra={input.r1s} \
            out={output.r1_norm} \
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
        r1s = os.path.join(ASSEMBLY, "{sample}_R1.singletons.fastq.gz"),
    output:
        contigs = os.path.join(ASSEMBLY, "{sample}", "{sample}.contigs.fa"),
        renamed = os.path.join(ASSEMBLY, "{sample}", "{sample}.rename.contigs.fa"),
        tar = os.path.join(ASSEMBLY,'{sample}.tar.zst')
    params:
        mh_dir=lambda w, output: os.path.split(output.contigs)[0],
        minlen=config['CONTIG_MINLENGTH']
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
        megahit -r {input.r1_norm},{input.r1s} \
            -o {params.mh_dir} --out-prefix {wildcards.sample} --min-contig-len {params.minlen} \
            --k-min 45 --k-max 225 --k-step 26 --min-count 2 -t {threads} &>> {log}
        sed 's/>/>{wildcards.sample}/' {output.contigs} > {output.renamed}
        tar cf - {params.mh_dir} | zstd -T8 -9 > {output.tar} 2> {log}
        rm {log}
        """

rule mapSampleAssemblyPairedReads:
    """Map the sample paired reads to the sample assembly"""
    input:
        r1 = os.path.join(ASSEMBLY,"{sample}_R1.unmapped.fastq.gz"),
        contigs = os.path.join(ASSEMBLY, "{sample}", "{sample}.contigs.fa"),
    output:
        temp(os.path.join(ASSEMBLY, '{sample}', '{sample}.se.bam'))
    conda:
        os.path.join('..', 'envs', 'minimap2.yaml')
    threads:
        BBToolsCPU
    resources:
        mem_mb = BBToolsMem
    log:
        os.path.join(STDERR, 'sampleAssemblyMapSe.{sample}.log')
    benchmark:
        os.path.join(BENCH, 'sampleAssemblyMapSe.{sample}.txt')
    shell:
        """
        {{
        minimap2 -t {threads} -ax sr {input.contigs} {input.r1} | \
            samtools sort -n -o {output[0]};
        }} 2> {log}
        rm {log}
        """

rule mapSampleAssemblyUnpairedReads:
    """Map the sample unpaired reads to the sample assembly"""
    input:
        r1s = os.path.join(ASSEMBLY,"{sample}_R1.singletons.fastq.gz"),
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
        {{
        minimap2 -t {threads} -ax sr {input.contigs} {input.r1s} | \
            samtools sort -n | \
            samtools fastq -f 4 > {output[0]};
        }} 2> {log}
        rm {log}
        """

rule pullPairedUnmappedReads:
    """Grab the paired unmapped reads (neither pair mapped)"""
    input:
        os.path.join(ASSEMBLY,'{sample}','{sample}.se.bam')
    output:
        r1 = temp(os.path.join(ASSEMBLY, '{sample}', '{sample}.assemblyUnmapped_R1.fastq')),
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
        samtools fastq -f 77 {input} > {output.r1} 2> {log}
        rm {log}
        """

rule pullPairedUnmappedReadsMateMapped:
    """Grab the paired unmapped reads (mate is mapped)"""
    input:
        os.path.join(ASSEMBLY,'{sample}','{sample}.se.bam')
    output:
        temp(os.path.join(ASSEMBLY, '{sample}', '{sample}.assemblyUnmapped.se.s.fastq'))
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
        samtools fastq -f5 -F8 {input[0]} > {output[0]} 2> {log}
        rm {log}
        """

rule poolR1Unmapped:
    """Concatenate the unmapped, paired R1 reads for all samples"""
    input:
        expand(os.path.join(ASSEMBLY, '{sample}', '{sample}.assemblyUnmapped_R1.fastq'), sample=SAMPLES)
    output:
        temp(os.path.join(ASSEMBLY, 'unmapRescue_R1.fastq'))
    shell:
        """cat {input} > {output}"""

rule poolUnpairedUnmapped:
    """Concatenate the unmapped, unpaired reads for all samples"""
    input:
        fq = expand(os.path.join(ASSEMBLY,'{sample}','{sample}.assemblyUnmapped.{sSe}.fastq'), sample=SAMPLES, sSe = ['s','se.s']),
    output:
        temp(os.path.join(ASSEMBLY, 'unmapRescue.s.fastq'))
    shell:
        """
        cat {input.fq} > {output}
        """

rule rescue_read_kmer_normalization:
    """Assembly step 01: Kmer normalization. Data reduction for assembly improvement"""
    input:
        r1 = os.path.join(ASSEMBLY, 'unmapRescue_R1.fastq'),
        s = os.path.join(ASSEMBLY, 'unmapRescue.s.fastq')
    output:
        r1_norm = temp(os.path.join(ASSEMBLY, 'unmapRescueNorm_R1.fastq'))
    benchmark:
        os.path.join(BENCH, "rescue_read_kmer_normalization.txt")
    log:
        os.path.join(STDERR, "rescue_read_kmer_normalization.log")
    resources:
        mem_mb = BBToolsMem,
        javaAlloc = int(0.9 * BBToolsMem)
    threads:
        BBToolsCPU
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbnorm.sh in={input.r1} \
            extra={input.s} \
            out={output.r1_norm} \
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
        s = os.path.join(ASSEMBLY, 'unmapRescue.s.fastq')
    output:
        contigs = os.path.join(ASSEMBLY, 'rescue', "rescue.contigs.fa"),
        renamed = os.path.join(ASSEMBLY, 'rescue', 'rescue.rename.contigs.fa')
    params:
        mh_dir=lambda w, output: os.path.split(output.contigs)[0],
        minlen=config['CONTIG_MINLENGTH']
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
        megahit -r {input.r1_norm},{input.s} \
            -o {params.mh_dir} --out-prefix rescue --min-contig-len {params.minlen} \
            --k-min 45 --k-max 225 --k-step 26 --min-count 2 -t {threads} &>> {log}
        sed 's/>/>rescue/' {output.contigs} > {output.renamed}
        rm {log}
        """

rule concatenate_contigs:
    """Assembly step 03: Concatenate individual assembly outputs (contigs) into a single file"""
    input:
        expand(os.path.join(ASSEMBLY, "{sample}", "{sample}.rename.contigs.fa"), sample=SAMPLES),
        os.path.join(ASSEMBLY,'rescue',"rescue.rename.contigs.fa")
    output:
        temp(os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "all_sample_contigs.fasta"))
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


rule coverage_calculations:
    """Assembly step 07: Calculate per sample contig coverage and extract unmapped reads"""
    input:
        r1 = os.path.join(ASSEMBLY, "{sample}_R1.all.fastq.gz"),
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
        javaAlloc = int(0.9 * BBToolsMem)
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