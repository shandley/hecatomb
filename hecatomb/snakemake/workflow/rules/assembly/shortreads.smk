rule cross_assembly:
    """Alternative to merged assembly; assemble everything together in one hit.

    Megahit: https://github.com/voutcn/megahit
    """
    input:
        targets.trimnami
    output:
        assembly=os.path.join(dir.out.results,"cross_assembly.fasta"),
        graph=os.path.join(dir.out.results,"cross_assembly_graph.gfa"),
        tmp=temp(os.path.join(dir.out.assembly,"crossAssembly","cross_assembly_graph.fastg")),
        tar=os.path.join(dir.out.assembly,"crossAssembly.tar.zst")
    params:
        r1p=targets.cross.r1,
        r2p=targets.cross.r2,
        rs=targets.cross.s,
        mh_dir=os.path.join(dir.out.assembly,"crossAssembly"),
        mh_int=os.path.join(dir.out.assembly,"crossAssembly","intermediate_contigs"),
        params=config.assembly.megahit,
        assembly=os.path.join(dir.out.assembly,"crossAssembly","final.contigs.fa"),
        graph=os.path.join(dir.out.assembly,"crossAssembly","assembly_graph.gfa"),
    benchmark:
        os.path.join(dir.out.bench,"cross_assembly.txt")
    log:
        os.path.join(dir.out.stderr,"cross_assembly.log")
    resources:
        mem_mb=config.resources.big.mem,
        time=config.resources.big.time
    threads:
        config.resources.big.cpu
    conda:
        os.path.join(dir.env,"megahit.yaml")
    shell:
        """
        if [ -d {params.mh_dir} ]
        then
            rm -rf {params.mh_dir}
        fi

        megahit \
            {params.r1p} \
            {params.r2p} \
            {params.rs} \
            -o {params.mh_dir} \
            -t {threads} \
            {params.params} \
            &> {log}

        kctg=$(ls -t {params.mh_int}/*.contigs.fa | grep -v final | head -1)

        kmax=$(head -1 $kctg | sed 's/>k\|_.*//g')

        megahit_toolkit \
            contig2fastg \
            $kmax \
            $kctg \
            > {output.tmp} \
            2> {log}

        Bandage reduce {output.tmp} {output.graph} 2> {log}

        cp {params.assembly} {output.assembly}

        tar cf - {params.mh_dir} \
            | zstd -T{threads} -9 \
            > {output.tar} \
            2> {log}
        rm {log}
        """


rule megahit_sample_paired:
    input:
        r1 = os.path.join(dir.out.trim, "{sample}" + config.args.hostStr + ".paired.R1.fastq.gz"),
        r2 = os.path.join(dir.out.trim, "{sample}" + config.args.hostStr + ".paired.R2.fastq.gz"),
        s  = os.path.join(dir.out.trim, "{sample}" + config.args.hostStr + ".paired.S.fastq.gz"),
    output:
        contigs = os.path.join(dir.out.assembly, "{sample}", "{sample}.contigs.fa"),
        renamed = os.path.join(dir.out.assembly, "{sample}", "{sample}.rename.contigs.fa"),
        tar = os.path.join(dir.out.assembly,"{sample}.tar.zst")
    params:
        mh_dir = lambda w, output: os.path.split(output.contigs)[0],
        params = config.assembly.megahit
    benchmark:
        os.path.join(dir.out.bench, "megahit_sample_paired.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "megahit_sample_paired.{sample}.log")
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env, "megahit.yaml")
    group:
        "assembly"
    shell:
        """
        if [ -d {params.mh_dir} ]
        then
            rm -rf {params.mh_dir}
        fi
        
        megahit \
            -1 {input.r1} \
            -2 {input.r2} \
            -r {input.s} \
            -o {params.mh_dir} \
            --out-prefix {wildcards.sample} \
            -t {threads} \
            {params.params} \
            &> {log}
        
        sed 's/>/>{wildcards.sample}/' {output.contigs} > {output.renamed}
        
        tar cf - {params.mh_dir} | zstd -T{threads} -9 > {output.tar} 2> {log}
        
        rm {log}
        """


rule megahit_sample_unpaired:
    input:
        os.path.join(dir.out.trim, "{sample}" + config.args.hostStr + ".single.fastq.gz"),
    output:
        contigs=os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.fa"),
        renamed=os.path.join(dir.out.assembly,"{sample}","{sample}.rename.contigs.fa"),
        tar=os.path.join(dir.out.assembly,"{sample}.tar.zst")
    params:
        mh_dir=lambda w, output: os.path.split(output.contigs)[0],
        params=config.assembly.megahit
    benchmark:
        os.path.join(dir.out.bench,"megahit_sample_unpaired.{sample}.txt")
    log:
        os.path.join(dir.out.stderr,"megahit_sample_unpaired.{sample}.log")
    resources:
        mem_mb=config.resources.med.mem,
        time=config.resources.med.time
    threads:
        config.resources.med.cpu
    conda:
        os.path.join(dir.env,"megahit.yaml")
    group:
        "assembly"
    shell:
        """
        if [ -d {params.mh_dir} ]
        then
            rm -rf {params.mh_dir}
        fi

        megahit \
            -r {input} \
            -o {params.mh_dir} \
            --out-prefix {wildcards.sample} \
            -t {threads} \
            {params.params} \
            &> {log}

        sed 's/>/>{wildcards.sample}/' {output.contigs} > {output.renamed}

        tar cf - {params.mh_dir} | zstd -T{threads} -9 > {output.tar} 2> {log}

        rm {log}
        """


rule minimap_sample_paired_contigs:
    """Map the sample paired reads to the sample assembly"""
    input:
        r1=os.path.join(dir.out.trim,"{sample}" + config.args.hostStr + ".paired.R1.fastq.gz"),
        r2=os.path.join(dir.out.trim,"{sample}" + config.args.hostStr + ".paired.R2.fastq.gz"),
        contigs = os.path.join(dir.out.assembly, "{sample}", "{sample}.contigs.fa"),
    output:
        temp(os.path.join(dir.out.assembly, "{sample}", "{sample}.pe.bam"))
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    log:
        os.path.join(dir.out.stderr, "minimap_sample_paired_contigs.{sample}.log")
    benchmark:
        os.path.join(dir.out.bench, "minimap_sample_paired_contigs.{sample}.txt")
    group:
        "assembly"
    shell:
        """
        {{
        minimap2 -t {threads} -ax sr {input.contigs} {input.r1} {input.r2} | \
            samtools sort -n -o {output};
        }} 2> {log} 
        rm {log}
        """


rule minimap_sample_paired_singletons_contigs:
    """Map the sample unpaired reads to the sample assembly"""
    input:
        s=os.path.join(dir.out.trim,"{sample}" + config.args.hostStr + ".paired.S.fastq.gz"),
        contigs = os.path.join(dir.out.assembly,"{sample}","{sample}.contigs.fa")
    output:
        temp(os.path.join(dir.out.assembly,"{sample}","{sample}.assemblyUnmapped.s.fastq"))
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    benchmark:
        os.path.join(dir.out.bench, "minimap_sample_paired_singletons_contigs.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "minimap_sample_paired_singletons_contigs.{sample}.log")
    group:
        "assembly"
    shell:
        """
        {{
        minimap2 -t {threads} -ax sr {input.contigs} {input.s} | \
            samtools sort -n | \
            samtools fastq -f 4 > {output};
        }} 2> {log}
        rm {log}
        """


rule minimap_sample_unpaired_contigs:
    """Map the sample paired reads to the sample assembly"""
    input:
        r1 = os.path.join(dir.out.trim, "{sample}" + config.args.hostStr + ".single.fastq.gz"),
        contigs = os.path.join(dir.out.assembly,"{sample}","{sample}.rename.contigs.fa"),
    output:
        temp(os.path.join(dir.out.assembly,"{sample}","{sample}.assemblyUnmapped.s.fastq"))
    conda:
        os.path.join(dir.env, "minimap2.yaml")
    threads:
        config.resources.med.cpu
    resources:
        mem_mb = config.resources.med.mem,
        time = config.resources.med.time
    log:
        os.path.join(dir.out.stderr, "minimap_sample_unpaired_contigs.{sample}.log")
    benchmark:
        os.path.join(dir.out.bench, "minimap_sample_unpaired_contigs.{sample}.txt")
    group:
        "assembly"
    shell:
        """
        minimap2 \
            -t {threads} \
            -ax sr \
            {input.contigs} \
            {input.r1} \
            2> {log} \
            | samtools fastq -f 4 > {output}
        """


rule samtools_fastq_paired:
    """Grab the paired unmapped reads (neither pair mapped)"""
    input:
        os.path.join(dir.out.assembly,"{sample}","{sample}.pe.bam")
    output:
        r1 = temp(os.path.join(dir.out.assembly, "{sample}", "{sample}.assemblyUnmapped_R1.fastq")),
        r2 = temp(os.path.join(dir.out.assembly, "{sample}", "{sample}.assemblyUnmapped_R2.fastq")),
    conda:
        os.path.join(dir.env,   "samtools.yaml")
    resources:
        time = config.resources.sml.time
    benchmark:
        os.path.join(dir.out.bench, "samtools_fastq_paired.{sample}.txt")
    log:
        os.path.join(dir.out.stderr, "samtools_fastq_paired.{sample}.log")
    group:
        "assembly"
    shell:
        """
        samtools fastq -f 77 {input} > {output.r1} 2> {log}
        samtools fastq -f 141 {input} > {output.r2} 2> {log}
        """


rule pool_paired_unmapped_R1:
    """Concatenate the unmapped, paired R1 reads for all samples"""
    input:
        targets.unmapped.r1
    output:
        os.path.join(dir.out.trim, "rescue" + config.args.hostStr + ".paired.R1.fastq.gz"),
    conda:
        os.path.join(dir.env, "pigz.yaml")
    threads:
        config.resources.med.cpu
    resources:
        time = config.resources.sml.time
    group:
        "assemblyRescue"
    shell:
        """cat {input} | pigz -c -p {threads} - > {output}"""


rule pool_paired_unmapped_R2:
    """Concatenate the unmapped, paired R2 reads for all samples"""
    input:
        targets.unmapped.r2
    output:
        os.path.join(dir.out.trim, "rescue" + config.args.hostStr + ".paired.R2.fastq.gz"),
    conda:
        os.path.join(dir.env, "pigz.yaml")
    threads:
        config.resources.med.cpu
    resources:
        time = config.resources.sml.time
    group:
        "assemblyRescue"
    shell:
        """cat {input} | pigz -c -p {threads} - > {output}"""


rule pool_unmapped_singletons:
    """Concatenate the unmapped, unpaired reads for all samples"""
    input:
        targets.unmapped.s
    output:
        os.path.join(dir.out.trim, "rescue" + config.args.hostStr + ".paired.S.fastq.gz"),
    conda:
        os.path.join(dir.env, "pigz.yaml")
    threads:
        config.resources.med.cpu
    resources:
        time = config.resources.sml.time
    group:
        "assemblyRescue"
    shell:
        """
        cat {input} | pigz -c -p {threads} - > {output}
        """


rule concatenate_contigs:
    """Assembly step 03: Concatenate individual assembly outputs (contigs) into a single file"""
    input:
        expand(os.path.join(dir.out.assembly, "{sample}", "{sample}.rename.contigs.fa"), sample=samples.names),
        os.path.join(dir.out.assembly, "rescue", "rescue.rename.contigs.fa")
    output:
        os.path.join(dir.out.assembly, "all_sample_contigs.fasta.gz")
    params:
        dirs = expand(os.path.join(dir.out.assembly,"{sample}"), sample=samples.names + ["rescue"]),
        compression= '-' + str(config.qc.compression)
    threads:
        config.resources.med.cpu
    resources:
        time = config.resources.sml.time
    conda:
        os.path.join(dir.env, "pigz.yaml")
    group:
        "assemblyRescue"
    shell:
        """
        cat {input} | pigz -p {threads} {params.compression} -c > {output}
        rm -rf {params.dirs}
        """