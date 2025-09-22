rule add_host_dl_refseq_viral:
    """Download the refseq viral genomic file needed"""
    output:
        config["virRefSeq"]
    params:
        url = config["dbvirRefSeq"]["url"],
        file = config["virRefSeq"]
    conda:
        os.path.join("..", "..", "envs", "krona_curl_zstd_pysam.yaml")
    container:
        config["hecatomb"]["container"]["krona_curl_zstd_pysam"]
    shell:
        "curl {params.url} -o {params.file}"


rule add_host_minimap_viral_refseq:
    """Map the viral refseq db to the host"""
    input:
        ref = config["args"]["hostFa"],
        vir = config["virRefSeq"]
    output:
        temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], config["args"]["hostName"] + ".bed"))
    params:
        config["addHost"]["minViralAlnLen"]
    conda:
        os.path.join("..", "..", "envs", "minimap2_samtools.yaml")
    container:
        config["hecatomb"]["container"]["minimap2_samtools"]
    resources:
        mem_mb = res["med"]["mem"],
        mem = str(res["med"]["mem"]) + "MB",
        time = res["med"]["time"]
    threads:
        res["med"]["cpu"]
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], 'minimap_viral_refseq.log')
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "minimap_viral_refseq.log")
    shell:
        """
        {{ minimap2 -t {threads} {input.ref} {input.vir}  \
            | awk '$11>{params}' \
            | cut -f6,8,9 \
            | sort -k1,1 -k2,2n -k3,3n; \
        }} > {output} 2> {log}
        """


rule add_host_mask_fasta:
    """Mask the host genome using bedtools"""
    input:
        fa = config["args"]["hostFa"],
        aln = os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], config["args"]["hostName"] + ".bed")
    output:
        mask = temp(os.path.join(config["hecatomb"]["args"]["output_paths"]["temp"], config["args"]["hostName"] + ".mask.bed")),
        fa = config["outFasta"]
    params:
        fa = os.path.join(config["hecatomb"]["args"]["database_paths"]["hostBase"],config["args"]["hostName"],"masked_ref.fa"),
    resources:
        mem_mb = res["med"]["mem"],
        mem = str(res["med"]["mem"]) + "MB",
        time = res["med"]["time"]
    threads:
        res["med"]["cpu"]
    conda:
        os.path.join("..", "..", "envs", "bbmap_bedtools_pigz_flye.yaml")
    container:
        config["hecatomb"]["container"]["bbmap_bedtools_pigz_flye"]
    log:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["log"], "mask_fasta.log")
    benchmark:
        os.path.join(config["hecatomb"]["args"]["output_paths"]["bench"], "mask_fasta.log")
    shell:
        """
        bedtools merge -i {input.aln} > {output.mask}
        bedtools maskfasta -fi {input.fa} -bed {output.mask} -fo {params.fa}
        pigz -p {threads} -1 {params.fa}
        """