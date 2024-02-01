
# LOAD CONFIG
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "immutable.yaml")
res = config["resources"]
config = config["hecatomb"]


# DIRECTORIES
include: os.path.join("rules", "preflight", "directories.smk")


# NEW HOST TARGET
config["outFasta"] = os.path.join(
    dir["dbs"]["hostBase"], config["args"]["hostName"], "masked_ref.fa.gz"
)
config["virRefSeq"] = os.path.join(dir["dbs"]["hostBase"], config["dbvirRefSeq"]["file"])


# TARGET RULE
rule all:
    input:
        config["outFasta"]


# ADD HOST RULES
rule dl_refseq_viral:
    """Download the refseq viral genomic file needed"""
    output:
        config["virRefSeq"]
    params:
        url = config["dbvirRefSeq"]["url"],
        file = config["virRefSeq"]
    conda:
        os.path.join(dir["env"], "curl.yaml")
    shell:
        "curl {params.url} -o {params.file}"


rule minimap_viral_refseq:
    """Map the viral refseq db to the host"""
    input:
        ref = config["args"]["hostFa"],
        vir = config["virRefSeq"]
    output:
        temp(os.path.join(dir["out"]["temp"], f'{config["args"]["hostName"]}.bed'))
    params:
        config["addHost"]["minViralAlnLen"]
    conda:
        os.path.join(dir["env"], "minimap2.yaml")
    resources:
        mem_mb = res["med"]["mem"],
        mem = str(res["med"]["mem"]) + "MB",
        time = res["med"]["time"]
    threads:
        res["med"]["cpu"]
    log:
        os.path.join(dir["out"]["stderr"], 'minimap_viral_refseq.log')
    benchmark:
        os.path.join(dir["out"]["bench"], "minimap_viral_refseq.log")
    shell:
        """
        {{ minimap2 -t {threads} {input.ref} {input.vir}  \
            | awk '$11>{params}' \
            | cut -f6,8,9 \
            | sort -k1,1 -k2,2n -k3,3n; \
        }} > {output} 2> {log}
        """


rule mask_fasta:
    """Mask the host genome using bedtools"""
    input:
        fa = config["args"]["hostFa"],
        aln = os.path.join(dir["out"]["temp"], f'{config["args"]["hostName"]}.bed')
    output:
        mask = temp(os.path.join(dir["out"]["temp"], f'{config["args"]["hostName"]}.mask.bed')),
        fa = config["outFasta"]
    params:
        fa = os.path.join(dir["dbs"]["hostBase"],config["args"]["hostName"],"masked_ref.fa"),
    resources:
        mem_mb = res["med"]["mem"],
        mem = str(res["med"]["mem"]) + "MB",
        time = res["med"]["time"]
    threads:
        res["med"]["cpu"]
    conda:
        os.path.join(dir["env"], "bedtools.yaml")
    log:
        os.path.join(dir["out"]["stderr"], "mask_fasta.log")
    benchmark:
        os.path.join(dir["out"]["bench"], "mask_fasta.log")
    shell:
        """
        bedtools merge -i {input.aln} > {output.mask}
        bedtools maskfasta -fi {input.fa} -bed {output.mask} -fo {params.fa}
        pigz -p {threads} -1 {params.fa}
        """
