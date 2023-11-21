
# LOAD CONFIG FILES
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "immutable.yaml")
resources = config["resources"]
config = config["hecatomb"]


# DIRECTORIES
include: os.path.join("rules", "preflight", "directories.smk")


# TARGET OUTPUTS
rule all:
    input:
        expand(os.path.join(dir["dbs"]["base"], "{file}"), file=config["dbs"]["files"]),
        expand(os.path.join(dir["dbs"]["base"], "{file}"), file=config["dbtax"]["files"]),


# DOWNLOAD RULES
rule download_db_file:
    """Download a Hecatomb-maintained DB file."""
    output:
        os.path.join(dir["dbs"]["base"], "{path}","{file}")
    wildcard_constraints:
        path="aa|nt|host|contaminants|tables"
    run:
        import urllib.request
        import urllib.parse
        import shutil
        dlUrl1 = urllib.parse.urljoin(config["dbs"]["mirror1"], os.path.join(wildcards.path, wildcards.file))
        dlUrl2 = urllib.parse.urljoin(config["dbs"]["mirror2"], os.path.join(wildcards.path, wildcards.file))
        try:
            with urllib.request.urlopen(dlUrl1) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)
        except:
            with urllib.request.urlopen(dlUrl2) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)


rule download_taxdump:
    """Download the current NCBI taxdump."""
    output:
        expand(os.path.join(dir["dbs"]["base"], "{file}"), file=config["dbtax"]["files"]),
        temp(config["dbtax"]["tar"])
    params:
        url = config["dbtax"]["url"],
        tar = config["dbtax"]["tar"],
        md5 = config["dbtax"]["md5"],
        dir = os.path.join(dir["dbs"]["base"], "tax", "taxonomy")
    conda:
        os.path.join(dir["env"], "curl.yaml")
    shell:
        """
        curl {params.url} -o {params.tar}
        curl {params.md5} | md5sum -c
        mkdir -p {params.dir}
        tar xvf {params.tar} -C {params.dir}
        """
