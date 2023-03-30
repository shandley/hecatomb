"""
Snakefile for downloading the databases. Only need to run this once.

Set a custom database install location in hecatomb/snakemake/config/config.yaml

Use the launcher to run:
hecatomb install

Rob Edwards, Feb 2020
Updated: Michael Roach, Q2 2021
"""

import attrmap as ap


# load db file config
configfile: os.path.join(workflow.basedir, "../", "config", "config.yaml")
configfile: os.path.join(workflow.basedir, "../", "config", "dbFiles.yaml")
config = ap.AttrMap(config)


# directories
include: "rules/preflight/directories.smk"


rule all:
    input:
        expand(os.path.join(dir.dbs.base, "{file}"), file=config.dbs.files),
        expand(os.path.join(dir.dbs.base, "{file}"), file=config.dbtax.files),


rule download_db_file:
    """Download a Hecatomb-maintained DB file."""
    output:
        os.path.join(dir.dbs.base, "{path}","{file}")
    wildcard_constraints:
        path="aa|nt|host|contaminants|tables"
    run:
        import urllib.request
        import urllib.parse
        import shutil
        dlUrl1 = urllib.parse.urljoin(config.dbs.mirror1, os.path.join(wildcards.path, wildcards.file))
        dlUrl2 = urllib.parse.urljoin(config.dbs.mirror2, os.path.join(wildcards.path, wildcards.file))
        try:
            with urllib.request.urlopen(dlUrl1) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)
        except:
            with urllib.request.urlopen(dlUrl2) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)


rule download_taxdump:
    """Download the current NCBI taxdump."""
    output:
        expand(os.path.join(dir.dbs.base, "{file}"), file=config.dbtax.files),
        temp(config.dbtax.tar)
    params:
        url = config.dbtax.url,
        tar = config.dbtax.tar,
        md5 = config.dbtax.md5,
        dir = os.path.join(dir.dbs.base, "tax", "taxonomy")
    conda:
        os.path.join(dir.env, "curl.yaml")
    shell:
        """
        curl {params.url} -o {params.tar}
        curl {params.md5} | md5sum -c
        mkdir -p {params.dir}
        tar xvf {params.tar} -C {params.dir}
        """
