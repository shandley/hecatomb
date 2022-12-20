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
        expand(os.path.join(dir.dbs.base, '{file}'), file=config.dbs.files)


rule download_db_file:
    """Generic rule to download a DB file."""
    output:
        os.path.join(dir.dbs.base, '{file}')
    run:
        import urllib.request
        import urllib.parse
        import shutil
        dlUrl1 = urllib.parse.urljoin(config.dbs.mirror1, wildcards.file)
        dlUrl2 = urllib.parse.urljoin(config.dbs.mirror2, wildcards.file)
        try:
            with urllib.request.urlopen(dlUrl1) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)
        except:
            with urllib.request.urlopen(dlUrl2) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)
