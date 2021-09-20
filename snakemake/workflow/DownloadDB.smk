"""
Snakefile for downloading the databases. Only need to run this once.

Set a custom database install location in hecatomb/snakemake/config/config.yaml

Use the launcher to run:
hecatomb install

Rob Edwards, Feb 2020
Updated: Michael Roach, Q2 2021
"""

import os
import sys


### LOAD DEFAULT CONFIG
configfile: os.path.join(workflow.basedir, '../', 'config', 'config.yaml')


# directories
include: "rules/00_directories.smk"


# Mirrors
mirror1 = config['mirror1']
mirror2 = config['mirror2']


# targets
allDbFiles = []
for f in config['dbFiles']:
    allDbFiles.append(os.path.join(DBDIR, f))


rule all:
    input:
        allDbFiles


# rule unpack_db_file:
#     """Generic rule to unpack included files."""
#     input:
#         os.path.join(DBDIR, 'contaminants', '{file}.zst')
#     output:
#         os.path.join(DBDIR, 'contaminants', '{file}')
#     wildcard_constraints:
#         file = '.*(?<!zst)'
#     shell:
#         """zstd -d {input} -o {output}"""


rule download_db_file:
    """Generic rule to download a DB file."""
    output:
        os.path.join(DBDIR, '{file}')
    # wildcard_constraints:
    #     file = '(?!contaminants).*'
    run:
        import urllib.request
        import urllib.parse
        import shutil
        dlUrl1 = urllib.parse.urljoin(mirror1, wildcards.file)
        dlUrl2 = urllib.parse.urljoin(mirror2, wildcards.file)
        try:
            with urllib.request.urlopen(dlUrl1) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)
        except:
            with urllib.request.urlopen(dlUrl2) as r, open(output[0],'wb') as o:
                shutil.copyfileobj(r,o)
