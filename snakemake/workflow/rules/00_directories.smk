"""
Database and output locations for Hecatomb

Ensures consistent variable names and file locations for the pipeline, the database download script,
and the addHost script.
"""

### DATABASE BASE DIRECTORY
if config['Databases'] is None:
    DBDIR = os.path.join(workflow.basedir,'../../databases')
else:
    DBDIR = config['Databases']

### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'hecatomb_out'
else:
    OUTPUT = config['Output']


### DATABASE SUBDIRs
BACPATH = os.path.join(DBDIR, "bac_giant_unique_species")
CONPATH = os.path.join(DBDIR, "contaminants")
BACBT2 = os.path.join(DBDIR, "bac_giant_unique_species", "bac_uniquespecies_giant.masked_Ns_removed")
TAX = os.path.join(DBDIR, "taxonomy")
TABLES = os.path.join(DBDIR, "tables")


### MMSEQS DBs
UNIVIRDB = os.path.join(DBDIR, "proteins", "virus_primary_aa")
UNIREF50VIR = os.path.join(DBDIR, "proteins", "virus_secondary_aa")
NCBIVIRDB = os.path.join(DBDIR, "nt", "virus_primary_nt")
POLYMICRODB = os.path.join(DBDIR, "nt", "virus_secondary_nt")


### OUTPUT DIRs
LOGS = 'logs'
RESULTS = os.path.join(OUTPUT, 'RESULTS')
WORKDIR = os.path.join(OUTPUT, 'PROCESSING')
TMPDIR = os.path.join(WORKDIR, 'TMP')
STDERR = os.path.join(OUTPUT,'STDERR')
BENCH = os.path.join(WORKDIR, 'BENCHMARKS')
QC = os.path.join(WORKDIR,'QC')
ASSEMBLY = os.path.join(WORKDIR,'ASSEMBLY')
STATS = os.path.join(WORKDIR,'STATS')

# MMSEQS OUTPUT DIRs
PRIMARY_AA_OUT = os.path.join(RESULTS, "MMSEQS_AA_PRIMARY")
SECONDARY_AA_OUT = os.path.join(RESULTS, "MMSEQS_AA_SECONDARY")
PRIMARY_NT_OUT = os.path.join(RESULTS, "MMSEQS_NT_PRIMARY")
SECONDARY_NT_OUT = os.path.join(RESULTS, "MMSEQS_NT_SECONDARY")


# DIRs TO INITIALIZE TODO: check and move to Hecatomb.smk if possible
for dir in [TMPDIR, WORKDIR, STDERR]:
    if not os.path.exists(dir):
        os.makedirs(dir, exist_ok=True)
