"""
Database and output locations for Hecatomb

Ensures consistent variable names and file locations for the pipeline, the database download script,
and the addHost script.
"""


### DATABASE BASE DIRECTORY
if config['Databases'] is None:
    DBDIR = os.path.join(workflow.basedir, '..', '..', 'databases')
else:
    DBDIR = config['Databases']


### OUTPUT DIRECTORY
if config['Output'] is None:
    OUTPUT = 'hecatomb_out'
else:
    OUTPUT = config['Output']


### DATABASE SUBDIRs
CONPATH = os.path.join(DBDIR, "contaminants")
TAX = os.path.join(DBDIR, "tax", "taxonomy")
TABLES = os.path.join(DBDIR, "tables")
HOSTPATH = os.path.join(DBDIR, "host")


### MMSEQS DBs
UNIVIRDB = os.path.join(DBDIR, "aa", "virus_primary_aa")
UNIREF50VIR = os.path.join(DBDIR, "aa", "virus_secondary_aa")
NCBIVIRDB = os.path.join(DBDIR, "nt", "virus_primary_nt")
POLYMICRODB = os.path.join(DBDIR, "nt", "virus_secondary_nt")


### OUTPUT DIRs
RESULTS = os.path.join(OUTPUT, 'RESULTS')
WORKDIR = os.path.join(OUTPUT, 'PROCESSING')
TMPDIR = os.path.join(WORKDIR, 'TMP')
STDERR = os.path.join(OUTPUT, 'STDERR')
BENCH = os.path.join(OUTPUT, 'BENCHMARKS')
SUMDIR = os.path.join('hecatomb_report')
ASSEMBLY = os.path.join(WORKDIR, 'ASSEMBLY')
MAPPING = os.path.join(WORKDIR, 'MAPPING')
STATS = os.path.join(WORKDIR, 'STATS')


# MMSEQS OUTPUT DIRs
PRIMARY_AA_OUT = os.path.join(RESULTS, "MMSEQS_AA_PRIMARY")
SECONDARY_AA_OUT = os.path.join(RESULTS, "MMSEQS_AA_SECONDARY")
PRIMARY_NT_OUT = os.path.join(RESULTS, "MMSEQS_NT_PRIMARY")
SECONDARY_NT_OUT = os.path.join(RESULTS, "MMSEQS_NT_SECONDARY")

