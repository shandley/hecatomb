
dir = dict()
dir["dbs"] = dict()
dir["out"] = dict()

### DATABASE LOCATION
try:
    assert(config["args"]["databases"]) is not None
    dir["dbs"]["base"] = config["args"]["databases"]
except (KeyError,AssertionError):
    try:
        assert(os.environ["HECATOMB_DB"]) is not None
        dir["dbs"]["base"] = os.environ["HECATOMB_DB"]
    except (KeyError, AssertionError):
        dir["dbs"]["base"] = os.path.join(workflow.basedir,"..","databases")


### OUTPUT LOCATION
try:
    assert(config["args"]["output"]) is not None
    dir["out"]["base"] = config["args"]["output"]
except (KeyError, AssertionError):
    dir["out"]["base"] = "hecatomb.out"


### WORKFLOW DIRs
dir["env"]     = os.path.join(workflow.basedir, "envs")
dir["scripts"] = os.path.join(workflow.basedir, "scripts")


### DATABASE DIRs
dir["dbs"]["contaminants"]  = os.path.join(dir["dbs"]["base"], "contaminants")
dir["dbs"]["taxonomy"]      = os.path.join(dir["dbs"]["base"], "tax", "taxonomy")
dir["dbs"]["tables"]        = os.path.join(dir["dbs"]["base"], "tables")
dir["dbs"]["hostBase"]      = os.path.join(dir["dbs"]["base"], "host")
dir["dbs"]["secondaryAA"]   = os.path.join(dir["dbs"]["base"], "aa", "virus_secondary_aa")
dir["dbs"]["secondaryNT"]   = os.path.join(dir["dbs"]["base"], "nt", "virus_secondary_nt")

if config["args"]["custom_aa"]:
    dir["dbs"]["primaryAA"] = config["args"]["custom_aa"]
else:
    dir["dbs"]["primaryAA"] = os.path.join(dir["dbs"]["base"], "aa", "virus_primary_aa", "sequenceDB")
if config["args"]["custom_nt"]:
    dir["dbs"]["primaryNT"] = config["args"]["custom_nt"]
else:
    dir["dbs"]["primaryNT"] = os.path.join(dir["dbs"]["base"], "nt", "virus_primary_nt", "sequenceDB")


### OUTPUT DIRs
dir["out"]["results"]         = os.path.join(dir["out"]["base"], "results")
dir["out"]["processing"]      = os.path.join(dir["out"]["base"], "processing")
if config["args"]["trim"]: # todo: something better
    dir["out"]["trim"]        = os.path.join(dir["out"]["base"], "trimnami", "results", config["args"]["trim"])
dir["out"]["temp"]            = os.path.join(dir["out"]["processing"], "temp")
dir["out"]["stderr"]          = os.path.join(dir["out"]["base"], "logs")
dir["out"]["bench"]           = os.path.join(dir["out"]["base"], "benchmarks")
dir["out"]["assembly"]        = os.path.join(dir["out"]["processing"], "assembly")
dir["out"]["mapping"]         = os.path.join(dir["out"]["processing"], "mapping")
dir["out"]["stats"]           = os.path.join(dir["out"]["processing"], "stats")
dir["out"]["primaryAA"]       = os.path.join(dir["out"]["processing"], "mmseqs_aa_primary")
dir["out"]["secondaryAA"]     = os.path.join(dir["out"]["processing"], "mmseqs_aa_secondary")
dir["out"]["primaryNT"]       = os.path.join(dir["out"]["processing"], "mmseqs_nt_primary")
dir["out"]["secondaryNT"]     = os.path.join(dir["out"]["processing"], "mmseqs_nt_secondary")
