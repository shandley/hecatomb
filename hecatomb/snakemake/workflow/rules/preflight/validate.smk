import glob
import re


# Concatenate Snakemake's own log file with the master log file
def copy_log():
    files = glob.glob(os.path.join(".snakemake", "log", "*.snakemake.log"))
    if not files:
        return None
    current_log = max(files, key=os.path.getmtime)
    shell("cat " + current_log + " >> " + config["args"]["log"])


# simulate db files
if config["args"]["simulate"]:
    for f in config["dbs"]["files"] + config["dbtax"]["files"]:
        dbFile = os.path.join(dir["dbs"]["base"], f)
        os.makedirs(os.path.dirname(dbFile), exist_ok=True)
        with open(dbFile, 'a'):
            os.utime(dbFile, None)


# Check for Database files
dbFail = False
for f in config["dbs"]["files"] + config["dbtax"]["files"]:
    dbFile = os.path.join(dir["dbs"]["base"], f)
    if not os.path.isfile(dbFile):
        dbFail = True
        sys.stderr.write(f"    ERROR: missing database file {dbFile}\n")
if dbFail:
    sys.stderr.write("\n"
        "    FATAL: One or more database files is missing.\n"
        "    Please run 'hecatomb install' to download the missing database files.\n"
        "\n")
    sys.exit(1)

# Cleanup old logfiles etc.
onstart:
    if os.path.isdir(dir["out"]["stderr"]):
        oldLogs = filter(re.compile(r'^(?!old_).*').match, os.listdir(dir["out"]["stderr"]))
        for logfile in oldLogs:
            os.rename(os.path.join(dir["out"]["stderr"], logfile), os.path.join(dir["out"]["stderr"], f'old_{logfile}'))

# Success message
onsuccess:
    copy_log()
    sys.stderr.write('\n\n    Hecatomb finished successfully!\n\n')

# Fail message and dump failed log outputs to a crash report file
onerror:
    copy_log()
    sys.stderr.write('\n\n    FATAL: Hecatomb encountered an error.\n\n')
    sys.stderr.write(f'Check the Hecatomb logs directory for command-related errors:\n\n{dir["out"]["stderr"]}\n\n')
    if config["args"]["profile"]:
        sys.stderr.write(
            'Also check your scheduler logs for sheduler-related errors. Your profile determins where these are saved'
            'but by default is usually logs/ in your working directory.\n\n'
        )
