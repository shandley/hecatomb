
def copy_log():
    try:
        assert (ap.utils.to_dict(config.args)["log"]) is not None
        shell("cat {log} >> " + config.args.log)
    except (KeyError, AssertionError):
        pass


# Check for Database files
dbFail = False
for f in config.dbs.files:
    dbFile = os.path.join(dir.dbs.base, f)
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
    if os.path.isdir(dir.out.stderr):
        oldLogs = filter(re.compile(r'^(?!old_).*.log').match, os.listdir(dir.out.stderr))
        for logfile in oldLogs:
            os.rename(os.path.join(dir.out.stderr, logfile), os.path.join(dir.out.stderr, f'old_{logfile}'))

# Success message
onsuccess:
    # copy_log()
    sys.stderr.write('\n\n    Hecatomb finished successfully!\n\n')

# Fail message and dump failed log outputs to a crash report file
onerror:
    # copy_log()
    sys.stderr.write('\n\n    FATAL: Hecatomb encountered an error.\n\n')
    logfiles = list(filter(re.compile(r'^(?!old_).*.log').match, os.listdir(dir.out.stderr)))
    if len(logfiles) > 0:
        sys.stderr.write('    Dumping all error logs to "hecatomb.crashreport.log"\n')
        with open('hecatomb.crashreport.log', 'w') as crashDump:
            for file in logfiles:
                if os.path.getsize(os.path.join(dir.out.stderr, file)) > 0:
                    rulename = re.sub(r'.log$', '', file)
                    crashDump.write(''.join(('\n','-'*10,f'\nstderr for rule {rulename}:\n','-'*10,'\n')))
                    with open(os.path.join(dir.out.stderr, file), 'r') as logfh:
                        for line in logfh:
                            crashDump.write(line)
