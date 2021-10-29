# Check for Database files
dbFail = False
for f in config['dbFiles']:
    dbFile = os.path.join(DBDIR, f)
    if not os.path.isfile(dbFile):
        dbFail = True
        sys.stderr.write("\n"
            f"    ERROR: missing database file {dbFile}"
            "\n")
if dbFail:
    sys.stderr.write("\n"
        "    FATAL: One or more database files is missing.\n"
        "    Please run 'hecatomb install' to download the missing database files.\n"
        "\n")
    sys.exit(1)

# Cleanup old logfiles etc.
onstart:
    sys.stderr.write('\n    Marking old log files.\n')
    oldLogs = filter(re.compile(r'^(?!old_).*.log').match, os.listdir(STDERR))
    for logfile in oldLogs:
        os.rename(os.path.join(STDERR, logfile), os.path.join(STDERR, f'old_{logfile}'))
    sys.stderr.write('    Beginning pipeline')

# Success message
onsuccess:
    sys.stderr.write('\n\n    Hecatomb finished successfully!\n\n')

# Fail message and dump failed log outputs to a crash report file
onerror:
    sys.stderr.write('\n    FATAL: Hecatomb encountered an error.')
    logfiles = filter(re.compile(r'^(?!old_).*.log').match, os.listdir(STDERR))
    if len(logfiles) > 0:
        sys.stderr.write('\n           Dumping all error logs to "hecatomb.errorLogs.txt"'))
        with open('hecatomb.crashreport.log', 'w') as crashDump:
            for file in logfiles:
                rulename = re.sub(r'.log$', '', file)
                crashDump.write('\n'
                                '-'*10
                                f'\nSTDERR for rule {rulename}:\n'
                                '-'*10
                                '\n')
                with open(os.path.join(STDERR,file),'r') as logfh:
                    for line in logfh:
                        crashDump.write(line)