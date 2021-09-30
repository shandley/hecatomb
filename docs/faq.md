# Hecatomb FAQ

## Where can I get help or support?

[Open an issue on github](https://github.com/shandley/hecatomb/issues).
There is no question too stupid.
The pipeline is still in the early days of development, so we are expecting that there will still be bugs in the code.

## The pipeline died, what went wrong?

A great many things may have gone wrong!
You can simply try rerunning the pipeline, and it will pick up where it left off. 

We try to make the logging as thorough as possible and have meaningful error messages where we can.
Look through the terminal output produced by Snakemake for the failed job.
It will appear in red text for the terminal, or you can look through the snakemake log in `.snakemake/logs/`.
Each rule has its own directory in `hecatomb_out/STDERR/` for log files, 
but the Snakemake error code will tell you exactly what file to look at.

If you're running with a profile on a HPC, the job may have failed for reasons relating to the scheduler.
For instance, the memory or time may have exceeded what was requested.
Look in your scheduler logs to see if this is the case.
In the example profile, the slurm logs are saved in `logs/{rule}/{jobid}`.
The rule name and jobid are supplied in the Snakemake error messages, making it easy to find the relevant logs.

## Hecatomb takes too long, what can I do?

First, try running Hecatomb with the `--fast` flag.
The MMSeqs steps are by far the most time consuming steps. 
The `--fast` flag will tell Hecatomb to use MMSeqs settings that are much much faster, but not quite as sensitive.
You should also configure your installation to utilise as many CPUs and as much memory as possible.
See [default resources config](https://hecatomb.readthedocs.io/en/latest/advanced/#default-resources) for more info.

## I've run the pipeline, now what?

Have a look at [the tutoria](#) which goes through some example plots and analyses.

