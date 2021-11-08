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

## I have too much data and am having memory problems

You can run your samples in batches, and even split samples into separate runs if you really need.
The `bigtable.tsv` files from different runs can be concatenated (just remember to remove the extra headers),
and your assemblies can be combined using Flye with the `--subassemblies` function.

If your sequencing depth is exceptionally high you could also subsample your FASTQ files 
(and run a couple at full coverage as a control).

Depending on which jobs are failing, you can reduce the cluster threshold (`CLUSTERID` in your `config.yaml`) to help 
reduce the size of your seqtable, and use more stringent e-value cutoffs for your MMSeqs searches 
(also in your `config.yaml`) to reduce the output file sizes.

## Hecatomb takes too long, what can I do?

First, try running Hecatomb with the `--fast` flag.
The MMSeqs steps are by far the most time consuming steps. 
The `--fast` flag will tell Hecatomb to use MMSeqs settings that are much much faster, but not quite as sensitive.
You should also configure your installation to utilise as many CPUs and as much memory as possible.
See [default resources config](https://hecatomb.readthedocs.io/en/latest/advanced/#default-resources) for more info.
If you don't need an assembly, you can skip those steps with `--skipAssembly` which will also save a bit of time.

## I've run the pipeline, now what?

Have a look at [the tutoria](#) which goes through some example plots and analyses.

