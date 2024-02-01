"""
Command line interface for installing and running Hecatomb
"""

import os
import click
import yaml
import glob

from snaketool_utils.cli_utils import (
    OrderedCommands,
    run_snakemake,
    copy_config,
    echo_click,
    msg_box,
)


def snake_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(snake_base("hecatomb.VERSION"), "r") as f:
        version = f.readline()
    return version


def print_citation():
    with open(snake_base("hecatomb.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def default_to_output(ctx, param, value):
    """Callback for click options; places value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="hecatomb.out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="hecatomb.config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/hecatomb.config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=32, show_default=True
        ),
        click.option(
            "--profile", help="Snakemake profile", default=None, show_default=False
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("snakemake", "conda")),
            help="Custom conda env directory",
            type=click.Path(),
            show_default=False,
        ),
        click.option(
            "--snake-default",
            multiple=True,
            default=[
                "--rerun-incomplete",
                "--printshellcmds",
                "--nolock",
                "--show-failed-logs",
                "--conda-frontend conda",
            ],
            help="Customise Snakemake runtime args",
            show_default=True,
        ),
        click.option(
            "--log",
            default="hecatomb.log",
            callback=default_to_output,
            hidden=True,
        ),
        click.option(
            "--simulate/--no-simulate",
            default=False,
            hidden=True,
            help="Create empty database files for simulated runs testing",
        ),
        click.option(
            "--system-config",
            default=snake_base(os.path.join("snakemake", "config", "config.yaml")),
            hidden=True,
        ),
        click.argument("snake_args", nargs=-1,),
    ]
    for option in reversed(options):
        func = option(func)
    return func


def run_test_options(func):
    """Common command line options for "run" and "test" """
    options = [
        click.option(
            "--trim",
            help="Trimming engine for trimnami",
            default="fastp",
            show_default=True,
            type=click.Choice(["fastp", "prinseq", "roundAB", "filtlong", "notrim", "cutadapt"]),
        ),
        click.option(
            "--fastqc/--no-fastqc",
            default=False,
            help="Generate fastqc reports",
            show_default=True,
        ),
        click.option(
            "--assembly",
            help="Assembly method: [cross]-assembly or [merged]-assembly",
            default="merged",
            show_default=True,
            type=click.Choice(["cross", "merged"]),
        ),
        click.option(
            "--custom-aa",
            help="Custom protein fasta for prefiltering",
            type=click.Path(readable=True, exists=True),
            default=None,
            show_default=False,
        ),
        click.option(
            "--custom-nt",
            help="Custom nucleotide fasta for prefiltering",
            type=click.Path(readable=True, exists=True),
            default=None,
            show_default=False,
        ),
        click.option(
            "--search",
            help="MMSeqs search speed settings",
            default="sensitive",
            type=click.Choice(["fast", "sensitive"]),
            show_default=True,
        ),
        click.option(
            "--host",
            help="Host genome name for filtering, or 'none' for no host removal.",
            default="human",
            show_default=True,
        ),
    ]
    for option in reversed(options):
        func = option(func)
    return func


@click.group(
    cls=OrderedCommands, context_settings=dict(help_option_names=["-h", "--help"])
)
@click.version_option(get_version(), "-v", "--version", is_flag=True)
def cli():
    """Viral metagenomics framework for short and longread sequencing.
    \b
    For more options, run:
    hecatomb command --help"""
    pass


def print_splash():
    click.echo(
        """
\b
██╗  ██╗███████╗ ██████╗ █████╗ ████████╗ ██████╗ ███╗   ███╗██████╗
██║  ██║██╔════╝██╔════╝██╔══██╗╚══██╔══╝██╔═══██╗████╗ ████║██╔══██╗
███████║█████╗  ██║     ███████║   ██║   ██║   ██║██╔████╔██║██████╔╝
██╔══██║██╔══╝  ██║     ██╔══██║   ██║   ██║   ██║██║╚██╔╝██║██╔══██╗
██║  ██║███████╗╚██████╗██║  ██║   ██║   ╚██████╔╝██║ ╚═╝ ██║██████╔╝
╚═╝  ╚═╝╚══════╝ ╚═════╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝     ╚═╝╚═════╝
"""
    )


def help_msg_epilog():
    return """
\b
CLUSTER EXECUTION:
hecatomb run ... --profile [profile]
\b
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           hecatomb run --reads [file/dir]
Specify threads:    hecatomb run ... --threads [threads]
Disable conda:      hecatomb run ... --no-use-conda 
Change defaults:    hecatomb run ... --snake-default="-k --nolock"
Add Snakemake args: hecatomb run ... --dry-run --keep-going --touch
Specify stages:     hecatomb run ... all print_stages
\b
AVAILABLE STAGES:
    all                 Run everything (default)
    preprocessing       Preprocessing steps only
    assembly            Assembly steps (+ preprocess)
    read_annotations    Read annotations (+ preprocess)
    contig_annotations  Contig annotations (+ preprocess,assemble)
    print_stages        List available stages
"""


@click.command(
    epilog=help_msg_epilog(),
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--reads",
    "reads",
    help="Input file/directory",
    type=str,
    default=None,
    required=True,
)
@run_test_options
@common_options
def run(**kwargs):
    """Run hecatomb"""

    merge_config = {
        "hecatomb":{
            "args": kwargs
        }
    }

    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(
            os.path.join("snakemake", "workflow", "Hecatomb.smk")
        ),
        merge_config=merge_config,
        **kwargs
    )


@click.command(
    epilog=help_msg_epilog(),
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@run_test_options
@common_options
def test(**kwargs):
    """Run the Hecatomb test dataset"""

    kwargs["reads"] = snake_base("test_data")

    merge_config = {
        "hecatomb": {
            "args": kwargs
        }
    }

    run_snakemake(
        snakefile_path=snake_base(
            os.path.join("snakemake", "workflow", "Hecatomb.smk")
        ),
        merge_config=merge_config,
        **kwargs
    )


@click.command(
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    )
)
@common_options
def config(**kwargs):
    """Copy the system default config file"""
    copy_config(kwargs["configfile"])


@click.command(
    epilog=help_msg_epilog(),
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@common_options
@run_test_options
def install(**kwargs):
    """Install the Hecatomb databases"""

    merge_config = {
        "hecatomb": {
            "args": kwargs
        }
    }
    run_snakemake(
        snakefile_path=snake_base(
            os.path.join("snakemake", "workflow", "DownloadDB.smk")
        ),
        merge_config=merge_config,
        **kwargs
    )


@click.command(
    epilog=help_msg_epilog(),
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--comb",
    multiple=True,
    required=True,
    show_default=False,
    help="Two or more Hecatomb output directories to combine. e.g. --comb dir1/ --comb dir2/ ...",
)
@common_options
def combine(**kwargs):
    """Combine multiple Hecatomb runs"""
    kwargs["combineRuns"] = list(kwargs["comb"])
    merge_config = {
        "hecatomb": {
            "args": kwargs
        }
    }
    run_snakemake(
        snakefile_path=snake_base(
            os.path.join("snakemake", "workflow", "combineOutputs.smk")
        ),
        merge_config=merge_config,
        **kwargs
    )


@click.command(
    epilog=help_msg_epilog(),
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option(
    "--host", help="Name for your host genome", show_default=False, required=True
)
@click.option(
    "--host-fa", help="Host genome fasta file", show_default=False, required=True
)
@common_options
def add_host(**kwargs):
    """Add a new host genome to use with Hecatomb"""
    kwargs["hostFa"] = kwargs["host_fa"]
    kwargs["hostName"] = kwargs["host"]
    merge_config = {
        "hecatomb": {
            "args": kwargs
        }
    }
    run_snakemake(
        snakefile_path=snake_base(os.path.join("snakemake", "workflow", "AddHost.smk")),
        merge_config=merge_config,
        **kwargs
    )


@click.command()
@common_options
def list_hosts(**kwargs):
    """List the available host genomes"""
    copy_config(kwargs["configfile"])
    with open(kwargs["configfile"], "r") as f:
        config = yaml.safe_load(f)
    dbdir = snake_base(os.path.join("snakemake", "databases"))
    try:
        if config["Databases"] is not None:
            dbdir = config["Databases"]
    except KeyError:
        pass
    host_path = os.path.join(dbdir, "host", "*")
    host_fastas = list([os.path.basename(x) for x in glob.glob(host_path)])
    try:
        host_fastas.remove("virus_shred.fasta.gz")
    except ValueError:
        pass
    msg_box("Available host genomes", "\n".join(host_fastas))


@click.command()
def citation(**kwargs):
    """Print the citation(s) for this tool"""
    print_citation()


cli.add_command(run)
cli.add_command(test)
cli.add_command(config)
cli.add_command(install)
cli.add_command(combine)
cli.add_command(add_host)
cli.add_command(list_hosts)
cli.add_command(citation)


def main():
    print_splash()
    cli()


if __name__ == "__main__":
    main()
