"""
Command line interface for installing and running Hecatomb
"""

import os
import click

from .cli_util import (
    snake_base,
    get_version,
    default_to_output,
    copy_config,
    run_snakemake,
    OrderedCommands,
    print_citation,
    run_list_hosts,
)


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
        click.argument("snake_args", nargs=-1),
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
    click.echo("""
    \b
    ██╗  ██╗███████╗ ██████╗ █████╗ ████████╗ ██████╗ ███╗   ███╗██████╗
    ██║  ██║██╔════╝██╔════╝██╔══██╗╚══██╔══╝██╔═══██╗████╗ ████║██╔══██╗
    ███████║█████╗  ██║     ███████║   ██║   ██║   ██║██╔████╔██║██████╔╝
    ██╔══██║██╔══╝  ██║     ██╔══██║   ██║   ██║   ██║██║╚██╔╝██║██╔══██╗
    ██║  ██║███████╗╚██████╗██║  ██║   ██║   ╚██████╔╝██║ ╚═╝ ██║██████╔╝
    ╚═╝  ╚═╝╚══════╝ ╚═════╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝     ╚═╝╚═════╝
    """)


def help_msg_epilog():
    return ("""
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
    all             Run everything (default)
    preprocess      Preprocessing steps only
    assemble        Assembly steps (+ preprocess)
    annotate        Read annotations (+ preprocess)
    ctg_annotate    Contig annotations (+ preprocess,assemble)
    print_stages    List available stages
""")


@click.command(
    epilog=help_msg_epilog(),
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option('--reads', 'reads', help='Input file/directory', type=str, default=None, required=True)
@click.option('--library', help='Library type', default='paired', show_default=True,
              type=click.Choice(['paired', 'single', 'longread', 'roundAB']))
@click.option('--assembly', help='Assembly method: [cross]-assembly or [co]-assembly', default='cross',
              show_default=True, type=click.Choice(['cross', 'co']))
@click.option('--search', help='MMSeqs search speed settings', default='sensitive',
              type=click.Choice(['fast', 'sensitive']), show_default=True)
@click.option('--host', help='Host genome name for filtering', default='human', show_default=True)
@common_options
def run(reads, library, assembly, search, host, output, log, **kwargs):
    """Run hecatomb"""

    merge_config = {
        'args': {
            'reads': reads,
            'output': output,
            'host': host,
            'library': library,
            'assembly': assembly,
            'search': search,
            'log': log
        }
    }

    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join('snakemake', 'workflow', 'Hecatomb.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(
    epilog=help_msg_epilog(),
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option('--library', help='Library type', default='paired', show_default=True,
              type=click.Choice(['paired', 'single', 'longread', 'roundAB']))
@click.option('--assembly', help='Assembly method: [cross]-assembly or [co]-assembly', default='cross',
              show_default=True, type=click.Choice(['cross', 'co']))
@click.option('--search', help='MMSeqs search speed settings', default='sensitive',
              type=click.Choice(['fast', 'sensitive']), show_default=True)
@click.option('--host', help='Host genome name for filtering', default='human', show_default=True)
@common_options
def test(library, assembly, search, host, output, log, **kwargs):
    """Run the Hecatomb test dataset"""

    reads = snake_base('test_data')

    merge_config = {
        'args': {
            'reads': reads,
            'output': output,
            'host': host,
            'library': library,
            'assembly': assembly,
            'search': search,
            'log': log
        }
    }

    run_snakemake(
        snakefile_path=snake_base(os.path.join('snakemake', 'workflow', 'Hecatomb.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@common_options
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile)


@click.command(epilog=help_msg_epilog(),
               context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@common_options
def install(output, log, **kwargs):
    """Install the Hecatomb databases"""

    merge_config = {
        'args': {
            'output': output,
            'log': log
        }
    }
    run_snakemake(
        snakefile_path=snake_base(os.path.join('snakemake', 'workflow', 'DownloadDB.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(epilog=help_msg_epilog(),
               context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--comb', multiple=True, required=True, show_default=False,
              help='Two or more Hecatomb output directories to combine. e.g. --comb dir1/ --comb dir2/ ...')
@common_options
def combine(comb, output, log, **kwargs):
    """Combine multiple Hecatomb runs"""

    merge_config = {
        'args': {
            'output': output,
            'combineRuns': list(comb),
            'log': log
        }
    }
    run_snakemake(
        snakefile_path=snake_base(os.path.join('snakemake', 'workflow', 'combineOutputs.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(epilog=help_msg_epilog(),
               context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@click.option('--host', help='Name for your host genome', show_default=False, required=True)
@click.option('--host-fa', help='Host genome fasta file', show_default=False, required=True)
@common_options
def add_host(host, host_fa, output, log, **kwargs):
    """Add a new host genome to use with Hecatomb"""

    merge_config = {
        'args': {
            'output': output,
            'hostFa': host_fa,
            'hostName': host,
            'log': log
        }
    }
    run_snakemake(
        snakefile_path=snake_base(os.path.join('snakemake', 'workflow', 'AddHost.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command()
@common_options
def list_hosts(configfile, **kwargs):
    """List the available host genomes"""

    run_list_hosts(configfile)


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
