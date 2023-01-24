"""
Entrypoint for {{cookiecutter.project_name}}

Check out the wiki for a detailed look at customising this file:
https://github.com/beardymcjohnface/Snaketool/wiki/Customising-your-Snaketool
"""

import os
import click

from .util import (
    snake_base,
    get_version,
    default_to_output,
    copy_config,
    run_snakemake,
    OrderedCommands,
    print_citation,
    run_list_hosts,
)


def print_splash():
    click.echo("""
\b
██╗  ██╗███████╗ ██████╗ █████╗ ████████╗ ██████╗ ███╗   ███╗██████╗
██║  ██║██╔════╝██╔════╝██╔══██╗╚══██╔══╝██╔═══██╗████╗ ████║██╔══██╗
███████║█████╗  ██║     ███████║   ██║   ██║   ██║██╔████╔██║██████╔╝
██╔══██║██╔══╝  ██║     ██╔══██║   ██║   ██║   ██║██║╚██╔╝██║██╔══██╗
██║  ██║███████╗╚██████╗██║  ██║   ██║   ╚██████╔╝██║ ╚═╝ ██║██████╔╝
╚═╝  ╚═╝╚══════╝ ╚═════╝╚═╝  ╚═╝   ╚═╝    ╚═════╝ ╚═╝     ╚═╝╚═════╝
""", err=True)


def common_options(func):
    """Common command line args
    Define common command line args here, and include them with the @common_options decorator below.
    """
    options = [
        click.option(
            "--output",
            help="Output directory",
            type=click.Path(dir_okay=True, writable=True, readable=True),
            default="{{cookiecutter.project_slug}}.out",
            show_default=True,
        ),
        click.option(
            "--configfile",
            default="config.yaml",
            show_default=False,
            callback=default_to_output,
            help="Custom config file [default: (outputDir)/config.yaml]",
        ),
        click.option(
            "--threads", help="Number of threads to use", default=1, show_default=True
        ),
        click.option(
            "--use-conda/--no-use-conda",
            default=True,
            help="Use conda for Snakemake rules",
            show_default=True,
        ),
        click.option(
            "--conda-prefix",
            default=snake_base(os.path.join("workflow", "conda")),
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
            default="{{cookiecutter.project_slug}}.log",
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
    """{{cookiecutter.project_description}}
    \b
    For more options, run:
    {{cookiecutter.project_slug}} command --help"""
    pass


help_msg_extra = """
\b
CLUSTER EXECUTION:
hecatomb run ... --profile [profile]
For information on Snakemake profiles see:
https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles
\b
RUN EXAMPLES:
Required:           hecatomb run --reads [file/dir]
Specify threads:    hecatomb run ... --threads [threads]
Disable conda:      hecatomb run ... --no-use-conda 
Change defaults:    hecatomb run ... --snake-default="-k --nolock"
Add Snakemake args: hecatomb run ... --dry-run --keep-going --touch
Specify stages:     hecatomb run ... all print_targets
Available stages:
    all                 Run everything (default)
    preprocessing       Preprocessing steps only
    assembly            Assembly steps (+ preprocessing)
    annotations         Read annotations (+ preprocessing)
    ctg_annotations     Contig annotations (+ preprocessing,assembly)
    print_stages        List available stages
"""


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option('--reads', 'reads', help='Input file/directory', type=str, default=None, required=True)
@click.option('--preprocess', help='Preprocessing method', default='paired', show_default=True,
              type=click.Choice(['paired', 'single', 'longread', 'roundAB']))
@click.option('--search', help='MMSeqs search speed settings', default='sensitive',
              type=click.Choice(['fast', 'sensitive']), show_default=True)
@click.option('--host', help='Host genome name for filtering', default='human', show_default=True)
@common_options
def run(reads, preprocess, search, host, output, log, **kwargs):
    """Run hecatomb"""

    merge_config = {
        'args': {
            'reads': reads,
            'output': output,
            'host': host,
            'preprocess': preprocess,
            'search': search,
            'log': log
        }
    }

    run_snakemake(
        # Full path to Snakefile
        snakefile_path=snake_base(os.path.join('', '../snakemake', 'workflow', 'Hecatomb.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(
    epilog=help_msg_extra,
    context_settings=dict(
        help_option_names=["-h", "--help"], ignore_unknown_options=True
    ),
)
@click.option('--preprocess', help='Preprocessing method', default='paired', show_default=True,
              type=click.Choice(['paired', 'single', 'longread', 'roundAB']))
@click.option('--search', help='MMSeqs search speed settings', default='sensitive',
              type=click.Choice(['fast', 'sensitive']), show_default=True)
@click.option('--host', help='Host genome name for filtering', default='human', show_default=True)
@common_options
def test(preprocess, search, host, output, log, **kwargs):
    """Run the Hecatomb test dataset"""

    reads = snake_base(os.path.join('', '../test_data'))

    merge_config = {
        'args': {
            'reads': reads,
            'output': output,
            'host': host,
            'preprocess': preprocess,
            'search': search,
            'log': log
        }
    }

    run_snakemake(
        snakefile_path=snake_base(os.path.join('', '../snakemake', 'workflow', 'Hecatomb.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
@common_options
def config(configfile, **kwargs):
    """Copy the system default config file"""
    copy_config(configfile)


@click.command(context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
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
        snakefile_path=snake_base(os.path.join('', '../snakemake', 'workflow', 'DownloadDB.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
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
        snakefile_path=snake_base(os.path.join('', '../snakemake', 'workflow', 'combineOutputs.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command(context_settings=dict(help_option_names=["-h", "--help"], ignore_unknown_options=True))
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
        snakefile_path=snake_base(os.path.join('', '../snakemake', 'workflow', 'AddHost.smk')),
        merge_config=merge_config,
        log=log,
        **kwargs
    )


@click.command()
@common_options
def list_hosts(configfile, **kwargs):
    """List the available host genomes"""

    run_list_hosts()


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
