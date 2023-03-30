"""MISC FUNCTIONS
You shouldn't need to tweak these much if at all
"""

import sys
import os
import subprocess
import yaml
import glob
import click
import collections.abc
from shutil import copyfile
from time import localtime, strftime


class OrderedCommands(click.Group):
    """This class will preserve the order of subcommands, which is useful when printing --help"""

    def list_commands(self, ctx: click.Context):
        return list(self.commands)


def snake_base(rel_path):
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), rel_path)


def get_version():
    with open(snake_base("hecatomb.VERSION"), "r") as f:
        version = f.readline()
    return version


def echo_click(_msg, log=None):
    click.echo(_msg, nl=False, err=True)
    if log:
        with open(log, "a") as lfh:
            lfh.write(_msg)


def print_citation():
    with open(snake_base("hecatomb.CITATION"), "r") as f:
        for line in f:
            echo_click(line)


def msg(err_message, log=None):
    tstamp = strftime("[%Y:%m:%d %H:%M:%S] ", localtime())
    echo_click(tstamp + err_message + "\n", log=log)


def msg_box(splash, errmsg=None, log=None):
    msg("-" * (len(splash) + 4), log=log)
    msg(f"| {splash} |", log=log)
    msg(("-" * (len(splash) + 4)), log=log)
    if errmsg:
        echo_click("\n" + errmsg + "\n", log=log)


def default_to_output(ctx, param, value):
    """Callback for --configfile; place value in output directory unless specified"""
    if param.default == value:
        return os.path.join(ctx.params["output"], value)
    return value


def read_config(file):
    with open(file, "r") as stream:
        _config = yaml.safe_load(stream)
    return _config


def recursive_merge_config(config, overwrite_config):
    def _update(d, u):
        for (key, value) in u.items():
            if isinstance(value, collections.abc.Mapping):
                d[key] = _update(d.get(key, {}), value)
            else:
                d[key] = value
        return d
    _update(config, overwrite_config)


def update_config(in_config=None, merge=None, output_config=None, log=None):
    """Update config with new values"""
    if output_config is None:
        output_config = in_config
    config = read_config(in_config)
    msg("Updating config file with new values", log=log)
    recursive_merge_config(config, merge)
    write_config(config, output_config, log=log)


def write_config(_config, file, log=None):
    msg(f"Writing config file to {file}", log=log)
    with open(file, "w") as stream:
        yaml.dump(_config, stream)


def copy_config(
    local_config,
    merge_config=None,
    system_config=snake_base(os.path.join("snakemake", "config", "config.yaml")),
    log=None,
):
    if not os.path.isfile(local_config):
        if len(os.path.dirname(local_config)) > 0:
            os.makedirs(os.path.dirname(local_config), exist_ok=True)
        msg(f"Copying system default config to {local_config}", log=log)

        if merge_config:
            update_config(
                in_config=system_config,
                merge=merge_config,
                output_config=local_config,
                log=log,
            )
        else:
            copyfile(system_config, local_config)
    else:
        msg(
            f"Config file {local_config} already exists. Using existing config file.",
            log=log,
        )


"""RUN A SNAKEFILE
Hopefully you shouldn't need to tweak this function at all.
- You must provide a Snakefile, all else is optional
- Highly recommend supplying a configfile and the default snakemake args"""


def run_snakemake(
    configfile=None,
    snakefile_path=None,
    merge_config=None,
    threads=1,
    use_conda=False,
    conda_prefix=None,
    snake_default=None,
    snake_args=None,
    log=None,
):
    """Run a Snakefile"""
    snake_command = ["snakemake", "-s", snakefile_path]

    # if using a configfile
    if configfile:
        # copy sys default config if not present
        copy_config(configfile, log=log)

        if merge_config:
            update_config(in_config=configfile, merge=merge_config, log=log)

        snake_command += ["--configfile", configfile]

        # display the runtime configuration
        snake_config = read_config(configfile)
        msg_box(
            "Runtime config",
            errmsg=yaml.dump(snake_config, Dumper=yaml.Dumper),
            log=log,
        )

    # add threads
    if "--profile" not in snake_args:
        snake_command += ["--jobs", threads]

    # add conda args if using conda
    if use_conda:
        snake_command += ["--use-conda"]
        if conda_prefix:
            snake_command += ["--conda-prefix", conda_prefix]

    # add snakemake default args
    if snake_default:
        snake_command += snake_default

    # add any additional snakemake commands
    if snake_args:
        snake_command += list(snake_args)

    # Run Snakemake!!!
    snake_command = " ".join(str(s) for s in snake_command)
    msg_box("Snakemake command", errmsg=snake_command, log=log)
    if not subprocess.run(snake_command, shell=True).returncode == 0:
        msg("ERROR: Snakemake failed", log=log)
        sys.exit(1)
    else:
        msg("Snakemake finished successfully", log=log)
    return 0


def run_list_hosts(configfile):
    copy_config(configfile)
    with open(configfile, 'r') as f:
        config = yaml.safe_load(f)
    dbdir = snake_base(os.path.join('snakemake', 'databases'))
    try:
        if config['Databases'] is not None:
            dbdir = config['Databases']
    except KeyError:
        pass
    host_path = os.path.join(dbdir, "host", "*")
    host_fastas = list([os.path.basename(x) for x in glob.glob(host_path)])
    try:
        host_fastas.remove('virus_shred.fasta.gz')
    except ValueError:
        pass
    msg_box('Available host genomes', '\n'.join(host_fastas))
