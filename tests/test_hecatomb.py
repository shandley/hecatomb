import subprocess
from pathlib import Path

import pytest


def exec_command(cmnd, stdout=subprocess.PIPE, stderr=subprocess.PIPE):
    """executes shell command and returns stdout if completes exit code 0
    Parameters
    ----------
    cmnd : str
      shell command to be executed
    stdout, stderr : streams
      Default value (PIPE) intercepts process output, setting to None
      blocks this."""

    proc = subprocess.Popen(cmnd, shell=True, stdout=stdout, stderr=stderr)
    out, err = proc.communicate()
    if proc.returncode != 0:
        raise RuntimeError(f"FAILED: {cmnd}\n{err}")
    return out.decode("utf8") if out is not None else None


def test_hecatomb_cli():
    exec_command("hecatomb -v")
    exec_command("hecatomb -h")
    exec_command("hecatomb run -h")
    exec_command("hecatomb config -h")


def test_koverage_commands(tmp_path):
    temp_dir = tmp_path / "tmp"
    temp_dir.mkdir()
    temp_file = temp_dir / "yeet.txt"
    temp_file.write_text("yeet")
    exec_command("hecatomb install -n")
    exec_command("hecatomb config")
    exec_command("hecatomb citation")
    simulate = ["hecatomb test --simulate -n "]
    hpc = ["--example-profile "]
    options = ["--assembly cross --fastqc --search fast --custom-aa ", str(temp_file), " --custom-nt ", str(temp_file)]
    exec_command(" ".join(simulate))
    exec_command(" ".join(simulate + hpc))
    exec_command(" ".join(simulate + options))
    exec_command(" ".join(simulate + hpc + options))
