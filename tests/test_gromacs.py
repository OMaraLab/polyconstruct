from polytop import *
import subprocess
import os


def test_run_gromacs():
    """
    Note you have to create a docker image of gromacs/gromacs.
    If you are running windows you need to install Docker Desktop first
    https://www.docker.com/products/docker-desktop/

    You then download the DOCKERFILE details for GROMACS by typing the following commands in a terminal:

    git clone https://github.com/bioexcel/gromacs-docker

    And then build the image by typing:
    cd gromacs-docker
    docker build .

    subsequent calls to "docker run --rm -it gromacs/gromacs" will use the image you just built
    and give you an interactive terminal, and delete the container when you exit.
    """
    docker_commands = [
        "docker",
        "run",
        "--rm",
        "-v",
        "./samples:/data",
        "-w",
        "/data",
        "gromacs/gromacs",
    ]
    # commands = docker_commands + ["gmx", "grompp", "-f", "md.mdp", "-c", "conf.gro", "-p", "topol.top", "-o", "md.tpr"]
    commands = docker_commands + ["gmx", "help", "commands"]
    process = subprocess.run(commands, capture_output=True, text=True)
    err = process.returncode
    assert err == 0, process.stderr


def test_samples_exist():
    # get current director
    dir = os.getcwd()
    assert os.path.exists(f"{dir}/tests/samples/")


def test_ARG_file_exists():
    # get current director
    dir = os.getcwd()
    assert os.path.exists(f"{dir}/tests/samples/arginine.itp")


def test_GLU_file_exists():
    # get current director
    dir = os.getcwd()
    assert os.path.exists(f"{dir}/tests/samples/glutamine.itp")
