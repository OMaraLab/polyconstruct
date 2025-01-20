import subprocess
from pathlib import Path

class Gromacs:
    """
    Runs GROMACS commands in a docker container.
    
    Note you have to create a docker image of gromacs/gromacs.
    If you are running windows you need to install Docker Desktop first
    https://www.docker.com/products/docker-desktop/

    You then download the DOCKERFILE details for GROMACS by typing the following commands 
    in a terminal to use the canonical Gromacs container:

    git clone https://github.com/bioexcel/gromacs-docker

    And then build the image by typing:
    cd gromacs-docker
    docker build .

    subsequent calls to "docker run --rm -it gromacs/gromacs" will use the image you just built
    and give you an interactive terminal, and delete the container when you exit.
    """    
    
    def __init__(self, data_dir: Path, output_dir: Path):
        self.data_dir = data_dir
        self.output_dir = output_dir

    def run(self, gromacs_commands: list[str], output_file_name: str):
        docker_commands = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{self.data_dir}:/data",
            "-v",
            f"{self.output_dir}:/output",
            "-w",
            "/data",
            "gromacs/gromacs",
        ]
        commands = docker_commands + gromacs_commands
        process = subprocess.run(commands, capture_output=True, text=True)
        err = process.returncode
        assert err == 0, process.stderr

        output_file = self.output_dir / output_file_name
        output_file.write_text(process.stdout)
        
    def run_simulation(self, mdp_file: str, gro_file: str, top_file: str, tpr_file: str, block: bool = True):
        """
        Runs a GROMACS simulation using the given mdp, gro, and top files.

        :param mdp_file: the name of the parameter file
        :type mdp_file: str
        :param gro_file: the name of the molecular structure file
        :type gro_file: str
        :param top_file: the name of the topology file
        :type top_file: str
        :param tpr_file: the name of the portable binary run input file
        :type tpr_file: str
        :param block: if True, then wait until the simulation has finished,
                otherwise run in the background, defaults to True
        :type block: bool, optional
        """
        gromacs_commands = ["gmx", "grompp", "-f", f"/data/{mdp_file}", "-c", f"/data/{gro_file}", "-p", f"/data/{top_file}", "-o", f"/output/{tpr_file}"]
        docker_commands = [
            "docker",
            "run",
            "--rm",
            "-v",
            f"{self.data_dir}:/data",
            "-v",
            f"{self.output_dir}:/output",
            "-w",
            "/data",
            "gromacs/gromacs",
        ]
        commands = docker_commands + gromacs_commands

        if block:
            process = subprocess.run(commands, capture_output=True, text=True)
            err = process.returncode
            assert err == 0, process.stderr
        else:
            subprocess.Popen(commands)